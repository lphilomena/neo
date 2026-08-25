from __future__ import annotations

import hashlib
import io
import json
import shutil
import subprocess
import tarfile
import urllib.parse
import urllib.request
from pathlib import Path, PurePosixPath
from typing import Any

from neoag.controlled_execution.io_utils import now_iso, write_json


DEFAULT_PUBLIC_ASSET_REPO = "open-neo/open-neo-public-assets"
DEFAULT_PUBLIC_ASSET_REVISION = "main"
ARCHIVE_PREFIX = "open_neo_public_assets_openrefs_"
SUPPLEMENT_PREFIXES = (
    "data/rna/rsem_reference/",
    "data/splice/splicemutr/r_library/BSgenome.Hsapiens.UCSC.hg38/",
)
REQUIRED_SENTINELS = (
    "data/ref",
    "data/normal",
    "data/easyfuse",
)


class PublicAssetSyncError(RuntimeError):
    pass


class _MultipartReader(io.RawIOBase):
    def __init__(self, parts: list[Path]) -> None:
        self._parts = iter(parts)
        self._handle: io.BufferedReader | None = None

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: bytearray) -> int:
        view = memoryview(buffer)
        total = 0
        while total < len(view):
            if self._handle is None:
                try:
                    self._handle = next(self._parts).open("rb")
                except StopIteration:
                    break
            count = self._handle.readinto(view[total:])
            if count:
                total += count
                continue
            self._handle.close()
            self._handle = None
        return total

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None
        super().close()


def _repo_url(repo_id: str, revision: str, path: str) -> str:
    encoded_path = "/".join(urllib.parse.quote(part, safe="") for part in path.split("/"))
    return f"https://huggingface.co/datasets/{repo_id}/resolve/{revision}/{encoded_path}"


def _api_json(url: str) -> Any:
    request = urllib.request.Request(url, headers={"User-Agent": "open-neo-public-assets/1"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return json.load(response)
    except (OSError, ValueError) as exc:
        raise PublicAssetSyncError(f"Unable to query public asset repository: {url}: {exc}") from exc


def _repository_state(repo_id: str, revision: str) -> tuple[str, list[dict[str, Any]]]:
    info = _api_json(f"https://huggingface.co/api/datasets/{repo_id}")
    tree_url = (
        f"https://huggingface.co/api/datasets/{repo_id}/tree/"
        f"{urllib.parse.quote(revision, safe='')}?recursive=true&expand=false"
    )
    tree = _api_json(tree_url)
    if not isinstance(tree, list):
        raise PublicAssetSyncError(f"Unexpected repository tree response for {repo_id}@{revision}")
    return str(info.get("sha") or revision), [item for item in tree if item.get("type") == "file"]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download(
    repo_id: str,
    revision: str,
    remote_path: str,
    destination: Path,
    *,
    expected_size: int = 0,
    expected_sha256: str = "",
) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file():
        size_ok = not expected_size or destination.stat().st_size == expected_size
        hash_ok = not expected_sha256 or _sha256(destination) == expected_sha256
        if size_ok and hash_ok:
            return "REUSED"
        destination.unlink()
    curl = shutil.which("curl")
    if not curl:
        raise PublicAssetSyncError("curl is required for resumable public asset downloads")
    partial = destination.with_name(destination.name + ".partial")
    command = [
        curl, "-fL", "--retry", "8", "--retry-delay", "5", "--connect-timeout", "30",
        "-C", "-", "-o", str(partial), _repo_url(repo_id, revision, remote_path),
    ]
    proc = subprocess.run(command, text=True, capture_output=True, check=False)
    if proc.returncode != 0:
        raise PublicAssetSyncError(
            f"Download failed for {remote_path}: {(proc.stderr or proc.stdout).strip()[-1200:]}"
        )
    partial.replace(destination)
    if expected_size and destination.stat().st_size != expected_size:
        raise PublicAssetSyncError(
            f"Size mismatch for {remote_path}: {destination.stat().st_size} != {expected_size}"
        )
    if expected_sha256 and _sha256(destination) != expected_sha256:
        raise PublicAssetSyncError(f"SHA256 mismatch for {remote_path}")
    return "DOWNLOADED"


def _parse_checksums(path: Path) -> dict[str, str]:
    checksums: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.strip().split(None, 1)
        if len(fields) == 2:
            checksums[fields[1].lstrip("*")] = fields[0].lower()
    return checksums


def _safe_tar_member(member: tarfile.TarInfo) -> None:
    path = PurePosixPath(member.name)
    if path.is_absolute() or ".." in path.parts:
        raise PublicAssetSyncError(f"Unsafe archive member: {member.name}")
    if member.issym() or member.islnk():
        target = PurePosixPath(member.linkname)
        if target.is_absolute() or ".." in target.parts:
            raise PublicAssetSyncError(f"Unsafe archive link: {member.name} -> {member.linkname}")
    if member.isdev() or member.isfifo():
        raise PublicAssetSyncError(f"Unsupported archive member: {member.name}")


def _extract_split_archive(parts: list[Path], destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    reader = _MultipartReader(parts)
    try:
        with tarfile.open(fileobj=io.BufferedReader(reader, 8 * 1024 * 1024), mode="r|gz") as archive:
            for member in archive:
                _safe_tar_member(member)
                archive.extract(member, path=destination, filter="data")
    except (OSError, tarfile.TarError) as exc:
        raise PublicAssetSyncError(f"Unable to extract public asset archive: {exc}") from exc
    finally:
        reader.close()


def _local_ready(asset_root: Path) -> bool:
    marker = asset_root / ".open_neo_public_assets.json"
    return marker.is_file() and all((asset_root / item).exists() for item in REQUIRED_SENTINELS)


def sync_public_assets(
    asset_root: str | Path,
    *,
    cache_dir: str | Path,
    repo_id: str = DEFAULT_PUBLIC_ASSET_REPO,
    revision: str = DEFAULT_PUBLIC_ASSET_REVISION,
    execute: bool = False,
) -> dict[str, Any]:
    """Synchronize redistributable fixed assets from the public Dataset.

    ``asset_root`` is the deployment reference root; archive members are
    restored beneath ``asset_root/data``. Licensed/restricted groups are never
    requested because they are absent from the public package by design.
    """
    root = Path(asset_root).expanduser().resolve()
    cache = Path(cache_dir).expanduser().resolve()
    marker = root / ".open_neo_public_assets.json"
    if not execute:
        return {
            "status": "REUSED" if _local_ready(root) else "PLANNED",
            "repo_id": repo_id,
            "revision": revision,
            "asset_root": str(root),
            "cache_dir": str(cache),
            "marker": str(marker),
        }

    root.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    commit_sha, files = _repository_state(repo_id, revision)
    file_map = {str(item["path"]): item for item in files}
    part_names = sorted(
        path for path in file_map
        if Path(path).name.startswith(ARCHIVE_PREFIX) and ".tar.gz.part-" in path
    )
    if not part_names or "SHA256SUMS" not in file_map:
        raise PublicAssetSyncError(f"Dataset {repo_id}@{revision} has no complete split archive")

    checksum_path = cache / "SHA256SUMS"
    _download(repo_id, revision, "SHA256SUMS", checksum_path,
              expected_size=int(file_map["SHA256SUMS"].get("size") or 0))
    checksums = _parse_checksums(checksum_path)
    parts: list[Path] = []
    downloaded = 0
    reused = 0
    for name in part_names:
        destination = cache / Path(name).name
        action = _download(
            repo_id, revision, name, destination,
            expected_size=int(file_map[name].get("size") or 0),
            expected_sha256=checksums.get(Path(name).name, ""),
        )
        downloaded += action == "DOWNLOADED"
        reused += action == "REUSED"
        parts.append(destination)

    previous: dict[str, Any] = {}
    if marker.is_file():
        try:
            previous = json.loads(marker.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            previous = {}
    archive_id = hashlib.sha256("\n".join(checksums.get(path.name, "") for path in parts).encode()).hexdigest()
    extracted = previous.get("archive_id") != archive_id or not _local_ready(root)
    if extracted:
        _extract_split_archive(parts, root / "data")

    supplement_count = 0
    for remote_path, item in sorted(file_map.items()):
        if not any(remote_path.startswith(prefix) for prefix in SUPPLEMENT_PREFIXES):
            continue
        relative = Path(remote_path).relative_to("data")
        action = _download(
            repo_id, revision, remote_path, root / "data" / relative,
            expected_size=int(item.get("size") or 0),
        )
        supplement_count += action == "DOWNLOADED"

    payload = {
        "schema_version": "open-neo-public-assets-v1",
        "repo_id": repo_id,
        "revision": revision,
        "commit_sha": commit_sha,
        "archive_id": archive_id,
        "archive_parts": len(parts),
        "downloaded_parts": downloaded,
        "reused_parts": reused,
        "archive_extracted": extracted,
        "supplement_files_downloaded": supplement_count,
        "asset_root": str(root),
        "cache_dir": str(cache),
        "restricted_assets_included": False,
        "updated_at": now_iso(),
    }
    write_json(marker, payload)
    return {"status": "PASS", "marker": str(marker), **payload}
