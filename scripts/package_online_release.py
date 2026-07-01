#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTDIR = ROOT / "work" / "releases"
VERSION = "v0.4.3"

INCLUDE_DIRS = [
    "assets",
    "bin",
    "conda",
    "conf",
    "docs",
    "modules",
    "profiles",
    "resources",
    "scripts",
    "src",
    "tests",
    "workflows",
]
INCLUDE_DATA_DIRS = [
    "data/fixtures",
    "data/fixtures_snv",
    "data/fixtures_sv",
    "data/improve",
]
INCLUDE_FILES = [
    ".gitignore",
    "CHANGELOG_V04_EVIDENCE_SAFETY_ESCAPE.md",
    "CHANGELOG_V041_APPM_CCF_IMMUNE_ESCAPE.md",
    "CHANGELOG_V042_P1_APPM_EXPLAINABILITY.md",
    "CHANGELOG_V043_CCF21.md",
    "CITATION.cff",
    "LICENSE",
    "NOTICE",
    "README.md",
    "README_zh.md",
    "README_PATCH_V043_CCF21.md",
    "RELEASE.md",
    "RELEASE_REFRESH_20260615.md",
    "RELEASE_V04_20260616.md",
    "pyproject.toml",
]

EXCLUDE_TOP_LEVEL = {
    ".git",
    ".nextflow",
    ".nextflow_user",
    ".pytest_cache",
    ".venv",
    ".venv.local",
    "conda_packs",
    "dist",
    "results",
    "tools",
    "work",
}
EXCLUDE_PARTS = {"__pycache__", ".pytest_cache"}
EXCLUDE_SUFFIXES = {
    ".pyc",
    ".pyo",
    ".docx",
    ".pptx",
    ".log",
    ".tmp",
}
EXCLUDE_NAMES = {
    "0",
    "FusionInspector.log",
    "tools.env.local.sh",
    "snp-pileup",
    "migrate_export.log",
    "wget_ad.log",
    "human_v102.tar.gz.ad",
    "neoantigen2-main.zip",
}
EXCLUDE_PREFIXES = (
    ".nextflow.log",
    ".nextflow.root_owned_backup_",
)
EXCLUDE_REL_PREFIXES = (
    "data/examples/",
    "data/external/",
    "data/ref/",
    "data/vep/",
    "resources/arriba_v2.5.1_official/",
    "conf/private/",
    "docs/release_audit/",
)
EXCLUDE_REL_SUFFIXES = (
    ".private.toml",
    ".local.toml",
    ".patient.toml",
)
EXCLUDE_NAME_TOKENS = (
    "chenxiaoliang",
    "hcc1395",
    "hcc1143",
    "ml150",
    "patient",
)


def is_allowed_conf(path: Path) -> bool:
    rel = path.relative_to(ROOT).as_posix()
    if not rel.startswith("conf/"):
        return True
    name = path.name
    return (
        name.endswith(".example")
        or name.endswith(".example.sh")
        or name.endswith(".example.toml")
        or name in {
            "demo.config",
            "local.config",
            "pip-constraints-tools.txt",
            "run.private.example.toml",
            "run.snv_wes.example.toml",
            "run.snv_wes_fixture.toml",
            "run.stub.toml",
            "run.sv_wgs_phase1.example.toml",
            "site.config.example",
            "snv_wes_demo.config",
            "sv_demo.config",
            "sv_wes_demo.config",
            "tools.config",
            "tools.env.local.example.sh",
            "tools.env.sh",
        }
    )


def is_allowed_doc(path: Path) -> bool:
    rel = path.relative_to(ROOT).as_posix()
    if not rel.startswith("docs/"):
        return True
    if rel.startswith("docs/release_audit/"):
        return False
    return path.suffix == ".md"


def is_excluded(path: Path) -> bool:
    rel = path.relative_to(ROOT)
    rel_posix = rel.as_posix()
    if rel.parts and rel.parts[0] in EXCLUDE_TOP_LEVEL:
        return True
    if set(rel.parts) & EXCLUDE_PARTS:
        return True
    if path.name in EXCLUDE_NAMES:
        return True
    if any(token in rel_posix.lower() for token in EXCLUDE_NAME_TOKENS):
        return True
    if path.name.startswith("._"):
        return True
    if any(path.name.startswith(prefix) for prefix in EXCLUDE_PREFIXES):
        return True
    if any(rel_posix.startswith(prefix) for prefix in EXCLUDE_REL_PREFIXES):
        return True
    if any(rel_posix.endswith(suffix) for suffix in EXCLUDE_REL_SUFFIXES):
        return True
    if any(path.name.endswith(marker) for marker in (".bak",)):
        return True
    if ".bak_" in path.name or ".backup_" in path.name:
        return True
    if path.suffix in EXCLUDE_SUFFIXES:
        return True
    if not is_allowed_conf(path):
        return True
    if not is_allowed_doc(path):
        return True
    return False


def iter_release_files() -> list[Path]:
    files: set[Path] = set()
    for rel_dir in INCLUDE_DIRS + INCLUDE_DATA_DIRS:
        base = ROOT / rel_dir
        if not base.exists():
            continue
        for path in base.rglob("*"):
            if path.is_file() and not is_excluded(path):
                files.add(path)
    for rel_file in INCLUDE_FILES:
        path = ROOT / rel_file
        if path.is_file() and not is_excluded(path):
            files.add(path)
    return sorted(files, key=lambda p: p.relative_to(ROOT).as_posix())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_manifest(path: Path, package_name: str, files: list[Path], archive_sha256: str, archive_size: int) -> None:
    top_level_counts: dict[str, int] = {}
    for file_path in files:
        rel = file_path.relative_to(ROOT).as_posix()
        top = rel.split("/", 1)[0]
        top_level_counts[top] = top_level_counts.get(top, 0) + 1
    manifest = {
        "schema_version": "neoag-online-release-manifest-v1",
        "project": "neoag_event_pipeline",
        "version": VERSION,
        "package": package_name,
        "created_at_utc": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
        "archive": {
            "size_bytes": archive_size,
            "sha256": archive_sha256,
            "sha256_file": f"{package_name}.sha256",
        },
        "file_count": len(files),
        "top_level_counts": dict(sorted(top_level_counts.items())),
        "included": {
            "source": ["src", "bin", "scripts"],
            "workflow": ["workflows", "modules"],
            "configuration": ["conf", "profiles", "conda"],
            "fixtures": INCLUDE_DATA_DIRS,
            "documentation": ["README.md", "README_zh.md", "docs/*.md", "RELEASE.md"],
        },
        "excluded": sorted(EXCLUDE_TOP_LEVEL | {"data/ref", "data/vep", "data/external", "data/examples", "patient/local scripts", "site/private configs"}),
        "verification_hint": [
            "python -m pip install -e '.[test]'",
            "pytest -q",
            "neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001",
            "bin/neoag-nextflow run workflows/main.nf --pvac_files data/fixtures/pvacseq_aggregated.tsv --outdir results/demo_nf",
        ],
        "interpretation_boundary": "Research triage and validation planning only; not clinical treatment recommendation software.",
    }
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build public online release package.")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR))
    parser.add_argument("--name", default=None)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now(dt.UTC).strftime("%Y%m%d")
    package_name = args.name or f"neoag_event_pipeline_v043_online_{stamp}.tar.gz"
    archive_path = outdir / package_name
    root_name = package_name.removesuffix(".tar.gz")
    files = iter_release_files()

    with tarfile.open(archive_path, "w:gz") as tar:
        for file_path in files:
            tar.add(file_path, arcname=str(Path(root_name) / file_path.relative_to(ROOT)))

    archive_sha256 = sha256(archive_path)
    sha_path = archive_path.with_suffix(archive_path.suffix + ".sha256")
    sha_path.write_text(f"{archive_sha256}  {archive_path.name}\n", encoding="utf-8")

    manifest_path = outdir / f"{root_name}.manifest.json"
    write_manifest(manifest_path, archive_path.name, files, archive_sha256, archive_path.stat().st_size)

    print(archive_path)
    print(f"files={len(files)} size={archive_path.stat().st_size} sha256={archive_sha256}")
    print(manifest_path)


if __name__ == "__main__":
    main()
