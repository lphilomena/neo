#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import platform
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run(cmd: list[str], *, env: dict[str, str] | None = None, timeout: int = 60) -> dict[str, object]:
    try:
        proc = subprocess.run(cmd, cwd=ROOT, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)
        return {"cmd": cmd, "returncode": proc.returncode, "stdout": proc.stdout.strip(), "stderr": proc.stderr.strip()}
    except Exception as exc:
        return {"cmd": cmd, "returncode": None, "stdout": "", "stderr": str(exc)}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def tool_env() -> dict[str, str]:
    cmd = "source conf/tools.env.sh >/dev/null 2>&1; env"
    proc = subprocess.run(["bash", "-lc", cmd], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    env = os.environ.copy()
    if proc.returncode == 0:
        for line in proc.stdout.splitlines():
            if "=" in line:
                key, value = line.split("=", 1)
                env[key] = value
    return env


def main() -> None:
    parser = argparse.ArgumentParser(description="Write a NeoAg release reproducibility manifest.")
    parser.add_argument("--out", default="release_manifest.json")
    parser.add_argument("--skip-pytest", action="store_true")
    parser.add_argument(
        "--release-dir",
        action="append",
        default=["dist", "work/releases"],
        help="Directory to scan for release tarballs; may be repeated.",
    )
    args = parser.parse_args()

    env = tool_env()
    manifest: dict[str, object] = {
        "created_at": datetime.datetime.now(datetime.UTC).isoformat().replace("+00:00", "Z"),
        "project_root": str(ROOT),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "git_head": run(["git", "rev-parse", "HEAD"]),
        "git_status_short": run(["git", "status", "--short"], timeout=120),
        "pytest": None if args.skip_pytest else run([str(ROOT / ".venv" / "bin" / "python"), "-m", "pytest", "-q"], env=env, timeout=300),
        "neoag_help": run([str(ROOT / "bin" / "neoag-v03"), "--help"], timeout=60),
        "java_version": run(["java", "-version"], env=env),
        "nextflow_version": run([str(ROOT / "bin" / "nextflow"), "-version"], env=env),
        "tool_check": run([str(ROOT / "bin" / "neoag-v03"), "check-tools"], env=env, timeout=180),
        "release_files": [],
        "known_optional_missing_tools": ["netmhcpan", "prime", "star_fusion", "fusioncatcher"],
    }
    files = []
    seen: set[Path] = set()
    for rel_dir in args.release_dir:
        release_dir = ROOT / rel_dir
        if not release_dir.is_dir():
            continue
        for path in sorted(release_dir.glob("*.tar.gz")):
            path = path.resolve()
            if path in seen:
                continue
            seen.add(path)
            files.append({"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    manifest["release_files"] = files
    out = ROOT / args.out
    out.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    print(out)


if __name__ == "__main__":
    main()
