#!/usr/bin/env python3
"""Validate the machine-readable NeoAg reference manifest.

The parser intentionally supports the small YAML subset used by
configs/references/reference_manifest.yaml so this check has no PyYAML runtime
requirement on a fresh machine.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
from pathlib import Path
from typing import Any


def parse_scalar(value: str) -> Any:
    value = value.strip()
    if value == "":
        return ""
    if (value.startswith(chr(34)) and value.endswith(chr(34))) or (value.startswith(chr(39)) and value.endswith(chr(39))):
        return value[1:-1]
    low = value.lower()
    if low in {"true", "yes"}:
        return True
    if low in {"false", "no"}:
        return False
    if low in {"null", "none"}:
        return None
    try:
        return int(value)
    except ValueError:
        return value


def strip_comment(line: str) -> str:
    in_single = False
    in_double = False
    out = []
    for ch in line:
        if ch == chr(39) and not in_double:
            in_single = not in_single
        elif ch == chr(34) and not in_single:
            in_double = not in_double
        elif ch == "#" and not in_single and not in_double:
            break
        out.append(ch)
    return "".join(out).rstrip()


def load_manifest(path: Path) -> dict[str, Any]:
    root: dict[str, Any] = {}
    current_section: str | None = None
    current_item: str | None = None
    for lineno, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = strip_comment(raw)
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip(" "))
        text = line.strip()
        if ":" not in text:
            raise ValueError(f"line {lineno}: expected key: value")
        key, value = text.split(":", 1)
        key = key.strip()
        value = value.strip()
        if indent == 0:
            if value == "":
                root[key] = {}
                current_section = key
                current_item = None
            else:
                root[key] = parse_scalar(value)
                current_section = None
                current_item = None
        elif indent == 2 and current_section:
            section = root.setdefault(current_section, {})
            if not isinstance(section, dict):
                raise ValueError(f"line {lineno}: section {current_section} is not a mapping")
            if value == "":
                section[key] = {}
                current_item = key
            else:
                section[key] = parse_scalar(value)
        elif indent == 4 and current_section and current_item:
            item = root[current_section][current_item]
            if not isinstance(item, dict):
                raise ValueError(f"line {lineno}: item {current_item} is not a mapping")
            item[key] = parse_scalar(value)
        else:
            raise ValueError(f"line {lineno}: unsupported indentation level {indent}")
    return root


def truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() not in {"", "0", "false", "no", "optional", "none"}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def marker_path(path: Path, marker: str) -> Path:
    if marker.startswith("/"):
        return Path(marker)
    if marker.startswith(".") and not path.is_dir():
        return Path(str(path) + marker)
    return path / marker


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate NeoAg reference_manifest.yaml")
    parser.add_argument("manifest", nargs="?", default="configs/references/reference_manifest.yaml")
    parser.add_argument("--strict", action="store_true", help="Treat optional missing references as failures")
    parser.add_argument("--vep-version", default=os.environ.get("NEOAG_VEP_CACHE_VERSION", ""))
    args = parser.parse_args()

    manifest_path = Path(args.manifest)
    if not manifest_path.is_file():
        print(f"[FAIL] reference manifest missing: {manifest_path}", file=sys.stderr)
        return 2

    try:
        manifest = load_manifest(manifest_path)
    except Exception as exc:
        print(f"[FAIL] could not parse {manifest_path}: {exc}", file=sys.stderr)
        return 2

    references = manifest.get("references")
    if not isinstance(references, dict) or not references:
        print("[FAIL] manifest has no references mapping", file=sys.stderr)
        return 2

    print(f"==> Reference manifest: {manifest_path}")
    print(f"manifest_version={manifest.get(chr(109)+chr(97)+chr(110)+chr(105)+chr(102)+chr(101)+chr(115)+chr(116)+chr(95)+chr(118)+chr(101)+chr(114)+chr(115)+chr(105)+chr(111)+chr(110), chr(117)+chr(110)+chr(115)+chr(101)+chr(116))} genome_build={manifest.get(chr(103)+chr(101)+chr(110)+chr(111)+chr(109)+chr(101)+chr(95)+chr(98)+chr(117)+chr(105)+chr(108)+chr(100), chr(117)+chr(110)+chr(115)+chr(101)+chr(116))}")

    failed = False
    warned = False

    def issue(level: str, message: str) -> None:
        nonlocal failed, warned
        print(f"[{level}] {message}")
        if level == "FAIL":
            failed = True
        elif level == "WARN":
            warned = True

    for name, entry in references.items():
        if not isinstance(entry, dict):
            issue("FAIL", f"{name}: entry is not a mapping")
            continue
        required = truthy(entry.get("required", True))
        level_missing = "FAIL" if required or args.strict else "WARN"
        raw_path = entry.get("path")
        if not raw_path:
            issue("FAIL" if required else "WARN", f"{name}: missing path field")
            continue
        path = Path(str(raw_path))
        exists = path.exists()
        if exists:
            issue("OK", f"{name}: path exists: {path}")
        else:
            issue(level_missing, f"{name}: path missing: {path}")
            continue

        marker = str(entry.get("marker", "") or "").strip()
        if marker and marker != "-":
            mpath = marker_path(path, marker)
            if mpath.exists():
                issue("OK", f"{name}: marker exists: {mpath}")
            else:
                issue(level_missing, f"{name}: marker missing: {mpath}")

        expected_sha = str(entry.get("sha256", "") or "").strip()
        if expected_sha and expected_sha != "-":
            if not path.is_file():
                issue("FAIL" if required else "WARN", f"{name}: sha256 is only supported for files: {path}")
            else:
                actual = sha256_file(path)
                if actual == expected_sha:
                    issue("OK", f"{name}: sha256 matched")
                else:
                    issue("FAIL" if required else "WARN", f"{name}: sha256 mismatch expected={expected_sha} actual={actual}")

        required_vep = str(entry.get("vep_version_required", "") or "").strip()
        if required_vep:
            version = str(entry.get("version", "") or "")
            if args.vep_version and args.vep_version != required_vep:
                issue("FAIL" if required else "WARN", f"{name}: tool VEP version {args.vep_version} != manifest requires {required_vep}")
            if required_vep not in version and required_vep not in path.name:
                issue("FAIL" if required else "WARN", f"{name}: VEP cache path/version does not include required version {required_vep}")

    if failed:
        print("==> Reference manifest verification failed")
        return 1
    if warned:
        print("==> Reference manifest verification passed with warnings")
    else:
        print("==> Reference manifest verification passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
