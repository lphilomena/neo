from __future__ import annotations

import argparse
import os
import re
from pathlib import Path
from typing import Any

from .io_utils import ensure_dir, markdown_table, now_iso, safe_rel, write_json, write_tsv

CACHE_PATTERNS = ["__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache", ".nextflow", "work", "results"]
DEFAULT_SKIP_DIRS = {".git", "__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache", ".nextflow", "work", "results"}
PRIVATE_PATH_PATTERNS = [
    re.compile(r"/mnt/[A-Za-z0-9_\-./]+"),
    re.compile(r"/home/[A-Za-z0-9_\-./]+"),
    re.compile(r"/root/[A-Za-z0-9_\-./]+"),
    re.compile(r"[A-Za-z]:\\\\"),
]
PATIENT_HINTS = ["patient", "患者", "肿瘤", "血液", "M1ML", "ML150", "chenxiaoliang", "陈"]
TEXT_SUFFIXES = {".py", ".md", ".txt", ".toml", ".yaml", ".yml", ".json", ".sh", ".nf", ".html", ".rst", ".cfg", ".ini"}


def scan_release_boundary(root: str | Path, *, max_file_bytes: int = 2_000_000, skip_dirs: set[str] | None = None) -> dict[str, Any]:
    base = Path(root)
    skip = DEFAULT_SKIP_DIRS if skip_dirs is None else set(skip_dirs)
    cache_hits: list[dict[str, str]] = []
    skipped_dirs: list[dict[str, str]] = []
    private_hits: list[dict[str, str]] = []
    patient_hits: list[dict[str, str]] = []

    for root_dir, dirs, files in os.walk(base):
        current = Path(root_dir)
        kept_dirs: list[str] = []
        for d in dirs:
            dp = current / d
            rel = safe_rel(dp, base)
            if d in CACHE_PATTERNS:
                cache_hits.append({"path": rel, "type": "cache_or_runtime_artifact"})
            if d in skip:
                skipped_dirs.append({"path": rel, "reason": "skipped_large_or_generated_dir"})
            else:
                kept_dirs.append(d)
        dirs[:] = kept_dirs

        for fname in files:
            p = current / fname
            rel = safe_rel(p, base)
            if p.name in CACHE_PATTERNS:
                cache_hits.append({"path": rel, "type": "cache_or_runtime_artifact"})
            if p.suffix.lower() not in TEXT_SUFFIXES:
                continue
            try:
                size = p.stat().st_size
            except OSError:
                continue
            if size > max_file_bytes:
                continue
            try:
                text = p.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                continue
            for pat in PRIVATE_PATH_PATTERNS:
                m = pat.search(text)
                if m:
                    private_hits.append({"path": rel, "match": m.group(0)[:200]})
                    break
            low = text.lower()
            if any(h.lower() in low for h in PATIENT_HINTS):
                patient_hits.append({"path": rel, "match": "patient_or_site_hint"})

    status = "PASS"
    if private_hits or patient_hits:
        status = "UNSAFE"
    elif cache_hits:
        status = "PARTIAL"
    return {
        "status": status,
        "root": str(base),
        "scanned_at": now_iso(),
        "cache_hits": cache_hits,
        "skipped_dirs": skipped_dirs,
        "private_path_hits": private_hits,
        "patient_hint_hits": patient_hits,
        "summary": {
            "cache_hits": len(cache_hits),
            "skipped_dirs": len(skipped_dirs),
            "private_path_hits": len(private_hits),
            "patient_hint_hits": len(patient_hits),
        },
    }


def write_release_audit(result: dict[str, Any], outdir: str | Path) -> dict[str, str]:
    od = ensure_dir(outdir)
    write_json(od / "release_audit.json", result)
    write_tsv(od / "release_cache_hits.tsv", result.get("cache_hits", []))
    write_tsv(od / "release_skipped_dirs.tsv", result.get("skipped_dirs", []))
    write_tsv(od / "release_private_path_hits.tsv", result.get("private_path_hits", []))
    write_tsv(od / "release_patient_hint_hits.tsv", result.get("patient_hint_hits", []))
    md = [
        "# NeoAg release boundary audit",
        "",
        f"Status: **{result.get('status')}**",
        "",
        "## Summary",
        markdown_table([result.get("summary", {})]),
        "",
        "## Private path hits",
        markdown_table(result.get("private_path_hits", []), max_rows=20),
        "",
        "## Patient/site hint hits",
        markdown_table(result.get("patient_hint_hits", []), max_rows=20),
        "",
        "## Cache/runtime artifacts",
        markdown_table(result.get("cache_hits", []), max_rows=20),
        "",
        "## Skipped large/generated directories",
        markdown_table(result.get("skipped_dirs", []), max_rows=20),
    ]
    (od / "release_audit.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    return {"release_audit_json": str(od / "release_audit.json"), "release_audit_md": str(od / "release_audit.md")}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="controlled-execution Phase 0 release boundary audit")
    ap.add_argument("--root", default=".")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--scan-generated-dirs", action="store_true", help="Recursively scan work/results/cache directories instead of recording and skipping them")
    args = ap.parse_args(argv)
    result = scan_release_boundary(args.root, skip_dirs=set() if args.scan_generated_dirs else None)
    outs = write_release_audit(result, args.outdir)
    print(f"Release audit status: {result['status']}")
    for k, v in outs.items():
        print(f"  {k}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
