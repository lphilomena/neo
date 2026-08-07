#!/usr/bin/env python3
"""Repair the duplicate SQANTI3 write_isoform_hits definitions."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sqanti-home", required=True)
    args = parser.parse_args()

    path = Path(args.sqanti_home) / "src/qc_output.py"
    original = path.read_text()
    updated = original

    old_finalizer = "def write_isoform_hits(outdir,prefix, isoforms_info):"
    new_finalizer = "def finalize_isoform_hits(outdir, prefix, isoforms_info):"
    if old_finalizer in updated:
        updated = updated.replace(old_finalizer, new_finalizer, 1)

    old_appender = "def write_isoform_hits(isoform_hits_name,data_list):\n"
    new_dispatch = (
        "def write_isoform_hits(*args):\n"
        "    # Classification writes one hit at a time; the pipeline finalizes all hits.\n"
        "    if len(args) == 3:\n"
        "        return finalize_isoform_hits(*args)\n"
        "    if len(args) != 2:\n"
        "        raise TypeError(f\"write_isoform_hits expected 2 or 3 arguments, got {len(args)}\")\n"
        "    isoform_hits_name, data_list = args\n"
    )
    if old_appender in updated:
        updated = updated.replace(old_appender, new_dispatch, 1)

    if "def finalize_isoform_hits(" not in updated or "if len(args) == 3:" not in updated:
        raise SystemExit("Unsupported SQANTI3 qc_output.py layout")

    changed = updated != original
    if changed:
        backup = path.with_suffix(path.suffix + ".neoag.bak")
        if not backup.exists():
            backup.write_text(original)
        path.write_text(updated)
    print(f"sqanti3_isoform_hits_api={'patched' if changed else 'already_ok'}")


if __name__ == "__main__":
    main()
