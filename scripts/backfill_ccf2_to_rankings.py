#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

CCF_FIELDS = ["ccf_estimate", "clonality_status", "ccf_confidence", "ccf_warning"]
TARGETS = [
    "scoring/ranked_peptides.tsv",
    "scoring/ranked_peptides.recommendation.tsv",
    "scoring/ranked_peptides.netmhcpan42.tsv",
]


def read_tsv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    if not path.is_file():
        return [], []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_ccf_index(ccf_path: Path) -> dict[str, dict[str, str]]:
    rows, _ = read_tsv(ccf_path)
    index: dict[str, dict[str, str]] = {}
    for row in rows:
        event_id = row.get("event_id", "")
        if not event_id:
            continue
        index[event_id] = {field: row.get(field, "") for field in CCF_FIELDS}
    return index


def backfill(path: Path, ccf_by_event: dict[str, dict[str, str]]) -> tuple[int, int]:
    rows, fields = read_tsv(path)
    if not rows:
        return 0, 0
    matched = 0
    for row in rows:
        ccf = ccf_by_event.get(row.get("event_id", ""))
        if not ccf:
            continue
        matched += 1
        row.update(ccf)
    out_fields = list(fields)
    for field in CCF_FIELDS:
        if field not in out_fields:
            if field == "ccf_estimate" and "ccf_multiplier" in out_fields:
                out_fields.insert(out_fields.index("ccf_multiplier"), field)
            elif field == "clonality_status" and "ccf_multiplier" in out_fields:
                out_fields.insert(out_fields.index("ccf_multiplier"), field)
            elif field == "ccf_confidence" and "ccf_multiplier" in out_fields:
                out_fields.insert(out_fields.index("ccf_multiplier") + 1, field)
            elif field == "ccf_warning" and "ccf_multiplier" in out_fields:
                out_fields.insert(out_fields.index("ccf_multiplier") + 1, field)
            else:
                out_fields.append(field)
    write_tsv(path, rows, out_fields)
    return len(rows), matched


def main() -> None:
    parser = argparse.ArgumentParser(description="Backfill CCF 2.1 explanatory fields into peptide ranking tables.")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--ccf", type=Path)
    parser.add_argument("--targets", nargs="*", default=TARGETS)
    args = parser.parse_args()

    outdir = args.outdir.resolve()
    ccf_path = args.ccf or outdir / "clonality" / "ccf_2.tsv"
    ccf_by_event = build_ccf_index(ccf_path)
    if not ccf_by_event:
        raise SystemExit(f"No CCF rows found: {ccf_path}")

    for target in args.targets:
        path = Path(target)
        if not path.is_absolute():
            path = outdir / path
        total, matched = backfill(path, ccf_by_event)
        print(f"ccf_backfill target={path} rows={total} matched={matched}")


if __name__ == "__main__":
    main()
