#!/usr/bin/env python3
"""Merge per-tissue normal-junction tables into a pan-tissue catalog."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path


def open_text(path: Path, mode: str):
    return gzip.open(path, mode, newline="") if path.suffix == ".gz" else path.open(mode, newline="")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    merged: dict[str, dict[str, object]] = {}
    for path in args.inputs:
        with open_text(path, "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = row["junction_id"]
                item = merged.setdefault(key, {
                    **row, "normal_samples": 0, "normal_reads": 0, "normal_total_reads": 0,
                    "normal_tissues": 0, "tissues": set(),
                })
                item["normal_samples"] = int(item["normal_samples"]) + int(row["normal_samples"])
                item["normal_reads"] = max(int(item["normal_reads"]), int(row["normal_reads"]))
                item["normal_total_reads"] = int(item["normal_total_reads"]) + int(row["normal_total_reads"])
                tissues = item["tissues"]
                assert isinstance(tissues, set)
                tissues.add(row["tissue"])
                item["normal_tissues"] = len(tissues)
                if row.get("annotated") == "1":
                    item["annotated"] = "1"
    fields = [
        "junction_id", "chromosome", "start", "end", "strand", "genome_build", "annotated",
        "normal_samples", "normal_reads", "normal_total_reads", "normal_tissues", "tissue",
        "source", "dataset",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.output, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for key in sorted(merged):
            row = merged[key]
            tissues = row.pop("tissues")
            assert isinstance(tissues, set)
            row["tissue"] = ",".join(sorted(tissues))
            writer.writerow({field: row.get(field, "") for field in fields})
    metadata = {
        "genome_build": "GRCh38", "source": "recount3_GTEx_v8",
        "input_files": [str(path) for path in args.inputs], "junctions_kept": len(merged),
        "output": str(args.output),
    }
    args.output.with_suffix(args.output.suffix + ".meta.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
