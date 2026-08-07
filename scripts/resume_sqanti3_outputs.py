#!/usr/bin/env python3
"""Finalize SQANTI3 outputs after the historical isoform-hits API failure."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path


HIT_FIELDS = [
    "Isoform", "Isoform_length", "Isoform_exon_number", "Hit", "Hit_length",
    "Hit_exon_number", "Match", "Diff_to_TSS", "Diff_to_TTS", "Matching_type",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", required=True)
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()

    outdir = Path(args.work) / "sqanti3"
    classification = outdir / f"{args.sample}_classification.txt"
    junctions = outdir / f"{args.sample}_junctions.txt"
    hits_tmp = outdir / f"{args.sample}_isoform_hits.txt_tmp"
    hits_final = outdir / f"{args.sample}_isoform_hits.txt"
    corrected_gtf = outdir / f"{args.sample}_corrected.gtf"

    for path in (classification, junctions, hits_tmp, corrected_gtf):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Required SQANTI3 output is missing or empty: {path}")

    primary: dict[str, set[str]] = {}
    category_counts: dict[str, int] = {}
    with classification.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            isoform = row["isoform"]
            transcripts = row.get("associated_transcript", "")
            primary[isoform] = {x for x in transcripts.replace(";", ",").split(",") if x and x != "novel"}
            category = row.get("structural_category", "NA")
            category_counts[category] = category_counts.get(category, 0) + 1

    rows = []
    with hits_tmp.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row["Matching_type"] = "primary" if row["Hit"] in primary.get(row["Isoform"], set()) else "secondary"
            rows.append(row)
    rows.sort(key=lambda row: row["Isoform"])

    staged = hits_final.with_suffix(hits_final.suffix + ".neoag.tmp")
    with staged.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HIT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(staged, hits_final)

    summary = outdir / f"{args.sample}.classification_summary.tsv"
    with summary.open("w") as handle:
        handle.write("metric\tvalue\n")
        handle.write(f"classified_isoforms\t{sum(category_counts.values())}\n")
        for category in sorted(category_counts):
            handle.write(f"structural_category:{category}\t{category_counts[category]}\n")

    print(f"classified_isoforms={sum(category_counts.values())}")
    print(f"isoform_hits={len(rows)}")
    print(f"hits_output={hits_final}")
    print(f"summary_output={summary}")


if __name__ == "__main__":
    main()
