#!/usr/bin/env python3
"""Convert recount3 exon-exon junction matrices to neoag normal-junction TSV."""

from __future__ import annotations

import argparse
from array import array
import csv
import gzip
import json
from pathlib import Path


FIELDS = [
    "junction_id", "chromosome", "start", "end", "strand", "genome_build",
    "annotated", "normal_samples", "normal_reads", "normal_total_reads",
    "tissue", "source", "dataset",
]


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("r")


def matrix_summary(path: Path) -> tuple[array, array, array, int]:
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            nrows, ncols, _ = map(int, line.split())
            break
        else:
            raise ValueError(f"Matrix Market dimensions missing: {path}")
        max_reads = array("Q", [0]) * nrows
        total_reads = array("Q", [0]) * nrows
        sample_count = array("I", [0]) * nrows
        for line in handle:
            if not line.strip() or line.startswith("%"):
                continue
            row, _column, raw_value = line.split()[:3]
            value = int(float(raw_value))
            if value <= 0:
                continue
            idx = int(row) - 1
            total_reads[idx] += value
            sample_count[idx] += 1
            if value > max_reads[idx]:
                max_reads[idx] = value
    return max_reads, total_reads, sample_count, ncols


def normalize_chromosome(value: str) -> str:
    chrom = value.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    if chrom in {"MT", "M"}:
        chrom = "M"
    return f"chr{chrom}"


def convert(rr: Path, mm: Path, output: Path, tissue: str, min_reads: int, min_samples: int) -> dict:
    max_reads, total_reads, sample_count, matrix_samples = matrix_summary(mm)
    output.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    seen = 0
    with open_text(rr) as source, gzip.open(output, "wt", newline="") as target:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(target, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for idx, row in enumerate(reader):
            seen += 1
            if idx >= len(max_reads):
                raise ValueError("RR has more rows than the Matrix Market matrix")
            if max_reads[idx] < min_reads or sample_count[idx] < min_samples:
                continue
            raw_chrom = row.get("chromosome", "")
            if raw_chrom.startswith("ERCC-"):
                continue
            chrom = normalize_chromosome(raw_chrom)
            start, end = int(row["start"]), int(row["end"])
            strand = row.get("strand", "?")
            writer.writerow({
                "junction_id": f"{chrom}:{start}-{end}:{strand}",
                "chromosome": chrom, "start": start, "end": end, "strand": strand,
                "genome_build": "GRCh38", "annotated": row.get("annotated", ""),
                "normal_samples": sample_count[idx], "normal_reads": max_reads[idx],
                "normal_total_reads": total_reads[idx], "tissue": tissue,
                "source": "recount3_GTEx_v8", "dataset": "GTEx_v8",
            })
            kept += 1
    if seen != len(max_reads):
        raise ValueError(f"RR/MM row mismatch: RR={seen}, MM={len(max_reads)}")
    result = {
        "tissue": tissue, "genome_build": "GRCh38", "matrix_samples": matrix_samples,
        "junction_rows": seen, "junctions_kept": kept, "min_reads": min_reads,
        "min_samples": min_samples, "rr": str(rr), "mm": str(mm), "output": str(output),
    }
    output.with_suffix(output.suffix + ".meta.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rr", type=Path, required=True)
    parser.add_argument("--mm", type=Path, required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-reads", type=int, default=1)
    parser.add_argument("--min-samples", type=int, default=1)
    args = parser.parse_args()
    print(json.dumps(convert(args.rr, args.mm, args.output, args.tissue, args.min_reads, args.min_samples), indent=2))


if __name__ == "__main__":
    main()
