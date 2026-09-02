#!/usr/bin/env python3
"""Export official SpliceMutr transcript output as NeoAg junction peptides."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from typing import TextIO


FIELDS = ["sample_id", "source_tool", "source_record_id", "source_junction_id", "chrom",
          "start", "end", "strand", "gene", "transcript_id", "event_type", "peptide",
          "hla_allele", "provided_rna_junction_reads", "rna_junction_reads", "frame_status",
          "normal_junction_status", "normal_background_reason", "source_file"]


def integer(value: str) -> int | None:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def matched_normal_junctions(path: Path | None, candidates: set[str]) -> set[str]:
    """Stream the normal catalog and retain only candidate exact matches."""
    if not path or not path.is_file() or not candidates:
        return set()
    matched: set[str] = set()
    with open_text(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            junction = row.get("canonical_junction_id") or row.get("junction_id") or ""
            if junction in candidates:
                matched.add(junction)
                if matched == candidates:
                    break
    return matched


def windows(protein: str, boundary: int, lengths: list[int]):
    protein = protein.rstrip("*").upper()
    for length in lengths:
        for start in range(max(0, boundary - length + 1), min(boundary, len(protein) - length) + 1):
            end = start + length
            if start < boundary < end:
                peptide = protein[start:end]
                if peptide and set(peptide) <= set("ACDEFGHIKLMNPQRSTVWY"):
                    yield peptide


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--peptide-lengths", default="8,9,10,11")
    parser.add_argument("--normal-junctions", type=Path)
    args = parser.parse_args()
    lengths = sorted({int(value) for value in args.peptide_lengths.split(",")})
    output: list[dict[str, str]] = []
    candidate_junctions: set[str] = set()
    with args.metadata.open(newline="", encoding="utf-8") as handle:
        for row_number, row in enumerate(csv.DictReader(handle, delimiter="\t"), 1):
            chrom = row.get("chr", ""); start = row.get("start", ""); end = row.get("end", "")
            strand = row.get("strand", ""); junction = row.get("juncs") or f"{chrom}:{start}:{end}:{strand}"
            protein = row.get("peptide", "")
            positions = [integer(value.lstrip("+")) for value in (row.get("pep_junc_loc") or "").split(",")]
            positions = [value for value in positions if value is not None and 0 < value < len(protein)]
            if not positions:
                continue
            candidate_junctions.add(junction)
            for boundary in positions:
                for peptide in windows(protein, boundary, lengths):
                    digest = hashlib.sha1(f"{junction}|{row.get('tx_id','')}|{peptide}".encode()).hexdigest()[:16]
                    output.append({"sample_id": args.sample_id, "source_tool": "SpliceMutr",
                        "source_record_id": f"SMR|{digest}", "source_junction_id": junction,
                        "chrom": chrom, "start": start, "end": end, "strand": strand,
                        "gene": row.get("gene", ""), "transcript_id": row.get("tx_id", ""),
                        "event_type": "Splice", "peptide": peptide, "hla_allele": "",
                        "provided_rna_junction_reads": "", "rna_junction_reads": "",
                        "frame_status": "IN_FRAME" if row.get("error") == "tx" else "ORF_REVIEW",
                        "normal_junction_status": "NOT_DETECTED_COVERAGE_UNASSESSED",
                        "normal_background_reason": "no exact normal match; coverage was not supplied",
                        "source_file": str(args.metadata.resolve())})
    normal_ids = matched_normal_junctions(args.normal_junctions, candidate_junctions)
    for row in output:
        if row["source_junction_id"] in normal_ids:
            row["normal_junction_status"] = "DETECTED_BROAD_NORMAL"
            row["normal_background_reason"] = "exact junction present in configured normal reference"
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader(); writer.writerows(output)
    manifest = {"status": "PASS", "source_rows": row_number if 'row_number' in locals() else 0,
                "candidate_rows": len(output), "peptide_lengths": lengths,
                "normal_background_policy": "absence without coverage remains UNASSESSED"}
    args.out.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
