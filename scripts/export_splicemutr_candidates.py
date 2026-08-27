#!/usr/bin/env python3
"""Export official SpliceMutr transcript output as NeoAg junction peptides."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


FIELDS = ["sample_id", "source_tool", "source_record_id", "source_junction_id", "chrom",
          "start", "end", "strand", "gene", "transcript_id", "event_type", "peptide",
          "hla_allele", "provided_rna_junction_reads", "rna_junction_reads", "frame_status",
          "normal_junction_status", "normal_background_reason", "source_file"]


def integer(value: str) -> int | None:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


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
    normal_ids: set[str] = set()
    if args.normal_junctions and args.normal_junctions.is_file():
        with args.normal_junctions.open(encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                normal_ids.add(row.get("canonical_junction_id") or row.get("junction_id") or "")

    output: list[dict[str, str]] = []
    with args.metadata.open(newline="", encoding="utf-8") as handle:
        for row_number, row in enumerate(csv.DictReader(handle, delimiter="\t"), 1):
            chrom = row.get("chr", ""); start = row.get("start", ""); end = row.get("end", "")
            strand = row.get("strand", ""); junction = row.get("juncs") or f"{chrom}:{start}:{end}:{strand}"
            protein = row.get("peptide", "")
            positions = [integer(value.lstrip("+")) for value in (row.get("pep_junc_loc") or "").split(",")]
            positions = [value for value in positions if value is not None and 0 < value < len(protein)]
            if not positions:
                continue
            normal_status = "DETECTED_BROAD_NORMAL" if junction in normal_ids else "NOT_DETECTED_COVERAGE_UNASSESSED"
            reason = "exact junction present in configured normal reference" if junction in normal_ids else "no exact normal match; coverage was not supplied"
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
                        "normal_junction_status": normal_status, "normal_background_reason": reason,
                        "source_file": str(args.metadata.resolve())})
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
