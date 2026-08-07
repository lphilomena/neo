#!/usr/bin/env python3
"""Annotate splice neoantigen candidates with SQANTI3 long-read support."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def truth(value: bool) -> str:
    return "true" if value else "false"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--sqanti-junctions", required=True, type=Path)
    parser.add_argument("--sqanti-classification", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    classes = {row["isoform"]: row for row in read_tsv(args.sqanti_classification)}
    support = defaultdict(list)
    for row in read_tsv(args.sqanti_junctions):
        chrom = row["chrom"]
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        # SQANTI3 stores intron bases; SNAF/SpliceMutr use flanking exon boundaries.
        key = (
            chrom,
            str(int(row["genomic_start_coord"]) - 1),
            str(int(row["genomic_end_coord"]) + 1),
            row["strand"],
        )
        support[key].append(row)

    with args.candidates.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        input_fields = list(reader.fieldnames or [])
        candidates = list(reader)

    extra_fields = [
        "longrna_support",
        "longrna_isoform_count",
        "longrna_isoforms",
        "longrna_structural_categories",
        "longrna_junction_categories",
        "longrna_max_unique_coverage",
        "longrna_max_multi_coverage",
        "longrna_all_canonical",
    ]
    matched = 0
    supported_pairs = set()
    for row in candidates:
        parts = row.get("junction_key", "").split(":")
        hits = support.get(tuple(parts), []) if len(parts) == 4 else []
        isoforms = sorted({hit["isoform"] for hit in hits})
        categories = sorted(
            {
                classes.get(isoform, {}).get("structural_category", "unknown")
                for isoform in isoforms
            }
        )
        junction_categories = sorted({hit.get("junction_category", "") for hit in hits})
        unique_cov = [int(float(hit.get("total_coverage_unique") or 0)) for hit in hits]
        multi_cov = [int(float(hit.get("total_coverage_multi") or 0)) for hit in hits]
        canonical = [hit.get("canonical", "").upper() == "CANONICAL" for hit in hits]
        row.update(
            {
                "longrna_support": truth(bool(hits)),
                "longrna_isoform_count": str(len(isoforms)),
                "longrna_isoforms": ";".join(isoforms),
                "longrna_structural_categories": ";".join(categories),
                "longrna_junction_categories": ";".join(junction_categories),
                "longrna_max_unique_coverage": str(max(unique_cov, default=0)),
                "longrna_max_multi_coverage": str(max(multi_cov, default=0)),
                "longrna_all_canonical": truth(bool(canonical) and all(canonical)),
            }
        )
        if hits:
            matched += 1
            supported_pairs.add((row.get("peptide", ""), row.get("hla", "")))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=input_fields + extra_fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(candidates)

    summary = {
        "status": "PASS",
        "coordinate_policy": "SQANTI3 intron start-1/end+1 to SNAF exon-boundary key",
        "candidate_rows": len(candidates),
        "longrna_supported_rows": matched,
        "longrna_supported_unique_peptide_hla_pairs": len(supported_pairs),
        "longrna_support_fraction": matched / len(candidates) if candidates else 0.0,
        "output": str(args.output),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
