#!/usr/bin/env python3
"""Migrate historical SNAF/SpliceMutr TSVs to normal-background P0 semantics.

The migration is intentionally conservative.  It never turns catalog
non-membership into a negative normal-tissue result and never upgrades evidence.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


NO_NORMAL_COHORT = "UNASSESSED_NO_COMPATIBLE_NORMAL_RNA_COHORT"
STATUS_MAP = {
    "ABSENT_GTEX_V11": "NOT_LISTED_IN_NORMAL_CATALOG",
    "SEEN_GTEX_V11": "SEEN_IN_NORMAL_CATALOG",
}


def is_splice_row(row: dict[str, str]) -> bool:
    text = " ".join(
        row.get(key, "")
        for key in ("event_type", "mutation_source", "source_tool", "source_generator", "peptide_consequence")
    ).casefold()
    return any(token in text for token in ("splice", "snaf", "splicemutr", "junction"))


def migrate(row: dict[str, str], kind: str, normal_catalog_label: str) -> dict[str, str]:
    out = dict(row)
    for field in ("normal_junction_status", "normal_junction_assessment_status"):
        old = out.get(field, "")
        if old in STATUS_MAP:
            out[field] = STATUS_MAP[old]
    if out.get("normal_junction_status") == "NOT_LISTED_IN_NORMAL_CATALOG":
        out["normal_junction_assessment_status"] = "NOT_LISTED_CATALOG_COVERAGE_UNASSESSED"
    out.setdefault("normal_catalog_label", normal_catalog_label)

    if kind in {"peptides", "events"} and is_splice_row(out):
        out["cohort_analysis_status"] = NO_NORMAL_COHORT
        out["tumor_specificity_status"] = NO_NORMAL_COHORT
        out["normal_safety_grade"] = "N1"
        out["splice_consensus_tier"] = "R3"
        out["safety_priority_cap"] = "R3"
        if out.get("priority_cap") in {"", "C_CAUTION"}:
            out["priority_cap"] = "R3"
        if out.get("rna_evidence_completeness") == "COMPLETE":
            out["rna_evidence_completeness"] = "PARTIAL_NO_COMPATIBLE_NORMAL_RNA_COHORT"

    if kind == "peptides" and is_splice_row(out):
        if out.get("contains_novel_aa", "").casefold() in {"yes", "true", "1"}:
            out["structural_novelty_status"] = "ALTERED_JUNCTION_SPANNING_SEQUENCE"
        if out.get("mutant_specificity_status") == "MT_SPECIFIC":
            out["mutant_specificity_status"] = "UNASSESSED"
            out["mutant_specificity_gate_status"] = "REVIEW_REQUIRED"
            out["mutant_specificity_priority_cap"] = "R3"
            out["mutant_specificity_reason"] = (
                "Historical structural splice novelty retained; tumor specificity is unassessed "
                "without a compatible normal RNA cohort."
            )
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--kind", choices=("candidates", "peptides", "events"), required=True)
    parser.add_argument("--normal-catalog-label", default="recount3_GTEx_v8_GRCh38")
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        source_fields = list(reader.fieldnames or [])
        migrated = [migrate(row, args.kind, args.normal_catalog_label) for row in reader]
    added_fields = sorted({key for row in migrated for key in row} - set(source_fields))
    fields = source_fields + added_fields
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(migrated)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
