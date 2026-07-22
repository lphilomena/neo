#!/usr/bin/env python3
"""Create consistency checks and a concise summary for a SNAF splice branch."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path: Path, source_path: Path, data: list[dict[str, str]]) -> None:
    with source_path.open(newline="", encoding="utf-8") as handle:
        fields = list(csv.DictReader(handle, delimiter="\t").fieldnames or [])
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in data)


def counts(items: list[dict[str, str]], field: str) -> dict[str, int]:
    return dict(sorted(Counter(row.get(field, "") or "BLANK" for row in items).items()))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--branch-run", required=True, type=Path)
    args = parser.parse_args()
    branch = args.branch_run.resolve()

    ranked_event_path = branch / "scoring" / "ranked_events.v03.tsv"
    ranked_peptide_path = branch / "scoring" / "ranked_peptides.v03.tsv"
    ranked_events = rows(ranked_event_path)
    ranked_peptides = rows(ranked_peptide_path)
    selected = rows(branch / "splice_snaf" / "snaf_candidates.selected.tsv")
    snaf_events = [row for row in ranked_events if row.get("mutation_source") == "SNAF"]
    snaf_peptides = [row for row in ranked_peptides if row.get("mutation_source") == "SNAF"]
    write_rows(
        branch / "splice_snaf" / "snaf_ranked_events.integrated.tsv",
        ranked_event_path,
        snaf_events,
    )
    write_rows(
        branch / "splice_snaf" / "snaf_ranked_peptides.integrated.tsv",
        ranked_peptide_path,
        snaf_peptides,
    )

    violations = {
        "snaf_priority_ab": sum(
            row.get("final_priority", "").startswith(("A", "B")) for row in snaf_peptides
        ),
        "snaf_ccf_numeric": sum(
            (row.get("ccf_estimate", "") or "").strip() not in {"", "NA", "N/A"}
            for row in snaf_peptides
        ),
        "snaf_ccf_wrong_status": sum(
            row.get("ccf_status") != "RNA_ONLY_UNRESOLVED" for row in snaf_peptides
        ),
        "selected_gtex_seen": sum(
            row.get("normal_junction_status") != "ABSENT_GTEX_V11" for row in selected
        ),
        "snaf_missing_junction_reads": sum(
            not (row.get("rna_junction_reads", "") or "").strip() for row in snaf_peptides
        ),
    }

    top_snaf = []
    for rank, row in enumerate(ranked_peptides, start=1):
        if row.get("mutation_source") != "SNAF":
            continue
        top_snaf.append(
            {
                "global_rank": rank,
                "gene": row.get("gene", ""),
                "event_id": row.get("event_id", ""),
                "peptide": row.get("peptide", ""),
                "wildtype_peptide": row.get("wildtype_peptide", ""),
                "hla_allele": row.get("hla_allele", ""),
                "junction_reads": row.get("rna_junction_reads", ""),
                "gene_tpm": row.get("gene_expression_tpm", ""),
                "transcript_tpm": row.get("transcript_expression_tpm", ""),
                "presentation_grade": row.get("presentation_evidence_grade", ""),
                "safety_status": row.get("safety_status", ""),
                "mutant_specificity_status": row.get("mutant_specificity_status", ""),
                "mutant_specificity_gate_status": row.get("mutant_specificity_gate_status", ""),
                "final_priority": row.get("final_priority", ""),
                "efficacy_score": row.get("efficacy_score", ""),
            }
        )
        if len(top_snaf) == 20:
            break

    summary = {
        "status": "PASS" if not any(violations.values()) else "REVIEW_REQUIRED",
        "branch_run": str(branch),
        "event_count": len(ranked_events),
        "peptide_hla_count": len(ranked_peptides),
        "event_type_counts": counts(ranked_events, "event_type"),
        "mutation_source_event_counts": counts(ranked_events, "mutation_source"),
        "priority_counts_all": counts(ranked_peptides, "final_priority"),
        "snaf": {
            "event_count": len(snaf_events),
            "peptide_hla_count": len(snaf_peptides),
            "with_gene_tpm": sum(bool(row.get("gene_expression_tpm", "").strip()) for row in snaf_peptides),
            "with_transcript_tpm": sum(
                bool(row.get("transcript_expression_tpm", "").strip()) for row in snaf_peptides
            ),
            "with_positive_expression_score": sum(
                float(row.get("l3_expression_score", "0") or 0) > 0 for row in snaf_peptides
            ),
            "priority_counts": counts(snaf_peptides, "final_priority"),
            "presentation_grade_counts": counts(snaf_peptides, "presentation_evidence_grade"),
            "safety_status_counts": counts(snaf_peptides, "safety_status"),
            "mutant_specificity_status_counts": counts(
                snaf_peptides, "mutant_specificity_status"
            ),
            "mutant_specificity_gate_counts": counts(
                snaf_peptides, "mutant_specificity_gate_status"
            ),
            "ccf_status_counts": counts(snaf_peptides, "ccf_status"),
            "normal_junction_assessment_counts": counts(
                snaf_peptides, "normal_junction_assessment_status"
            ),
        },
        "violations": violations,
        "top_snaf_peptide_hla": top_snaf,
        "interpretation_boundary": (
            "SNAF splice candidates are RNA-derived computational candidates. Gene and transcript "
            "TPM come from the matched WTS RSEM/Gencode expression tables; junction reads remain a "
            "separate event-level signal. RNA-only events do not have resolved DNA CCF. Candidates "
            "remain capped at C_CAUTION until orthogonal junction, frame, and peptide validation."
        ),
    }

    qc_dir = branch / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    json_path = qc_dir / "snaf_splice_branch_qc.json"
    json_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# SNAF splice branch summary",
        "",
        f"- Status: **{summary['status']}**",
        f"- Combined events: {len(ranked_events)}",
        f"- Combined peptide-HLA rows: {len(ranked_peptides)}",
        f"- SNAF events: {len(snaf_events)}",
        f"- SNAF peptide-HLA rows: {len(snaf_peptides)}",
        f"- SNAF rows with gene TPM: {summary['snaf']['with_gene_tpm']}",
        f"- SNAF rows with transcript TPM: {summary['snaf']['with_transcript_tpm']}",
        f"- SNAF rows with positive expression score: {summary['snaf']['with_positive_expression_score']}",
        f"- SNAF priorities: {summary['snaf']['priority_counts']}",
        f"- SNAF presentation grades: {summary['snaf']['presentation_grade_counts']}",
        f"- SNAF safety states: {summary['snaf']['safety_status_counts']}",
        f"- SNAF mutant-specificity states: {summary['snaf']['mutant_specificity_status_counts']}",
        f"- SNAF mutant-specificity gates: {summary['snaf']['mutant_specificity_gate_counts']}",
        f"- Consistency violations: {violations}",
        "",
        "## Top SNAF candidates",
        "",
        "| Global rank | Gene | MT peptide | WT peptide | HLA | Junction reads | Gene TPM | Transcript TPM | Presentation | MT specificity | Safety | Priority |",
        "|---:|---|---|---|---|---:|---:|---:|---|---|---|---|",
    ]
    for row in top_snaf:
        lines.append(
            "| {global_rank} | {gene} | {peptide} | {wildtype_peptide} | {hla_allele} | {junction_reads} | "
            "{gene_tpm} | {transcript_tpm} | "
            "{presentation_grade} | {mutant_specificity_status} | {safety_status} | {final_priority} |".format(**row)
        )
    lines.extend(["", "## Interpretation boundary", "", summary["interpretation_boundary"], ""])
    md_path = qc_dir / "snaf_splice_branch_summary.md"
    md_path.write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
