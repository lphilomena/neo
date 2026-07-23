from __future__ import annotations

from collections.abc import Iterable
import json
import math
from pathlib import Path
from typing import Any

from .utils import read_tsv, write_tsv


EVIDENCE_SOURCE_PRECEDENCE_VERSION = "1.0"

IDENTITY_FIELDS = {
    "peptide_id", "event_id", "sample_id", "peptide", "mutant_peptide",
    "wildtype_peptide", "hla_allele", "mhc_class", "crosses_junction",
    "contains_novel_aa", "peptide_consequence", "mutation_positions_in_peptide",
}
EVENT_FIELDS = {
    "event_type", "mutation_source", "gene", "gene_id", "transcript_id",
    "transcript", "chrom", "pos", "ref", "alt", "protein_change", "hgvsp",
    "consequence", "tumor_vaf", "tumor_depth", "tumor_alt_count", "phase_group_id",
    "haplotype_status", "phase_support_reads", "phase_total_informative_reads",
    "phase_confidence", "component_event_ids", "combined_protein_change",
}
PRESENTATION_FIELDS = {
    "source_tool", "binding_rank", "el_rank", "presentation_score",
    "immunogenicity_score", "wildtype_binding_rank", "netmhcpan_mt_ic50",
    "netmhcpan_mt_rank_ba", "netmhcpan_mt_rank_el", "netmhcpan_wt_ic50",
    "netmhcpan_wt_rank_ba", "netmhcpan_wt_rank_el", "netmhcpan_ba_rank",
    "netmhcpan_el_rank", "netmhcstabpan_score", "netmhcstabpan_rank",
    "mhcflurry_affinity_percentile", "mhcflurry_processing_score",
    "mhcflurry_presentation_score", "mhcflurry_wt_affinity_percentile",
    "mhcflurry_wt_processing_score", "mhcflurry_wt_presentation_score",
    "binding_evidence_score", "presentation_evidence_score",
    "presentation_evidence_grade", "iedb_immunogenicity_score",
    "immunogenicity_resolved", "prime_score", "prime_rank", "prime_wt_score",
    "prime_wt_rank", "bigmhc_im_score", "bigmhc_im_wt_score", "deepimmuno_score",
    "immunogenicity_composite_score", "immunogenicity_source",
    "presentation_gate_status", "presentation_gate_reason", "presentation_gate_multiplier",
}
EXPRESSION_FIELDS = {
    "gene_expression_tpm", "transcript_expression_tpm", "expression_tpm",
    "expression_source", "expression_evidence_status", "rna_evidence_completeness",
    "rna_evidence_score",
}
RNA_FIELDS = {
    "rna_support_status", "rna_vaf", "rna_alt_reads", "rna_depth", "rna_vaf_source",
    "rna_junction_reads", "rna_junction_source", "rna_frame_status", "junction_reads",
    "junction_source", "junction_status",
}
CCF_FIELDS = {
    "raw_ccf", "ccf_estimate", "ccf_best", "ccf_status", "clonality_status",
    "ccf_confidence", "ccf_warning", "ccf_method", "ccf_resolution",
    "ccf_resolution_reason", "ccf_multiplier", "purity", "total_cn", "minor_cn",
    "multiplicity",
}
APPM_FIELDS = {
    "appm_multiplier", "appm_multiplier_reason", "appm_integrity_status",
    "appm_evidence_completeness", "appm_review_required", "appm_action",
    "appm_call_confidence", "restricting_locus_expression_status",
}
ESCAPE_FIELDS = {
    "escape_status", "escape_flag", "escape_reason", "resistance_risk",
    "escape_action", "escape_multiplier", "restricting_hla_lost", "hla_loh_status",
    "hla_loh_alleles",
}
SAFETY_FIELDS = {
    "safety_tier", "safety_status", "safety_reason", "safety_multiplier",
    "review_required", "reference_proteome_exact_match", "normal_ligand_tissue",
    "mutation_anchor_only", "normal_tissue_max_tpm", "normal_tissue_max_tissue",
    "critical_tissue_max_tpm", "critical_tissue_name", "normal_hspc_tpm",
    "normal_hspc_unit", "normal_expression_status", "normal_hspc_status",
    "reference_proteome_status", "normal_ligandome_status", "anchor_assessment_status",
    "normal_junction_assessment_status", "safety_evidence_completeness",
    "safety_missing_layers", "safety_priority_cap",
}
EVENT_SAFETY_FIELDS = {f"event_{field}" for field in SAFETY_FIELDS}
LEGACY_RANKING_FIELDS = {
    "priority_cap", "wes_confidence_tier", "l3_event_confidence_score",
    "l3_expression_score", "l3_clonality_score", "l3_tumor_specificity_score",
    "l3_hla_binding_score", "l3_hla_presentation_score", "l3_rna_support_score",
    "l3_rna_junction_support_score", "l3_normal_tissue_safety_score",
    "l3_apm_integrity_score", "l3_immunogenicity_score", "immunology_composite_score",
    "efficacy_score", "final_priority", "recommended_use",
}
VALIDATION_FIELDS = {
    "validation_mode", "recommended_assay", "validation_priority",
    "validation_rationale", "wt_control_required", "minigene_design",
}

# Public registry used by audits and downstream code.
AUTHORITATIVE_FIELDS = {
    "annotated_peptides": IDENTITY_FIELDS | PRESENTATION_FIELDS,
    "raw_peptides": IDENTITY_FIELDS,
    "raw_events": EVENT_FIELDS,
    "presentation_evidence": PRESENTATION_FIELDS,
    "expression_evidence": EXPRESSION_FIELDS,
    "rna_junction_evidence": RNA_FIELDS,
    "ccf_2": CCF_FIELDS,
    "appm_peptide_modifiers": APPM_FIELDS,
    "peptide_escape_flags": ESCAPE_FIELDS,
    "peptide_safety": SAFETY_FIELDS,
    "event_safety": EVENT_SAFETY_FIELDS,
    "ranked_peptides": LEGACY_RANKING_FIELDS,
    "validation_plan": VALIDATION_FIELDS,
}

FIELD_SOURCE_PRECEDENCE: dict[str, tuple[str, ...]] = {}
for fields, sources in (
    (IDENTITY_FIELDS, ("annotated_peptides", "raw_peptides", "ranked_peptides")),
    (EVENT_FIELDS, ("raw_events", "annotated_peptides", "ranked_peptides")),
    (PRESENTATION_FIELDS, ("presentation_evidence", "annotated_peptides", "raw_peptides", "ranked_peptides")),
    (EXPRESSION_FIELDS, ("expression_evidence", "ranked_peptides", "raw_events")),
    (RNA_FIELDS, ("rna_junction_evidence", "raw_events", "raw_peptides", "ranked_peptides")),
    (CCF_FIELDS, ("ccf_2", "ranked_peptides", "raw_events")),
    (APPM_FIELDS, ("appm_peptide_modifiers", "ranked_peptides")),
    (ESCAPE_FIELDS, ("peptide_escape_flags", "ranked_peptides")),
    (SAFETY_FIELDS, ("peptide_safety", "event_safety", "ranked_peptides")),
    (EVENT_SAFETY_FIELDS, ("event_safety",)),
    (LEGACY_RANKING_FIELDS, ("ranked_peptides",)),
    (VALIDATION_FIELDS, ("validation_plan", "ranked_peptides")),
):
    for field in fields:
        FIELD_SOURCE_PRECEDENCE[field] = sources

PEPTIDE_SOURCES = (
    "ranked_peptides", "raw_peptides", "presentation_evidence",
    "appm_peptide_modifiers", "peptide_safety", "peptide_escape_flags",
    "validation_plan",
)
EVENT_SOURCES = (
    "raw_events", "ccf_2", "expression_evidence", "rna_junction_evidence",
    "event_safety",
)
SOURCE_ORDER = (
    "annotated_peptides", "raw_peptides", "raw_events", "presentation_evidence",
    "expression_evidence", "rna_junction_evidence", "ccf_2",
    "appm_peptide_modifiers", "peptide_escape_flags", "peptide_safety",
    "event_safety", "ranked_peptides", "validation_plan",
)
CONFLICT_FIELDS = (
    "peptide_id", "event_id", "gene", "field", "selected_source", "selected_value",
    "other_source", "other_value", "precedence_version", "conflict_type",
)


def _read_optional(path: str | Path | None) -> list[dict[str, str]]:
    if not path:
        return []
    source = Path(path)
    return read_tsv(source) if source.is_file() else []


def _field_order(rows: Iterable[dict[str, str]]) -> list[str]:
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for field in row:
            if field not in seen:
                seen.add(field)
                fields.append(field)
    return fields


def _index(rows: list[dict[str, str]], key: str) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        value = str(row.get(key) or "")
        if not value:
            continue
        current = result.setdefault(value, {})
        for field, field_value in row.items():
            if field_value not in (None, "") and not current.get(field):
                current[field] = str(field_value)
    return result


def _event_index(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    result = _index(rows, "event_id")
    for row in rows:
        event_id = str(row.get("event_id") or "")
        if not event_id:
            continue
        current = result[event_id]
        for field in ("junction_reads", "rna_junction_reads", "rna_alt_reads", "rna_depth"):
            value = str(row.get(field) or "")
            if not value:
                continue
            try:
                if float(value) > float(current.get(field) or "-inf"):
                    current[field] = value
            except ValueError:
                if not current.get(field):
                    current[field] = value
    return result


def _base_rows(
    annotated_peptides: str | Path | None,
    ranked_peptides: str | Path,
) -> tuple[list[dict[str, str]], str]:
    annotated = _read_optional(annotated_peptides)
    if annotated:
        return annotated, "annotated_peptides"
    ranked = _read_optional(ranked_peptides)
    if not ranked:
        raise ValueError("Comprehensive evidence requires ranked_peptides or annotated_peptides")
    return ranked, "ranked_peptides"


def _equivalent(left: str, right: str) -> bool:
    a = str(left).strip()
    b = str(right).strip()
    if a == b or a.casefold() == b.casefold():
        return True
    boolean = {"true": "1", "yes": "1", "false": "0", "no": "0"}
    if boolean.get(a.casefold(), a.casefold()) == boolean.get(b.casefold(), b.casefold()):
        return True
    try:
        x = float(a)
        y = float(b)
        return math.isclose(x, y, rel_tol=2e-3, abs_tol=1e-10)
    except ValueError:
        return False


def _select_value(
    field: str,
    candidates: list[tuple[str, str]],
    base_source: str,
) -> tuple[str, str]:
    by_source = {source: value for source, value in candidates}
    for source in FIELD_SOURCE_PRECEDENCE.get(field, ()):
        if source in by_source:
            return source, by_source[source]
    if base_source in by_source:
        return base_source, by_source[base_source]
    for source in SOURCE_ORDER:
        if source in by_source:
            return source, by_source[source]
    return candidates[0]


def build_comprehensive_peptide_evidence(
    *,
    output_tsv: str | Path,
    ranked_peptides: str | Path,
    annotated_peptides: str | Path | None = None,
    raw_peptides: str | Path | None = None,
    raw_events: str | Path | None = None,
    presentation_evidence: str | Path | None = None,
    appm_peptide_modifiers: str | Path | None = None,
    ccf_2: str | Path | None = None,
    expression_evidence: str | Path | None = None,
    rna_junction_evidence: str | Path | None = None,
    peptide_safety: str | Path | None = None,
    event_safety: str | Path | None = None,
    peptide_escape_flags: str | Path | None = None,
    validation_plan: str | Path | None = None,
    conflicts_tsv: str | Path | None = None,
) -> dict[str, Any]:
    base_rows, base_name = _base_rows(annotated_peptides, ranked_peptides)
    source_paths = {
        "ranked_peptides": ranked_peptides,
        "raw_peptides": raw_peptides,
        "raw_events": raw_events,
        "presentation_evidence": presentation_evidence,
        "appm_peptide_modifiers": appm_peptide_modifiers,
        "ccf_2": ccf_2,
        "expression_evidence": expression_evidence,
        "rna_junction_evidence": rna_junction_evidence,
        "peptide_safety": peptide_safety,
        "event_safety": event_safety,
        "peptide_escape_flags": peptide_escape_flags,
        "validation_plan": validation_plan,
    }
    source_rows = {name: _read_optional(path) for name, path in source_paths.items()}
    peptide_indexes = {name: _index(source_rows[name], "peptide_id") for name in PEPTIDE_SOURCES}
    event_indexes = {name: _event_index(source_rows[name]) for name in EVENT_SOURCES}

    output_rows: list[dict[str, str]] = []
    conflict_rows: list[dict[str, str]] = []
    source_counts: dict[str, int] = {name: 0 for name in (base_name, *source_paths)}
    for base in base_rows:
        peptide_id = str(base.get("peptide_id") or "")
        event_id = str(base.get("event_id") or "")
        matches: list[tuple[str, dict[str, str]]] = [(base_name, base)]
        ranked = peptide_indexes["ranked_peptides"].get(peptide_id, {})
        if not event_id and ranked.get("event_id"):
            event_id = ranked["event_id"]
        for name in PEPTIDE_SOURCES:
            evidence = peptide_indexes[name].get(peptide_id)
            if evidence and name != base_name:
                matches.append((name, evidence))
        for name in EVENT_SOURCES:
            evidence = event_indexes[name].get(event_id)
            if evidence:
                if name == "event_safety":
                    evidence = {
                        (f"event_{field}" if field in SAFETY_FIELDS else field): value
                        for field, value in evidence.items()
                    }
                matches.append((name, evidence))

        evidence_sources: list[str] = []
        candidates_by_field: dict[str, list[tuple[str, str]]] = {}
        for source, evidence in matches:
            if source not in evidence_sources:
                evidence_sources.append(source)
                source_counts[source] = source_counts.get(source, 0) + 1
            for field, value in evidence.items():
                text = str(value or "").strip()
                if text:
                    candidates_by_field.setdefault(field, []).append((source, text))

        row: dict[str, str] = {}
        conflict_fields: set[str] = set()
        conflict_details: list[dict[str, str]] = []
        selected_overrides: list[str] = []
        for field, candidates in candidates_by_field.items():
            selected_source, selected_value = _select_value(field, candidates, base_name)
            row[field] = selected_value
            if selected_source != base_name:
                selected_overrides.append(f"{field}:{selected_source}")
            allowed_sources = set(FIELD_SOURCE_PRECEDENCE.get(field, ()))
            conflict_candidates = [
                (source, value) for source, value in candidates
                if not allowed_sources or source in allowed_sources
            ]
            if field.endswith("_source") or field in {"source_tool", "comprehensive_evidence_sources"}:
                conflict_candidates = []
            seen_other: set[tuple[str, str]] = set()
            for other_source, other_value in conflict_candidates:
                if other_source == selected_source or _equivalent(selected_value, other_value):
                    continue
                marker = (other_source, other_value)
                if marker in seen_other:
                    continue
                seen_other.add(marker)
                conflict_fields.add(field)
                detail = {
                    "peptide_id": peptide_id,
                    "event_id": event_id,
                    "gene": str(row.get("gene") or base.get("gene") or ""),
                    "field": field,
                    "selected_source": selected_source,
                    "selected_value": selected_value,
                    "other_source": other_source,
                    "other_value": other_value,
                    "precedence_version": EVIDENCE_SOURCE_PRECEDENCE_VERSION,
                    "conflict_type": "NONEMPTY_SOURCE_DISAGREEMENT",
                }
                conflict_rows.append(detail)
                conflict_details.append(detail)

        if event_id and not row.get("event_id"):
            row["event_id"] = event_id
        row["comprehensive_evidence_sources"] = ",".join(evidence_sources)
        row["comprehensive_evidence_source_count"] = str(len(evidence_sources))
        row["comprehensive_evidence_status"] = (
            "COMPLETE" if "presentation_evidence" in evidence_sources and "ranked_peptides" in evidence_sources else "PARTIAL"
        )
        row["evidence_source_precedence_version"] = EVIDENCE_SOURCE_PRECEDENCE_VERSION
        row["evidence_conflict_fields"] = ",".join(sorted(conflict_fields))
        row["evidence_conflict_count"] = str(len(conflict_details))
        row["evidence_conflict_details"] = json.dumps(conflict_details, ensure_ascii=True, separators=(",", ":")) if conflict_details else "[]"
        row["evidence_selected_source_overrides"] = ";".join(sorted(selected_overrides))
        output_rows.append(row)

    fields = _field_order(base_rows)
    for source_name in SOURCE_ORDER:
        for field in _field_order(source_rows.get(source_name, [])):
            if field not in fields:
                fields.append(field)
    for field in (
        "comprehensive_evidence_sources", "comprehensive_evidence_source_count",
        "comprehensive_evidence_status", "evidence_source_precedence_version",
        "evidence_conflict_fields", "evidence_conflict_count", "evidence_conflict_details",
        "evidence_selected_source_overrides",
    ):
        if field not in fields:
            fields.append(field)
    for field in _field_order(output_rows):
        if field not in fields:
            fields.append(field)

    output_path = Path(output_tsv)
    conflict_path = Path(conflicts_tsv) if conflicts_tsv else output_path.with_name("evidence_conflicts.tsv")
    write_tsv(output_path, output_rows, fields)
    write_tsv(conflict_path, conflict_rows, CONFLICT_FIELDS)
    return {
        "output_tsv": str(output_path),
        "conflicts_tsv": str(conflict_path),
        "rows": len(output_rows),
        "columns": len(fields),
        "conflicts": len(conflict_rows),
        "conflict_peptides": sum(bool(row["evidence_conflict_fields"]) for row in output_rows),
        "precedence_version": EVIDENCE_SOURCE_PRECEDENCE_VERSION,
        "base_source": base_name,
        "source_counts": source_counts,
    }
