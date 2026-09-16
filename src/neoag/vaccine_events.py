from __future__ import annotations

import re
from collections import defaultdict
from typing import Any, Mapping


def _text(row: Mapping[str, Any], *fields: str, default: str = "") -> str:
    for field in fields:
        value = str(row.get(field) or "").strip()
        if value:
            return value
    return default


def _integer(value: Any, default: int = 10**9) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def event_member_ids(row: Mapping[str, Any]) -> list[str]:
    values: list[str] = []
    for field in ("event_id", "source_event_id"):
        value = str(row.get(field) or "").strip()
        if value and value not in values:
            values.append(value)
    for value in re.split(r"[;,]", str(row.get("member_event_ids") or "")):
        value = value.strip()
        if value and value not in values:
            values.append(value)
    return values


def _event_kind(event_type: str) -> str:
    text = str(event_type or "").lower()
    if "fusion" in text:
        return "FUSION"
    if "splice" in text or "junction" in text:
        return "SPLICE"
    if any(token in text for token in ("frameshift", "indel", "insertion", "deletion")):
        return "FRAMESHIFT"
    if any(token in text for token in ("dna_sv", "structural", "bnd")):
        return "DNA_SV"
    if any(token in text for token in ("snv", "missense", "substitution")):
        return "MISSENSE"
    return "OTHER"


def _construct_strategy(event_type: str) -> str:
    return {
        "MISSENSE": "MUTATION_CENTERED_LONG_SEQUENCE_OR_MINIGENE",
        "FRAMESHIFT": "NOVEL_TAIL_LONG_SEQUENCE_OR_MINIGENE",
        "FUSION": "JUNCTION_CENTERED_LONG_SEQUENCE_OR_MINIGENE",
        "SPLICE": "JUNCTION_CENTERED_LONG_SEQUENCE_OR_MINIGENE",
        "DNA_SV": "RNA_CONFIRMED_ORF_REQUIRED_BEFORE_MINIGENE",
    }.get(_event_kind(event_type), "EVENT_CONTEXT_SEQUENCE_REVIEW_REQUIRED")


def build_vaccine_event_tables(
    events: list[dict[str, str]],
    peptides: list[dict[str, str]],
    review_rows: list[dict[str, str]] | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Build one vaccine-selection row per event and nested peptide-HLA evidence.

    Identical peptide-HLA predictions are collapsed within an event while their
    source peptide IDs and evidence-row count remain auditable.
    """
    reviews: dict[str, dict[str, str]] = {}
    for row in review_rows or []:
        for key in (str(row.get("event_group_id") or ""), *event_member_ids(row)):
            if key:
                reviews[key] = row

    summaries: list[dict[str, str]] = []
    details: list[dict[str, str]] = []
    for fallback_rank, event in enumerate(events, 1):
        event_rank = _text(event, "event_evidence_rank", "vaccine_event_rank", default=str(fallback_rank))
        event_group_id = _text(event, "event_group_id", "event_id")
        member_ids = event_member_ids(event)
        member_set = set(member_ids)
        matched = [row for row in peptides if _text(row, "event_id") in member_set]
        matched.sort(key=lambda row: _integer(row.get("evidence_rank")))

        grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for peptide in matched:
            sequence = _text(peptide, "peptide", "mutant_peptide", "mt_peptide").upper()
            hla = re.sub(r"\s+", "", _text(peptide, "hla_allele", "hla", "allele", "restricting_hla").upper())
            if sequence:
                grouped[(sequence, hla)].append(peptide)

        representatives: list[dict[str, str]] = []
        for epitope_rank, ((sequence, hla), members) in enumerate(grouped.items(), 1):
            representative = members[0]
            representatives.append(representative)
            detail = dict(representative)
            detail.update({
                "vaccine_event_rank": event_rank,
                "event_group_id": event_group_id,
                "vaccine_selection_unit": "MUTATION_OR_BIOLOGICAL_EVENT",
                "neoepitope_role": "SUPPORTING_PEPTIDE_HLA_EVIDENCE",
                "neoepitope_rank_within_event": str(epitope_rank),
                "peptide": sequence,
                "hla_allele": hla,
                "duplicate_evidence_row_count": str(len(members)),
                "member_peptide_ids": ";".join(sorted({_text(row, "peptide_id") for row in members if _text(row, "peptide_id")})),
            })
            details.append(detail)

        review = reviews.get(event_group_id) or next((reviews[key] for key in member_ids if key in reviews), {})
        best = representatives[0] if representatives else {}
        event_type = _text(event, "event_type", "biological_event_track", "evidence_track", default=_text(best, "event_type"))
        peptide_values = {_text(row, "peptide", "mutant_peptide", "mt_peptide").upper() for row in representatives}
        hla_values = {re.sub(r"\s+", "", _text(row, "hla_allele", "hla", "allele", "restricting_hla")) for row in representatives}
        hla_values.discard("")
        pair_labels = [f"{sequence}/{hla}" if hla else sequence for sequence, hla in grouped]
        summaries.append({
            "vaccine_event_rank": event_rank,
            "vaccine_selection_unit": "MUTATION_OR_BIOLOGICAL_EVENT",
            "event_group_id": event_group_id,
            "event_id": _text(event, "event_id", default=_text(best, "event_id")),
            "member_event_ids": ";".join(member_ids),
            "gene": _text(event, "gene", default=_text(best, "gene")),
            "event_type": event_type,
            "protein_change": _text(event, "protein_change", "hgvsp", default=_text(best, "protein_change", "hgvsp")),
            "transcript_id": _text(event, "transcript_id", default=_text(best, "transcript_id", "fusion_transcript_id")),
            "orf_id": _text(event, "orf_id", default=_text(best, "orf_id")),
            "best_evidence_grade": _text(event, "best_evidence_grade", default=_text(best, "evidence_grade")),
            "review_status": str(review.get("review_status") or "EVIDENCE_RANKED_NOT_SELECTED"),
            "experiment_priority": str(review.get("experiment_priority") or "UNASSESSED_FOR_EXPERIMENT"),
            "clonality_state": _text(event, "clonality_state", default=_text(best, "clonality_state", "clonality_status")),
            "ccf_estimate": _text(event, "ccf_estimate", default=_text(best, "ccf_estimate", "ccf_best")),
            "rna_support_state": _text(event, "rna_support_state", default=_text(best, "rna_support_state", "rna_support_status")),
            "raw_peptide_hla_row_count": str(len(matched)),
            "unique_neoepitope_hla_count": str(len(grouped)),
            "duplicate_evidence_rows_collapsed": str(len(matched) - len(grouped)),
            "unique_peptide_count": str(len(peptide_values)),
            "hla_coverage_count": str(len(hla_values)),
            "hla_alleles": ";".join(sorted(hla_values)),
            "neoepitope_hla_pairs": ";".join(pair_labels),
            "best_peptide_id": _text(event, "best_peptide_id", default=_text(best, "peptide_id")),
            "best_peptide": _text(event, "best_peptide", default=_text(best, "peptide", "mutant_peptide", "mt_peptide")),
            "best_hla_allele": _text(event, "best_hla_allele", default=_text(best, "hla_allele", "hla", "allele")),
            "construct_strategy": _construct_strategy(event_type),
            "construct_sequence_status": "DESIGN_REQUIRED_FROM_CONFIRMED_TRANSCRIPT_OR_ORF",
            "recommended_validation": str(review.get("recommended_validation") or ""),
            "selection_note": "Select and review at event level; peptide-HLA rows are supporting evidence, not independent vaccine targets.",
        })
    return summaries, details
