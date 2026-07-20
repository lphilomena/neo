from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

from .utils import read_tsv, write_tsv


PEPTIDE_SOURCES = (
    "ranked_peptides",
    "raw_peptides",
    "presentation_evidence",
    "appm_peptide_modifiers",
    "peptide_safety",
    "peptide_escape_flags",
    "validation_plan",
)

EVENT_SOURCES = (
    "raw_events",
    "ccf_2",
    "expression_evidence",
    "rna_junction_evidence",
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
        for field in ("junction_reads", "rna_alt_reads", "rna_depth"):
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


def _merge_missing(target: dict[str, str], source: dict[str, str]) -> None:
    for field, value in source.items():
        if value not in (None, "") and target.get(field) in (None, ""):
            target[field] = str(value)


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
    peptide_escape_flags: str | Path | None = None,
    validation_plan: str | Path | None = None,
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
        "peptide_escape_flags": peptide_escape_flags,
        "validation_plan": validation_plan,
    }
    source_rows = {name: _read_optional(path) for name, path in source_paths.items()}
    peptide_indexes = {
        name: _index(source_rows[name], "peptide_id")
        for name in PEPTIDE_SOURCES
    }
    event_indexes = {
        name: _event_index(source_rows[name])
        for name in EVENT_SOURCES
    }

    output_rows: list[dict[str, str]] = []
    source_counts: dict[str, int] = {name: 0 for name in (base_name, *source_paths)}
    for base in base_rows:
        row = {field: str(value or "") for field, value in base.items()}
        peptide_id = str(row.get("peptide_id") or "")
        event_id = str(row.get("event_id") or "")
        evidence_sources = [base_name]
        source_counts[base_name] = source_counts.get(base_name, 0) + 1

        ranked = peptide_indexes["ranked_peptides"].get(peptide_id, {})
        if not event_id and ranked.get("event_id"):
            event_id = ranked["event_id"]
            row["event_id"] = event_id

        for name in PEPTIDE_SOURCES:
            evidence = peptide_indexes[name].get(peptide_id)
            if not evidence:
                continue
            _merge_missing(row, evidence)
            if name not in evidence_sources:
                evidence_sources.append(name)
            source_counts[name] += 1

        for name in EVENT_SOURCES:
            evidence = event_indexes[name].get(event_id)
            if not evidence:
                continue
            _merge_missing(row, evidence)
            if name not in evidence_sources:
                evidence_sources.append(name)
            source_counts[name] += 1

        row["comprehensive_evidence_sources"] = ",".join(evidence_sources)
        row["comprehensive_evidence_source_count"] = str(len(evidence_sources))
        row["comprehensive_evidence_status"] = (
            "COMPLETE"
            if "presentation_evidence" in evidence_sources
            and "ranked_peptides" in evidence_sources
            else "PARTIAL"
        )
        output_rows.append(row)

    fields = _field_order(base_rows)
    for source_name in PEPTIDE_SOURCES:
        for field in _field_order(source_rows[source_name]):
            if field not in fields:
                fields.append(field)
    for source_name in EVENT_SOURCES:
        for field in _field_order(source_rows[source_name]):
            if field not in fields:
                fields.append(field)
    for field in (
        "comprehensive_evidence_sources",
        "comprehensive_evidence_source_count",
        "comprehensive_evidence_status",
    ):
        if field not in fields:
            fields.append(field)

    write_tsv(output_tsv, output_rows, fields)
    return {
        "output_tsv": str(output_tsv),
        "rows": len(output_rows),
        "columns": len(fields),
        "base_source": base_name,
        "source_counts": source_counts,
    }
