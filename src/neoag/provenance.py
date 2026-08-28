"""Deterministic provenance-preserving entity merge for NeoAg v0.4.4.

The v0.4.3 production runner discarded every duplicate row after the first one.
That behaviour hid independent caller support and made later auditing impossible.
This module keeps one entity row for backward compatibility, while materialising
one provenance row per input record and one conflict row per incompatible field.

Important safety rule
---------------------
Numeric evidence may be aggregated only after callers have established a valid
entity key.  For splice evidence that key must be the canonical junction ID (or
an event/peptide key containing it); this module never performs gene-level or
nearest-locus matching.
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from typing import Iterable, Mapping

from .utils import MISSING, to_float

PROVENANCE_FIELDS = [
    "entity_type",
    "merge_key",
    "source_row_index",
    "stage_name",
    "stage_source",
    "source_tool",
    "source_file",
    "source_row_number",
    "source_record_id",
    "source_row_sha256",
]

CONFLICT_FIELDS = [
    "entity_type",
    "merge_key",
    "field",
    "selected_value",
    "observed_values",
    "source_tools",
    "conflict_type",
]

# Values in these fields are semantic sets rather than scalar values.
_LIST_FIELDS = {
    "source_tools",
    "source_records",
    "rna_junction_source",
    "component_event_ids",
    "cancer_gene_sources",
}

# These fields may legitimately vary between exact measurements of the same
# entity. The maximum is retained for backward-compatible ranking, while all
# measurements remain visible in the provenance/conflict tables.
_NUMERIC_MAX_FIELDS = {
    "rna_junction_reads",
    "provided_rna_junction_reads",
    "rna_alt_reads",
    "rna_depth",
    "tumor_depth",
    "tumor_alt_count",
    "phase_support_reads",
    "phase_total_informative_reads",
}

# Multiple source files/row numbers are retained in the entity row as a
# semicolon-separated audit hint. Full one-row-per-source detail is in the
# provenance table.
_AUDIT_LIST_FIELDS = {"source_file", "source_row_number"}

# These fields describe how a source represented an already identical entity.
# Differences are retained in long-form provenance and are not biological
# conflicts after the canonical merge key has been established.
_NON_CONFLICT_SOURCE_VARIATION_FIELDS = {
    "event_name",
    "pos",
    "source_junction_id",
    "junction_resolution_status",
    "junction_resolution_reason",
    "junction_match_status",
    "junction_match_method",
    "junction_support_reason",
}


def row_sha256(row: Mapping[str, object]) -> str:
    payload = json.dumps(
        {
            str(key): "" if value is None else str(value)
            for key, value in sorted(row.items(), key=lambda item: str(item[0]))
            if not str(key).startswith("_provenance_")
        },
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _present(value: object) -> bool:
    return str(value or "").strip() not in MISSING


def _tokens(value: object) -> list[str]:
    raw = str(value or "").replace(",", ";")
    return [
        part.strip()
        for part in raw.split(";")
        if part.strip() and part.strip() not in MISSING
    ]


def _source_tools(row: Mapping[str, object]) -> list[str]:
    """Return explicit tool/stage provenance, falling back to legacy ``source``.

    ``source`` is often a descriptive pipeline label (for example
    ``splice_source:RegTools``), not an independent tool.  Mixing it with
    explicit ``source_tool(s)`` inflated caller counts in v0.4.3.
    """

    values: list[str] = []
    explicit = (
        row.get("source_tool"),
        row.get("source_tools"),
        row.get("_provenance_stage_source"),
    )
    for candidate in explicit:
        for token in _tokens(candidate):
            if token not in values:
                values.append(token)
    if not values:
        for token in _tokens(row.get("source")):
            if token not in values:
                values.append(token)
    return values


def _source_record_ids(row: Mapping[str, object]) -> list[str]:
    """Return explicit source records; use entity IDs only as a legacy fallback."""

    values: list[str] = []
    for candidate in (row.get("source_record_id"), row.get("source_records")):
        for token in _tokens(candidate):
            if token not in values:
                values.append(token)
    if not values:
        for candidate in (row.get("event_id"), row.get("peptide_id")):
            for token in _tokens(candidate):
                if token not in values:
                    values.append(token)
    return values


def _merge_key(row: Mapping[str, object], key_fields: tuple[str, ...]) -> tuple[str, ...]:
    key = tuple(str(row.get(field) or "").strip() for field in key_fields)
    if any(key):
        return key
    return (f"UNKEYED:{row_sha256(row)}",)


def merge_rows_preserving_provenance(
    rows: Iterable[Mapping[str, object]],
    fields: list[str],
    key_fields: tuple[str, ...],
    *,
    entity_type: str,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    """Merge duplicate entities without discarding source records.

    Returns ``(merged_rows, provenance_rows, conflict_rows)``.

    ``key_fields`` must identify a biologically equivalent entity. For splice
    rows this means a canonical junction-based event ID, never a gene name.
    """

    groups: dict[tuple[str, ...], list[dict[str, object]]] = defaultdict(list)
    provenance: list[dict[str, str]] = []

    for index, raw in enumerate(rows, 1):
        row = dict(raw)
        key = _merge_key(row, key_fields)
        groups[key].append(row)
        merge_key = "|".join(key)
        source_tools = _source_tools(row)
        record_ids = _source_record_ids(row)
        provenance.append(
            {
                "entity_type": entity_type,
                "merge_key": merge_key,
                "source_row_index": str(index),
                "stage_name": str(row.get("_provenance_stage_name") or ""),
                "stage_source": str(row.get("_provenance_stage_source") or ""),
                "source_tool": ";".join(source_tools),
                "source_file": str(
                    row.get("_provenance_source_file")
                    or row.get("source_file")
                    or ""
                ),
                "source_row_number": str(row.get("source_row_number") or ""),
                "source_record_id": ";".join(record_ids),
                "source_row_sha256": row_sha256(row),
            }
        )

    merged: list[dict[str, str]] = []
    conflicts: list[dict[str, str]] = []

    for key in sorted(groups, key=lambda value: tuple(str(part) for part in value)):
        group = groups[key]
        merge_key = "|".join(key)
        out: dict[str, str] = {field: "" for field in fields}
        tools = sorted({token for row in group for token in _source_tools(row)})
        record_ids = sorted({token for row in group for token in _source_record_ids(row)})
        biological_conflict = False
        numeric_variation = False

        for field in fields:
            if field in {
                "source_tools",
                "source_records",
                "provenance_record_count",
                "evidence_conflict_status",
            }:
                continue

            values = [
                str(row.get(field) or "").strip()
                for row in group
                if _present(row.get(field))
            ]
            unique: list[str] = []
            for value in values:
                if value not in unique:
                    unique.append(value)
            if not unique:
                continue

            if field in _LIST_FIELDS or field in _AUDIT_LIST_FIELDS:
                out[field] = ";".join(
                    sorted({token for value in unique for token in _tokens(value)})
                )
                continue

            if field in _NUMERIC_MAX_FIELDS:
                numeric = [(to_float(value, float("-inf")), value) for value in unique]
                selected = max(numeric, key=lambda item: item[0])[1]
                out[field] = selected
                if len(unique) > 1:
                    numeric_variation = True
                    conflicts.append(
                        {
                            "entity_type": entity_type,
                            "merge_key": merge_key,
                            "field": field,
                            "selected_value": selected,
                            "observed_values": ";".join(unique),
                            "source_tools": ";".join(tools),
                            "conflict_type": "NUMERIC_EVIDENCE_VARIATION",
                        }
                    )
                continue

            # Preserve deterministic source precedence. Source-specific
            # representation differences remain visible in provenance but do
            # not make a canonical entity biologically discordant.
            out[field] = unique[0]
            if (
                len(unique) > 1
                and field not in {"source", "source_tool"}
                and field not in _NON_CONFLICT_SOURCE_VARIATION_FIELDS
            ):
                biological_conflict = True
                conflicts.append(
                    {
                        "entity_type": entity_type,
                        "merge_key": merge_key,
                        "field": field,
                        "selected_value": unique[0],
                        "observed_values": ";".join(unique),
                        "source_tools": ";".join(tools),
                        "conflict_type": "FIELD_CONFLICT_FIRST_NONEMPTY_RETAINED",
                    }
                )

        if "source_tool" in out and not out["source_tool"] and tools:
            out["source_tool"] = tools[0]
        if "source" in out and not out["source"] and tools:
            out["source"] = tools[0]
        if "source_tools" in out:
            out["source_tools"] = ";".join(tools)
        if "source_records" in out:
            out["source_records"] = ";".join(record_ids)
        if "source_record_id" in out and not out["source_record_id"] and record_ids:
            out["source_record_id"] = record_ids[0]
        if "provenance_record_count" in out:
            out["provenance_record_count"] = str(len(group))
        if "evidence_conflict_status" in out:
            if biological_conflict:
                out["evidence_conflict_status"] = "FIELD_CONFLICT"
            elif numeric_variation:
                out["evidence_conflict_status"] = "NUMERIC_EVIDENCE_VARIATION"
            else:
                out["evidence_conflict_status"] = "NONE"

        merged.append(out)

    return merged, provenance, conflicts
