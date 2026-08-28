"""Strict adapter for independently validated high-order splice evidence."""
from __future__ import annotations

from pathlib import Path

from neoag.splice.identifiers import stable_id

from .base import clean, get, read_delimited, row_hash, source_record_id


ALLOWED_GROUPS = {"LONG_READ", "DNA_CAUSAL", "PROTEIN_VALIDATION", "LIGANDOME"}
SUPPORTED_STATES = {"SUPPORTED", "CONFIRMED", "PASS", "DETECTED", "VALIDATED"}
ENTITY_TABLES = {
    "JUNCTION": ("junctions", "junction_id"),
    "SPLICE_EVENT": ("events", "splice_event_id"),
    "TRANSCRIPT_HYPOTHESIS": ("transcripts", "transcript_hypothesis_id"),
    "ORF": ("orfs", "orf_id"),
    "PEPTIDE_ORIGIN": ("peptide_origins", "origin_peptide_id"),
}


def parse_high_order_evidence(
    path: str | Path,
    *,
    sample_id: str,
    entity_bundle: dict[str, list[dict[str, str]]],
    strict: bool = False,
) -> dict[str, list[dict[str, str]]]:
    """Accept evidence only through an exact, existing formal entity ID."""
    p = Path(path)
    registries = {
        entity_type: {row.get(id_field, "") for row in entity_bundle.get(table, []) if row.get(id_field)}
        for entity_type, (table, id_field) in ENTITY_TABLES.items()
    }
    evidence: list[dict[str, str]] = []
    conflicts: list[dict[str, str]] = []
    for row_no, row in enumerate(read_delimited(p), start=2):
        entity_type = get(row, "entity_type").upper()
        entity_id = get(row, "entity_id")
        group = get(row, "evidence_group").upper()
        status = get(row, "evidence_status", "status").upper()
        tool = get(row, "source_tool", "tool")
        version = get(row, "source_tool_version", "tool_version")
        record_id = get(row, "source_record_id") or source_record_id(tool or "HighOrderEvidence", p, row_no, row)
        reasons: list[str] = []
        if entity_type not in registries:
            reasons.append("unsupported entity_type")
        elif entity_id not in registries[entity_type]:
            reasons.append("entity_id does not exactly match the formal layer")
        if group not in ALLOWED_GROUPS:
            reasons.append("unsupported evidence_group")
        if status not in SUPPORTED_STATES:
            reasons.append("evidence status is not an accepted positive state")
        if not tool:
            reasons.append("source_tool is missing")
        if strict and clean(version).upper() in {"", "UNASSESSED", "UNKNOWN", "NA", "N/A"}:
            reasons.append("source_tool_version is not locked")
        if reasons:
            conflict = {
                "entity_type": entity_type or "HIGH_ORDER_EVIDENCE", "entity_id": entity_id or record_id,
                "sample_id": sample_id, "conflict_type": "HIGH_ORDER_EVIDENCE_UNRESOLVED",
                "field_name": "entity_id/evidence_group/status", "observed_values": row_hash(row),
                "source_tools": tool or "UNASSESSED", "source_record_ids": record_id,
                "severity": "ERROR" if strict else "WARNING", "resolution_status": "UNRESOLVED",
                "resolution_reason": "; ".join(reasons),
            }
            conflict["conflict_id"] = stable_id("CFL", conflict)
            conflicts.append(conflict)
            continue
        evidence_row = {
            "entity_type": entity_type, "entity_id": entity_id, "sample_id": sample_id,
            "evidence_group": group, "evidence_type": get(row, "evidence_type", default=f"{group}_VALIDATION"),
            "source_tool": tool, "source_tool_version": version,
            "source_assay_id": get(row, "source_assay_id", "assay_id") or "ASSAY_UNRESOLVED",
            "source_file": str(p), "source_row_number": str(row_no), "source_record_id": record_id,
            "provided_value": status, "verified_value": entity_id,
            "resolution_status": "RESOLVED_EXACT",
            "resolution_reason": "High-order evidence was linked through an exact formal entity identifier.",
            "raw_payload_sha256": row_hash(row),
        }
        evidence_row["evidence_id"] = stable_id("EVD", evidence_row)
        evidence.append(evidence_row)
    return {"tool_evidence": evidence, "conflicts": conflicts}
