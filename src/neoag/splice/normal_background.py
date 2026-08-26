"""Coverage-aware normal-background registration."""
from __future__ import annotations

from pathlib import Path

from neoag.splice.coordinates import iter_junction_records
from neoag.splice.identifiers import stable_id

from .adapters.base import as_float_text, as_int, get, read_delimited, source_record_id


_ADEQUATE_COVERAGE = {"ADEQUATE", "PASS", "SUFFICIENT", "LOCUS_COVERAGE_ADEQUATE"}
_LOW_COVERAGE = {"LOW", "LOW_COVERAGE", "INSUFFICIENT", "FAIL", "LOCUS_COVERAGE_LOW"}
_NEGATIVE = {"NOT_DETECTED", "ABSENT", "NEGATIVE"}
_DETECTED = {"DETECTED", "PRESENT", "POSITIVE"}
_LOW_LEVEL = {"LOW_LEVEL", "LOW_LEVEL_DETECTED", "TRACE", "WEAK_POSITIVE"}

NORMAL_DETECTED_ASSESSMENTS = {
    "NORMAL_DETECTED", "DETECTED_MATCHED_NORMAL", "DETECTED_CRITICAL_TISSUE",
    "DETECTED_BROAD_NORMAL", "LOW_LEVEL_NONCRITICAL_NORMAL",
}


def assessment_is_detected(row: dict[str, str]) -> bool:
    return row.get("assessment_status", "") in NORMAL_DETECTED_ASSESSMENTS


def assessment_is_critical(row: dict[str, str]) -> bool:
    return row.get("assessment_status") == "DETECTED_CRITICAL_TISSUE" or (
        row.get("assessment_status") == "NORMAL_DETECTED"
        and _is_true(row.get("critical_tissue", ""))
    )


def _is_true(value: str | bool) -> bool:
    return str(value).strip().casefold() in {"1", "true", "yes", "y"}


def _assessment_status(
    *,
    detection: str,
    coverage: str,
    source_type: str,
    critical_tissue: str | bool,
) -> tuple[str, str]:
    detection = str(detection or "UNASSESSED").upper()
    coverage = str(coverage or "UNASSESSED").upper()
    source_type_upper = str(source_type or "").upper()
    detected = detection in _DETECTED | _LOW_LEVEL
    if detected and _is_true(critical_tissue):
        return "DETECTED_CRITICAL_TISSUE", "Normal signal detected in a critical tissue."
    if detected and source_type_upper in {"MATCHED_NORMAL", "PATIENT_MATCHED_NORMAL", "MATCHED_NORMAL_RNA"}:
        return "DETECTED_MATCHED_NORMAL", "Signal detected in the patient-matched normal sample."
    if detection in _LOW_LEVEL:
        return "LOW_LEVEL_NONCRITICAL_NORMAL", "Low-level normal signal detected outside a classified critical tissue."
    if detected:
        return "DETECTED_BROAD_NORMAL", "Signal detected in a normal tissue, cohort, or panel."
    if detection in _NEGATIVE and coverage in _ADEQUATE_COVERAGE:
        return "NOT_DETECTED_ADEQUATE_COVERAGE", "No signal detected with explicit adequate locus coverage."
    if detection in _NEGATIVE and coverage in _LOW_COVERAGE:
        return "NOT_DETECTED_LOW_COVERAGE", "No signal detected, but locus coverage was insufficient."
    return "UNASSESSED", "Normal absence cannot be assessed without an explicit adequate-coverage result."


def parse_normal_junctions(
    path: str | Path,
    *,
    sample_id: str,
    genome_build: str = "GRCh38",
    coordinate_system: str = "auto",
    source_name: str = "NormalPanel",
    source_type: str = "PROTOCOL_MATCHED_NORMAL_PANEL",
    tissue: str = "",
    critical_tissue: bool = False,
    strict: bool = False,
    allowed_junction_ids: set[str] | None = None,
) -> dict[str, list[dict[str, str]]]:
    p = Path(path)
    rows: list[dict[str, str]] = []
    conflicts: list[dict[str, str]] = []
    evidence: list[dict[str, str]] = []
    scanned = 0
    retained = 0
    for record in iter_junction_records(
        p, sample_id=sample_id, source_tool=source_name,
        genome_build=genome_build, coordinate_system=coordinate_system, strict=strict,
    ):
        scanned += 1
        if record.junction is None:
            # In targeted mode, unresolved background records cannot match an
            # exact candidate and need not create millions of conflicts.
            if allowed_junction_ids is not None:
                continue
            conflicts.append({
                "entity_type": "NORMAL_BACKGROUND", "entity_id": record.source_record_id,
                "sample_id": sample_id, "conflict_type": "NORMAL_JUNCTION_UNRESOLVED",
                "field_name": "coordinates", "observed_values": str(record.row),
                "source_tools": source_name, "source_record_ids": record.source_record_id,
                "severity": "ERROR" if strict else "WARNING", "resolution_status": "UNRESOLVED",
                "resolution_reason": record.coordinate_warning or "No canonical junction.",
            })
            continue
        jid = record.junction.junction_id
        if allowed_junction_ids is not None and jid not in allowed_junction_ids:
            continue
        retained += 1
        detected = record.total_split_reads > 0
        detection = "DETECTED" if detected else "NOT_DETECTED"
        coverage = "LOCUS_COVERAGE_UNASSESSED"
        assessment, reason = _assessment_status(
            detection=detection, coverage=coverage, source_type=source_type, critical_tissue=critical_tissue,
        )
        background_id = stable_id("NBG", source_name, jid, tissue, source_type)
        rows.append({
            "normal_background_id": background_id, "splice_event_id": "", "junction_id": jid,
            "origin_peptide_id": "", "sample_id": sample_id, "normal_source": source_name,
            "normal_source_type": source_type, "normal_tissue": tissue,
            "critical_tissue": "true" if critical_tissue else "false",
            "detection_status": detection,
            "coverage_status": coverage,
            "junction_reads": str(record.total_split_reads), "sample_prevalence": "",
            "kmer_prevalence": "", "assessment_status": assessment,
            "assessment_reason": reason,
            "source_file": str(p), "source_record_id": record.source_record_id,
            "evidence_conflict_status": "NONE",
        })
        evidence.append({
            "entity_type": "JUNCTION", "entity_id": jid, "sample_id": sample_id,
            "evidence_group": "NORMAL_BACKGROUND", "evidence_type": "NORMAL_JUNCTION_DETECTION",
            "source_tool": source_name, "source_tool_version": "UNASSESSED", "source_file": str(p),
            "source_row_number": str(record.source_row_number), "source_record_id": record.source_record_id,
            "provided_value": str(record.total_split_reads), "verified_value": str(record.total_split_reads),
            "resolution_status": "RESOLVED_EXACT", "resolution_reason": rows[-1]["assessment_reason"],
            "raw_payload_sha256": record.record_sha256,
        })
    return {
        "normal_background": rows,
        "tool_evidence": evidence,
        "conflicts": conflicts,
        "manifest": [{
            "adapter": "normal_junctions",
            "input_path": str(p),
            "rows_scanned": str(scanned),
            "rows_retained": str(retained),
            "filter_policy": "EXACT_CANONICAL_JUNCTION" if allowed_junction_ids is not None else "UNFILTERED",
            "target_junction_count": str(len(allowed_junction_ids or set())),
        }],
    }


def parse_normal_coverage(
    path: str | Path,
    *,
    sample_id: str,
    source_name: str = "NormalCoverage",
) -> dict[str, list[dict[str, str]]]:
    """Parse an explicit coverage-aware normal assessment table.

    Required identity is one of ``junction_id``, ``splice_event_id`` or
    ``origin_peptide_id``.  A negative assessment is accepted only when
    ``coverage_status`` explicitly denotes adequate coverage.
    """
    p = Path(path)
    rows: list[dict[str, str]] = []
    conflicts: list[dict[str, str]] = []
    for row_no, row in enumerate(read_delimited(p), start=2):
        rid = source_record_id(source_name, p, row_no, row)
        jid = get(row, "junction_id", "canonical_junction_id")
        event_id = get(row, "splice_event_id", "event_id")
        por = get(row, "origin_peptide_id")
        if not any((jid, event_id, por)):
            conflicts.append({
                "entity_type": "NORMAL_BACKGROUND", "entity_id": rid, "sample_id": sample_id,
                "conflict_type": "NORMAL_BACKGROUND_IDENTITY_MISSING", "field_name": "entity_id",
                "observed_values": "", "source_tools": source_name, "source_record_ids": rid,
                "severity": "ERROR", "resolution_status": "UNRESOLVED",
                "resolution_reason": "No exact junction/event/origin identifier was supplied.",
            })
            continue
        detection = get(row, "detection_status", "status", default="UNASSESSED").upper()
        coverage = get(row, "coverage_status", "coverage", default="UNASSESSED").upper()
        source_type = get(row, "normal_source_type", "source_type", default="EXPLICIT_COVERAGE_TABLE")
        critical = get(row, "critical_tissue", default="false")
        assessment, reason = _assessment_status(
            detection=detection, coverage=coverage, source_type=source_type, critical_tissue=critical,
        )
        rows.append({
            "normal_background_id": stable_id("NBG", source_name, jid, event_id, por, rid),
            "splice_event_id": event_id, "junction_id": jid, "origin_peptide_id": por,
            "sample_id": sample_id, "normal_source": source_name,
            "normal_source_type": source_type,
            "normal_tissue": get(row, "normal_tissue", "tissue"),
            "critical_tissue": critical,
            "detection_status": detection, "coverage_status": coverage,
            "junction_reads": str(as_int(get(row, "junction_reads", "reads"), 0)),
            "sample_prevalence": as_float_text(get(row, "sample_prevalence", "prevalence")),
            "kmer_prevalence": as_float_text(get(row, "kmer_prevalence")),
            "assessment_status": assessment, "assessment_reason": reason,
            "source_file": str(p), "source_record_id": rid, "evidence_conflict_status": "NONE",
        })
    return {"normal_background": rows, "conflicts": conflicts}
