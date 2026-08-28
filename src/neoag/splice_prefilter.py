"""Auditable splice prefilter applied before expensive presentation models."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any, Callable, Mapping

from .schemas import PEPTIDE_FIELDS
from .utils import read_tsv, write_tsv


Decision = tuple[bool | None, str]

SPLICE_EVENT_TYPES = {"splice", "splice_junction", "junction"}
PASS_TOKENS = {"PASS", "PASSED", "RESOLVED", "EXACT", "EXACT_MATCH", "NOT_DETECTED", "ABSENT", "NOVEL", "UNANNOTATED", "IN_FRAME", "VALID", "FALSE", "NO"}
FAIL_TOKENS = {"FAIL", "FAILED", "UNRESOLVED", "CONFLICT", "DETECTED", "PRESENT", "KNOWN_NORMAL", "ANNOTATED_NORMAL", "OUT_OF_FRAME", "INVALID", "TRUE", "YES"}


def _text(rows: list[Mapping[str, Any]], *fields: str) -> str:
    for row in rows:
        for field in fields:
            value = str(row.get(field, "") or "").strip()
            if value and value.upper() not in {"NA", "N/A", "NONE", "UNASSESSED", "UNKNOWN", "NOT_AVAILABLE"}:
                return value
    return ""


def _number(rows: list[Mapping[str, Any]], *fields: str) -> float | None:
    value = _text(rows, *fields)
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _status(rows: list[Mapping[str, Any]], fields: tuple[str, ...], *, pass_tokens: set[str], fail_tokens: set[str]) -> Decision:
    value = _text(rows, *fields)
    if not value:
        return None, "not_recorded"
    upper = value.upper()
    def matches(token: str) -> bool:
        return upper == token if token in {"NO", "YES", "TRUE", "FALSE", "PASS", "FAIL"} else token in upper
    negated_pass_tokens = {"NOT_DETECTED", "NOT_FOUND", "NOT_PREDICTED", "ABSENT", "NEGATIVE", "FALSE", "NO", "ESCAPES_NMD", "UNANNOTATED"}
    if any(matches(token) for token in pass_tokens & negated_pass_tokens):
        return True, value
    if any(matches(token) for token in fail_tokens):
        return False, value
    if any(matches(token) for token in pass_tokens):
        return True, value
    return None, value


def _aggregate_exclusion_status(rows: list[Mapping[str, Any]], fields: tuple[str, ...], label: str) -> Decision:
    """Conservatively combine peptide-level exclusion evidence for an event.

    One exact normal match is sufficient to fail the event.  Absence is only
    accepted when every recorded value is an assessed negative; mixed or
    missing records remain unassessed.
    """
    values: list[str] = []
    for row in rows:
        for field in fields:
            value = str(row.get(field, "") or "").strip()
            if value and value.upper() not in {"NA", "N/A", "NONE", "UNKNOWN", "NOT_AVAILABLE"}:
                values.append(value)
    if not values:
        return None, f"{label}_not_recorded"
    uppers = [value.upper() for value in values]
    if any(
        any(token in value for token in ("DETECTED", "EXACT_MATCH", "FOUND", "PRESENT", "POSITIVE", "FAIL"))
        and not any(token in value for token in ("NOT_DETECTED", "NOT_FOUND"))
        for value in uppers
    ):
        return False, f"{label}=normal_exact_match_detected"
    assessed_negative = [
        value for value in uppers
        if any(token in value for token in ("NOT_DETECTED", "NOT_FOUND", "ABSENT", "NEGATIVE", "PASS"))
    ]
    unassessed = [value for value in uppers if "UNASSESSED" in value or "PARTIAL" in value]
    if assessed_negative and not unassessed and len(assessed_negative) == len(uppers):
        return True, f"{label}=all_recorded_peptides_negative"
    return None, f"{label}=mixed_or_incomplete_evidence"


def _minimum(rows: list[Mapping[str, Any]], fields: tuple[str, ...], threshold: float, label: str) -> Decision:
    value = _number(rows, *fields)
    if value is None:
        return None, f"{label}_not_recorded"
    return value >= threshold, f"{label}={value:g}; threshold={threshold:g}"


def _normal_absence_with_coverage(
    rows: list[Mapping[str, Any]],
    status_fields: tuple[str, ...],
    coverage_fields: tuple[str, ...],
    threshold: float,
    label: str,
) -> Decision:
    status = _text(rows, *status_fields)
    coverage = _number(rows, *coverage_fields)
    upper = status.upper()
    detected = any(token in upper for token in (
        "SEEN_IN_NORMAL_CATALOG", "DETECTED_BROAD_NORMAL", "DETECTED_CRITICAL_TISSUE",
        "DETECTED_MATCHED_NORMAL", "SUPPORTED_IN_NORMAL", "EXACT_MATCH", "PRESENT", "POSITIVE",
    )) and "NOT_DETECTED" not in upper
    if detected:
        return False, f"{label}={status}; normal_coverage={coverage}"
    adequate = any(token in upper for token in (
        "NOT_DETECTED_ADEQUATE_COVERAGE", "ADEQUATE_COVERAGE_NEGATIVE", "EXCLUDED_WITH_ADEQUATE_COVERAGE",
    ))
    absent = any(token in upper for token in ("NOT_DETECTED", "ABSENT", "NEGATIVE", "NOT_FOUND"))
    if adequate or (absent and coverage is not None and coverage >= threshold):
        return True, f"{label}={status}; normal_coverage={coverage}; threshold={threshold:g}"
    if "NOT_LISTED_IN_NORMAL_CATALOG" in upper:
        return None, f"{label}=NOT_LISTED_IN_NORMAL_CATALOG; presence-only catalog non-membership is not a coverage-qualified negative"
    if absent:
        return None, f"{label}={status}; normal_coverage_unassessed"
    if status:
        return None, f"{label}={status}; normal_coverage={coverage}"
    return None, f"{label}_not_recorded"


def _is_splice(row: Mapping[str, Any]) -> bool:
    event_type = str(row.get("event_type", "") or "").strip().lower()
    consequence = str(row.get("peptide_consequence", "") or "").strip().lower()
    if event_type:
        return event_type in SPLICE_EVENT_TYPES
    return consequence == "splice_junction"


def _event_key(row: Mapping[str, Any], index: int) -> str:
    return str(row.get("event_id") or row.get("splice_event_id") or f"splice-row-{index}")


def _stage_definitions(profile: Mapping[str, Any]) -> list[tuple[str, str, Callable[[list[Mapping[str, Any]]], Decision]]]:
    gates = profile.get("gates", {}) if isinstance(profile, Mapping) else {}
    min_unique = float(gates.get("min_splice_unique_reads", gates.get("min_rna_junction_reads", 3)))
    min_total = float(gates.get("min_splice_total_coverage", min_unique))
    min_psi = float(gates.get("min_splice_psi", 0.05))
    min_normal_coverage = float(gates.get("min_splice_normal_coverage", 10))
    return [
        (
            "ALIGNMENT_COORDINATE_QC",
            "exact/resolved junction coordinates and no mapping conflict",
            lambda rows: _status(
                rows,
                ("splice_alignment_qc_status", "junction_resolution_status", "junction_match_status"),
                pass_tokens={"PASS", "RESOLVED", "EXACT", "PRIMARY_SOURCE_RECORD"},
                fail_tokens={"FAIL", "UNRESOLVED", "UNSTRANDED", "CONFLICT", "AMBIGUOUS"},
            ),
        ),
        (
            "FORMAL_JUNCTION_READ_QC",
            "complete junction read QC passed (reads, unique starts, anchor, MAPQ, multimapping, PSI and caller filter)",
            lambda rows: _status(
                rows,
                ("junction_read_qc_status",),
                pass_tokens={"PASS"},
                fail_tokens={"FAIL", "INCOMPLETE"},
            ),
        ),
        (
            "UNIQUE_JUNCTION_READS",
            f"explicit unique junction reads >= {min_unique:g}; total/caller reads are not substituted",
            lambda rows: _minimum(rows, ("unique_junction_reads", "unique_split_reads", "junction_unique_reads"), min_unique, "unique_junction_reads"),
        ),
        (
            "TOTAL_JUNCTION_COVERAGE",
            f"explicit total junction coverage >= {min_total:g}",
            lambda rows: _minimum(rows, ("junction_total_coverage", "total_junction_coverage", "total_split_reads", "splice_total_coverage"), min_total, "junction_total_coverage"),
        ),
        (
            "PSI",
            f"tumor splice PSI >= {min_psi:g}",
            lambda rows: _minimum(rows, ("splice_psi", "junction_psi", "tumor_psi", "psi"), min_psi, "splice_psi"),
        ),
        (
            "MATCHED_NORMAL_JUNCTION",
            f"junction absent from matched/adjacent normal sample with coverage >= {min_normal_coverage:g}",
            lambda rows: _normal_absence_with_coverage(
                rows,
                ("matched_normal_junction_status", "patient_normal_junction_status"),
                ("matched_normal_junction_coverage", "matched_normal_depth", "normal_junction_depth"),
                min_normal_coverage,
                "matched_normal_junction",
            ),
        ),
        (
            "NORMAL_COHORT_JUNCTION",
            "junction absent from a coverage-aware normal cohort; catalog non-membership alone is unassessed",
            lambda rows: _normal_absence_with_coverage(
                rows,
                ("normal_cohort_junction_status", "normal_junction_assessment_status", "normal_junction_status"),
                ("normal_cohort_junction_coverage", "normal_junction_coverage", "normal_junction_depth"),
                min_normal_coverage,
                "normal_cohort_junction",
            ),
        ),
        (
            "ANNOTATED_NORMAL_ISOFORM",
            "not a known/annotated normal isoform",
            lambda rows: _status(
                rows,
                ("annotated_normal_isoform_status", "splice_annotation_status", "known_junction"),
                pass_tokens={"NOVEL", "UNANNOTATED", "ABSENT", "NOT_DETECTED", "FALSE", "NO"},
                fail_tokens={"KNOWN_NORMAL", "ANNOTATED_NORMAL", "KNOWN", "PRESENT", "TRUE", "YES"},
            ),
        ),
        (
            "CREDIBLE_ORF",
            "credible translated ORF/frame",
            lambda rows: _status(
                rows,
                ("splice_orf_status", "orf_evidence_grade", "rna_frame_status"),
                pass_tokens={"PASS", "VALID", "IN_FRAME", "TRANSLATED", "O1", "O2"},
                fail_tokens={"FAIL", "INVALID", "OUT_OF_FRAME", "UNTRANSLATED", "O4"},
            ),
        ),
        (
            "NMD",
            "not predicted to undergo nonsense-mediated decay",
            lambda rows: _status(
                rows,
                ("splice_nmd_status", "nmd_status", "predicted_nmd"),
                pass_tokens={"NOT_PREDICTED", "ESCAPES_NMD", "PASS", "FALSE", "NO"},
                fail_tokens={"NMD", "PREDICTED_NMD", "FAIL", "TRUE", "YES"},
            ),
        ),
        (
            "JUNCTION_SPANNING_PEPTIDE",
            "peptide contains residues from both sides of the abnormal junction",
            lambda rows: _status(
                rows,
                ("crosses_junction", "junction_spanning_status"),
                pass_tokens={"YES", "TRUE", "PASS", "JUNCTION_SPANNING"},
                fail_tokens={"NO", "FALSE", "FAIL", "ONE_SIDED"},
            ),
        ),
        (
            "NORMAL_PROTEOME_EXCLUSION",
            "full peptide not detected in the configured normal proteome",
            lambda rows: _aggregate_exclusion_status(
                rows,
                ("normal_proteome_exact_match_status", "reference_proteome_status", "reference_proteome_exact_match"),
                "normal_proteome",
            ),
        ),
        (
            "NORMAL_TRANSCRIPTOME_EXCLUSION",
            "full junction peptide/transcript not detected in the configured normal transcriptome",
            lambda rows: _aggregate_exclusion_status(
                rows,
                ("normal_transcriptome_exact_match_status", "normal_transcriptome_status"),
                "normal_transcriptome",
            ),
        ),
    ]


def prefilter_splice_peptides(
    raw_peptides: str | Path,
    profile: Mapping[str, Any],
    outdir: str | Path,
    raw_events: str | Path | None = None,
) -> dict[str, Any]:
    """Filter splice events before presentation while retaining all rows for audit."""
    raw_path = Path(raw_peptides)
    output_dir = Path(outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(raw_path)
    non_splice: list[dict[str, str]] = []
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    event_evidence: dict[str, list[dict[str, str]]] = defaultdict(list)
    event_order: list[str] = []
    if raw_events and Path(raw_events).is_file():
        for index, row in enumerate(read_tsv(raw_events)):
            if not _is_splice(row):
                continue
            key = _event_key(row, index)
            if key not in event_evidence:
                event_order.append(key)
            event_evidence[key].append(row)
    for index, row in enumerate(rows):
        if not _is_splice(row):
            non_splice.append(row)
            continue
        key = _event_key(row, index)
        if key not in grouped and key not in event_evidence:
            event_order.append(key)
        grouped[key].append(row)

    archive_path = output_dir / "raw_peptides.before_splice_prefilter.tsv"
    write_tsv(archive_path, rows, PEPTIDE_FIELDS)
    stages = _stage_definitions(profile)
    decisions: list[dict[str, str]] = []
    stage_results: dict[str, dict[str, Decision]] = {name: {} for name, _, _ in stages}
    status_by_event: dict[str, str] = {}
    score_by_event: dict[str, tuple[float, float, float, float]] = {}
    for event_id in event_order:
        event_rows = event_evidence.get(event_id, []) + grouped.get(event_id, [])
        failures: list[str] = []
        missing: list[str] = []
        details: list[str] = []
        for stage, _, evaluator in stages:
            passed, reason = evaluator(event_rows)
            stage_results[stage][event_id] = (passed, reason)
            details.append(f"{stage}={reason}")
            if passed is False:
                failures.append(stage)
            elif passed is None:
                missing.append(stage)
        status = "REJECT" if failures else "PASS" if not missing else "REVIEW"
        status_by_event[event_id] = status
        score_by_event[event_id] = (
            _number(event_rows, "unique_junction_reads", "unique_split_reads", "junction_unique_reads") or -1.0,
            _number(event_rows, "junction_total_coverage", "total_junction_coverage", "total_split_reads") or -1.0,
            _number(event_rows, "splice_psi", "junction_psi", "tumor_psi", "psi") or -1.0,
            _number(event_rows, "rna_junction_reads", "junction_reads", "provided_rna_junction_reads") or -1.0,
        )
        decisions.append({
            "event_id": event_id,
            "prefilter_status": status,
            "splice_candidate_pool": (
                "FORMAL_SPLICE_CANDIDATE" if status == "PASS"
                else "REJECTED_SPLICE_CANDIDATE" if status == "REJECT"
                else "EXPLORATION_EVIDENCE_INCOMPLETE"
            ),
            "splice_formal_gate_pass": "yes" if status == "PASS" else "no",
            "failed_stages": ";".join(failures),
            "unassessed_stages": ";".join(missing),
            "stage_details": ";".join(details),
        })

    gates = profile.get("gates", {}) if isinstance(profile, Mapping) else {}
    review_cap = int(gates.get("max_splice_review_prediction_events", 200))
    max_pairs = int(gates.get("max_splice_peptide_hla_per_event", 20))
    predictable_events = [event_id for event_id in event_order if grouped.get(event_id)]
    pass_events = [event_id for event_id in predictable_events if status_by_event[event_id] == "PASS"]
    review_events = sorted(
        (event_id for event_id in predictable_events if status_by_event[event_id] == "REVIEW"),
        key=lambda event_id: tuple(-value for value in score_by_event[event_id]) + (event_id,),
    )[:review_cap]
    selected_events = set(pass_events + review_events)
    selected_splice: list[dict[str, str]] = []
    for event_id in event_order:
        if event_id not in selected_events:
            continue
        for row in grouped.get(event_id, [])[:max_pairs]:
            enriched = dict(row)
            enriched["splice_prefilter_status"] = status_by_event[event_id]
            enriched["splice_candidate_pool"] = (
                "FORMAL_SPLICE_CANDIDATE"
                if status_by_event[event_id] == "PASS"
                else "EXPLORATION_EVIDENCE_INCOMPLETE"
            )
            enriched["splice_formal_gate_pass"] = "yes" if status_by_event[event_id] == "PASS" else "no"
            selected_splice.append(enriched)
    write_tsv(raw_path, non_splice + selected_splice, PEPTIDE_FIELDS)

    funnel_rows: list[dict[str, str]] = [{
        "stage": "RAW_SPLICE_EVENTS_IN_UNIFIED_INPUT",
        "entered_events": str(len(event_order)),
        "assessed_events": str(len(event_order)),
        "passed_events": str(len(event_order)),
        "failed_events": "0",
        "unassessed_events": "0",
        "possible_remaining_range": str(len(event_order)),
        "criterion": "unique splice event IDs entering the unified peptide input; upstream aligner raw count may be larger",
    }]
    active = set(event_order)
    for stage, criterion, _ in stages:
        results = stage_results[stage]
        passed = {event_id for event_id in active if results[event_id][0] is True}
        failed = {event_id for event_id in active if results[event_id][0] is False}
        unassessed = active - passed - failed
        funnel_rows.append({
            "stage": stage,
            "entered_events": str(len(active)),
            "assessed_events": str(len(passed) + len(failed)),
            "passed_events": str(len(passed)),
            "failed_events": str(len(failed)),
            "unassessed_events": str(len(unassessed)),
            "possible_remaining_range": f"{len(passed)}-{len(passed) + len(unassessed)}",
            "criterion": criterion,
        })
        active = passed | unassessed
    funnel_rows.append({
        "stage": "FORMAL_GATE_COMPLETE",
        "entered_events": str(len(event_order)),
        "assessed_events": str(len(pass_events) + sum(status_by_event[event_id] == "REJECT" for event_id in event_order)),
        "passed_events": str(len(pass_events)),
        "failed_events": str(sum(status_by_event[event_id] == "REJECT" for event_id in event_order)),
        "unassessed_events": str(sum(status_by_event[event_id] == "REVIEW" for event_id in event_order)),
        "possible_remaining_range": str(len(pass_events)),
        "criterion": "all upstream biological gates passed; REVIEW events remain exploration-only and are not formal HLA-presentation passes",
    })
    funnel_rows.append({
        "stage": "SELECTED_FOR_PRESENTATION",
        "entered_events": str(len(active)),
        "assessed_events": str(len(pass_events) + len(review_events)),
        "passed_events": str(len(selected_events)),
        "failed_events": str(max(0, len(active) - len(selected_events))),
        "unassessed_events": "0",
        "possible_remaining_range": str(len(selected_events)),
        "criterion": f"all formal PASS plus top {review_cap} REVIEW events for exploratory prediction only; at most {max_pairs} peptide-HLA rows per event",
    })
    funnel_path = output_dir / "splice_prefilter_funnel.tsv"
    decision_path = output_dir / "splice_prefilter_decisions.tsv"
    write_tsv(funnel_path, funnel_rows, list(funnel_rows[0]))
    write_tsv(decision_path, decisions, [
        "event_id", "prefilter_status", "splice_candidate_pool", "splice_formal_gate_pass",
        "failed_stages", "unassessed_stages", "stage_details",
    ])
    return {
        "raw_splice_events": len(event_order),
        "pass_events": len(pass_events),
        "review_events_selected": len(review_events),
        "selected_events": len(selected_events),
        "selected_splice_peptide_hla_rows": len(selected_splice),
        "funnel": str(funnel_path),
        "decisions": str(decision_path),
        "archive": str(archive_path),
    }
