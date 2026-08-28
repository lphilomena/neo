from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from neoag.controlled_execution.io_utils import read_tsv, write_json, write_tsv


REQUIRED_ARTIFACTS = (
    "run_manifest",
    "consensus_events",
    "consensus_peptides",
    "weighted_baseline",
    "all_tool_results",
    "validation_plan",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _truthy(value: Any) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y"}


def _safe_int(value: Any) -> int:
    try:
        return int(float(str(value or "0").strip()))
    except (TypeError, ValueError):
        return 0


def _safe_float(value: Any) -> float | None:
    try:
        text = str(value or "").strip()
        return float(text) if text else None
    except (TypeError, ValueError):
        return None


def _event_track(row: dict[str, str]) -> str:
    value = str(
        row.get("event_type")
        or row.get("biological_event_track")
        or row.get("evidence_track")
        or ""
    ).strip().upper().replace("-", "_")
    aliases = {"MISSENSE": "SNV", "INSERTION": "INDEL", "DELETION": "INDEL"}
    return aliases.get(value, value or "UNASSESSED")


def _manifest_hash_records(value: Any) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    if isinstance(value, dict):
        path = value.get("path")
        digest = value.get("sha256")
        if isinstance(path, str) and isinstance(digest, str) and digest not in {"", "-", "NA", "not_computed_large_file"}:
            records.append((path, digest))
        for nested in value.values():
            records.extend(_manifest_hash_records(nested))
    elif isinstance(value, list):
        for nested in value:
            records.extend(_manifest_hash_records(nested))
    return records


def _add(checks: list[dict[str, str]], check: str, status: str, severity: str, detail: str, artifact: str = "") -> None:
    checks.append({"check": check, "status": status, "severity": severity, "artifact": artifact, "detail": detail})


def audit_review_inputs(artifacts: dict[str, str], outdir: str | Path) -> dict[str, Any]:
    od = Path(outdir)
    od.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, str]] = []

    missing = [key for key in REQUIRED_ARTIFACTS if not artifacts.get(key) or not Path(artifacts[key]).is_file()]
    for key in REQUIRED_ARTIFACTS:
        path = artifacts.get(key, "")
        _add(checks, f"required_artifact:{key}", "PASS" if key not in missing else "FAIL", "ERROR", path or "missing", path)
    if "consensus_events" in missing:
        overall = "NEEDS_RANKING"
    elif missing:
        overall = "BLOCKED"
    else:
        overall = "PASS"

    manifest: dict[str, Any] = {}
    if artifacts.get("run_manifest") and Path(artifacts["run_manifest"]).is_file():
        try:
            manifest = json.loads(Path(artifacts["run_manifest"]).read_text(encoding="utf-8"))
            algorithm = str(manifest.get("algorithm") or manifest.get("pipeline_version") or "")
            status = str(manifest.get("status") or "")
            consensus_ok = "pareto" in algorithm.lower() or "evidence" in algorithm.lower() or bool(manifest.get("outputs"))
            _add(checks, "evidence_consensus_completed", "PASS" if consensus_ok else "WARN", "WARNING", f"algorithm={algorithm or 'UNASSESSED'}; status={status or 'UNASSESSED'}", artifacts["run_manifest"])
        except (OSError, json.JSONDecodeError) as exc:
            _add(checks, "run_manifest_parse", "FAIL", "ERROR", str(exc), artifacts["run_manifest"])
            overall = "BLOCKED"

    auxiliary_manifests: list[dict[str, Any]] = []
    for key in ("all_tool_results_manifest", "provenance"):
        path = artifacts.get(key)
        if not path or not Path(path).is_file():
            continue
        try:
            payload = json.loads(Path(path).read_text(encoding="utf-8"))
            if isinstance(payload, dict):
                auxiliary_manifests.append(payload)
        except (OSError, json.JSONDecodeError) as exc:
            _add(checks, f"manifest_parse:{key}", "WARN", "WARNING", str(exc), path)
    hash_records = _manifest_hash_records(manifest)
    for auxiliary in auxiliary_manifests:
        hash_records.extend(_manifest_hash_records(auxiliary))
    if not hash_records:
        _add(checks, "input_reference_hash_consistency", "UNASSESSED", "WARNING", "No verifiable path+sha256 records in run_manifest")
    else:
        mismatches = []
        assessed = 0
        for raw_path, expected in hash_records:
            path = Path(raw_path)
            if not path.is_file():
                continue
            assessed += 1
            actual = _sha256(path)
            if actual.lower() != expected.lower():
                mismatches.append(f"{path}:{expected}!={actual}")
        _add(checks, "input_reference_hash_consistency", "FAIL" if mismatches else ("PASS" if assessed else "UNASSESSED"), "ERROR" if mismatches else "WARNING", ";".join(mismatches) if mismatches else f"verified={assessed}")
        if mismatches:
            overall = "BLOCKED"

    table_rows: dict[str, list[dict[str, str]]] = {}
    for key in ("consensus_events", "consensus_peptides", "all_tool_results"):
        path = artifacts.get(key)
        if path and Path(path).is_file():
            _, table_rows[key] = read_tsv(path)

    run_ids = set()
    for rows in table_rows.values():
        run_ids.update(str(row.get("run_id") or "").strip() for row in rows if str(row.get("run_id") or "").strip())
    manifest_run = str(manifest.get("run_id") or "").strip()
    if manifest_run:
        run_ids.add(manifest_run)
    _add(checks, "single_run_id", "PASS" if len(run_ids) == 1 else ("UNASSESSED" if not run_ids else "FAIL"), "ERROR" if len(run_ids) > 1 else "WARNING", ",".join(sorted(run_ids)) or "run_id not embedded in artifacts")
    if len(run_ids) > 1:
        overall = "BLOCKED"

    peptides = table_rows.get("consensus_peptides", [])
    events = table_rows.get("consensus_events", [])
    peptide_ids = {str(row.get("peptide_id") or "") for row in peptides}
    event_ids = {str(row.get("event_id") or "") for row in peptides}
    mapping_errors = []
    for event in events:
        if str(event.get("event_id") or "") and str(event.get("event_id")) not in event_ids:
            mapping_errors.append(f"event:{event.get('event_id')}")
        for index in (1, 2):
            peptide_id = str(event.get(f"representative_{index}_peptide_id") or "")
            if peptide_id and peptide_id not in peptide_ids:
                mapping_errors.append(f"peptide:{peptide_id}")
    _add(checks, "event_peptide_mapping", "FAIL" if mapping_errors else "PASS", "ERROR", ";".join(mapping_errors[:100]) if mapping_errors else f"events={len(events)}; peptides={len(peptides)}")
    if mapping_errors:
        overall = "BLOCKED"

    promoted_hard_fail = [
        row for row in peptides
        if (_truthy(row.get("hard_failure")) or str(row.get("hard_failure_codes") or "").strip())
        and str(row.get("evidence_grade") or "").upper() in {"R1", "R2"}
    ]
    promoted_hard_fail.extend(
        row for row in events
        if _safe_int(row.get("hard_failure_peptide_count")) > 0
        and str(row.get("best_evidence_grade") or "").upper() in {"R1", "R2"}
    )
    _add(checks, "hard_fail_not_promoted", "FAIL" if promoted_hard_fail else "PASS", "ERROR", f"violations={len(promoted_hard_fail)}")
    if promoted_hard_fail:
        overall = "BLOCKED"

    grade_counts: dict[str, int] = {}
    for row in peptides:
        grade = str(row.get("evidence_grade") or "UNASSESSED").strip().upper()
        grade_counts[grade] = grade_counts.get(grade, 0) + 1
    all_r4 = len(peptides) >= 20 and grade_counts.get("R4", 0) == len(peptides)
    _add(
        checks,
        "evidence_grade_distribution",
        "WARN" if all_r4 else "PASS",
        "WARNING",
        ";".join(f"{key}={grade_counts[key]}" for key in sorted(grade_counts))
        + ("; all candidates are R4: verify core evidence coverage before review" if all_r4 else ""),
    )

    net_assessed = 0
    mhc_assessed = 0
    for row in peptides:
        net_rank = _safe_float(row.get("netmhcpan_mt_rank_el") or row.get("netmhcpan_el_rank"))
        mhc_score = _safe_float(row.get("mhcflurry_presentation_score"))
        if net_rank is not None and net_rank < 99:
            net_assessed += 1
        if mhc_score is not None and mhc_score > 0:
            mhc_assessed += 1
    core_missing = len(peptides) >= 20 and (net_assessed == 0 or mhc_assessed == 0)
    _add(
        checks,
        "core_presentation_coverage",
        "FAIL" if core_missing else "PASS",
        "ERROR" if core_missing else "WARNING",
        f"peptides={len(peptides)}; netmhcpan_assessed={net_assessed}; mhcflurry_assessed={mhc_assessed}",
    )
    if core_missing:
        overall = "BLOCKED"

    stabpan_eligible: set[tuple[str, str]] = set()
    stabpan_covered: set[tuple[str, str]] = set()
    for row in peptides:
        peptide = str(row.get("peptide") or "").strip().upper()
        hla = str(row.get("hla_allele") or "").strip().upper()
        if 8 <= len(peptide) <= 11 and hla.startswith(("HLA-A", "HLA-B", "HLA-C")):
            key = (peptide, hla)
            stabpan_eligible.add(key)
            rank = _safe_float(row.get("netmhcstabpan_rank"))
            score = _safe_float(row.get("netmhcstabpan_score"))
            if (rank is not None and rank < 99) or score is not None:
                stabpan_covered.add(key)
    if not stabpan_eligible:
        stabpan_status = "NOT_APPLICABLE"
    elif not stabpan_covered:
        stabpan_status = "UNASSESSED"
    elif len(stabpan_covered) < len(stabpan_eligible):
        stabpan_status = "WARN"
    else:
        stabpan_status = "PASS"
    _add(
        checks,
        "netmhcstabpan_applicable_coverage",
        stabpan_status,
        "WARNING",
        f"covered_unique={len(stabpan_covered)}; eligible_unique={len(stabpan_eligible)}; "
        f"missing_unique={len(stabpan_eligible - stabpan_covered)}",
    )

    tracks = sorted({_event_track(row) for row in events if _event_track(row) != "UNASSESSED"})
    _add(
        checks,
        "applicable_event_tracks",
        "PASS" if tracks else "UNASSESSED",
        "WARNING",
        "tracks=" + (",".join(tracks) if tracks else "none"),
    )

    consensus_path = Path(artifacts.get("consensus_events") or "")
    search_roots = [consensus_path.parent, consensus_path.parent.parent, consensus_path.parent.parent.parent]
    funnel = next(
        (
            root / "parsed/splice_prefilter_funnel.tsv"
            for root in search_roots
            if (root / "parsed/splice_prefilter_funnel.tsv").is_file()
        ),
        None,
    )
    if "SPLICE" not in tracks:
        _add(checks, "splice_prefilter_funnel", "NOT_APPLICABLE", "INFO", "No Splice events in formal event ranking")
    elif funnel is None:
        _add(checks, "splice_prefilter_funnel", "WARN", "WARNING", "Splice is applicable but splice_prefilter_funnel.tsv was not found")
    else:
        _, funnel_rows = read_tsv(funnel)
        unassessed = sum(_safe_int(row.get("unassessed")) for row in funnel_rows)
        _add(
            checks,
            "splice_prefilter_funnel",
            "PASS" if unassessed == 0 else "WARN",
            "WARNING",
            f"stages={len(funnel_rows)}; unassessed={unassessed}; REVIEW candidates are not strict funnel PASS",
            str(funnel),
        )

    missing_as_negative = []
    for row in peptides:
        rna_state = str(row.get("rna_support_state") or "").upper()
        rna_status = str(row.get("rna_support_status") or "").upper()
        if "UNASSESSED" in rna_state and any(token in rna_status for token in ("NOT_DETECTED", "NEGATIVE")):
            missing_as_negative.append(str(row.get("peptide_id") or ""))
        missing_safety = str(row.get("safety_missing_layers") or row.get("event_safety_missing_layers") or "").strip()
        safety = str(row.get("safety_state") or row.get("safety_status") or "").upper()
        if missing_safety and safety in {"PASS", "SAFETY_PASS"}:
            missing_as_negative.append(str(row.get("peptide_id") or ""))
    _add(checks, "missing_evidence_not_negative", "FAIL" if missing_as_negative else "PASS", "ERROR", f"violations={len(set(missing_as_negative))}")
    if missing_as_negative:
        overall = "BLOCKED"

    if overall == "PASS" and any(row["status"] in {"WARN", "UNASSESSED"} for row in checks):
        overall = "PARTIAL"
    blocking = [row for row in checks if row["status"] == "FAIL"]
    write_tsv(od / "review_integrity_checks.tsv", checks)
    write_tsv(od / "review_blocking_issues.tsv", blocking)
    summary = {"status": overall, "required_artifacts": list(REQUIRED_ARTIFACTS), "missing_artifacts": missing, "checks": checks, "blocking_issues": blocking}
    write_json(od / "review_integrity.json", summary)
    return summary
