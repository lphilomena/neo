"""Auditable read-level quality gates for canonical splice junctions.

The formal splice consensus must not interpret a non-zero caller count as
validated RNA support.  This module converts caller rows into an explicit QC
record.  Missing required measurements remain missing and therefore produce
``INCOMPLETE`` rather than being silently imputed as passing.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from neoag.splice.identifiers import stable_id


@dataclass(frozen=True)
class JunctionReadQCThresholds:
    min_unique_split_reads: int = 3
    min_unique_fragment_starts: int = 2
    min_overhang: int = 10
    min_mapping_quality: float = 20.0
    max_multimapping_fraction: float = 0.20
    min_tumor_psi: float = 0.05
    policy: str = "complete"

    def __post_init__(self) -> None:
        if self.policy not in {"complete", "reads_only"}:
            raise ValueError("junction QC policy must be 'complete' or 'reads_only'")
        if self.min_unique_split_reads < 1:
            raise ValueError("min_unique_split_reads must be >= 1")
        if self.min_unique_fragment_starts < 1:
            raise ValueError("min_unique_fragment_starts must be >= 1")
        if self.min_overhang < 1:
            raise ValueError("min_overhang must be >= 1")
        if self.min_mapping_quality < 0:
            raise ValueError("min_mapping_quality must be >= 0")
        if not 0 <= self.max_multimapping_fraction <= 1:
            raise ValueError("max_multimapping_fraction must be in [0, 1]")
        if not 0 <= self.min_tumor_psi <= 1:
            raise ValueError("min_tumor_psi must be in [0, 1]")


def _text(row: Mapping[str, Any], *names: str) -> str:
    folded = {str(k).strip().casefold(): v for k, v in row.items()}
    for name in names:
        value = folded.get(name.casefold())
        if value is not None and str(value).strip() != "":
            return str(value).strip()
    return ""


def _number(value: Any) -> float | None:
    raw = str(value or "").strip().replace(",", "")
    if not raw or raw.upper() in {"NA", "N/A", "NONE", "UNASSESSED", "UNKNOWN", "."}:
        return None
    try:
        return float(raw)
    except (TypeError, ValueError):
        return None


def _format(value: float | None) -> str:
    if value is None:
        return ""
    return str(int(value)) if value.is_integer() else f"{value:.12g}"


def _bed12_overhang(row: Mapping[str, Any]) -> float | None:
    sizes = _text(row, "blockSizes", "block_sizes")
    if not sizes:
        return None
    parsed = [_number(token) for token in sizes.strip(",").split(",")]
    parsed = [value for value in parsed if value is not None]
    return min(parsed[0], parsed[-1]) if len(parsed) >= 2 else None


def build_junction_read_qc_row(
    *,
    row: Mapping[str, Any],
    junction_id: str,
    sample_id: str,
    source_tool: str,
    source_tool_version: str,
    source_assay_id: str,
    source_file: str,
    source_record_id: str,
    resolution_status: str,
    unique_split_reads: int,
    multi_split_reads: int,
    total_split_reads: int,
    thresholds: JunctionReadQCThresholds,
) -> dict[str, str]:
    unique_starts = _number(_text(
        row, "unique_fragment_starts", "unique_start_positions", "distinct_fragment_starts",
        "distinct_read_starts", "unique_molecule_starts", "umi_families",
    ))
    overhang = _number(_text(row, "max_overhang", "overhang", "anchor", "max_anchor"))
    if overhang is None:
        overhang = _bed12_overhang(row)
    mapq = _number(_text(
        row, "median_mapping_quality", "median_mapq", "min_mapping_quality", "min_mapq", "mapq",
    ))
    multimapping_fraction = _number(_text(
        row, "multimapping_fraction", "multi_mapping_fraction", "multimap_fraction",
    ))
    if multimapping_fraction is None and total_split_reads > 0:
        multimapping_fraction = multi_split_reads / total_split_reads
    tumor_psi = _number(_text(row, "tumor_psi", "junction_psi", "splice_psi", "psi"))
    caller_filter = _text(row, "caller_filter", "filter", "qc_status", "warnings", "warning")
    caller_folded = caller_filter.upper()
    caller_failed = bool(caller_folded) and any(
        token in caller_folded for token in ("FAIL", "REJECT", "BLACKLIST", "ARTIFACT")
    )

    failed: list[str] = []
    missing: list[str] = []
    if resolution_status != "RESOLVED_EXACT":
        failed.append("EXACT_STRANDED_JUNCTION_REQUIRED")
    if unique_split_reads < thresholds.min_unique_split_reads:
        failed.append("UNIQUE_SPLIT_READS_BELOW_THRESHOLD")
    if caller_failed:
        failed.append("CALLER_FILTER_FAILED")

    if thresholds.policy == "complete":
        measurements = (
            ("UNIQUE_FRAGMENT_STARTS", unique_starts, thresholds.min_unique_fragment_starts, "min"),
            ("MAX_OVERHANG", overhang, thresholds.min_overhang, "min"),
            ("MAPPING_QUALITY", mapq, thresholds.min_mapping_quality, "min"),
            ("MULTIMAPPING_FRACTION", multimapping_fraction, thresholds.max_multimapping_fraction, "max"),
            ("TUMOR_PSI", tumor_psi, thresholds.min_tumor_psi, "min"),
        )
        for name, value, threshold, direction in measurements:
            if value is None:
                missing.append(name)
            elif direction == "min" and value < threshold:
                failed.append(f"{name}_BELOW_THRESHOLD")
            elif direction == "max" and value > threshold:
                failed.append(f"{name}_ABOVE_THRESHOLD")

    if failed:
        status = "FAIL"
    elif missing:
        status = "INCOMPLETE"
    else:
        status = "PASS"

    return {
        "junction_read_qc_id": stable_id("JQC", junction_id, source_tool, source_assay_id, source_record_id),
        "junction_id": junction_id,
        "sample_id": sample_id,
        "source_tool": source_tool,
        "source_tool_version": source_tool_version,
        "source_assay_id": source_assay_id or "RNA_ASSAY_UNRESOLVED",
        "source_file": source_file,
        "source_record_id": source_record_id,
        "resolution_status": resolution_status,
        "unique_split_reads": str(unique_split_reads),
        "total_split_reads": str(total_split_reads),
        "unique_fragment_starts": _format(unique_starts),
        "max_overhang": _format(overhang),
        "mapping_quality": _format(mapq),
        "multimapping_fraction": _format(multimapping_fraction),
        "tumor_psi": _format(tumor_psi),
        "caller_filter_status": caller_filter or "NOT_REPORTED",
        "qc_policy": thresholds.policy,
        "min_unique_split_reads": str(thresholds.min_unique_split_reads),
        "min_unique_fragment_starts": str(thresholds.min_unique_fragment_starts),
        "min_overhang": str(thresholds.min_overhang),
        "min_mapping_quality": f"{thresholds.min_mapping_quality:g}",
        "max_multimapping_fraction": f"{thresholds.max_multimapping_fraction:g}",
        "min_tumor_psi": f"{thresholds.min_tumor_psi:g}",
        "qc_status": status,
        "failed_checks": ";".join(failed),
        "missing_checks": ";".join(missing),
        "qc_reason": (
            "All required read-level junction checks passed."
            if status == "PASS" else
            f"failed={','.join(failed) or 'NONE'}; missing={','.join(missing) or 'NONE'}"
        ),
    }


def passing_junctions(tables: Mapping[str, list[dict[str, str]]]) -> set[str]:
    """Return exact junctions with at least one explicit PASS QC record."""
    return {
        row.get("junction_id", "")
        for row in tables.get("junction_read_qc", [])
        if row.get("junction_id") and row.get("qc_status") == "PASS"
        and row.get("resolution_status") == "RESOLVED_EXACT"
    }
