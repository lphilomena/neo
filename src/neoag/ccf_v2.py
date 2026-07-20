"""CCF 2.1: event-type-aware clonality evidence layer.

This module upgrades CCF 2.0 while preserving the legacy columns consumed by
score. It remains a computational triage layer: SNV/InDel events with
read counts, purity and allele-specific copy number have the highest confidence;
SV/WES-SV and RNA-only events are explicitly marked approximate or unresolved.

P0/P1 additions in v0.4.3-style CCF 2.1:
  * input QC sidecar: ccf_input_qc.tsv
  * multiplicity candidates/confidence/ambiguity flags
  * clonality confidence and probability-style coarse summaries
  * external clonality adapter: PyClone-VI / PhylogicNDT-like TSVs
  * SVclone-like adapter for WGS SV CCF evidence
  * conflict sidecar between internal and external estimates
  * event-type policy: SNV/InDel, SV, WES-SV, RNA-only and CNV-derived events
    keep distinct ccf_method/ccf_confidence semantics.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Mapping

from .utils import clamp, first, read_tsv, to_float, write_tsv

CCF_FIELDS = [
    # Legacy / stable identifiers
    "event_id", "sample_id", "gene", "event_type", "mutation_source", "chrom", "pos",
    "tumor_alt_count", "tumor_depth", "tumor_vaf", "vaf_ci_low", "vaf_ci_high",
    "purity", "purity_source", "purity_confidence", "ploidy",
    "purity_consensus_status", "purity_range", "purity_n_tools", "purity_tool_values",
    "purity_recommendation_file",
    "total_cn", "major_cn", "minor_cn", "loh_status", "cnv_confidence",
    "multiplicity_candidates", "multiplicity_best", "multiplicity_confidence", "multiplicity_ambiguity",
    "ccf_best", "ccf_min", "ccf_max", "ccf_ci_low", "ccf_ci_high",
    "probability_ccf_gt_0_8", "probability_ccf_lt_0_3",
    # Backward-compatible columns consumed by older score code
    "ccf_estimate", "ccf_status", "clonality_status", "clonality_confidence", "clonality_multiplier",
    "ccf_method", "ccf_confidence", "ccf_warning", "mutation_multiplicity_assumption",
    "total_copy_number",
    # External evidence integration
    "external_clonality_tool", "external_cluster_id", "external_ccf", "external_ccf_low", "external_ccf_high",
    "external_clonality_status", "external_assignment_probability", "svclone_ccf", "svclone_copy_number",
    "ccf_resolution", "ccf_resolution_reason",
]

CCF_INPUT_QC_FIELDS = [
    "sample_id", "purity", "purity_source", "purity_confidence", "ploidy", "cnv_source",
    "purity_consensus_status", "purity_range", "purity_n_tools", "purity_tool_values",
    "purity_recommendation_file",
    "cnv_quality_status", "cnv_segment_count", "fraction_genome_cna", "fraction_genome_loh",
    "major_minor_cn_available", "subclonal_cna_detected", "ccf_ready_status", "ccf_ready_reason",
]

CCF_CONFLICT_FIELDS = [
    "event_id", "sample_id", "internal_ccf", "external_ccf", "internal_status", "external_status",
    "internal_confidence", "external_tool", "conflict_type", "conflict_reason", "resolution",
]

CCF_CLUSTER_FIELDS = [
    "cluster_id", "sample_id", "cluster_ccf", "cluster_ccf_low", "cluster_ccf_high", "n_events",
    "cluster_status", "external_tool", "cluster_assignment_probability",
]


def _path_exists(path: str | Path | None) -> bool:
    return bool(path and Path(path).exists())


def _clean_status(x: str) -> str:
    return str(x or "").strip().lower().replace(" ", "_")


def load_purity(path: str | Path | None) -> tuple[float | None, float | None]:
    """Backward-compatible purity loader."""
    info = load_purity_info(path)
    return info["purity"], info["ploidy"]


def _load_purity_recommendation(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    status = str(payload.get("status") or "NO_PURITY").upper()
    purity = to_float(payload.get("recommended_purity", ""), -1.0)
    tool_values = payload.get("tool_values") if isinstance(payload.get("tool_values"), dict) else {}
    n_tools = int(to_float(payload.get("n_tools", len(tool_values)), len(tool_values)))
    if status == "CONCORDANT":
        confidence = "high"
    elif status in {"SINGLE_TOOL", "MODERATE_DISCORDANCE"}:
        confidence = "medium"
    else:
        confidence = "low"
    source = "multi_tool_median" if n_tools > 1 else next(iter(tool_values), "purity_recommendation")
    return {
        "purity": purity if purity > 0 else None,
        "ploidy": None,
        "purity_source": source,
        "purity_confidence": confidence,
        "purity_consensus_status": status,
        "purity_range": str(payload.get("range") or ""),
        "purity_n_tools": n_tools,
        "purity_tool_values": json.dumps(tool_values, sort_keys=True, separators=(",", ":")),
        "purity_recommendation_file": str(path),
        "reason": f"purity_consensus_{status.lower()}",
    }


def _with_purity_consensus_defaults(info: dict[str, Any]) -> dict[str, Any]:
    return {
        **info,
        "purity_consensus_status": info.get("purity_consensus_status", "NOT_ASSESSED"),
        "purity_range": info.get("purity_range", ""),
        "purity_n_tools": info.get("purity_n_tools", ""),
        "purity_tool_values": info.get("purity_tool_values", ""),
        "purity_recommendation_file": info.get("purity_recommendation_file", ""),
    }


def load_purity_info(path: str | Path | None) -> dict[str, Any]:
    if not path or not Path(path).exists():
        return _with_purity_consensus_defaults({
            "purity": None,
            "ploidy": None,
            "purity_source": "not_provided",
            "purity_confidence": "low",
            "reason": "missing_purity",
        })
    purity_path = Path(path)
    if purity_path.suffix.lower() == ".json":
        return _load_purity_recommendation(purity_path)
    rows = read_tsv(purity_path)
    if not rows:
        return _with_purity_consensus_defaults({
            "purity": None,
            "ploidy": None,
            "purity_source": "provided_empty",
            "purity_confidence": "low",
            "reason": "empty_purity_file",
        })
    r = rows[0]
    purity = to_float(first(r, ["purity", "tumor_purity", "cellularity", "purity_estimate"], ""), -1.0)
    ploidy = to_float(first(r, ["ploidy", "tumor_ploidy"], ""), -1.0)
    source = first(r, ["source", "tool", "purity_source", "caller"], "provided")
    conf = first(r, ["confidence", "purity_confidence", "quality"], "")
    if not conf:
        if purity <= 0:
            conf = "low"
        elif purity < 0.20:
            conf = "low"
        elif purity < 0.35:
            conf = "medium"
        else:
            conf = "high"
    return _with_purity_consensus_defaults({
        "purity": purity if purity > 0 else None,
        "ploidy": ploidy if ploidy > 0 else None,
        "purity_source": source,
        "purity_confidence": conf,
        "reason": "purity_loaded" if purity > 0 else "invalid_purity",
    })


def _row_chrom(r: Mapping[str, Any]) -> str:
    return first(r, ["chrom", "chromosome", "chr", "Chromosome"], "")


def _row_start(r: Mapping[str, Any]) -> int:
    return int(to_float(first(r, ["start", "Start", "loc.start", "seg_start", "chromStart"], "0"), 0))


def _row_end(r: Mapping[str, Any]) -> int:
    return int(to_float(first(r, ["end", "End", "loc.end", "seg_end", "chromEnd"], "0"), 0))


def _copy_status_from_minor(minor_cn: float | None, total_cn: float | None) -> str:
    if total_cn is None:
        return "unknown"
    if total_cn <= 0.2:
        return "homozygous_deletion"
    if minor_cn is not None and minor_cn <= 0.2 and total_cn >= 1.0:
        return "loh"
    if total_cn < 1.5:
        return "copy_loss"
    if total_cn > 2.5:
        return "copy_gain"
    return "neutral"


def _major_minor_available(row: Mapping[str, Any]) -> bool:
    major_raw = first(row, ["major_cn", "major_copy_number", "major", "tcn.em", "cf.em"], "")
    minor_raw = first(row, ["minor_cn", "minor_copy_number", "minor", "lcn.em", "lcn"], "")
    return bool(major_raw and minor_raw)


def find_cn(chrom: str, pos: str | int | float, cnv_rows: list[dict[str, str]]) -> dict[str, Any]:
    p = int(to_float(pos, -1))
    best = {
        "total_cn": 2.0,
        "major_cn": "",
        "minor_cn": "",
        "loh_status": "unknown",
        "matched": False,
        "cnv_confidence": "low",
        "cnv_source": "default_diploid",
    }
    if p < 0:
        return best
    for r in cnv_rows:
        if _row_chrom(r) != chrom:
            continue
        start, end = _row_start(r), _row_end(r)
        if start <= p <= end:
            total = to_float(first(r, ["total_cn", "copy_number", "total_copy_number", "tcn", "cn", "CNt"], "2"), 2.0)
            major_raw = first(r, ["major_cn", "major_copy_number", "major", "tcn.em", "cf.em"], "")
            minor_raw = first(r, ["minor_cn", "minor_copy_number", "minor", "lcn.em", "lcn"], "")
            major = to_float(major_raw, math.nan) if major_raw else math.nan
            minor = to_float(minor_raw, math.nan) if minor_raw else math.nan
            loh = first(r, ["loh_status", "LOH", "loh", "cnloh", "status", "call"], "")
            if not loh:
                loh = _copy_status_from_minor(None if math.isnan(minor) else minor, total)
            conf = first(r, ["cnv_confidence", "confidence", "seg_confidence", "quality"], "")
            if not conf:
                conf = "high" if not math.isnan(major) and not math.isnan(minor) else "medium"
            return {
                "total_cn": max(0.1, total),
                "major_cn": "" if math.isnan(major) else f"{major:.4f}",
                "minor_cn": "" if math.isnan(minor) else f"{minor:.4f}",
                "loh_status": loh,
                "matched": True,
                "cnv_confidence": conf,
                "cnv_source": first(r, ["source", "tool", "caller"], "segment"),
            }
    return best


def estimate_ccf(vaf: float, purity: float, total_cn: float, multiplicity: float) -> float:
    if purity <= 0 or multiplicity <= 0:
        return 0.0
    denom = purity * multiplicity
    numerator = vaf * (purity * total_cn + 2.0 * (1.0 - purity))
    return clamp(numerator / denom, 0.0, 1.5)


def enumerate_ccf(vaf: float, purity: float, total_cn: float) -> dict[str, Any]:
    max_m = max(1, int(round(max(total_cn, 1.0))))
    vals: list[tuple[int, float]] = []
    for m in range(1, max_m + 1):
        c = estimate_ccf(vaf, purity, total_cn, float(m))
        if 0 <= c <= 1.5:
            vals.append((m, c))
    if not vals:
        c = estimate_ccf(vaf, purity, total_cn, 1.0)
        return {
            "best_m": 1,
            "best": c,
            "min": c,
            "max": c,
            "candidates": "1",
            "multiplicity_confidence": "low",
            "multiplicity_ambiguity": "no_valid_multiplicity",
        }
    valid_under = [(m, c) for m, c in vals if c <= 1.0]
    if valid_under:
        # Conservative default: smallest multiplicity that remains biologically valid.
        best_m, best_c = valid_under[0]
    else:
        best_m, best_c = min(vals, key=lambda x: abs(x[1] - 1.0))

    near_valid = [(m, c) for m, c in vals if 0.70 <= c <= 1.10]
    if len(vals) == 1:
        mult_conf = "high"
        ambiguity = "none"
    elif len(near_valid) > 1:
        mult_conf = "ambiguous"
        ambiguity = "multiple_plausible_multiplicities"
    elif len(valid_under) == 1:
        mult_conf = "medium"
        ambiguity = "single_biologically_valid_multiplicity"
    else:
        mult_conf = "low"
        ambiguity = "weak_multiplicity_resolution"

    return {
        "best_m": best_m,
        "best": best_c,
        "min": min(c for _, c in vals),
        "max": max(c for _, c in vals),
        "candidates": ";".join(str(m) for m, _ in vals),
        "multiplicity_confidence": mult_conf,
        "multiplicity_ambiguity": ambiguity,
    }


def _vaf_from_event(e: Mapping[str, Any]) -> tuple[float, int | None, int | None]:
    vaf = to_float(first(e, ["tumor_vaf", "vaf", "sv_vaf_like"], ""), -1.0)
    alt = int(to_float(first(e, ["tumor_alt_count", "alt_count", "tumor_alt_support", "sv_alt_support"], ""), -1))
    depth = int(to_float(first(e, ["tumor_depth", "depth", "local_depth", "tumor_local_depth"], ""), -1))
    if vaf < 0 and alt >= 0 and depth > 0:
        vaf = alt / depth
    return max(0.0, vaf if vaf >= 0 else 0.0), (alt if alt >= 0 else None), (depth if depth > 0 else None)


def _ci_for_vaf(vaf: float, depth: int | None) -> tuple[float | None, float | None]:
    if not depth or depth <= 0:
        return None, None
    # Wilson interval would be better; normal approximation is transparent and stable for this triage layer.
    se = math.sqrt(max(vaf * (1 - vaf), 0.0) / max(depth, 1))
    return max(0.0, vaf - 1.96 * se), min(1.0, vaf + 1.96 * se)


def _event_method(e: Mapping[str, Any]) -> str:
    source = " ".join([
        str(e.get("mutation_source") or ""),
        str(e.get("event_type") or ""),
        str(e.get("peptide_consequence") or ""),
        str(e.get("source") or ""),
        str(e.get("evidence_scope") or ""),
    ]).lower()
    if "rna_only" in source or ("rna" in source and not e.get("tumor_vaf") and not e.get("tumor_alt_count")):
        return "RNA_ONLY_UNRESOLVED"
    if "wes" in source and "sv" in source:
        return "WES_SV_CAPTURE_LIMITED_APPROX"
    if any(x in source for x in ["sv", "bnd", "structural"]):
        return "SV_BREAKPOINT_APPROX"
    if any(x in source for x in ["fusion", "junction", "splice"]):
        return "JUNCTION_APPROX"
    if any(x in source for x in ["snv", "indel", "missense", "frameshift"]):
        return "SNV_INDEL_COPY_NUMBER_AWARE"
    if any(x in source for x in ["cnv", "copy", "exon"]):
        return "COPY_NUMBER_AWARE_APPROX"
    return "COPY_NUMBER_AWARE_APPROX"


def _probability_from_interval(ccf_best: float | None, lo: float | None, hi: float | None, threshold: float, direction: str) -> str:
    if ccf_best is None:
        return ""
    if lo is None or hi is None or hi <= lo:
        if direction == "gt":
            return "1.0000" if ccf_best > threshold else "0.0000"
        return "1.0000" if ccf_best < threshold else "0.0000"
    if direction == "gt":
        if lo >= threshold:
            p = 1.0
        elif hi <= threshold:
            p = 0.0
        else:
            p = (hi - threshold) / (hi - lo)
    else:
        if hi <= threshold:
            p = 1.0
        elif lo >= threshold:
            p = 0.0
        else:
            p = (threshold - lo) / (hi - lo)
    return f"{clamp(p):.4f}"


def clonality_status(ccf: float | None, confidence: str, profile: Mapping[str, Any]) -> str:
    if ccf is None:
        return "unresolved"
    cfg = profile.get("ccf_lite", {}) if profile else {}
    clonal = float(cfg.get("clonal_threshold", 0.80))
    sub = float(cfg.get("subclonal_threshold", 0.30))
    if ccf >= clonal:
        return "clonal_like"
    if ccf >= sub:
        return "subclonal_like"
    return "low_frequency_subclonal"


def clonality_multiplier(status: str, profile: Mapping[str, Any]) -> float:
    cfg = profile.get("ccf_lite", {}) if profile else {}
    if status == "clonal_like":
        return 1.0
    if status == "subclonal_like":
        return float(cfg.get("subclonal_multiplier", 0.75))
    if status == "low_frequency_subclonal":
        return float(cfg.get("low_frequency_multiplier", 0.45))
    return float(cfg.get("missing_ccf_multiplier", 0.65))


def _clonality_confidence(
    *,
    ccf_confidence: str,
    multiplicity_confidence: str,
    ccf_ci_low: float | None,
    ccf_ci_high: float | None,
    method: str,
) -> str:
    if method == "RNA_ONLY_UNRESOLVED":
        return "unresolved"
    if ccf_confidence == "low" or multiplicity_confidence in {"low", "ambiguous"}:
        return "low"
    if ccf_ci_low is not None and ccf_ci_high is not None and (ccf_ci_high - ccf_ci_low) > 0.65:
        return "low"
    if ccf_confidence == "high" and multiplicity_confidence == "high":
        return "high"
    return "medium"


def _cnv_qc(cnv_rows: list[dict[str, str]]) -> dict[str, str]:
    if not cnv_rows:
        return {
            "cnv_source": "not_provided",
            "cnv_quality_status": "missing",
            "cnv_segment_count": "0",
            "fraction_genome_cna": "",
            "fraction_genome_loh": "",
            "major_minor_cn_available": "no",
            "subclonal_cna_detected": "unknown",
        }
    n = len(cnv_rows)
    with_mm = sum(1 for r in cnv_rows if _major_minor_available(r))
    cna = 0
    loh = 0
    subclonal = False
    total_len = 0
    cna_len = 0
    loh_len = 0
    for r in cnv_rows:
        start, end = _row_start(r), _row_end(r)
        length = max(0, end - start + 1)
        total_len += length
        total_cn = to_float(first(r, ["total_cn", "copy_number", "total_copy_number", "tcn", "cn", "CNt"], "2"), 2.0)
        minor_raw = first(r, ["minor_cn", "minor_copy_number", "minor", "lcn.em", "lcn"], "")
        minor = to_float(minor_raw, math.nan) if minor_raw else math.nan
        loh_status = _clean_status(first(r, ["loh_status", "LOH", "loh", "cnloh", "status", "call"], ""))
        if abs(total_cn - 2.0) > 0.35:
            cna += 1; cna_len += length
        if loh_status in {"loh", "loss", "copy_neutral_loh", "cnloh"} or (not math.isnan(minor) and minor <= 0.2 and total_cn >= 1.0):
            loh += 1; loh_len += length
        if _clean_status(first(r, ["subclonal", "subclonal_status", "cellular_fraction", "cf"], "")) in {"yes", "true", "subclonal"}:
            subclonal = True
    qual = "allele_specific" if with_mm / max(n, 1) >= 0.5 else "total_cn_only"
    return {
        "cnv_source": first(cnv_rows[0], ["source", "tool", "caller"], "provided"),
        "cnv_quality_status": qual,
        "cnv_segment_count": str(n),
        "fraction_genome_cna": "" if total_len <= 0 else f"{cna_len / total_len:.4f}",
        "fraction_genome_loh": "" if total_len <= 0 else f"{loh_len / total_len:.4f}",
        "major_minor_cn_available": "yes" if with_mm > 0 else "no",
        "subclonal_cna_detected": "yes" if subclonal else "no",
    }


def _ccf_ready_status(purity_info: Mapping[str, Any], cnv_rows: list[dict[str, str]]) -> tuple[str, str]:
    reasons = []
    purity = purity_info.get("purity")
    if purity is None:
        reasons.append("missing_purity")
    elif purity < 0.20:
        reasons.append("low_purity")
    consensus_status = str(purity_info.get("purity_consensus_status", "")).upper()
    if consensus_status == "STRONG_DISCORDANCE":
        reasons.append("purity_strong_discordance")
    elif consensus_status == "MODERATE_DISCORDANCE":
        reasons.append("purity_moderate_discordance")
    if not cnv_rows:
        reasons.append("missing_cnv_segments")
    elif not any(_major_minor_available(r) for r in cnv_rows):
        reasons.append("major_minor_cn_missing")
    if not reasons:
        return "ready", "purity_cnv_available"
    if set(reasons).issubset({"major_minor_cn_missing", "purity_moderate_discordance"}):
        return "approximate", ";".join(reasons)
    return "limited", ";".join(reasons)


def build_input_qc(sample_id: str, purity_info: Mapping[str, Any], cnv_rows: list[dict[str, str]]) -> dict[str, str]:
    q = _cnv_qc(cnv_rows)
    ready, reason = _ccf_ready_status(purity_info, cnv_rows)
    return {
        "sample_id": sample_id,
        "purity": "" if purity_info.get("purity") is None else f"{float(purity_info['purity']):.4f}",
        "purity_source": str(purity_info.get("purity_source", "")),
        "purity_confidence": str(purity_info.get("purity_confidence", "low")),
        "ploidy": "" if purity_info.get("ploidy") is None else f"{float(purity_info['ploidy']):.4f}",
        "purity_consensus_status": str(purity_info.get("purity_consensus_status", "")),
        "purity_range": str(purity_info.get("purity_range", "")),
        "purity_n_tools": str(purity_info.get("purity_n_tools", "")),
        "purity_tool_values": str(purity_info.get("purity_tool_values", "")),
        "purity_recommendation_file": str(purity_info.get("purity_recommendation_file", "")),
        **q,
        "ccf_ready_status": ready,
        "ccf_ready_reason": reason,
    }


def _load_external_clonality(path: str | Path | None) -> dict[str, dict[str, str]]:
    if not _path_exists(path):
        return {}
    out: dict[str, dict[str, str]] = {}
    for r in read_tsv(path):
        eid = first(r, ["event_id", "mutation_id", "variant_id", "neoag_event_id"], "")
        if not eid:
            continue
        ccf = first(r, ["cellular_prevalence", "ccf", "ccf_best", "cluster_ccf", "cellular_prevalence_mean"], "")
        out[eid] = {
            "event_id": eid,
            "tool": first(r, ["tool", "source", "external_tool"], "external_clonality"),
            "cluster_id": first(r, ["cluster_id", "clone_id", "cluster"], ""),
            "external_ccf": ccf,
            "external_ccf_low": first(r, ["cellular_prevalence_low", "ccf_low", "ccf_ci_low", "cluster_ccf_low"], ""),
            "external_ccf_high": first(r, ["cellular_prevalence_high", "ccf_high", "ccf_ci_high", "cluster_ccf_high"], ""),
            "external_status": first(r, ["clonal_status", "clonality_status", "status", "cluster_status"], ""),
            "assignment_probability": first(r, ["cluster_assignment_probability", "assignment_probability", "probability", "posterior"], ""),
        }
    return out


def _load_svclone(path: str | Path | None) -> dict[str, dict[str, str]]:
    if not _path_exists(path):
        return {}
    out: dict[str, dict[str, str]] = {}
    for r in read_tsv(path):
        eid = first(r, ["event_id", "sv_event_id", "sv_id", "variant_id"], "")
        if not eid:
            continue
        out[eid] = {
            "event_id": eid,
            "tool": first(r, ["tool", "source"], "SVclone"),
            "svclone_ccf": first(r, ["sv_ccf", "ccf", "cellular_prevalence", "sv_ccf_best"], ""),
            "svclone_ccf_low": first(r, ["sv_ccf_low", "ccf_low", "ccf_ci_low"], ""),
            "svclone_ccf_high": first(r, ["sv_ccf_high", "ccf_high", "ccf_ci_high"], ""),
            "svclone_copy_number": first(r, ["sv_copy_number", "copy_number", "sv_cn"], ""),
            "svclone_confidence": first(r, ["confidence", "svclone_confidence", "quality"], ""),
        }
    return out


def _external_status_from_ccf(x: float) -> str:
    if x >= 0.85:
        return "clonal_like"
    if x >= 0.25:
        return "subclonal_like"
    if x > 0:
        return "low_frequency_subclonal"
    return "unresolved"


def _resolve_external(row: dict[str, str], external: Mapping[str, Mapping[str, str]], svclone: Mapping[str, Mapping[str, str]]) -> dict[str, str]:
    eid = row["event_id"]
    ext = external.get(eid, {})
    sv = svclone.get(eid, {})
    if sv and row.get("ccf_method") in {"SV_BREAKPOINT_APPROX", "WES_SV_CAPTURE_LIMITED_APPROX"}:
        ccf = sv.get("svclone_ccf", "")
        return {
            "external_clonality_tool": sv.get("tool", "SVclone"),
            "external_cluster_id": "",
            "external_ccf": ccf,
            "external_ccf_low": sv.get("svclone_ccf_low", ""),
            "external_ccf_high": sv.get("svclone_ccf_high", ""),
            "external_clonality_status": _external_status_from_ccf(to_float(ccf, 0.0)) if ccf else "",
            "external_assignment_probability": "",
            "svclone_ccf": ccf,
            "svclone_copy_number": sv.get("svclone_copy_number", ""),
        }
    if ext:
        ccf = ext.get("external_ccf", "")
        return {
            "external_clonality_tool": ext.get("tool", "external_clonality"),
            "external_cluster_id": ext.get("cluster_id", ""),
            "external_ccf": ccf,
            "external_ccf_low": ext.get("external_ccf_low", ""),
            "external_ccf_high": ext.get("external_ccf_high", ""),
            "external_clonality_status": ext.get("external_status", "") or (_external_status_from_ccf(to_float(ccf, 0.0)) if ccf else ""),
            "external_assignment_probability": ext.get("assignment_probability", ""),
            "svclone_ccf": "",
            "svclone_copy_number": "",
        }
    return {
        "external_clonality_tool": "",
        "external_cluster_id": "",
        "external_ccf": "",
        "external_ccf_low": "",
        "external_ccf_high": "",
        "external_clonality_status": "",
        "external_assignment_probability": "",
        "svclone_ccf": "",
        "svclone_copy_number": "",
    }


def _resolve_ccf_resolution(row: dict[str, str]) -> tuple[str, str]:
    if row.get("external_ccf"):
        if row.get("external_clonality_tool") == "SVclone":
            return "external_svclone_preferred", "svclone_available_for_sv_event"
        internal = to_float(row.get("ccf_best", ""), -1.0)
        external = to_float(row.get("external_ccf", ""), -1.0)
        if internal >= 0 and external >= 0 and abs(internal - external) > 0.35:
            return "external_conflict_review", "internal_external_ccf_difference_gt_0.35"
        return "external_supported", "external_clonality_available"
    if row.get("ccf_method") == "RNA_ONLY_UNRESOLVED":
        return "unresolved", "rna_only_no_dna_ccf"
    return "internal", "internal_copy_number_aware_estimate"


def _conflict_row(row: Mapping[str, str]) -> dict[str, str] | None:
    if not row.get("external_ccf"):
        return None
    internal = to_float(row.get("ccf_best", ""), -1.0)
    external = to_float(row.get("external_ccf", ""), -1.0)
    if internal < 0 or external < 0:
        return None
    internal_status = row.get("clonality_status", "")
    external_status = row.get("external_clonality_status", "")
    diff = abs(internal - external)
    if diff <= 0.35 and (not external_status or external_status == internal_status):
        return None
    return {
        "event_id": row.get("event_id", ""),
        "sample_id": row.get("sample_id", ""),
        "internal_ccf": row.get("ccf_best", ""),
        "external_ccf": row.get("external_ccf", ""),
        "internal_status": internal_status,
        "external_status": external_status,
        "internal_confidence": row.get("ccf_confidence", ""),
        "external_tool": row.get("external_clonality_tool", ""),
        "conflict_type": "ccf_delta" if diff > 0.35 else "status_discordance",
        "conflict_reason": f"internal_external_delta={diff:.4f}",
        "resolution": row.get("ccf_resolution", "review"),
    }


def _cluster_rows(external: Mapping[str, Mapping[str, str]]) -> list[dict[str, str]]:
    clusters: dict[tuple[str, str], dict[str, Any]] = {}
    for ext in external.values():
        cid = ext.get("cluster_id") or "unclustered"
        tool = ext.get("tool") or "external_clonality"
        key = (tool, cid)
        ccf = to_float(ext.get("external_ccf", ""), -1.0)
        if key not in clusters:
            clusters[key] = {"ccfs": [], "n": 0, "prob": [], "lo": [], "hi": []}
        clusters[key]["n"] += 1
        if ccf >= 0:
            clusters[key]["ccfs"].append(ccf)
        prob = to_float(ext.get("assignment_probability", ""), -1.0)
        if prob >= 0:
            clusters[key]["prob"].append(prob)
        lo = to_float(ext.get("external_ccf_low", ""), -1.0)
        hi = to_float(ext.get("external_ccf_high", ""), -1.0)
        if lo >= 0:
            clusters[key]["lo"].append(lo)
        if hi >= 0:
            clusters[key]["hi"].append(hi)
    rows = []
    for (tool, cid), vals in sorted(clusters.items()):
        ccfs = vals["ccfs"] or [0.0]
        ccf = sum(ccfs) / len(ccfs)
        rows.append({
            "cluster_id": cid,
            "sample_id": "",
            "cluster_ccf": f"{ccf:.4f}",
            "cluster_ccf_low": f"{min(vals['lo']):.4f}" if vals["lo"] else "",
            "cluster_ccf_high": f"{max(vals['hi']):.4f}" if vals["hi"] else "",
            "n_events": str(vals["n"]),
            "cluster_status": _external_status_from_ccf(ccf),
            "external_tool": tool,
            "cluster_assignment_probability": f"{sum(vals['prob']) / len(vals['prob']):.4f}" if vals["prob"] else "",
        })
    return rows


def _base_confidence(
    *,
    method: str,
    cn: Mapping[str, Any],
    purity_info: Mapping[str, Any],
    depth: int | None,
    multiplicity_confidence: str,
) -> tuple[str, list[str]]:
    warning: list[str] = []
    conf = "medium"
    purity = purity_info.get("purity")
    consensus_status = str(purity_info.get("purity_consensus_status", "")).upper()
    if consensus_status == "STRONG_DISCORDANCE":
        return "low", ["purity_strong_discordance"]
    if consensus_status == "MODERATE_DISCORDANCE":
        warning.append("purity_moderate_discordance")
    if purity is None:
        return "low", ["missing_purity"]
    if purity < 0.20:
        conf = "low"; warning.append("low_purity")
    elif purity_info.get("purity_confidence") == "high" and cn.get("matched") and cn.get("cnv_confidence") == "high":
        conf = "high"
    if consensus_status == "MODERATE_DISCORDANCE" and conf == "high":
        conf = "medium"
    if method == "RNA_ONLY_UNRESOLVED":
        return "low", ["rna_only_no_dna_ccf"]
    if method == "WES_SV_CAPTURE_LIMITED_APPROX":
        conf = "low"; warning.append("wes_sv_capture_limited")
    elif method in {"SV_BREAKPOINT_APPROX", "JUNCTION_APPROX"}:
        conf = "low" if conf != "high" else "medium"; warning.append("breakpoint_support_approx")
    if not cn.get("matched"):
        conf = "low" if conf == "medium" else conf
        warning.append("default_copy_number_2")
    elif cn.get("major_cn", "") == "" or cn.get("minor_cn", "") == "":
        conf = "medium" if conf == "high" else conf
        warning.append("major_minor_cn_missing")
    if depth is not None and depth < 20:
        conf = "low"; warning.append("low_depth")
    if multiplicity_confidence in {"low", "ambiguous"}:
        conf = "low" if conf == "medium" else ("medium" if conf == "high" else conf)
        warning.append(f"multiplicity_{multiplicity_confidence}")
    return conf, warning


def build_ccf_2(
    events_tsv: str | Path,
    purity_tsv: str | Path | None,
    cnv_tsv: str | Path | None,
    profile: Mapping[str, Any],
    out: str | Path,
    *,
    external_clonality_tsv: str | Path | None = None,
    svclone_tsv: str | Path | None = None,
    sidecar_dir: str | Path | None = None,
    input_qc_out: str | Path | None = None,
    conflicts_out: str | Path | None = None,
    clusters_out: str | Path | None = None,
) -> list[dict[str, str]]:
    events = read_tsv(events_tsv)
    sample_id = events[0].get("sample_id", "") if events else ""
    purity_info = load_purity_info(purity_tsv)
    purity = purity_info.get("purity")
    ploidy = purity_info.get("ploidy")
    cnv = read_tsv(cnv_tsv) if cnv_tsv and Path(cnv_tsv).exists() else []
    external = _load_external_clonality(external_clonality_tsv)
    svclone = _load_svclone(svclone_tsv)

    out_path = Path(out)
    sd = Path(sidecar_dir) if sidecar_dir else out_path.parent
    input_qc_path = Path(input_qc_out) if input_qc_out else sd / "ccf_input_qc.tsv"
    conflicts_path = Path(conflicts_out) if conflicts_out else sd / "ccf_conflicts.tsv"
    clusters_path = Path(clusters_out) if clusters_out else sd / "ccf_cluster.tsv"

    qc_row = build_input_qc(sample_id, purity_info, cnv)
    rows: list[dict[str, str]] = []
    for e in events:
        event_id = e.get("event_id", "")
        chrom = e.get("chrom", "")
        pos = e.get("pos", "")
        vaf, alt, depth = _vaf_from_event(e)
        vaf_lo, vaf_hi = _ci_for_vaf(vaf, depth)
        cn = find_cn(chrom, pos, cnv)
        total_cn = float(cn["total_cn"])
        method = _event_method(e)
        ccf_best = ccf_min = ccf_max = None
        best_m = ""
        ci_low = ci_high = None
        mult_candidates = ""
        mult_conf = "unresolved"
        mult_ambiguity = "unresolved"
        status = "unresolved"
        confidence = "low"
        warning: list[str] = []
        if purity is None or purity <= 0:
            status = "vaf_only_unresolved"
            warning.append("missing_purity")
        elif method == "RNA_ONLY_UNRESOLVED":
            status = "unresolved"
            warning.append("rna_only_no_dna_ccf")
        else:
            enum = enumerate_ccf(vaf, float(purity), total_cn)
            ccf_best, ccf_min, ccf_max = enum["best"], enum["min"], enum["max"]
            best_m = str(enum["best_m"])
            mult_candidates = str(enum["candidates"])
            mult_conf = str(enum["multiplicity_confidence"])
            mult_ambiguity = str(enum["multiplicity_ambiguity"])
            if vaf_lo is not None and purity and float(purity) > 0:
                ci_low = estimate_ccf(vaf_lo, float(purity), total_cn, float(enum["best_m"]))
                ci_high = estimate_ccf(vaf_hi, float(purity), total_cn, float(enum["best_m"]))
            confidence, warning = _base_confidence(
                method=method,
                cn=cn,
                purity_info=purity_info,
                depth=depth,
                multiplicity_confidence=mult_conf,
            )
            status = clonality_status(ccf_best, confidence, profile)
        clon_conf = _clonality_confidence(
            ccf_confidence=confidence,
            multiplicity_confidence=mult_conf,
            ccf_ci_low=ci_low,
            ccf_ci_high=ci_high,
            method=method,
        )
        mult = clonality_multiplier(status, profile)
        row: dict[str, str] = {
            "event_id": event_id,
            "sample_id": e.get("sample_id", ""),
            "gene": e.get("gene", ""),
            "event_type": e.get("event_type", ""),
            "mutation_source": e.get("mutation_source", ""),
            "chrom": chrom,
            "pos": str(pos),
            "tumor_alt_count": "" if alt is None else str(alt),
            "tumor_depth": "" if depth is None else str(depth),
            "tumor_vaf": f"{vaf:.4f}",
            "vaf_ci_low": "" if vaf_lo is None else f"{vaf_lo:.4f}",
            "vaf_ci_high": "" if vaf_hi is None else f"{vaf_hi:.4f}",
            "purity": "" if purity is None else f"{float(purity):.4f}",
            "purity_source": str(purity_info.get("purity_source", "")),
            "purity_confidence": str(purity_info.get("purity_confidence", "low")),
            "ploidy": "" if ploidy is None else f"{float(ploidy):.4f}",
            "purity_consensus_status": str(purity_info.get("purity_consensus_status", "")),
            "purity_range": str(purity_info.get("purity_range", "")),
            "purity_n_tools": str(purity_info.get("purity_n_tools", "")),
            "purity_tool_values": str(purity_info.get("purity_tool_values", "")),
            "purity_recommendation_file": str(purity_info.get("purity_recommendation_file", "")),
            "total_cn": f"{total_cn:.4f}",
            "major_cn": cn.get("major_cn", ""),
            "minor_cn": cn.get("minor_cn", ""),
            "loh_status": cn.get("loh_status", "unknown"),
            "cnv_confidence": cn.get("cnv_confidence", "low"),
            "multiplicity_candidates": mult_candidates,
            "multiplicity_best": best_m,
            "multiplicity_confidence": mult_conf,
            "multiplicity_ambiguity": mult_ambiguity,
            "ccf_best": "" if ccf_best is None else f"{ccf_best:.4f}",
            "ccf_min": "" if ccf_min is None else f"{ccf_min:.4f}",
            "ccf_max": "" if ccf_max is None else f"{ccf_max:.4f}",
            "ccf_ci_low": "" if ci_low is None else f"{ci_low:.4f}",
            "ccf_ci_high": "" if ci_high is None else f"{ci_high:.4f}",
            "probability_ccf_gt_0_8": _probability_from_interval(ccf_best, ci_low, ci_high, 0.8, "gt"),
            "probability_ccf_lt_0_3": _probability_from_interval(ccf_best, ci_low, ci_high, 0.3, "lt"),
            "ccf_estimate": "" if ccf_best is None else f"{ccf_best:.4f}",
            "ccf_status": status,
            "clonality_status": status,
            "clonality_confidence": clon_conf,
            "clonality_multiplier": f"{mult:.4f}",
            "ccf_method": method,
            "ccf_confidence": confidence,
            "ccf_warning": ";".join(dict.fromkeys(warning)) if warning else "copy_number_aware_ccf",
            "mutation_multiplicity_assumption": best_m or "unresolved",
            "total_copy_number": f"{total_cn:.4f}",
        }
        row.update(_resolve_external(row, external, svclone))
        resolution, reason = _resolve_ccf_resolution(row)
        row["ccf_resolution"] = resolution
        row["ccf_resolution_reason"] = reason
        rows.append(row)

    write_tsv(out_path, rows, CCF_FIELDS)
    write_tsv(input_qc_path, [qc_row], CCF_INPUT_QC_FIELDS)
    conflicts = [r for r in (_conflict_row(row) for row in rows) if r]
    write_tsv(conflicts_path, conflicts, CCF_CONFLICT_FIELDS)
    write_tsv(clusters_path, _cluster_rows(external), CCF_CLUSTER_FIELDS)
    return rows
