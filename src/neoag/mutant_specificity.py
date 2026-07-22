from __future__ import annotations

from typing import Any, Mapping

from .utils import to_float


def _number(row: Mapping[str, Any], *keys: str) -> float | None:
    for key in keys:
        raw = str(row.get(key, "")).strip()
        if raw and raw.upper() not in {"NA", "N/A", ".", "NONE"}:
            try:
                return float(raw)
            except ValueError:
                continue
    return None


def _mutation_geometry(mt: str, wt: str) -> tuple[str, str, str]:
    mt = str(mt or "").strip().upper()
    wt = str(wt or "").strip().upper()
    if not mt or not wt or len(mt) != len(wt):
        return "", "unknown", "unknown"
    positions = [i + 1 for i, (m, w) in enumerate(zip(mt, wt)) if m != w]
    if not positions:
        return "", "no", "no"
    anchors = {2, len(mt)}
    anchor_only = all(pos in anchors for pos in positions)
    tcr_facing = any(pos not in anchors for pos in positions)
    return ",".join(str(pos) for pos in positions), "yes" if anchor_only else "no", "yes" if tcr_facing else "no"


def evaluate_mutant_specificity(
    peptide: Mapping[str, Any],
    presentation: Mapping[str, Any] | None,
    profile: Mapping[str, Any],
) -> dict[str, Any]:
    """Evaluate whether presentation evidence is mutant-specific rather than merely strong."""
    cfg = profile.get("mutant_specificity", {})
    pres = presentation or {}
    mt_peptide = str(peptide.get("peptide", "")).strip().upper()
    wt_peptide = str(peptide.get("wildtype_peptide", "")).strip().upper()
    positions, anchor_only, tcr_facing = _mutation_geometry(mt_peptide, wt_peptide)

    result: dict[str, Any] = {
        "mutation_positions_in_peptide": positions,
        "mutation_anchor_only": anchor_only,
        "mutation_tcr_facing": tcr_facing,
        "mutant_specificity_status": "UNASSESSED",
        "mutant_specificity_gate_status": "UNASSESSED",
        "mutant_specificity_reason": "wildtype_peptide_or_prediction_missing",
        "mutant_specificity_multiplier": "1.0000",
        "mutant_specificity_priority_cap": "",
    }
    if not wt_peptide:
        return result
    if mt_peptide == wt_peptide:
        result.update({
            "mutant_specificity_status": "NON_MUTANT_SEQUENCE",
            "mutant_specificity_gate_status": "FAIL",
            "mutant_specificity_reason": "mutant_and_wildtype_peptides_identical",
            "mutant_specificity_multiplier": f"{float(cfg.get('non_mutant_multiplier', 0.0)):.4f}",
            "mutant_specificity_priority_cap": str(cfg.get("non_mutant_priority_cap", "D")),
        })
        return result

    mt_el = _number(pres, "netmhcpan_mt_rank_el", "netmhcpan_el_rank")
    if mt_el is None:
        mt_el = _number(peptide, "netmhcpan_mt_rank_el", "netmhcpan_el_rank", "el_rank")
    wt_el = _number(pres, "netmhcpan_wt_rank_el")
    if wt_el is None:
        wt_el = _number(peptide, "netmhcpan_wt_rank_el")
    mt_mhcflurry = _number(pres, "mhcflurry_presentation_score")
    if mt_mhcflurry is None:
        mt_mhcflurry = _number(peptide, "mhcflurry_presentation_score")
    wt_mhcflurry = _number(pres, "mhcflurry_wt_presentation_score")
    if wt_mhcflurry is None:
        wt_mhcflurry = _number(peptide, "mhcflurry_wt_presentation_score")
    mt_prime = _number(pres, "prime_score")
    if mt_prime is None:
        mt_prime = _number(peptide, "prime_score")
    wt_prime = _number(pres, "prime_wt_score")
    if wt_prime is None:
        wt_prime = _number(peptide, "prime_wt_score")
    mt_bigmhc = _number(pres, "bigmhc_im_score")
    if mt_bigmhc is None:
        mt_bigmhc = _number(peptide, "bigmhc_im_score")
    wt_bigmhc = _number(pres, "bigmhc_im_wt_score")
    if wt_bigmhc is None:
        wt_bigmhc = _number(peptide, "bigmhc_im_wt_score")

    if mt_el is not None:
        result["mt_wt_el_rank_difference"] = f"{(wt_el - mt_el):.6f}" if wt_el is not None else ""
        result["agretopicity_el"] = f"{(wt_el / max(mt_el, 1e-6)):.6f}" if wt_el is not None else ""
    if mt_mhcflurry is not None:
        result["mhcflurry_mt_wt_presentation_difference"] = f"{(mt_mhcflurry - wt_mhcflurry):.6f}" if wt_mhcflurry is not None else ""
    if mt_prime is not None:
        result["prime_mt_wt_score_difference"] = f"{(mt_prime - wt_prime):.6f}" if wt_prime is not None else ""
    if mt_bigmhc is not None:
        result["bigmhc_mt_wt_score_difference"] = f"{(mt_bigmhc - wt_bigmhc):.6f}" if wt_bigmhc is not None else ""

    if mt_el is None or wt_el is None:
        return result

    delta = wt_el - mt_el  # Positive means the mutant has the lower (better) EL rank.
    ratio = wt_el / max(mt_el, 1e-6)
    near_equal_abs = float(cfg.get("near_equal_el_rank_difference", 0.01))
    positive_ratio = float(cfg.get("positive_agretopicity_ratio", 2.0))
    positive_delta = float(cfg.get("positive_el_rank_difference", 0.10))

    if wt_el <= mt_el:
        status = "MT_WT_SIMILAR" if abs(delta) <= near_equal_abs else "WT_BETTER"
    elif abs(delta) <= near_equal_abs:
        status = "MT_WT_SIMILAR"
    elif ratio >= positive_ratio and delta >= positive_delta:
        status = "MT_SPECIFIC"
    else:
        status = "MARGINAL_MT_ADVANTAGE"

    if status == "MT_SPECIFIC":
        gate, reason, mult, cap = "PASS", "mutant_el_rank_clearly_better_than_wildtype", 1.0, ""
    else:
        gate = "CAUTION"
        reason = {
            "WT_BETTER": "wildtype_el_rank_better_than_mutant",
            "MT_WT_SIMILAR": "mutant_and_wildtype_el_ranks_similar",
            "MARGINAL_MT_ADVANTAGE": "mutant_el_rank_advantage_below_specificity_threshold",
        }[status]
        mult = float(cfg.get({
            "WT_BETTER": "wt_better_multiplier",
            "MT_WT_SIMILAR": "similar_multiplier",
            "MARGINAL_MT_ADVANTAGE": "marginal_multiplier",
        }[status], {"WT_BETTER": 0.55, "MT_WT_SIMILAR": 0.70, "MARGINAL_MT_ADVANTAGE": 0.85}[status]))
        cap = str(cfg.get("caution_priority_cap", "C_CAUTION"))

    result.update({
        "mutant_specificity_status": status,
        "mutant_specificity_gate_status": gate,
        "mutant_specificity_reason": reason,
        "mutant_specificity_multiplier": f"{mult:.4f}",
        "mutant_specificity_priority_cap": cap,
    })
    return result
