from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping


@dataclass(frozen=True)
class EvidenceState:
    """Availability/consistency state, separate from biological favorability."""

    name: str
    score_modifier: float


STATES = {
    "SUPPORTED": EvidenceState("SUPPORTED", 1.0),
    "PARTIAL": EvidenceState("PARTIAL", 0.75),
    "MISSING": EvidenceState("MISSING", 0.0),
    "CONFLICT": EvidenceState("CONFLICT", 0.5),
}

MISSING_TOKENS = ("UNASSESSED", "MISSING", "NOT_AVAILABLE", "NOT ASSESSED", "UNRESOLVED", "NO_EVIDENCE")
CONFLICT_TOKENS = ("CONFLICT", "DISCORDANT", "INCONSISTENT")
PARTIAL_TOKENS = ("PARTIAL", "CAUTION", "REVIEW", "LOW_CONFIDENCE", "LOW EVIDENCE")
SUPPORTED_TOKENS = ("PASS", "SUPPORTED", "ASSESSED", "INTACT", "COMPLETE", "HIGH_CONFIDENCE")


def evidence_state(value: str | None, *, numeric_present: bool = False) -> str:
    """Classify availability; FAIL/REJECT are assessed negative evidence."""

    text = str(value or "").strip().upper()
    if any(token in text for token in CONFLICT_TOKENS):
        return "CONFLICT"
    if any(token in text for token in MISSING_TOKENS):
        return "MISSING"
    if any(token in text for token in PARTIAL_TOKENS):
        return "PARTIAL"
    if text in {"A", "B"}:
        return "SUPPORTED"
    if text == "C":
        return "PARTIAL"
    if text in {"D", "FAIL", "REJECT", "ESCAPE_REJECT"}:
        return "SUPPORTED"
    if any(token in text for token in SUPPORTED_TOKENS):
        return "SUPPORTED"
    if text:
        return "PARTIAL"
    return "SUPPORTED" if numeric_present else "MISSING"


def _section(rules: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    value = rules.get(name, {}) if isinstance(rules, Mapping) else {}
    return value if isinstance(value, Mapping) else {}


def _float(row: Mapping[str, Any], *fields: str) -> float | None:
    for field in fields:
        raw = str(row.get(field, "")).strip()
        if not raw or raw.upper() in {"NA", "N/A", "NONE", "NULL", ".", "NAN", "UNASSESSED", "UNRESOLVED"}:
            continue
        try:
            return float(raw)
        except ValueError:
            continue
    return None


def _text(row: Mapping[str, Any], *fields: str) -> str:
    return " ".join(str(row.get(field, "")).strip().upper() for field in fields if str(row.get(field, "")).strip())


def _result(
    state: str,
    grade: int,
    reason: str,
    *,
    assessed: bool = True,
    hard_fail: bool = False,
    hard_code: str = "",
    conflict: bool = False,
) -> dict[str, Any]:
    return {
        "state": state,
        "grade": max(0, min(3, int(grade))),
        "assessed": bool(assessed),
        "hard_fail": bool(hard_fail),
        "hard_code": hard_code,
        "conflict": bool(conflict),
        "reason": reason,
    }


def event_track(row: Mapping[str, Any]) -> str:
    text = _text(row, "event_type", "mutation_source", "peptide_consequence")
    if "FUSION" in text:
        return "FUSION"
    if "SPLICE" in text or "JUNCTION" in text:
        return "SPLICE"
    if any(token in text for token in ("SV", "BND", "STRUCTURAL")):
        return "DNA_SV"
    if any(token in text for token in ("FRAMESHIFT", "FRAME_SHIFT", "FS_VARIANT")):
        return "FRAMESHIFT"
    if any(token in text for token in ("SNV", "INDEL", "MISSENSE", "INFRAME", "IN_FRAME")):
        return "MISSENSE"
    return "OTHER"


def derive_event_authenticity(row: Mapping[str, Any], rules: Mapping[str, Any]) -> dict[str, Any]:
    track = event_track(row)
    cross = _text(row, "cross_platform_status", "comparison_status")
    phasing = _text(row, "haplotype_status", "phase_confidence")
    filter_status = _text(row, "filter_status")
    normal_junction = _text(row, "normal_junction_assessment_status", "normal_junction_status")
    if "SOURCE_PASS_NOT_REPRODUCED_BY_PILEUP" in cross or any(token in cross for token in CONFLICT_TOKENS):
        return _result("EVENT_CONFLICT", 0, f"discordant event evidence: {cross}", conflict=True)
    normal_junction_assessed = normal_junction and not any(token in normal_junction for token in ("UNASSESSED", "MISSING", "NOT_AVAILABLE", "NOT_DETECTED"))
    if "NORMAL_SUPPORT_REVIEW" in cross or normal_junction_assessed and any(token in normal_junction for token in ("DETECTED", "EXACT_MATCH", "SUPPORTED_IN_NORMAL")):
        if normal_junction_assessed:
            return _result("EVENT_ARTIFACT_RISK", 0, f"normal junction detected: {normal_junction}", hard_fail=True, hard_code="HARD_NORMAL_JUNCTION")
        return _result("EVENT_ARTIFACT_RISK", 0, f"normal-background review: {cross or normal_junction}")
    if any(token in filter_status for token in ("FAIL", "REJECT", "ARTIFACT")):
        return _result("EVENT_ARTIFACT_RISK", 0, f"filter_status={filter_status}", hard_fail=True, hard_code="HARD_EVENT_ARTIFACT")
    normal_alt = _float(row, "normal_alt_count")
    normal_vaf = _float(row, "normal_alt_vaf")
    normal_status = _text(row, "matched_normal_status", "normal_support_status")
    if (normal_alt is not None and normal_alt >= 2 and normal_vaf is not None and normal_vaf >= 0.01) or "ALT_SUPPORTED" in normal_status:
        return _result("EVENT_ARTIFACT_RISK", 0, f"matched-normal ALT support: count={normal_alt}; VAF={normal_vaf}; status={normal_status}", hard_fail=True, hard_code="HARD_MATCHED_NORMAL_SUPPORT")
    if "CROSS_PLATFORM_PASS_CONCORDANT" in cross:
        return _result("EVENT_CONFIRMED", 3, cross)
    if "ALT_PRESENT_BELOW_PASS_OR_CALLER_DIFFERENCE" in cross:
        return _result("EVENT_STRONG", 2, cross)
    if "COVERED_NO_ALT_SAMPLE_OR_ASSAY_DIFFERENCE" in cross:
        return _result("EVENT_SAMPLE_SPECIFIC", 1, cross)
    if track in {"FUSION", "SPLICE"}:
        reads = _float(row, "rna_junction_reads", "junction_reads")
        minimum = float(_section(rules, "rna").get("junction_min_reads", 3))
        strong = float(_section(rules, "rna").get("junction_strong_reads", 10))
        frame = _text(row, "rna_frame_status", "frame", "in_frame")
        tools = str(row.get("tools_detected", "")).strip()
        try:
            caller_count = int(float(tools))
        except ValueError:
            caller_count = len([value for value in tools.replace(";", ",").split(",") if value.strip()])
        if reads is not None and reads >= strong and caller_count >= 2 and any(token in frame for token in ("IN_FRAME", "INFRAME", "FRAME_OK")):
            return _result("EVENT_CONFIRMED", 3, f"junction reads={reads:g}; callers={caller_count}; frame={frame}")
        if reads is not None and reads >= minimum:
            return _result("EVENT_STRONG", 2, f"junction reads={reads:g}; callers={caller_count}; frame={frame or 'unassessed'}")
        if reads is not None:
            return _result("EVENT_PARTIAL", 1, f"junction reads={reads:g} below {minimum:g}")
    if "PASS" in cross or "CONCORDANT" in cross or "HIGH_CONFIDENCE" in phasing:
        return _result("EVENT_STRONG", 2, f"cross-platform/phasing support: {cross or phasing}")
    if cross or phasing:
        return _result("EVENT_PARTIAL", 1, f"event evidence requires review: {cross or phasing}")
    return _result("EVENT_UNASSESSED", 0, "event authenticity evidence unavailable", assessed=False)


def derive_rna_support(row: Mapping[str, Any], rules: Mapping[str, Any]) -> dict[str, Any]:
    cfg = _section(rules, "rna")
    track = event_track(row)
    status = _text(row, "rna_support_status", "expression_evidence_status", "rna_evidence_completeness")
    if "RNA_ALT_SUPPORTED" in status or "RNA_JUNCTION_SUPPORTED" in status:
        return _result("RNA_CONFIRMED", 3, status)
    if "RNA_ALT_NOT_DETECTED" in status or "RNA_JUNCTION_NOT_DETECTED" in status:
        return _result("RNA_NEGATIVE", 0, status)
    if any(token in status for token in CONFLICT_TOKENS):
        return _result("RNA_LOW_SUPPORT", 1, status, conflict=True)
    if track in {"FUSION", "SPLICE"}:
        reads = _float(row, "rna_junction_reads", "junction_reads")
        minimum = float(cfg.get("junction_min_reads", 3))
        strong = float(cfg.get("junction_strong_reads", 10))
        if reads is None:
            gene_tpm = _float(row, "gene_expression_tpm", "expression_tpm")
            if gene_tpm is not None:
                return _result("GENE_EXPRESSION_ONLY", 1, f"gene TPM={gene_tpm:g}; junction unassessed")
            return _result("RNA_UNASSESSED", 0, "junction reads unavailable", assessed=False)
        if reads >= strong:
            return _result("RNA_CONFIRMED", 3, f"junction reads={reads:g} >= {strong:g}")
        if reads >= minimum:
            return _result("RNA_CONFIRMED", 3, f"junction reads={reads:g} >= {minimum:g}")
        if reads == 0:
            return _result("RNA_NEGATIVE", 0, "junction reads=0")
        return _result("RNA_LOW_SUPPORT", 1, f"junction reads={reads:g} below {minimum:g}")
    alt_reads = _float(row, "rna_alt_reads")
    vaf = _float(row, "rna_vaf")
    minimum_reads = float(cfg.get("snv_min_alt_reads", 3))
    minimum_vaf = float(cfg.get("snv_min_vaf", 0.02))
    if alt_reads is None and vaf is None:
        gene_tpm = _float(row, "gene_expression_tpm", "expression_tpm")
        if gene_tpm is not None:
            return _result("GENE_EXPRESSION_ONLY", 1, f"gene TPM={gene_tpm:g}; mutant allele unassessed")
        return _result("RNA_UNASSESSED", 0, "RNA alt reads/VAF unavailable", assessed=False)
    if (alt_reads or 0) >= minimum_reads and (vaf or 0) >= minimum_vaf:
        grade = 3 if (alt_reads or 0) >= 2 * minimum_reads else 2
        return _result("RNA_CONFIRMED", grade, f"RNA alt reads={alt_reads}; VAF={vaf}")
    if (alt_reads or 0) == 0:
        return _result("RNA_NEGATIVE", 0, f"RNA alt reads={alt_reads}; VAF={vaf}")
    return _result("RNA_LOW_SUPPORT", 1, f"RNA support below provisional thresholds: reads={alt_reads}; VAF={vaf}")


def derive_presentation_consensus(row: Mapping[str, Any], rules: Mapping[str, Any]) -> dict[str, Any]:
    cfg = _section(rules, "presentation")
    gate = _text(row, "presentation_gate_status")
    if any(token in gate for token in ("FAIL", "REJECT", "BLOCKED")):
        return _result("PRESENTATION_WEAK", 0, f"presentation gate={gate}")
    el = _float(row, "netmhcpan_mt_rank_el", "netmhcpan_el_rank", "el_rank", "binding_rank")
    net_support = el is not None and el <= float(cfg.get("netmhcpan_el_supported", 2.0))
    net_strong = el is not None and el <= float(cfg.get("netmhcpan_el_strong", 0.5))
    mhcf = _float(row, "mhcflurry_presentation_score")
    mhcf_support = mhcf is not None and mhcf >= float(cfg.get("mhcflurry_presentation_supported", 0.50))
    mhcf_strong = mhcf is not None and mhcf >= float(cfg.get("mhcflurry_presentation_strong", 0.70))
    stability = _float(row, "netmhcstabpan_rank")
    stability_support = stability is not None and stability <= float(cfg.get("netmhcstabpan_supported", 1.4))
    prime = _float(row, "prime_rank")
    bigmhc = _float(row, "bigmhc_im_score")
    deep = _float(row, "deepimmuno_score")
    immunogenicity_support = (
        (prime is not None and prime <= float(cfg.get("prime_rank_supported", 5.0)))
        or (bigmhc is not None and bigmhc >= float(cfg.get("bigmhc_im_supported", 0.50)))
        or (deep is not None and deep >= float(cfg.get("deepimmuno_supported", 0.50)))
    )
    auxiliary = f"stability={'support' if stability_support else 'no/NA'}; immunogenicity_group={'support' if immunogenicity_support else 'no/NA'}"
    if el is not None and mhcf is not None:
        if net_support and mhcf_support:
            strength = "strong" if net_strong and mhcf_strong else "supported"
            return _result("PRESENTATION_CONSISTENT_STRONG", 3, f"NetMHCpan+MHCflurry {strength}; {auxiliary}")
        if net_support != mhcf_support:
            return _result("PRESENTATION_DISCORDANT", 1, f"NetMHCpan={net_support}; MHCflurry={mhcf_support}; {auxiliary}", conflict=True)
        return _result("PRESENTATION_WEAK", 0, f"both core groups below threshold; {auxiliary}")
    if el is not None or mhcf is not None:
        support = net_support if el is not None else mhcf_support
        return _result("PRESENTATION_SINGLE_TOOL" if support else "PRESENTATION_WEAK", 2 if support else 0, f"single core group; NetMHCpan={net_support if el is not None else 'NA'}; MHCflurry={mhcf_support if mhcf is not None else 'NA'}; {auxiliary}")
    if stability is not None or prime is not None or bigmhc is not None or deep is not None:
        return _result("PRESENTATION_UNASSESSED", 0, f"core groups unavailable; auxiliary evidence cannot replace them; {auxiliary}", assessed=False)
    return _result("PRESENTATION_UNASSESSED", 0, "presentation evidence unavailable", assessed=False)


def derive_mutant_specificity(row: Mapping[str, Any], rules: Mapping[str, Any]) -> dict[str, Any]:
    del rules
    gate = _text(row, "mutant_specificity_gate_status", "mutant_specificity_status")
    mutant = str(row.get("peptide") or row.get("mutant_peptide") or "").strip().upper()
    wildtype = str(row.get("wildtype_peptide") or "").strip().upper()
    novel = str(row.get("contains_novel_aa", "")).strip().lower()
    crosses = str(row.get("crosses_junction", "")).strip().lower()
    if mutant and wildtype and mutant == wildtype:
        return _result("NON_MUTANT_SEQUENCE", 0, "mutant peptide equals wild-type peptide", hard_fail=True, hard_code="HARD_NON_MUTANT_SEQUENCE")
    if novel in {"false", "no", "0"} and crosses in {"false", "no", "0"} and not str(row.get("mutation_positions_in_peptide", "")).strip():
        return _result("NON_MUTANT_SEQUENCE", 0, "peptide contains neither mutation nor novel junction", hard_fail=True, hard_code="HARD_NON_MUTANT_SEQUENCE")
    if "NON_MUTANT_SEQUENCE" in gate:
        return _result("NON_MUTANT_SEQUENCE", 0, gate, hard_fail=True, hard_code="HARD_NON_MUTANT_SEQUENCE")
    if "WT_BETTER" in gate:
        return _result("WT_BETTER", 0, gate)
    if "MT_WT_SIMILAR" in gate:
        return _result("MT_WT_SIMILAR", 1, gate)
    if "MARGINAL_MT_ADVANTAGE" in gate:
        return _result("MARGINAL_MT_ADVANTAGE", 2, gate)
    if "MT_SPECIFIC" in gate or "PASS" in gate or "MT_BETTER" in gate:
        return _result("MT_SPECIFIC", 3, gate)
    if "UNASSESSED" in gate or not gate:
        return _result("UNASSESSED", 0, "mutant specificity unassessed", assessed=False)
    if any(token in gate for token in ("CAUTION", "REVIEW", "EQUIVALENT")):
        return _result("MARGINAL_MT_ADVANTAGE", 2, gate)
    agretopicity = _float(row, "agretopicity_el")
    difference = _float(row, "mt_wt_el_rank_difference")
    if agretopicity is not None or difference is not None:
        if (agretopicity or 0) > 1 or (difference or 0) > 0:
            return _result("MARGINAL_MT_ADVANTAGE", 2, f"agretopicity={agretopicity}; EL difference={difference}")
        return _result("MT_WT_SIMILAR", 1, f"no positive MT-vs-WT differential: agretopicity={agretopicity}; difference={difference}")
    return _result("UNASSESSED", 0, "MT-vs-WT evidence unavailable", assessed=False)


def derive_clonality_state(row: Mapping[str, Any], rules: Mapping[str, Any]) -> dict[str, Any]:
    del rules
    resolution = _text(row, "ccf_resolution", "ccf_status", "clonality_status", "ccf_confidence")
    ccf = _float(row, "ccf_estimate", "ccf_best", "l3_clonality_score")
    if any(token in resolution for token in CONFLICT_TOKENS):
        return _result("CONFLICT", 0, resolution, conflict=True)
    if ccf is None or "UNRESOLVED" in resolution:
        return _result("UNASSESSED", 0, resolution or "CCF unavailable", assessed=False)
    if ccf >= 0.8 and "LOW" not in resolution:
        return _result("CLONAL", 3, f"CCF={ccf:g}; {resolution}")
    if ccf >= 0.5:
        return _result("SUPPORTED", 2, f"CCF={ccf:g}; {resolution}")
    return _result("SUBCLONAL", 1, f"CCF={ccf:g}; {resolution}")


def derive_hla_appm_state(row: Mapping[str, Any], rules: Mapping[str, Any]) -> dict[str, Any]:
    del rules
    status = _text(row, "appm_integrity_status", "appm_evidence_completeness", "escape_status")
    lost = _text(row, "restricting_hla_lost")
    if any(token in status for token in CONFLICT_TOKENS):
        return _result("CONFLICT", 0, status, conflict=True)
    if lost in {"YES", "TRUE", "1"}:
        return _result("RESTRICTING_HLA_LOST", 0, f"restricting HLA confirmed lost: {lost}", hard_fail=True, hard_code="HARD_RESTRICTING_HLA_LOST")
    action = _text(row, "appm_action")
    if "ESCAPE_REJECT" in status or any(token in action for token in ("REJECT", "BLOCK")):
        return _result("MAJOR_APPM_DEFECT", 0, f"major APPM defect: {status}; action={action}")
    loh_status = _text(row, "hla_loh_status")
    if not lost and (not loh_status or "UNASSESSED" in loh_status):
        return _result("HLA_LOH_UNASSESSED", 1, "restricting-HLA LOH not assessed", assessed=False)
    value = _float(row, "l3_apm_integrity_score", "appm_multiplier", "escape_multiplier")
    if value is None and not status:
        return _result("HLA_LOH_UNASSESSED", 0, "HLA/APM evidence unavailable", assessed=False)
    if (value is not None and value >= 0.8) and not any(token in status for token in PARTIAL_TOKENS):
        return _result("HLA_APPM_RETAINED", 3, f"HLA/APM score={value:g}; {status}")
    if value is not None and value >= 0.5:
        return _result("HLA_APPM_RETAINED", 2, f"HLA/APM score={value:g}; {status}")
    return _result("HLA_APPM_CAUTION", 1, f"HLA/APM score={value}; {status}")


def derive_safety_state(row: Mapping[str, Any], rules: Mapping[str, Any]) -> dict[str, Any]:
    del rules
    status = _text(
        row, "safety_status", "safety_tier", "safety_evidence_completeness",
        "event_safety_status", "event_safety_tier", "event_safety_evidence_completeness",
    )
    exact = str(row.get("reference_proteome_exact_match", "")).strip().lower() in {"true", "yes", "1"}
    if exact:
        return _result("SAFETY_REJECT", 0, "reference proteome exact match", hard_fail=True, hard_code="HARD_REFERENCE_PROTEOME_MATCH")
    normal_junction = _text(row, "normal_junction_assessment_status", "event_normal_junction_assessment_status")
    if not any(token in normal_junction for token in ("NOT_DETECTED", "UNASSESSED", "MISSING")) and any(token in normal_junction for token in ("DETECTED", "EXACT_MATCH", "SUPPORTED_IN_NORMAL")):
        return _result("SAFETY_REJECT", 0, f"normal junction clearly detected: {normal_junction}", hard_fail=True, hard_code="HARD_NORMAL_JUNCTION")
    if any(token in status for token in ("UNSAFE", "FAIL", "REJECT", "BLOCKED")):
        return _result("SAFETY_REJECT", 0, status, hard_fail=True, hard_code="HARD_SAFETY_REJECT")
    if not status:
        return _result("SAFETY_PARTIAL", 1, "safety evidence unavailable", assessed=False)
    missing_layers = _text(row, "safety_missing_layers")
    reference_statuses = _text(
        row, "reference_proteome_status", "normal_ligandome_status",
        "normal_junction_assessment_status",
    )
    if missing_layers or any(token in reference_statuses for token in MISSING_TOKENS):
        return _result("SAFETY_PARTIAL", 1, f"missing safety layers={missing_layers or reference_statuses}")
    if any(token in status for token in ("HIGH_RISK", "CRITICAL_TISSUE")):
        return _result("SAFETY_HIGH_RISK", 0, status)
    if any(token in status for token in ("CAUTION", "REVIEW")):
        return _result("SAFETY_REVIEW", 1, status)
    if "PASS" in status and not any(token in status for token in PARTIAL_TOKENS):
        return _result("SAFETY_PASS", 3, status)
    if any(token in status for token in PARTIAL_TOKENS):
        return _result("SAFETY_PARTIAL", 1, status)
    return _result("SAFETY_REVIEW", 2, status)


def derive_evidence_completeness(
    row: Mapping[str, Any],
    track: str,
    rules: Mapping[str, Any],
    states: Mapping[str, Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    cfg = _section(rules, "completeness")
    states = states or {
        "event_authenticity": derive_event_authenticity(row, rules),
        "rna_support": derive_rna_support(row, rules),
        "presentation_consensus": derive_presentation_consensus(row, rules),
        "mutant_specificity": derive_mutant_specificity(row, rules),
        "clonality": derive_clonality_state(row, rules),
        "hla_appm": derive_hla_appm_state(row, rules),
        "safety": derive_safety_state(row, rules),
    }
    required = list(states)
    if track in {"FUSION", "SPLICE"}:
        required = [name for name in required if name != "clonality"]
    explicit_novel = (
        str(row.get("contains_novel_aa", "")).strip().lower() in {"true", "yes", "1"}
        or str(row.get("crosses_junction", "")).strip().lower() in {"true", "yes", "1"}
    )
    if explicit_novel and not states["mutant_specificity"]["assessed"]:
        required = [name for name in required if name != "mutant_specificity"]
    assessed = [name for name in required if states[name]["assessed"]]
    ratio = len(assessed) / len(required) if required else 0.0
    blockers: list[str] = []
    if not bool(cfg.get("allow_r1_with_hla_loh_unassessed", False)) and not states["hla_appm"]["assessed"]:
        blockers.append("hla_appm_unassessed")
    if not bool(cfg.get("allow_r1_with_safety_partial", False)) and states["safety"]["grade"] < 3:
        blockers.append("safety_not_complete")
    if track not in {"FUSION", "SPLICE"} and not bool(cfg.get("allow_r1_with_ccf_low_confidence", False)) and states["clonality"]["grade"] < 2:
        blockers.append("clonality_low_or_unassessed")
    if ratio >= 0.90 and not blockers:
        grade = 3
    elif ratio >= 0.70:
        grade = 2
    elif ratio >= 0.40:
        grade = 1
    else:
        grade = 0
    state = "COMPLETE" if grade == 3 else "PARTIAL" if grade else "LOW"
    return _result(state, grade, f"track={track}; assessed={len(assessed)}/{len(required)}; blockers={','.join(blockers) or 'none'}", assessed=bool(assessed))


DERIVERS = {
    "event_authenticity": derive_event_authenticity,
    "rna_support": derive_rna_support,
    "presentation_consensus": derive_presentation_consensus,
    "mutant_specificity": derive_mutant_specificity,
    "clonality": derive_clonality_state,
    "hla_appm": derive_hla_appm_state,
    "safety": derive_safety_state,
}


def derive_all_states(row: Mapping[str, Any], rules: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    track = event_track(row)
    result = {name: function(row, rules) for name, function in DERIVERS.items()}
    result["evidence_completeness"] = derive_evidence_completeness(row, track, rules, result)
    return result
