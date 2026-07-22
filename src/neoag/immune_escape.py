"""Immune Escape 2.0 evidence layer.

Detects HLA LOH, B2M/APM defects, IFNG-JAK-STAT defects and projects them to
peptide-level flags. The output is immune-escape mechanism evidence, not a
clinical resistance diagnosis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .adapters.peptide_input import hla_allele_class, normalize_hla_allele
from .appm_v2 import (
    IFNG_GENES,
    MHC_I_GENES,
    MHC_II_GENES,
    build_appm_2,
    load_cnv_gene_status,
    load_expression,
    load_hla_loh,
    load_variants,
)
from .utils import first, read_tsv, to_float, write_tsv

KILLING_GENES = {"CASP8", "FAS", "FADD", "TNFRSF10A", "TNFRSF10B", "BAX", "BAK1", "SERPINB9"}

IMMUNE_ESCAPE_EVENT_FIELDS = [
    "sample_id", "mechanism", "gene_or_hla", "pathway", "alteration_type", "biallelic_status",
    "functional_status", "escape_event_source", "escape_event_source_event_ids",
    "escape_event_ccf_best", "escape_event_ccf_min", "escape_event_ccf_max",
    "escape_event_ccf_confidence", "escape_event_clonality", "affected_candidate_count",
    "affected_top_candidate_count", "affected_hla_alleles", "affected_event_ids",
    "affected_peptide_ids", "risk_level", "evidence_level", "reason",
]

IMMUNE_ESCAPE_SUMMARY_FIELDS = [
    "sample_id", "therapy_context", "mhc_i_escape_status", "mhc_ii_escape_status", "ifng_response_status",
    "cytotoxic_killing_resistance_status", "hla_loh_status", "lost_hla_alleles", "b2m_biallelic_loss",
    "lost_hla_i_alleles", "lost_hla_ii_alleles", "unclassified_lost_hla_alleles",
    "jak1_biallelic_loss", "jak2_biallelic_loss", "tap_defect", "nlrc5_defect", "ciita_defect",
    "n_peptides_affected_by_hla_loh", "n_top_peptides_affected_by_hla_loh",
    "n_mhc_i_peptides_affected_by_b2m", "n_top_mhc_i_peptides_affected_by_b2m",
    "escape_burden_summary", "overall_immune_escape_risk", "mechanism_summary",
    "evidence_completeness", "interpretation",
]

PEPTIDE_ESCAPE_FIELDS = [
    "peptide_id", "event_id", "peptide", "hla_allele", "mhc_class", "therapy_context", "restricting_hla_lost",
    "lost_hla_alleles", "b2m_status", "hla_class_i_global_status", "jak_stat_status", "tap_processing_status",
    "nlrc5_status", "ciita_status", "escape_status", "escape_reason", "escape_multiplier", "priority_cap",
]


def _path_exists(path: str | Path | None) -> bool:
    return bool(path and Path(path).exists())


def _read_map(path: str | Path | None, key: str) -> dict[str, dict[str, str]]:
    if not _path_exists(path):
        return {}
    return {r.get(key, ""): r for r in read_tsv(path) if r.get(key)}


def _load_lost_hla_with_conf(path: str | Path | None) -> dict[str, dict[str, str]]:
    if not _path_exists(path):
        return {}
    out: dict[str, dict[str, str]] = {}
    for r in read_tsv(path):
        allele = first(r, ["hla_allele", "allele", "HLA", "LossAllele", "loss_allele"], "")
        status = first(r, ["loh_status", "LOH", "status", "loss", "Loss"], "")
        if allele and (not status or status.lower() in {"loh", "loss", "lost", "yes", "true", "1"}):
            h = normalize_hla_allele(allele)
            out[h] = {
                "hla_allele": h,
                "confidence": first(r, ["confidence", "loh_confidence", "Pval_unique", "pval"], ""),
                "method": first(r, ["method", "tool", "source"], ""),
                "raw_status": status or "loh",
            }
    return out


def _therapy_context(profile: Mapping[str, Any] | None, override: str | None = None) -> str:
    if override:
        return override
    return str((profile or {}).get("immune_escape", {}).get("therapy_context", "vaccine")).lower()


def _strict_lost_hla_policy(ctx: str) -> bool:
    return ctx in {"vaccine", "tcr_target", "tcr-t", "cell_therapy"}


def _risk_rank(level: str) -> int:
    return {"LOW": 0, "MEDIUM": 1, "HIGH": 2, "INCONCLUSIVE": -1}.get(str(level).upper(), 0)


def _gene_status_from_appm(appm_gene_status: str | Path | None) -> dict[str, dict[str, str]]:
    return _read_map(appm_gene_status, "gene")


def _build_or_load_appm(
    *,
    sample_id: str,
    workdir: Path,
    appm_gene_status: str | Path | None,
    appm_pathway_status: str | Path | None,
    vep_tsv: str | Path | None,
    expression_tsv: str | Path | None,
    cnv_tsv: str | Path | None,
    hla_loh_tsv: str | Path | None,
    raw_peptides: str | Path,
    profile: Mapping[str, Any] | None,
) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]], dict[str, str]]:
    if _path_exists(appm_gene_status):
        gene_map = _gene_status_from_appm(appm_gene_status)
        pathway_map = _read_map(appm_pathway_status, "pathway") if _path_exists(appm_pathway_status) else {}
        summary_path = Path(appm_gene_status).parent / "appm_summary.tsv"
        summary_rows = read_tsv(summary_path) if summary_path.exists() else []
        summary = summary_rows[0] if summary_rows else {}
        return gene_map, pathway_map, summary
    tmp = workdir / "appm_2_input_resolved"
    paths = build_appm_2(
        sample_id=sample_id,
        outdir=tmp,
        vep_tsv=vep_tsv,
        expression_tsv=expression_tsv,
        hla_loh_tsv=hla_loh_tsv,
        cnv_tsv=cnv_tsv,
        raw_peptides=raw_peptides,
        profile=profile,
    )
    return _gene_status_from_appm(paths["appm_gene_status"]), _read_map(paths["appm_pathway_status"], "pathway"), read_tsv(paths["appm_summary"])[0]


def _biallelic(gene_map: Mapping[str, Mapping[str, str]], gene: str) -> bool:
    return gene_map.get(gene, {}).get("biallelic_status") == "BIALLELIC_LOSS"


def _functional(gene_map: Mapping[str, Mapping[str, str]], gene: str) -> str:
    return gene_map.get(gene, {}).get("functional_status", "intact")


def _event_ccf_map(ccf_tsv: str | Path | None) -> dict[str, dict[str, str]]:
    if not _path_exists(ccf_tsv):
        return {}
    return {r.get("event_id", ""): r for r in read_tsv(ccf_tsv) if r.get("event_id")}


def _extract_source_event_ids(gene_status: Mapping[str, str]) -> list[str]:
    raw = gene_status.get("source_event_ids") or gene_status.get("source_event_id") or gene_status.get("event_id") or gene_status.get("variant_id") or ""
    return [x.strip() for x in str(raw).replace(",", ";").split(";") if x.strip()]


def _resolve_escape_ccf(source_ids: list[str], ccf_by_event: Mapping[str, Mapping[str, str]]) -> dict[str, str]:
    rows = [ccf_by_event[eid] for eid in source_ids if eid in ccf_by_event]
    if not rows:
        return {
            "escape_event_ccf_best": "", "escape_event_ccf_min": "", "escape_event_ccf_max": "",
            "escape_event_ccf_confidence": "not_mapped", "escape_event_clonality": "unresolved",
        }
    best_row = max(rows, key=lambda r: to_float(r.get("ccf_best", r.get("ccf_estimate", "0")), 0.0))
    best = to_float(best_row.get("ccf_best", best_row.get("ccf_estimate", "0")), 0.0)
    if best >= 0.85:
        clonality = "clonal"
    elif best >= 0.25:
        clonality = "subclonal"
    elif best > 0:
        clonality = "low_frequency"
    else:
        clonality = "unresolved"
    return {
        "escape_event_ccf_best": best_row.get("ccf_best", best_row.get("ccf_estimate", "")),
        "escape_event_ccf_min": best_row.get("ccf_min", ""),
        "escape_event_ccf_max": best_row.get("ccf_max", ""),
        "escape_event_ccf_confidence": best_row.get("ccf_confidence", "mapped"),
        "escape_event_clonality": clonality,
    }


def _escape_event_row(
    sample_id: str, mechanism: str, gene_or_hla: str, pathway: str, alteration: str,
    gene_status: Mapping[str, str] | None = None, risk: str = "MEDIUM", reason: str = "",
    ccf_by_event: Mapping[str, Mapping[str, str]] | None = None,
) -> dict[str, str]:
    gene_status = gene_status or {}
    ccf_by_event = ccf_by_event or {}
    source_ids = _extract_source_event_ids(gene_status)
    ccf = _resolve_escape_ccf(source_ids, ccf_by_event)
    return {
        "sample_id": sample_id,
        "mechanism": mechanism,
        "gene_or_hla": gene_or_hla,
        "pathway": pathway,
        "alteration_type": alteration,
        "biallelic_status": gene_status.get("biallelic_status", ""),
        "functional_status": gene_status.get("functional_status", ""),
        "escape_event_source": "appm_gene_status" if gene_status else "hla_loh_tool",
        "escape_event_source_event_ids": ";".join(source_ids),
        **ccf,
        "affected_candidate_count": "0",
        "affected_top_candidate_count": "0",
        "affected_hla_alleles": "",
        "affected_event_ids": "",
        "affected_peptide_ids": "",
        "risk_level": risk,
        "evidence_level": gene_status.get("evidence_level", "MEDIUM"),
        "reason": reason or gene_status.get("reason", alteration),
    }


def _mhc_class_i(mhc: str) -> bool:
    return str(mhc or "I").upper() in {"I", "MHC-I", "CLASSI"}


def build_immune_escape_evidence(
    *,
    sample_id: str,
    raw_peptides: str | Path,
    outdir: str | Path,
    vep_tsv: str | Path | None = None,
    expression_tsv: str | Path | None = None,
    cnv_tsv: str | Path | None = None,
    hla_loh_tsv: str | Path | None = None,
    profile: Mapping[str, Any] | None = None,
    appm_gene_status: str | Path | None = None,
    appm_pathway_status: str | Path | None = None,
    ccf_tsv: str | Path | None = None,
    therapy_context: str | None = None,
    ranked_peptides: str | Path | None = None,
    top_priority_threshold: str = "B_CAUTION",
) -> dict[str, str]:
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    ctx = _therapy_context(profile, therapy_context)
    lost_hla_info = _load_lost_hla_with_conf(hla_loh_tsv)
    lost_hla = set(lost_hla_info)
    lost_hla_i = {hla for hla in lost_hla if hla_allele_class(hla) == "I"}
    lost_hla_ii = {hla for hla in lost_hla if hla_allele_class(hla) == "II"}
    lost_hla_unclassified = lost_hla - lost_hla_i - lost_hla_ii
    gene_map, pathway_map, appm_summary = _build_or_load_appm(
        sample_id=sample_id,
        workdir=out,
        appm_gene_status=appm_gene_status,
        appm_pathway_status=appm_pathway_status,
        vep_tsv=vep_tsv,
        expression_tsv=expression_tsv,
        cnv_tsv=cnv_tsv,
        hla_loh_tsv=hla_loh_tsv,
        raw_peptides=raw_peptides,
        profile=profile,
    )
    ccf_by_event = _event_ccf_map(ccf_tsv)

    b2m = _biallelic(gene_map, "B2M")
    jak1 = _biallelic(gene_map, "JAK1")
    jak2 = _biallelic(gene_map, "JAK2")
    tap_def = _biallelic(gene_map, "TAP1") or _biallelic(gene_map, "TAP2")
    nlrc5_def = _functional(gene_map, "NLRC5") in {"defective", "caution"}
    ciita_def = _functional(gene_map, "CIITA") in {"defective", "caution"}

    events: list[dict[str, str]] = []
    for hla in sorted(lost_hla_i):
        events.append(_escape_event_row(sample_id, "HLA_ALLELE_SPECIFIC_LOSS", hla, "HLA-I", "HLA_LOH", risk="HIGH" if len(lost_hla_i) >= 2 else "MEDIUM", reason="restricting_HLA_class_I_allele_loss", ccf_by_event=ccf_by_event))
    for hla in sorted(lost_hla_ii):
        events.append(_escape_event_row(sample_id, "HLA_CLASS_II_ALLELE_SPECIFIC_LOSS", hla, "HLA-II", "HLA_LOH", risk="MEDIUM", reason="HLA_class_II_background_loss", ccf_by_event=ccf_by_event))
    if b2m:
        events.append(_escape_event_row(sample_id, "MHC_I_GLOBAL_LOSS", "B2M", "MHC-I", "BIALLELIC_LOSS", gene_map.get("B2M"), risk="HIGH", ccf_by_event=ccf_by_event))
    if tap_def:
        gene = "TAP1" if _biallelic(gene_map, "TAP1") else "TAP2"
        events.append(_escape_event_row(sample_id, "ANTIGEN_PROCESSING_DEFECT", gene, "MHC-I", "BIALLELIC_LOSS", gene_map.get(gene), risk="HIGH", ccf_by_event=ccf_by_event))
    if nlrc5_def:
        events.append(_escape_event_row(sample_id, "MHC_I_TRANSCRIPTION_CAUTION", "NLRC5", "MHC-I", "DEFECT_OR_LOW", gene_map.get("NLRC5"), risk="MEDIUM", ccf_by_event=ccf_by_event))
    if jak1 or jak2:
        gene = "JAK1" if jak1 else "JAK2"
        events.append(_escape_event_row(sample_id, "IFNG_RESPONSE_DEFECT", gene, "IFNG-JAK-STAT", "BIALLELIC_LOF", gene_map.get(gene), risk="HIGH", ccf_by_event=ccf_by_event))
    if ciita_def:
        events.append(_escape_event_row(sample_id, "MHC_II_PRESENTATION_DEFECT", "CIITA", "MHC-II", "DEFECT_OR_LOW", gene_map.get("CIITA"), risk="MEDIUM", ccf_by_event=ccf_by_event))

    mechanisms = [e["mechanism"] for e in events]
    appm_unassessed = appm_summary.get("appm_overall_status") == "UNASSESSED" or appm_summary.get("appm_evidence_completeness") == "UNASSESSED"
    no_direct_escape_inputs = not any([vep_tsv, expression_tsv, cnv_tsv, hla_loh_tsv])
    if b2m or jak1 or jak2:
        overall = "HIGH"
    elif events:
        overall = "MEDIUM"
    elif not events and appm_unassessed:
        overall = "INCONCLUSIVE"
    elif no_direct_escape_inputs and not appm_gene_status:
        overall = "INCONCLUSIVE"
    else:
        overall = "LOW"

    summary = {
        "sample_id": sample_id,
        "therapy_context": ctx,
        "mhc_i_escape_status": "HIGH" if b2m else ("MEDIUM" if lost_hla_i or tap_def or nlrc5_def else "LOW"),
        "mhc_ii_escape_status": "MEDIUM" if ciita_def or lost_hla_ii else "LOW",
        "ifng_response_status": "HIGH" if jak1 or jak2 else "LOW",
        "cytotoxic_killing_resistance_status": "LOW",
        "hla_loh_status": "LOH_DETECTED" if lost_hla else ("NOT_ASSESSED" if not hla_loh_tsv else "NO_LOH_DETECTED"),
        "lost_hla_alleles": ",".join(sorted(lost_hla)),
        "lost_hla_i_alleles": ",".join(sorted(lost_hla_i)),
        "lost_hla_ii_alleles": ",".join(sorted(lost_hla_ii)),
        "unclassified_lost_hla_alleles": ",".join(sorted(lost_hla_unclassified)),
        "b2m_biallelic_loss": "yes" if b2m else "no",
        "jak1_biallelic_loss": "yes" if jak1 else "no",
        "jak2_biallelic_loss": "yes" if jak2 else "no",
        "tap_defect": "yes" if tap_def else "no",
        "nlrc5_defect": "yes" if nlrc5_def else "no",
        "ciita_defect": "yes" if ciita_def else "no",
        "overall_immune_escape_risk": overall,
        "mechanism_summary": ";".join(mechanisms) if mechanisms else "no_major_signal",
        "evidence_completeness": "unassessed" if (not events and appm_unassessed) else ("partial" if any([vep_tsv, expression_tsv, cnv_tsv, hla_loh_tsv, appm_gene_status]) else "unassessed"),
        "interpretation": "immune_escape_evidence_not_clinical_resistance_diagnosis",
    }

    flags: list[dict[str, str]] = []
    strict_lost = _strict_lost_hla_policy(ctx)
    for p in read_tsv(raw_peptides):
        peptide_id = p.get("peptide_id", "")
        event_id = p.get("event_id", "")
        hla = normalize_hla_allele(p.get("hla_allele", ""))
        mhc_i = _mhc_class_i(p.get("mhc_class", "I"))
        restricting_hla_lost = hla in (lost_hla_i if mhc_i else lost_hla_ii)
        reasons: list[str] = []
        status = "ESCAPE_PASS"; mult = 1.0; cap = ""
        if restricting_hla_lost:
            if strict_lost:
                status = "LOST_RESTRICTING_HLA"; mult = 0.0; cap = "D"
            else:
                status = "LOST_RESTRICTING_HLA_RETAINED_FOR_REVIEW"; mult = min(mult, 0.35); cap = "C_CAUTION"
            reasons.append("restricting_hla_lost")
        if mhc_i and b2m:
            status = "GLOBAL_MHC_I_LOSS"; mult = 0.0; cap = "D"; reasons.append("b2m_biallelic_loss")
        elif mhc_i and tap_def:
            if mult > 0:
                status = "APM_PROCESSING_DEFECT" if status == "ESCAPE_PASS" else status
                mult = min(mult, 0.35); cap = cap or "C"
            reasons.append("tap_processing_defect")
        elif mhc_i and nlrc5_def:
            if mult > 0:
                status = "APM_TRANSCRIPTION_CAUTION" if status == "ESCAPE_PASS" else status
                mult = min(mult, 0.70); cap = cap or "B_CAUTION"
            reasons.append("nlrc5_caution")
        if (not mhc_i) and ciita_def:
            status = "MHC_II_PRESENTATION_DEFECT" if status == "ESCAPE_PASS" else status
            mult = min(mult, 0.40); cap = cap or "C"; reasons.append("ciita_or_mhc_ii_defect")
        if jak1 or jak2:
            if mult > 0:
                status = "IFNG_RESPONSE_DEFECT_CAUTION" if status == "ESCAPE_PASS" else status
                mult = min(mult, 0.60); cap = cap or "B_CAUTION"
            reasons.append("jak_stat_defect")
        # Escape event clonality context if the peptide's source event has CCF.
        ccf = ccf_by_event.get(event_id, {})
        if ccf.get("ccf_status") and ccf.get("ccf_status") not in {"clonal_like", ""} and status not in {"ESCAPE_PASS", "GLOBAL_MHC_I_LOSS", "LOST_RESTRICTING_HLA"}:
            reasons.append(f"escape_or_source_ccf_{ccf.get('ccf_status')}")
        if not reasons:
            reasons.append("no_major_signal")
        flags.append({
            "peptide_id": peptide_id,
            "event_id": event_id,
            "peptide": p.get("peptide", ""),
            "hla_allele": p.get("hla_allele", ""),
            "mhc_class": p.get("mhc_class", ""),
            "therapy_context": ctx,
            "restricting_hla_lost": "yes" if restricting_hla_lost else "no",
            "lost_hla_alleles": ",".join(sorted(lost_hla)),
            "b2m_status": "BIALLELIC_LOSS" if b2m else gene_map.get("B2M", {}).get("biallelic_status", "NO_EVIDENCE"),
            "hla_class_i_global_status": "GLOBAL_MHC_I_LOSS" if b2m else "NO_GLOBAL_LOSS_SIGNAL",
            "jak_stat_status": "DEFECT" if jak1 or jak2 else "NO_HIGH_RISK_SIGNAL",
            "tap_processing_status": "DEFECT" if tap_def else "NO_HIGH_RISK_SIGNAL",
            "nlrc5_status": "DEFECT_OR_LOW" if nlrc5_def else "NO_HIGH_RISK_SIGNAL",
            "ciita_status": "DEFECT_OR_LOW" if ciita_def else "NO_HIGH_RISK_SIGNAL",
            "escape_status": status,
            "escape_reason": ";".join(reasons),
            "escape_multiplier": f"{mult:.4f}",
            "priority_cap": cap,
        })


    # v0.4.2 P1: annotate mechanism burden against candidate peptide space.
    ranked_map = {r.get("peptide_id", ""): r for r in read_tsv(ranked_peptides)} if _path_exists(ranked_peptides) else {}
    priority_order = {"A": 0, "B": 1, "B_CAUTION": 2, "C": 3, "C_CAUTION": 4, "D": 5, "REJECT": 6, "": 99}
    top_cut = priority_order.get(str(top_priority_threshold or "B_CAUTION"), 2)

    def _is_top(pid: str) -> bool:
        if not ranked_map:
            return False
        pr = ranked_map.get(pid, {}).get("final_priority", ranked_map.get(pid, {}).get("priority", ""))
        return priority_order.get(pr, 99) <= top_cut

    def _affects_event(ev: Mapping[str, str], flag: Mapping[str, str]) -> bool:
        mech = ev.get("mechanism", "")
        hla = normalize_hla_allele(flag.get("hla_allele", ""))
        mhc_i_flag = _mhc_class_i(flag.get("mhc_class", "I"))
        if mech == "HLA_ALLELE_SPECIFIC_LOSS":
            return mhc_i_flag and hla == normalize_hla_allele(ev.get("gene_or_hla", ""))
        if mech == "HLA_CLASS_II_ALLELE_SPECIFIC_LOSS":
            return (not mhc_i_flag) and hla == normalize_hla_allele(ev.get("gene_or_hla", ""))
        if mech == "MHC_I_GLOBAL_LOSS":
            return mhc_i_flag
        if mech in {"ANTIGEN_PROCESSING_DEFECT", "MHC_I_TRANSCRIPTION_CAUTION"}:
            return mhc_i_flag
        if mech == "MHC_II_PRESENTATION_DEFECT":
            return not mhc_i_flag
        if mech == "IFNG_RESPONSE_DEFECT":
            return True
        return False

    annotated_events: list[dict[str, str]] = []
    for ev in events:
        affected = [f for f in flags if _affects_event(ev, f)]
        top_affected = [f for f in affected if _is_top(f.get("peptide_id", ""))]
        ev = dict(ev)
        ev.update({
            "affected_candidate_count": str(len(affected)),
            "affected_top_candidate_count": str(len(top_affected)),
            "affected_hla_alleles": ";".join(sorted({f.get("hla_allele", "") for f in affected if f.get("hla_allele")})),
            "affected_event_ids": ";".join(sorted({f.get("event_id", "") for f in affected if f.get("event_id")})),
            "affected_peptide_ids": ";".join(f.get("peptide_id", "") for f in affected[:50] if f.get("peptide_id")),
        })
        annotated_events.append(ev)
    events = annotated_events

    hla_loh_affected = [f for f in flags if f.get("restricting_hla_lost") == "yes"]
    hla_loh_top = [f for f in hla_loh_affected if _is_top(f.get("peptide_id", ""))]
    b2m_affected = [f for f in flags if f.get("hla_class_i_global_status") == "GLOBAL_MHC_I_LOSS"]
    b2m_top = [f for f in b2m_affected if _is_top(f.get("peptide_id", ""))]
    material_events = [ev for ev in events if int(ev.get("affected_candidate_count", "0") or 0) > 0]
    completeness_status = str(appm_summary.get("appm_evidence_completeness", "")).upper()
    if b2m or jak1 or jak2 or hla_loh_top:
        overall = "HIGH"
    elif material_events:
        overall = "MEDIUM"
    elif events or completeness_status == "LOW":
        overall = "REVIEW_REQUIRED"
    elif appm_unassessed:
        overall = "INCONCLUSIVE"
    else:
        overall = "LOW"
    summary.update({
        "n_peptides_affected_by_hla_loh": str(len(hla_loh_affected)),
        "n_top_peptides_affected_by_hla_loh": str(len(hla_loh_top)),
        "n_mhc_i_peptides_affected_by_b2m": str(len(b2m_affected)),
        "n_top_mhc_i_peptides_affected_by_b2m": str(len(b2m_top)),
        "escape_burden_summary": f"hla_loh={len(hla_loh_affected)};b2m_mhc_i={len(b2m_affected)};top_hla_loh={len(hla_loh_top)};top_b2m={len(b2m_top)}",
        "overall_immune_escape_risk": overall,
    })

    paths = {
        "immune_escape_events": str(out / "immune_escape_events.tsv"),
        "immune_escape_summary": str(out / "immune_escape_summary.tsv"),
        "peptide_escape_flags": str(out / "peptide_escape_flags.tsv"),
    }
    write_tsv(paths["immune_escape_events"], events, IMMUNE_ESCAPE_EVENT_FIELDS)
    write_tsv(paths["immune_escape_summary"], [summary], IMMUNE_ESCAPE_SUMMARY_FIELDS)
    write_tsv(paths["peptide_escape_flags"], flags, PEPTIDE_ESCAPE_FIELDS)
    return paths


def load_peptide_escape_flags(path: str | Path | None) -> dict[str, dict[str, str]]:
    if not path or not Path(path).exists():
        return {}
    return {r.get("peptide_id", ""): r for r in read_tsv(path) if r.get("peptide_id")}
