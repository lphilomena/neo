"""Immune escape / HLA LOH resistance module (independent from APPM-lite expression penalties)."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .adapters.peptide_input import normalize_hla_allele
from .evidence_provenance import ProvenanceRecord, attach_provenance, provenance_derived
from .schemas import (
    IMMUNE_ESCAPE_EVENT_FIELDS,
    IMMUNE_ESCAPE_SUMMARY_FIELDS,
    PEPTIDE_ESCAPE_FLAG_FIELDS,
)
from .utils import clamp, first, read_tsv, safe_id, to_float, write_tsv

RESISTANCE_RISKS = frozenset({"HIGH", "MEDIUM", "LOW", "INCONCLUSIVE"})

DAMAGING_LOF = frozenset({
    "frameshift_variant",
    "stop_gained",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "stop_lost",
})

APM_GENES = ["TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5", "TAPBP"]
IFNG_GENES = ["JAK1", "JAK2", "STAT1", "IFNGR1", "IFNGR2"]
MHC_I_GENES = ["B2M", "HLA-A", "HLA-B", "HLA-C", *APM_GENES, *IFNG_GENES]
MHC_II_GENES = ["HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "CIITA", "RFX5", "RFXANK", "RFXAP"]
HLA_CLASS_I = ["HLA-A", "HLA-B", "HLA-C"]


def _norm_risk(level: str) -> str:
    u = str(level or "").strip().upper()
    return u if u in RESISTANCE_RISKS else "INCONCLUSIVE"


def _max_risk(*levels: str) -> str:
    order = {"INCONCLUSIVE": 0, "LOW": 1, "MEDIUM": 2, "HIGH": 3}
    best = "INCONCLUSIVE"
    for lv in levels:
        lv = _norm_risk(lv)
        if order[lv] > order[best]:
            best = lv
    return best


def load_expression(path: str | Path | None) -> dict[str, float]:
    if not path or not Path(path).exists():
        return {}
    out: dict[str, float] = {}
    for row in read_tsv(path):
        gene = first(row, ["gene", "Gene", "gene_name"], "")
        if gene:
            out[gene] = to_float(first(row, ["TPM", "tpm", "expression"], "0"), 0.0)
    return out


def load_variants(path: str | Path | None) -> dict[str, list[str]]:
    if not path or not Path(path).exists():
        return {}
    out: dict[str, list[str]] = {}
    for row in read_tsv(path):
        gene = first(row, ["gene", "Gene", "SYMBOL"], "")
        cons = first(row, ["consequence", "Consequence"], "")
        if gene and cons:
            out.setdefault(gene, []).append(cons)
    return out


def load_loh_alleles(path: str | Path | None) -> set[str]:
    if not path or not Path(path).exists():
        return set()
    alleles: set[str] = set()
    for row in read_tsv(path):
        allele = first(row, ["hla_allele", "allele", "HLA"], "")
        status = first(row, ["loh_status", "LOH", "status"], "")
        if allele and str(status).lower() in {"loh", "loss", "lost", "yes", "true"}:
            alleles.add(normalize_hla_allele(allele))
    return alleles


def load_hla_typing(path: str | Path | None = None, alleles: list[str] | None = None) -> list[str]:
    typed: list[str] = []
    if alleles:
        typed.extend(normalize_hla_allele(a) for a in alleles if str(a).strip())
    if path and Path(path).exists():
        for row in read_tsv(path):
            for key in ("hla_allele", "allele", "HLA", "allele_name"):
                val = first(row, [key], "")
                if val:
                    typed.append(normalize_hla_allele(val))
    seen: set[str] = set()
    out: list[str] = []
    for a in typed:
        if a and a not in seen:
            seen.add(a)
            out.append(a)
    return out


def load_cnv_homdel_flag(path: str | Path | None, *, homdel_threshold: float) -> bool:
    if not path or not Path(path).exists():
        return False
    for row in read_tsv(path):
        cn = to_float(first(row, ["total_cn", "copy_number"], "2"), 2.0)
        if cn <= homdel_threshold:
            return True
    return False


def _damaging_lof_hits(consequences: list[str]) -> list[str]:
    return [c for c in consequences if c in DAMAGING_LOF]


def _gene_defect_status(
    gene: str,
    variants: dict[str, list[str]],
    expression: dict[str, float],
    *,
    expression_assessed: bool,
    low_tpm: float,
    cnv_homdel: bool,
) -> tuple[str, str, list[str]]:
    """Return (loss_status, mechanism, evidence_parts)."""
    cons = variants.get(gene, [])
    lof = _damaging_lof_hits(cons)
    evidence: list[str] = []
    if len(lof) >= 2:
        evidence.append(f"biallelic_lof:{','.join(lof[:2])}")
        return "biallelic_loss", "biallelic_lof", evidence
    if lof:
        evidence.append(f"lof:{lof[0]}")
        if cnv_homdel:
            evidence.append("cnv_homdel_support")
            return "biallelic_loss", "lof_plus_cnv_homdel", evidence
        return "monoallelic_loss", "monoallelic_lof", evidence
    if expression_assessed:
        tpm = expression.get(gene)
        if tpm is not None and tpm < low_tpm:
            evidence.append(f"low_expression:{tpm:.4f}")
            return "expression_defect", "low_expression", evidence
    return "intact", "none", evidence


def _pathway_for_gene(gene: str) -> str:
    if gene in HLA_CLASS_I:
        return "HLA-I"
    if gene in APM_GENES:
        return "APM"
    if gene in IFNG_GENES:
        return "IFNG"
    if gene in MHC_II_GENES:
        return "MHC-II"
    if gene == "B2M":
        return "MHC-I"
    return "other"


def _risk_for_gene_defect(gene: str, loss_status: str, mechanism: str) -> str:
    if loss_status == "biallelic_loss":
        if gene == "B2M" or gene in IFNG_GENES[:3]:
            return "HIGH"
        if gene in APM_GENES or gene in HLA_CLASS_I:
            return "HIGH"
        if gene == "CIITA":
            return "HIGH"
        return "MEDIUM"
    if loss_status == "monoallelic_loss":
        if gene in ("B2M", "JAK1", "JAK2", "CIITA"):
            return "MEDIUM"
        if gene in APM_GENES:
            return "MEDIUM"
        return "LOW"
    if loss_status == "expression_defect":
        if gene in ("NLRC5", "CIITA"):
            return "MEDIUM"
        return "LOW"
    return "INCONCLUSIVE"


def build_immune_escape_events(
    sample_id: str,
    variants: dict[str, list[str]],
    expression: dict[str, float],
    lost_hla: set[str],
    hla_typing: list[str],
    *,
    expression_assessed: bool,
    cnv_homdel: bool,
    profile: Mapping[str, Any],
) -> list[dict[str, str]]:
    cfg = profile.get("immune_escape", {}) or {}
    low_tpm = float(profile.get("appm", {}).get("low_expression_tpm", 1.0))
    rows: list[dict[str, str]] = []

    for gene in sorted(set(MHC_I_GENES + MHC_II_GENES)):
        loss_status, mechanism, evidence = _gene_defect_status(
            gene, variants, expression,
            expression_assessed=expression_assessed,
            low_tpm=low_tpm,
            cnv_homdel=cnv_homdel,
        )

        if loss_status == "intact" and not evidence:
            continue

        risk = _risk_for_gene_defect(gene, loss_status, mechanism)
        if loss_status == "hla_loh":
            risk = _max_risk(risk, "MEDIUM")

        action = "none"
        if gene == "B2M" and loss_status == "biallelic_loss":
            action = "reject_mhc_i"
        elif gene in IFNG_GENES[:3] and loss_status == "biallelic_loss":
            action = "ifng_escape"
        elif gene in APM_GENES and loss_status in {"biallelic_loss", "monoallelic_loss", "expression_defect"}:
            action = "apm_penalty"
        elif gene == "CIITA" and loss_status in {"biallelic_loss", "monoallelic_loss", "expression_defect"}:
            action = "mhc_ii_penalty"

        rows.append({
            "event_id": safe_id(f"{sample_id}_{gene}_{mechanism}"),
            "sample_id": sample_id,
            "gene": gene,
            "pathway": _pathway_for_gene(gene),
            "mechanism": mechanism,
            "gene_status": loss_status,
            "evidence": ";".join(evidence) if evidence else "none",
            "resistance_risk": risk,
            "peptide_action": action,
        })

    if lost_hla and not any(r["mechanism"] == "hla_allele_loh" for r in rows):
        rows.append({
            "event_id": safe_id(f"{sample_id}_HLA_LOH"),
            "sample_id": sample_id,
            "gene": "HLA-I",
            "pathway": "HLA-I",
            "mechanism": "hla_allele_loh",
            "gene_status": "hla_loh",
            "evidence": "lost_alleles:" + ",".join(sorted(lost_hla)),
            "resistance_risk": "MEDIUM",
            "peptide_action": "cap_restricting_hla",
        })

    if hla_typing and not variants and not lost_hla and not expression_assessed:
        rows.append({
            "event_id": safe_id(f"{sample_id}_INCONCLUSIVE"),
            "sample_id": sample_id,
            "gene": "NA",
            "pathway": "NA",
            "mechanism": "insufficient_evidence",
            "gene_status": "inconclusive",
            "evidence": "hla_typing_only",
            "resistance_risk": "INCONCLUSIVE",
            "peptide_action": "none",
        })

    return rows


def build_immune_escape_summary(
    sample_id: str,
    events: list[dict[str, str]],
    lost_hla: set[str],
    hla_typing: list[str],
    *,
    expression_assessed: bool,
    profile: Mapping[str, Any],
) -> dict[str, str]:
    cfg = profile.get("immune_escape", {}) or {}
    penalties = profile.get("appm_penalty", {})

    b2m_events = [e for e in events if e.get("gene") == "B2M"]
    jak_events = [e for e in events if e.get("gene") in {"JAK1", "JAK2"}]
    apm_events = [e for e in events if e.get("gene") in APM_GENES]
    ciita_events = [e for e in events if e.get("gene") == "CIITA"]

    b2m_status = b2m_events[0]["gene_status"] if b2m_events else "intact"
    jak_status = "biallelic_loss" if any(e["gene_status"] == "biallelic_loss" for e in jak_events) else (
        "monoallelic_loss" if jak_events else "intact"
    )
    apm_status = "defect" if apm_events else "intact"
    ciita_status = ciita_events[0]["gene_status"] if ciita_events else "intact"

    mhc_i_risk = _max_risk(
        *(e["resistance_risk"] for e in events if e.get("pathway") in {"MHC-I", "HLA-I", "APM", "IFNG"}),
        "MEDIUM" if lost_hla else "INCONCLUSIVE",
    )
    mhc_ii_risk = _max_risk(*(e["resistance_risk"] for e in events if e.get("pathway") == "MHC-II"), "INCONCLUSIVE")
    ifng_risk = _max_risk(*(e["resistance_risk"] for e in events if e.get("pathway") == "IFNG"), "INCONCLUSIVE")
    apm_risk = _max_risk(*(e["resistance_risk"] for e in events if e.get("pathway") == "APM"), "INCONCLUSIVE")

    overall = _max_risk(mhc_i_risk, mhc_ii_risk, ifng_risk, apm_risk)

    mhc_i_score = 1.0
    if b2m_status == "biallelic_loss":
        mhc_i_score -= float(penalties.get("b2m_damaging", 0.60))
    elif b2m_status == "monoallelic_loss":
        mhc_i_score -= float(penalties.get("b2m_damaging", 0.60)) * 0.5
    if apm_status == "defect":
        mhc_i_score -= float(penalties.get("tap1_tap2_damaging", 0.35))
    if jak_status == "biallelic_loss":
        mhc_i_score -= float(penalties.get("ifng_pathway_damaging", 0.25))
    elif jak_status == "monoallelic_loss":
        mhc_i_score -= float(penalties.get("ifng_pathway_damaging", 0.25)) * 0.5
    if lost_hla:
        mhc_i_score -= float(penalties.get("hla_allele_loh", 0.50)) * min(1.0, len(lost_hla) * 0.25)

    mhc_ii_score = 1.0
    if ciita_status in {"biallelic_loss", "monoallelic_loss", "expression_defect"}:
        mhc_ii_score -= float(penalties.get("ciita_low_expression", 0.35))

    if not events and not lost_hla and not expression_assessed:
        overall = "INCONCLUSIVE"

    return {
        "sample_id": sample_id,
        "resistance_risk": overall,
        "mhc_i_resistance_risk": mhc_i_risk,
        "mhc_ii_resistance_risk": mhc_ii_risk,
        "ifng_resistance_risk": ifng_risk,
        "apm_resistance_risk": apm_risk,
        "b2m_gene_status": b2m_status,
        "jak_pathway_status": jak_status,
        "tap_apm_status": apm_status,
        "ciita_status": ciita_status,
        "hla_loh_alleles": ",".join(sorted(lost_hla)),
        "hla_typing_alleles": ",".join(hla_typing),
        "mhc_i_integrity_score": f"{clamp(mhc_i_score):.4f}",
        "mhc_ii_integrity_score": f"{clamp(mhc_ii_score):.4f}",
        "expression_assessment_status": "assessed" if expression_assessed else "unassessed",
        "evidence_event_count": str(len(events)),
    }


def build_peptide_escape_flags(
    sample_id: str,
    peptides: list[dict[str, str]],
    summary: dict[str, str],
    lost_hla: set[str],
    profile: Mapping[str, Any],
) -> list[dict[str, str]]:
    cfg = profile.get("immune_escape", {}) or {}
    b2m_reject = bool(cfg.get("b2m_biallelic_reject", True))
    b2m_cap = float(cfg.get("b2m_biallelic_multiplier", 0.25))
    hla_loh_cap = float(cfg.get("hla_loh_multiplier", profile.get("appm_penalty", {}).get("hla_allele_loh", 0.50)))
    hla_loh_mult = clamp(1.0 - hla_loh_cap)
    mhc_ii_cap = float(cfg.get("mhc_ii_defect_multiplier", 0.50))
    ifng_cap = float(cfg.get("ifng_defect_multiplier", 0.65))

    b2m_status = summary.get("b2m_gene_status", "intact")
    jak_status = summary.get("jak_pathway_status", "intact")
    ciita_status = summary.get("ciita_status", "intact")
    apm_status = summary.get("tap_apm_status", "intact")
    sample_risk = summary.get("resistance_risk", "INCONCLUSIVE")

    rows: list[dict[str, str]] = []
    for pep in peptides:
        pid = pep.get("peptide_id", "")
        hla = normalize_hla_allele(pep.get("hla_allele", ""))
        mhc = str(pep.get("mhc_class") or "").upper()
        is_mhc_ii = mhc in {"II", "MHC-II", "CLASSII"}

        reasons: list[str] = []
        risk = "INCONCLUSIVE"
        action = "none"
        mult = 1.0

        if not is_mhc_ii and b2m_status == "biallelic_loss":
            reasons.append("b2m_biallelic_loss")
            risk = "HIGH"
            if b2m_reject:
                action = "reject"
                mult = b2m_cap
            else:
                action = "cap"
                mult = min(mult, b2m_cap)
        elif not is_mhc_ii and hla and hla in lost_hla:
            reasons.append("hla_allele_loh")
            risk = _max_risk(risk, "MEDIUM")
            action = "cap"
            mult = min(mult, hla_loh_mult)
        elif not is_mhc_ii and apm_status == "defect":
            reasons.append("apm_defect")
            risk = _max_risk(risk, "MEDIUM")
            action = "penalty"
            mult = min(mult, float(cfg.get("apm_defect_multiplier", 0.70)))

        if jak_status == "biallelic_loss":
            reasons.append("jak_biallelic_lof")
            risk = _max_risk(risk, "HIGH")
            action = action if action != "none" else "cap"
            mult = min(mult, ifng_cap)

        if is_mhc_ii and ciita_status in {"biallelic_loss", "monoallelic_loss", "expression_defect"}:
            reasons.append("ciita_defect")
            risk = _max_risk(risk, "MEDIUM")
            action = "cap"
            mult = min(mult, mhc_ii_cap)

        if not reasons:
            risk = sample_risk if sample_risk != "INCONCLUSIVE" else "LOW"

        rows.append({
            "peptide_id": pid,
            "sample_id": sample_id,
            "peptide": pep.get("peptide", ""),
            "hla_allele": pep.get("hla_allele", ""),
            "mhc_class": pep.get("mhc_class", ""),
            "escape_flag": "yes" if reasons else "no",
            "escape_reason": ";".join(reasons) if reasons else "none",
            "resistance_risk": _norm_risk(risk),
            "escape_action": action,
            "escape_multiplier": f"{clamp(mult):.4f}",
        })
    return rows


def build_immune_escape_resistance(
    sample_id: str,
    raw_peptides: str | Path,
    profile: Mapping[str, Any],
    outdir: str | Path,
    *,
    vep_tsv: str | Path | None = None,
    expression_tsv: str | Path | None = None,
    hla_loh_tsv: str | Path | None = None,
    cnv_tsv: str | Path | None = None,
    hla_typing_tsv: str | Path | None = None,
    hla_alleles: list[str] | None = None,
    provenance: ProvenanceRecord | None = None,
) -> tuple[list[dict[str, str]], dict[str, str], list[dict[str, str]]]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cfg = profile.get("immune_escape", {}) or {}
    homdel_threshold = float(cfg.get("homdel_cn_threshold", 0.5))

    variants = load_variants(vep_tsv)
    expression = load_expression(expression_tsv)
    expression_assessed = bool(expression_tsv and Path(expression_tsv).exists() and expression)
    lost_hla = load_loh_alleles(hla_loh_tsv)
    hla_typing = load_hla_typing(hla_typing_tsv, hla_alleles)
    cnv_homdel = load_cnv_homdel_flag(cnv_tsv, homdel_threshold=homdel_threshold)

    events = build_immune_escape_events(
        sample_id, variants, expression, lost_hla, hla_typing,
        expression_assessed=expression_assessed,
        cnv_homdel=cnv_homdel,
        profile=profile,
    )
    summary = build_immune_escape_summary(
        sample_id, events, lost_hla, hla_typing,
        expression_assessed=expression_assessed,
        profile=profile,
    )
    peptides = read_tsv(raw_peptides)
    pep_flags = build_peptide_escape_flags(sample_id, peptides, summary, lost_hla, profile)

    upstream = (
        f"vep:{vep_tsv};expression:{expression_tsv};hla_loh:{hla_loh_tsv};"
        f"cnv:{cnv_tsv};hla_typing:{hla_typing_tsv}"
    )
    prov = provenance or provenance_derived("immune_escape", outdir, upstream=upstream)

    write_tsv(outdir / "immune_escape_events.tsv", attach_provenance(events, prov), IMMUNE_ESCAPE_EVENT_FIELDS)
    write_tsv(outdir / "immune_escape_summary.tsv", attach_provenance([summary], prov), IMMUNE_ESCAPE_SUMMARY_FIELDS)
    write_tsv(outdir / "peptide_escape_flags.tsv", attach_provenance(pep_flags, prov), PEPTIDE_ESCAPE_FLAG_FIELDS)

    return events, summary, pep_flags
