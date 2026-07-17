from __future__ import annotations
from pathlib import Path
from typing import Any, Mapping

from .utils import read_tsv, write_tsv, to_float, norm_tpm, clamp
from .config import source_priority
from .safety import safety_multiplier, apply_event_safety, apply_peptide_safety, load_normal_expression, load_normal_hla_ligands
from .gates import evaluate_presentation_gate
from .schemas import EVENT_FIELDS, PEPTIDE_FIELDS
from .model_layers import enrich_event_layers, enrich_peptide_layers, compute_l3_dimension_scores

from .immunogenicity_composite import has_resolved_immunogenicity, resolve_immunogenicity_score
from .peptide_safety_gate import load_peptide_safety
from .immune_escape import load_peptide_escape_flags


def map_by(rows, key):
    return {r.get(key,""): r for r in rows if r.get(key,"")}



def load_optional_map(path, key):
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    return map_by(read_tsv(p), key)


def resolve_appm_peptide_modifiers_tsv(
    appm_peptide_modifiers_tsv: str | Path | None = None,
    *,
    appm_summary_tsv: str | Path | None = None,
    appm_dir: str | Path | None = None,
) -> str | None:
    """Resolve APPM 2.0 peptide modifier sidecar next to appm_summary or appm outdir."""
    candidates: list[Path] = []
    if appm_peptide_modifiers_tsv:
        candidates.append(Path(appm_peptide_modifiers_tsv))
    if appm_summary_tsv:
        candidates.append(Path(appm_summary_tsv).parent / "appm_peptide_modifiers.tsv")
    if appm_dir:
        candidates.append(Path(appm_dir) / "appm_peptide_modifiers.tsv")
    for path in candidates:
        if path.is_file():
            return str(path)
    return None

_PRIORITY_ORDER = {"D":0, "C":1, "C_CAUTION":1, "B_CAUTION":2, "B":3, "A":4}

def apply_priority_cap_value(priority_value: str, cap: str) -> str:
    cap = str(cap or "").strip().upper()
    if not cap:
        return priority_value
    pri = str(priority_value or "D").strip().upper()
    if cap in {"A", "NONE"}:
        return priority_value
    if _PRIORITY_ORDER.get(pri, 0) <= _PRIORITY_ORDER.get(cap, 0):
        return priority_value
    return cap

def merge_appm_modifier(p: dict[str, Any], appm_row: Mapping[str, Any] | None) -> dict[str, Any]:
    if not appm_row:
        return p
    for k in [
        "appm_multiplier", "appm_multiplier_reason", "appm_integrity_status",
        "appm_evidence_completeness", "appm_review_required", "appm_action", "priority_cap",
    ]:
        val = appm_row.get(k)
        if val in {None, ""}:
            continue
        if k == "priority_cap" and p.get("priority_cap"):
            p["priority_cap"] = apply_priority_cap_value(p.get("priority_cap"), str(val))
        else:
            p[k] = val
    return p

def merge_safety_gate(peptide: dict, gate: Mapping[str, Any] | None) -> dict:
    if not gate:
        return peptide
    for k in ["safety_tier", "safety_status", "safety_reason", "safety_multiplier"]:
        if gate.get(k) not in {None, ""}:
            peptide[k] = gate.get(k)
    return peptide

def merge_escape_flags(peptide: dict, flags: Mapping[str, Any] | None) -> dict:
    if not flags:
        peptide.setdefault("escape_multiplier", "1.0000")
        return peptide
    for k in ["escape_status", "escape_severity", "escape_flag", "escape_reason", "resistance_risk", "escape_action", "escape_multiplier", "priority_cap"]:
        val = flags.get(k)
        if val not in {None, ""}:
            if k == "priority_cap" and peptide.get("priority_cap"):
                # keep the stricter cap later via apply_priority_cap_value
                peptide["priority_cap"] = apply_priority_cap_value(peptide.get("priority_cap"), str(val))
            else:
                peptide[k] = val
    # `escape_severity` is a normalized 3-tier field ("ESCAPE_REJECT" /
    # "ESCAPE_CAUTION" / "ESCAPE_PASS") derived directly from the escape
    # multiplier in immune_escape.py, unlike the free-form `escape_status`
    # string (e.g. "GLOBAL_MHC_I_LOSS") this used to (incorrectly) compare
    # against literal "ESCAPE_REJECT"/"ESCAPE_CAUTION" values that
    # immune_escape.py never actually produced (scoring audit fix #5).
    if flags.get("escape_severity") == "ESCAPE_REJECT":
        peptide["safety_status"] = "FAIL"
        peptide["safety_reason"] = ";".join(x for x in [peptide.get("safety_reason"), flags.get("escape_reason")] if x)
    elif flags.get("escape_severity") == "ESCAPE_CAUTION" and peptide.get("safety_status") != "FAIL":
        peptide["safety_status"] = "CAUTION"
        peptide["safety_reason"] = ";".join(x for x in [peptide.get("safety_reason"), flags.get("escape_reason")] if x)
    return peptide

def is_placeholder_immunogenicity(peptide: Mapping[str, Any]) -> bool:
    raw = str(peptide.get("immunogenicity_score", "")).strip()
    if not raw:
        return True
    val = to_float(raw, 0.5)
    return abs(val - 0.5) < 1e-6 or val <= 0.0


def normalize_immunogenicity_score(
    peptide: Mapping[str, Any],
    presentation: Mapping[str, Any] | None = None,
    profile: Mapping[str, Any] | None = None,
) -> float:
    """Map PRIME/BigMHC_IM composite (or legacy pVAC/IEDB) to a 0–1 score."""
    if profile is not None:
        composite, _ = resolve_immunogenicity_score(peptide, presentation, profile)
        if composite > 0:
            return composite
    pres = presentation or {}
    raw = peptide.get("immunogenicity_score") or pres.get("immunogenicity_score", "")
    if str(raw).strip() not in {"", "0.5", "0.5000"}:
        val = to_float(raw, -1.0)
        if val >= 0.0:
            return clamp(val)
    return 0.0


def effective_v03_weights(profile: Mapping[str, Any], peptide: Mapping[str, Any], presentation: Mapping[str, Any] | None = None) -> dict[str, float]:
    """Redistribute immunogenicity weight when only the pVAC placeholder is available."""
    w = {k: float(v) for k, v in profile.get("v03_score_weights", {}).items()}
    imm_w = float(w.get("immunogenicity", 0.15))
    has_real = has_resolved_immunogenicity(peptide, presentation, profile)

    if not has_real:
        pres_w = float(w.get("presentation_evidence", 0.25))
        bind_w = float(w.get("binding_evidence", 0.20))
        denom = pres_w + bind_w
        if denom > 0 and imm_w > 0:
            w["presentation_evidence"] = pres_w + imm_w * (pres_w / denom)
            w["binding_evidence"] = bind_w + imm_w * (bind_w / denom)
        w["immunogenicity"] = 0.0
    return w



def _priority_rank(label: str) -> int:
    order = {"D": 0, "C": 1, "C_CAUTION": 2, "B_CAUTION": 3, "B": 4, "A": 5}
    return order.get(str(label or ""), 5)

def apply_priority_cap(priority_label: str, cap: str | None) -> str:
    cap = str(cap or "").strip()
    if not cap:
        return priority_label
    order = ["D", "C", "C_CAUTION", "B_CAUTION", "B", "A"]
    if cap not in order:
        return priority_label
    return priority_label if _priority_rank(priority_label) <= _priority_rank(cap) else cap

def merge_peptide_safety(p: dict[str, Any], safety_row: Mapping[str, Any] | None) -> dict[str, Any]:
    if not safety_row:
        return p
    for k in [
        "safety_tier", "safety_status", "safety_reason", "safety_multiplier", "review_required",
        "reference_proteome_exact_match", "normal_ligand_tissue", "mutation_anchor_only",
    ]:
        if k in safety_row and str(safety_row.get(k, "")).strip() != "":
            p[k] = safety_row.get(k, "")
    return p

def merge_peptide_escape(p: dict[str, Any], escape_row: Mapping[str, Any] | None) -> dict[str, Any]:
    if not escape_row:
        return p
    for k in ["escape_status", "escape_reason", "escape_multiplier", "restricting_hla_lost", "priority_cap"]:
        if k == "restricting_hla_lost":
            p[k] = escape_row.get("restricting_hla_lost", "")
        elif k == "priority_cap":
            p[k] = escape_row.get("priority_cap") or p.get("priority_cap", "")
        elif k in escape_row:
            p[k] = escape_row.get(k, "")
    return p

def score_event(e, profile):
    w = profile.get("event_weights", {})
    base = (
        float(w.get("event_confidence",0.2))*clamp(to_float(e.get("event_confidence"),0)) +
        float(w.get("event_expression",0.2))*norm_tpm(e.get("event_expression"), float(profile.get("scoring",{}).get("high_expression_tpm",20))) +
        float(w.get("driver_relevance",0.2))*clamp(to_float(e.get("driver_relevance"),0)) +
        float(w.get("clonality",0.15))*clamp(to_float(e.get("clonality"),0)) +
        float(w.get("persistence",0.15))*clamp(to_float(e.get("persistence"),0)) +
        float(w.get("tumor_specificity",0.1))*clamp(to_float(e.get("tumor_specificity"),0))
    )
    base *= source_priority(profile, e.get("event_type", "Other"), e.get("mutation_source"))
    base *= safety_multiplier(e.get("safety_status",""))
    base *= to_float(e.get("clonality_multiplier"), 1.0)
    e["event_score"] = f"{clamp(base,0,2):.4f}"
    return e

def appm_multiplier(p, summary, profile=None):
    if str(p.get("appm_multiplier", "")).strip():
        return clamp(to_float(p.get("appm_multiplier"), 1.0))
    if (p.get("mhc_class") or "").upper() in {"II","MHC-II","CLASSII"}:
        base = to_float(summary.get("mhc_ii_integrity_score"), 1.0)
    else:
        base = to_float(summary.get("mhc_i_integrity_score"), 1.0)
    loh_str = str(summary.get("hla_loh_alleles", "")).strip()
    if loh_str and profile:
        from .adapters.peptide_input import normalize_hla_allele
        loh_set = {normalize_hla_allele(x) for x in loh_str.split(",") if x.strip()}
        hla = normalize_hla_allele(p.get("hla_allele", ""))
        if hla in loh_set:
            pen = float(profile.get("appm_penalty", {}).get("hla_allele_loh", 0.50))
            base = clamp(base - pen)
    return clamp(base)

def priority(safety, score):
    if safety == "FAIL": return "D"
    if safety == "CAUTION": return "B_CAUTION" if score >= 0.55 else "C_CAUTION"
    if score >= 0.75: return "A"
    if score >= 0.55: return "B"
    if score >= 0.35: return "C"
    return "D"

def recommended(p, appm, ccf):
    if p.get("safety_status") == "FAIL":
        return "Do not advance; safety gate failed"
    notes = []
    if p.get("safety_status") == "CAUTION": notes.append("requires focused safety validation")
    if p.get("presentation_gate_status") == "FAIL":
        notes.append(f"presentation gate: {p.get('presentation_gate_reason','')}")
    if appm < 0.65: notes.append("APPM caution")
    if ccf < 0.75: notes.append("clonality/persistence caution")
    from .validation_design import classify_validation_mode, recommended_assay_text
    mode = classify_validation_mode(p)
    if mode not in {"do_not_advance", "safety_caution"}:
        notes.append(recommended_assay_text(mode))
    return "; ".join(notes)


def compute_peptide_efficacy(
    peptide: Mapping[str, Any],
    event: Mapping[str, Any],
    presentation: Mapping[str, Any],
    profile: Mapping[str, Any],
    *,
    appm: float = 1.0,
    ccf: float = 1.0,
) -> dict[str, Any]:
    """Shared peptide scoring used by v0.3 pipeline and IMPROVE benchmark."""
    gate = evaluate_presentation_gate(peptide, event, presentation, profile)
    gate_mult = to_float(gate["presentation_gate_multiplier"], 1.0)
    w = effective_v03_weights(profile, peptide, presentation)
    immuno = normalize_immunogenicity_score(peptide, presentation, profile)
    if w.get("immunogenicity", 0.0) <= 0:
        immuno = 0.0
    _, immuno_meta = resolve_immunogenicity_score(peptide, presentation, profile)

    high_tpm = float(profile.get("scoring", {}).get("high_expression_tpm", 20))
    l3 = compute_l3_dimension_scores(
        peptide, event, presentation, profile,
        appm=appm, ccf=ccf, immuno=immuno, high_expression_tpm=high_tpm,
    )
    composite = to_float(l3.get("immunology_composite_score"), 0.0)
    safety = peptide.get("safety_status", event.get("safety_status", "PASS"))
    score = composite * safety_multiplier(safety) * appm * gate_mult * ccf
    score *= to_float(peptide.get("escape_multiplier"), 1.0)
    return {
        "efficacy_score": f"{clamp(score):.4f}",
        "immunogenicity_resolved": "yes" if w.get("immunogenicity", 0.0) > 0 else "no",
        **immuno_meta,
        **gate,
        **l3,
    }


def score_peptide(p, e, profile, pres, summary):
    for k in [
        "netmhcpan_ba_rank", "netmhcpan_el_rank", "netmhcstabpan_score", "netmhcstabpan_rank",
        "mhcflurry_affinity_percentile", "mhcflurry_processing_score", "mhcflurry_presentation_score",
        "binding_evidence_score", "presentation_evidence_score", "presentation_evidence_grade",
        "prime_score", "prime_rank", "bigmhc_im_score", "deepimmuno_score",
        "iedb_immunogenicity_score", "immunogenicity_composite_score", "immunogenicity_source",
    ]:
        p[k] = pres.get(k, p.get(k, ""))

    appm = appm_multiplier(p, summary, profile)
    ccf = to_float(e.get("clonality_multiplier"), 1.0)
    scored = compute_peptide_efficacy(p, e, pres, profile, appm=appm, ccf=ccf)
    p.update(scored)
    escape_mult = to_float(p.get("escape_multiplier"), 1.0)
    # escape_multiplier is already applied inside compute_peptide_efficacy;
    # do not apply it a second time here.
    p["appm_multiplier"] = f"{appm:.4f}"
    p.setdefault("appm_multiplier_reason", "sample_level_appm_summary")
    p.setdefault("appm_integrity_status", summary.get("mhc_ii_integrity_status" if (p.get("mhc_class") or "").upper() in {"II", "MHC-II", "CLASSII"} else "mhc_i_integrity_status", ""))
    p.setdefault("appm_evidence_completeness", summary.get("appm_evidence_completeness", ""))
    p.setdefault("appm_review_required", "yes" if summary.get("appm_evidence_completeness") in {"LOW", "UNASSESSED"} else "no")
    p["ccf_multiplier"] = f"{ccf:.4f}"
    p["escape_multiplier"] = f"{escape_mult:.4f}"
    pri = priority(p.get("safety_status"), to_float(p["efficacy_score"], 0.0))
    p["final_priority"] = apply_priority_cap_value(pri, p.get("priority_cap", ""))
    p["recommended_use"] = recommended(p, appm, ccf)
    return p

def filter_by_enabled_sources(events: list[dict], profile: Mapping[str, Any]) -> list[dict]:
    if not profile.get("enforce_enabled_sources", False):
        return events
    enabled = profile.get("enabled_sources", {}).get("types", [])
    if not enabled:
        return events
    allowed = {str(x) for x in enabled}
    legacy = {str(x) for x in enabled}
    out = []
    for e in events:
        ms = str(e.get("mutation_source") or "").strip()
        et = str(e.get("event_type") or "").strip()
        if ms in allowed or et in legacy:
            out.append(e)
    return out


def score_v03(raw_events, raw_peptides, presentation_evidence, appm_summary_tsv, ccf_lite_tsv, normal_expression_tsv, normal_hla_ligands_tsv, profile, out_events, out_peptides, peptide_safety_tsv=None, peptide_escape_flags_tsv=None, appm_peptide_modifiers_tsv=None):
    events = [enrich_event_layers(e) for e in read_tsv(raw_events)]
    events = filter_by_enabled_sources(events, profile)
    peptides = read_tsv(raw_peptides)
    pres_map = map_by(read_tsv(presentation_evidence), "peptide_id")
    ccf_map = map_by(read_tsv(ccf_lite_tsv), "event_id") if ccf_lite_tsv else {}
    appm_rows = read_tsv(appm_summary_tsv) if appm_summary_tsv and Path(appm_summary_tsv).exists() else []
    summary = appm_rows[0] if appm_rows else {"mhc_i_integrity_score":"1.0", "mhc_ii_integrity_score":"1.0"}
    norm_expr = load_normal_expression(normal_expression_tsv)
    norm_lig = load_normal_hla_ligands(normal_hla_ligands_tsv)
    safety_map = load_optional_map(peptide_safety_tsv, "peptide_id")
    escape_map = load_optional_map(peptide_escape_flags_tsv, "peptide_id")
    appm_modifier_map = load_optional_map(appm_peptide_modifiers_tsv, "peptide_id")
    event_map, scored_e = {}, []
    for e in events:
        if e.get("event_id") in ccf_map:
            c = ccf_map[e["event_id"]]
            e["ccf_estimate"] = c.get("ccf_estimate","")
            e["ccf_status"] = c.get("ccf_status","")
            e["clonality_multiplier"] = c.get("clonality_multiplier","1.0")
        else:
            e["clonality_multiplier"] = e.get("clonality_multiplier","1.0")
        e["appm_mhc_i_integrity"] = summary.get("mhc_i_integrity_score","1.0")
        e["appm_mhc_ii_integrity"] = summary.get("mhc_ii_integrity_score","1.0")
        e = apply_event_safety(e, profile, norm_expr)
        e = score_event(e, profile)
        event_map[e["event_id"]] = e
        scored_e.append(e)
    scored_p = []
    for p in peptides:
        e = event_map.get(p.get("event_id"))
        if not e:
            continue
        p = enrich_peptide_layers(p, e)
        p = apply_peptide_safety(p, e, profile, norm_lig)
        p = merge_safety_gate(p, safety_map.get(p.get("peptide_id", "")))
        p = merge_appm_modifier(p, appm_modifier_map.get(p.get("peptide_id", "")))
        p = merge_escape_flags(p, escape_map.get(p.get("peptide_id", "")))
        p = score_peptide(p, e, profile, pres_map.get(p.get("peptide_id",""), {}), summary)
        scored_p.append(p)
    scored_e.sort(key=lambda r: to_float(r.get("event_score"),0), reverse=True)
    scored_p.sort(key=lambda r: to_float(r.get("efficacy_score"),0), reverse=True)
    write_tsv(out_events, scored_e, EVENT_FIELDS)
    write_tsv(out_peptides, scored_p, PEPTIDE_FIELDS)
    return scored_e, scored_p
