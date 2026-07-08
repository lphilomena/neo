from __future__ import annotations
from .utils import read_tsv, to_float, truthy

def load_normal_expression(path):
    if not path: return {}
    d = {}
    for r in read_tsv(path):
        gene = (r.get("gene") or r.get("Gene") or "").strip()
        if gene:
            d[gene] = {
                "normal_tissue_max_tpm": to_float(r.get("normal_tissue_max_tpm"), 0.0),
                "normal_hspc_tpm": to_float(r.get("normal_hspc_tpm"), 0.0),
                "critical_tissue_hit": 1.0 if truthy(r.get("critical_tissue_hit")) else 0.0,
            }
    return d

def load_normal_hla_ligands(path):
    if not path: return set()
    return { (r.get("peptide") or r.get("Peptide") or "").strip().upper() for r in read_tsv(path) if (r.get("peptide") or r.get("Peptide") or "").strip() }

def _split_event_genes(gene_value):
    raw = str(gene_value or "")
    for sep in ("::", ";", ",", "|"):
        raw = raw.replace(sep, " ")
    return [g.strip() for g in raw.split() if g.strip()]

def _normal_expr_for_event_gene(gene_value, normal_expr):
    genes = _split_event_genes(gene_value)
    if not genes:
        return {}
    nt = 0.0; nh = 0.0; crit = 0.0
    found = False
    for g in genes:
        row = normal_expr.get(g, {})
        if row:
            found = True
            nt = max(nt, to_float(row.get("normal_tissue_max_tpm"), 0.0))
            nh = max(nh, to_float(row.get("normal_hspc_tpm"), 0.0))
            crit = max(crit, to_float(row.get("critical_tissue_hit"), 0.0))
    return {"normal_tissue_max_tpm": nt, "normal_hspc_tpm": nh, "critical_tissue_hit": crit} if found else {}

def apply_event_safety(e, profile, normal_expr):
    safety = profile.get("safety", {})
    expr = _normal_expr_for_event_gene(e.get("gene",""), normal_expr)
    nt = to_float(e.get("normal_tissue_max_tpm"), expr.get("normal_tissue_max_tpm", 0.0))
    nh = to_float(e.get("normal_hspc_tpm"), expr.get("normal_hspc_tpm", 0.0))
    crit = bool(to_float(e.get("critical_tissue_hit"), expr.get("critical_tissue_hit", 0.0)))
    status, reasons = "PASS", []
    if crit and nt >= float(safety.get("hard_exclusion_critical_tissue_tpm", 10.0)):
        status = "FAIL"; reasons.append("critical_tissue_expression")
    elif crit and nt >= float(safety.get("caution_critical_tissue_tpm", 1.0)):
        status = "CAUTION"; reasons.append("critical_tissue_expression")
    if nh >= float(safety.get("normal_hspc_fail_tpm", 8.0)):
        status = "FAIL"; reasons.append("normal_HSPC_expression")
    elif nh >= float(safety.get("normal_hspc_caution_tpm", 2.0)) and status != "FAIL":
        status = "CAUTION"; reasons.append("normal_HSPC_expression")
    e.update({
        "normal_tissue_max_tpm": f"{nt:.4f}", "normal_hspc_tpm": f"{nh:.4f}",
        "critical_tissue_hit": "yes" if crit else "no",
        "safety_status": status, "safety_reason": ";".join(reasons) if reasons else "no_major_signal"
    })
    return e

def apply_peptide_safety(p, e, profile, normal_ligands):
    safety = profile.get("safety", {})
    status = e.get("safety_status", "PASS")
    reasons = [] if e.get("safety_reason") in {"", "no_major_signal"} else [e.get("safety_reason")]
    if p.get("peptide", "").upper() in normal_ligands:
        status = "FAIL" if str(safety.get("normal_hla_ligand_overlap_action", "CAUTION")).upper() == "FAIL" else ("CAUTION" if status != "FAIL" else status)
        p["normal_hla_ligand_overlap"] = "yes"; reasons.append("normal_HLA_ligand_overlap")
    if to_float(p.get("wildtype_binding_rank"), 99) <= float(safety.get("wildtype_strong_binding_rank", 0.5)) and to_float(p.get("self_similarity_score"), 0) >= float(safety.get("self_similarity_caution", 0.85)):
        if status != "FAIL": status = "CAUTION"
        reasons.append("wildtype_strong_binding_and_high_self_similarity")
    if p.get("event_type") == "TAA":
        if status != "FAIL": status = "CAUTION"
        reasons.append("TAA_not_true_neoantigen")
    p["safety_status"] = status
    p["safety_reason"] = ";".join([x for x in reasons if x]) if reasons else "no_major_signal"
    return p

def safety_multiplier(status):
    s = (status or "").upper()
    if s == "PASS": return 1.0
    if s == "CAUTION": return 0.45
    if s == "FAIL": return 0.0
    return 0.25
