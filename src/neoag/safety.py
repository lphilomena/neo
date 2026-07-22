from __future__ import annotations
from .utils import read_tsv, to_float, truthy

SAFETY_SEVERITY = {
    "PASS": 0, "SAFETY_PASS": 0,
    "SAFETY_PARTIAL": 1, "PARTIAL": 1,
    "CAUTION": 2, "SAFETY_REVIEW": 2, "SAFETY_HIGH_RISK": 2,
    "FAIL": 3, "SAFETY_REJECT": 3,
}


def worst_safety_status(*statuses):
    normalized = [str(status or "").strip().upper() for status in statuses if str(status or "").strip()]
    if not normalized:
        return "SAFETY_PARTIAL"
    return max(normalized, key=lambda status: SAFETY_SEVERITY.get(status, 1))


def combine_safety_reasons(*reasons):
    parts = []
    for reason in reasons:
        for item in str(reason or "").split(";"):
            item = item.strip()
            if item and item != "no_major_signal" and item not in parts:
                parts.append(item)
    return ";".join(parts) if parts else "no_major_signal"

def load_normal_expression(path):
    if not path: return {}
    d = {}
    for r in read_tsv(path):
        gene = (r.get("gene") or r.get("Gene") or "").strip()
        ensembl_ids = (r.get("ensembl_gene_id") or r.get("gene_id") or "").split(";")
        aliases = [value for value in [gene, *(item.strip() for item in ensembl_ids)] if value]
        if not aliases:
            continue
        tissue_tpm = to_float(r.get("normal_tissue_max_tpm"), 0.0)
        critical_tpm = to_float(r.get("critical_tissue_max_tpm"), tissue_tpm)
        hspc = to_float(r.get("normal_hspc_tpm"), 0.0)
        expression_status = (r.get("normal_expression_status") or "ASSESSED").upper()
        hspc_status = (r.get("normal_hspc_status") or "ASSESSED").upper()
        current = next((d[alias] for alias in aliases if alias in d), None)
        if current is None:
            current = {
                "normal_tissue_max_tpm": 0.0,
                "normal_tissue_max_tissue": "",
                "critical_tissue_max_tpm": 0.0,
                "critical_tissue_name": "",
                "normal_hspc_tpm": 0.0,
                "normal_hspc_unit": (r.get("normal_hspc_unit") or "TPM").strip(),
                "critical_tissue_hit": 0.0,
                "normal_expression_status": "UNASSESSED",
                "normal_hspc_status": "UNASSESSED",
            }
        for alias in aliases:
            d[alias] = current
        if tissue_tpm >= current["normal_tissue_max_tpm"]:
            current["normal_tissue_max_tpm"] = tissue_tpm
            current["normal_tissue_max_tissue"] = r.get("normal_tissue_max_tissue", "")
        if critical_tpm >= current["critical_tissue_max_tpm"]:
            current["critical_tissue_max_tpm"] = critical_tpm
            current["critical_tissue_name"] = r.get("critical_tissue_name", "")
        if hspc >= current["normal_hspc_tpm"]:
            current["normal_hspc_tpm"] = hspc
            current["normal_hspc_unit"] = (r.get("normal_hspc_unit") or "TPM").strip()
        current["critical_tissue_hit"] = max(
            current["critical_tissue_hit"],
            1.0 if truthy(r.get("critical_tissue_hit")) else 0.0,
        )
        if expression_status == "ASSESSED":
            current["normal_expression_status"] = "ASSESSED"
        if hspc_status == "ASSESSED":
            current["normal_hspc_status"] = "ASSESSED"
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
    nt = 0.0; nh = 0.0; crit = 0.0; critical_tpm = 0.0
    nt_tissue = ""; critical_tissue = ""; hspc_unit = ""; found = False
    expr_assessed = True; hspc_assessed = True
    for g in genes:
        row = normal_expr.get(g, {})
        if row:
            found = True
            row_nt = to_float(row.get("normal_tissue_max_tpm"), 0.0)
            if row_nt >= nt:
                nt = row_nt; nt_tissue = row.get("normal_tissue_max_tissue", "")
            row_hspc = to_float(row.get("normal_hspc_tpm"), 0.0)
            if row_hspc >= nh:
                nh = row_hspc; hspc_unit = row.get("normal_hspc_unit", "TPM")
            row_critical = to_float(row.get("critical_tissue_max_tpm"), row_nt)
            if row_critical >= critical_tpm:
                critical_tpm = row_critical; critical_tissue = row.get("critical_tissue_name", "")
            crit = max(crit, to_float(row.get("critical_tissue_hit"), 0.0))
            expr_assessed = expr_assessed and row.get("normal_expression_status", "ASSESSED") == "ASSESSED"
            hspc_assessed = hspc_assessed and row.get("normal_hspc_status", "ASSESSED") == "ASSESSED"
        else:
            expr_assessed = False; hspc_assessed = False
    return {
        "normal_tissue_max_tpm": nt, "normal_tissue_max_tissue": nt_tissue,
        "critical_tissue_max_tpm": critical_tpm, "critical_tissue_name": critical_tissue,
        "normal_hspc_tpm": nh, "normal_hspc_unit": hspc_unit or "TPM", "critical_tissue_hit": crit,
        "normal_expression_status": "ASSESSED" if found and expr_assessed else "UNASSESSED",
        "normal_hspc_status": "ASSESSED" if found and hspc_assessed else "UNASSESSED",
    } if found else {}

def apply_event_safety(e, profile, normal_expr):
    safety = profile.get("safety", {})
    expr = _normal_expr_for_event_gene(e.get("gene",""), normal_expr)
    nt = to_float(e.get("normal_tissue_max_tpm"), expr.get("normal_tissue_max_tpm", 0.0))
    critical_tpm = to_float(e.get("critical_tissue_max_tpm"), expr.get("critical_tissue_max_tpm", nt))
    nh = to_float(e.get("normal_hspc_tpm"), expr.get("normal_hspc_tpm", 0.0))
    crit = bool(to_float(e.get("critical_tissue_hit"), expr.get("critical_tissue_hit", 0.0)))
    status, reasons = "PASS", []
    expr_assessed = bool(expr) and expr.get("normal_expression_status") == "ASSESSED"
    hspc_assessed = bool(expr) and expr.get("normal_hspc_status") == "ASSESSED"
    hspc_unit = str(expr.get("normal_hspc_unit") or "TPM")
    hspc_is_ncpm = "NCPM" in hspc_unit.upper()
    if not expr_assessed:
        status = "SAFETY_PARTIAL"
        reasons.append("normal_expression_unassessed")
    if not hspc_assessed:
        status = worst_safety_status(status, "SAFETY_PARTIAL")
        reasons.append("normal_hspc_unassessed")
    critical_action = str(safety.get("critical_tissue_expression_action", "FAIL")).upper()
    if expr_assessed and crit and critical_tpm >= float(safety.get("hard_exclusion_critical_tissue_tpm", 10.0)):
        status = "FAIL" if critical_action in {"FAIL", "REJECT"} else worst_safety_status(status, "CAUTION")
        reasons.append("critical_tissue_expression")
    elif expr_assessed and crit and critical_tpm >= float(safety.get("caution_critical_tissue_tpm", 1.0)):
        status = worst_safety_status(status, "CAUTION"); reasons.append("critical_tissue_expression")
    fail_key = "normal_hspc_fail_ncpm" if hspc_is_ncpm else "normal_hspc_fail_tpm"
    caution_key = "normal_hspc_caution_ncpm" if hspc_is_ncpm else "normal_hspc_caution_tpm"
    hspc_fail = float(safety.get(fail_key, safety.get("normal_hspc_fail_tpm", 8.0)))
    hspc_caution = float(safety.get(caution_key, safety.get("normal_hspc_caution_tpm", 2.0)))
    hspc_action = str(safety.get("normal_hspc_expression_action", "FAIL")).upper()
    if hspc_assessed and nh >= hspc_fail:
        status = "FAIL" if hspc_action in {"FAIL", "REJECT"} else worst_safety_status(status, "CAUTION")
        reasons.append("normal_HSPC_expression")
    elif hspc_assessed and nh >= hspc_caution and status != "FAIL":
        status = worst_safety_status(status, "CAUTION"); reasons.append("normal_HSPC_expression")
    e.update({
        "normal_tissue_max_tpm": f"{nt:.4f}", "normal_hspc_tpm": f"{nh:.4f}",
        "normal_hspc_unit": hspc_unit,
        "normal_tissue_max_tissue": expr.get("normal_tissue_max_tissue", ""),
        "critical_tissue_max_tpm": f"{critical_tpm:.4f}",
        "critical_tissue_name": expr.get("critical_tissue_name", ""),
        "critical_tissue_hit": "yes" if crit else "no",
        "normal_expression_status": "ASSESSED" if expr_assessed else "UNASSESSED",
        "normal_hspc_status": "ASSESSED" if hspc_assessed else "UNASSESSED",
        "safety_status": status, "safety_reason": combine_safety_reasons(*reasons)
    })
    return e

def apply_peptide_safety(p, e, profile, normal_ligands):
    safety = profile.get("safety", {})
    status = e.get("safety_status", "PASS")
    reasons = [] if e.get("safety_reason") in {"", "no_major_signal"} else [e.get("safety_reason")]
    if p.get("peptide", "").upper() in normal_ligands:
        status = "FAIL" if str(safety.get("normal_hla_ligand_overlap_action", "CAUTION")).upper() == "FAIL" else worst_safety_status(status, "CAUTION")
        p["normal_hla_ligand_overlap"] = "yes"; reasons.append("normal_HLA_ligand_overlap")
    if to_float(p.get("wildtype_binding_rank"), 99) <= float(safety.get("wildtype_strong_binding_rank", 0.5)) and to_float(p.get("self_similarity_score"), 0) >= float(safety.get("self_similarity_caution", 0.85)):
        status = worst_safety_status(status, "CAUTION")
        reasons.append("wildtype_strong_binding_and_high_self_similarity")
    if p.get("event_type") == "TAA":
        status = worst_safety_status(status, "CAUTION")
        reasons.append("TAA_not_true_neoantigen")
    p["safety_status"] = status
    p["safety_reason"] = combine_safety_reasons(*reasons)
    return p

def safety_multiplier(status, profile=None):
    cfg = (profile or {}).get("safety", {})
    s = (status or "").upper()
    if s == "PASS": return 1.0
    if s == "CAUTION": return float(cfg.get("caution_multiplier", 0.45))
    if s in {"SAFETY_PARTIAL", "PARTIAL"}: return float(cfg.get("partial_multiplier", 0.75))
    if s == "FAIL": return 0.0
    return 0.25
