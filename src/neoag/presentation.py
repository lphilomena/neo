from __future__ import annotations
from pathlib import Path
from .utils import read_tsv, write_tsv, to_float, norm_rank, clamp, safe_id
from .schemas import PRESENTATION_FIELDS
from .evidence_provenance import ProvenanceRecord, ProvenanceRegistry, provenance_derived, attach_provenance
from .adapters.peptide_input import normalize_hla_allele

def by_key(rows):
    result = {}
    for row in rows:
        peptide = str(row.get("peptide") or "").strip().upper()
        hla = normalize_hla_allele(str(row.get("hla_allele") or ""))
        # A caller-local peptide_hla_key may refer to an older peptide ID.
        # Exact sequence + normalized HLA is the canonical biological key.
        key = safe_id(f"{peptide}_{hla}") if peptide and hla else str(row.get("peptide_hla_key") or "")
        if key:
            result[key] = row
    return result

def grade(binding, presentation, complete):
    if complete <= 0:
        return "MISSING"
    if binding >= 0.85 and presentation >= 0.85:
        return "A"
    if binding >= 0.75 and presentation >= 0.65:
        return "B"
    if binding >= 0.65:
        return "C_BINDING_ONLY"
    return "D_WEAK"

def build_presentation_evidence(
    raw_peptides,
    netmhcpan,
    mhcflurry,
    profile,
    out=None,
    netmhcstabpan=None,
    netchop=None,
    provenance_registry: ProvenanceRegistry | None = None,
):
    peptides = read_tsv(raw_peptides)
    net = by_key(read_tsv(netmhcpan)) if netmhcpan else {}
    mhc = by_key(read_tsv(mhcflurry)) if mhcflurry else {}
    stab = by_key(read_tsv(netmhcstabpan)) if netmhcstabpan else {}
    chop_rows = read_tsv(netchop) if netchop else []
    chop_by_id = {row.get("peptide_id", ""): row for row in chop_rows if row.get("peptide_id")}
    chop_by_peptide = {row.get("peptide", ""): row for row in chop_rows if row.get("peptide")}
    w = profile.get("presentation_weights", {})
    w_ba = float(w.get("netmhcpan_ba", 0.25))
    w_el = float(w.get("netmhcpan_el", 0.35))
    w_mhcf = float(w.get("mhcflurry_presentation", 0.30))
    w_proc = float(w.get("mhcflurry_processing", 0.10))
    w_chop = float(w.get("netchop_processing", 0.10))
    rows = []
    for p in peptides:
        hla = normalize_hla_allele(str(p.get("hla_allele") or ""))
        key = safe_id(f"{str(p.get('peptide','')).strip().upper()}_{hla}")
        wt_peptide = p.get("wildtype_peptide", "")
        wt_key = safe_id(f"{str(wt_peptide).strip().upper()}_{hla}") if wt_peptide else ""
        n = net.get(key, {})
        m = mhc.get(key, {})
        wt_n = net.get(wt_key, {}) if wt_key else {}
        wt_m = mhc.get(wt_key, {}) if wt_key else {}
        s = stab.get(key, {})
        c = chop_by_id.get(p.get("peptide_id", ""), {}) or chop_by_peptide.get(p.get("peptide", ""), {})
        ba = n.get("netmhcpan_ba_rank", "")
        el = n.get("netmhcpan_el_rank", "")
        pct = m.get("mhcflurry_affinity_percentile", "")
        proc = m.get("mhcflurry_processing_score", "")
        pres = m.get("mhcflurry_presentation_score", "")
        def evidence_value(key, default=""):
            return n.get(key) or m.get(key) or p.get(key) or default
        ba_s = norm_rank(ba) if ba != "" else None
        el_s = norm_rank(el) if el != "" else None
        pct_s = norm_rank(pct) if pct != "" else None
        proc_s = clamp(to_float(proc, -1)) if proc != "" else None
        pres_s = clamp(to_float(pres, -1)) if pres != "" else None
        chop_raw = c.get("netchop_31d_cterm_score") or c.get("netchop_31d_max_score", "")
        chop_s = clamp(to_float(chop_raw, -1)) if chop_raw != "" else None
        binding_parts = [x for x in [ba_s, pct_s] if x is not None]
        binding = max(binding_parts) if binding_parts else norm_rank(p.get("binding_rank", 99))
        num = den = 0.0
        is_mhci = str(p.get("mhc_class") or "I").upper() not in {"II", "MHC-II", "CLASSII"}
        weighted_evidence = [(ba_s,w_ba),(el_s,w_el),(pres_s,w_mhcf),(proc_s,w_proc)]
        if is_mhci:
            weighted_evidence.append((chop_s, w_chop))
        for val, wt in weighted_evidence:
            if val is not None:
                num += val * wt; den += wt
        if den:
            presentation = num / den
            complete = min(1.0, den / sum(wt for _, wt in weighted_evidence))
        else:
            presentation = clamp(to_float(p.get("presentation_score"), 0.0))
            complete = 0.25 if p.get("presentation_score") else 0.0
        rows.append({
            "peptide_id": p.get("peptide_id",""),
            "event_id": p.get("event_id",""),
            "sample_id": p.get("sample_id",""),
            "peptide": p.get("peptide",""),
            # Retained in memory so downstream immunogenicity merges can look up
            # the matched WT peptide. The canonical presentation schema stores
            # the resulting WT scores rather than this helper field.
            "wildtype_peptide": wt_peptide,
            "hla_allele": p.get("hla_allele",""),
            "mhc_class": p.get("mhc_class",""),
            "netmhcpan_ba_rank": str(to_float(ba, 99.0)),
            "netmhcpan_el_rank": str(to_float(el, 99.0)),
            "netmhcpan_mt_rank_ba": evidence_value("netmhcpan_mt_rank_ba", ba),
            "netmhcpan_mt_rank_el": evidence_value("netmhcpan_mt_rank_el", el),
            "netmhcpan_wt_rank_ba": wt_n.get("netmhcpan_ba_rank") or evidence_value("netmhcpan_wt_rank_ba"),
            "netmhcpan_wt_rank_el": wt_n.get("netmhcpan_el_rank") or evidence_value("netmhcpan_wt_rank_el"),
            "netmhcstabpan_score": str(to_float(s.get("netmhcstabpan_score"), 0.0)) if s else "",
            "netmhcstabpan_rank": str(to_float(s.get("netmhcstabpan_rank"), 99.0)) if s else "",
            "netchop_31d_max_score": c.get("netchop_31d_max_score", ""),
            "netchop_31d_mean_score": c.get("netchop_31d_mean_score", ""),
            "netchop_31d_cterm_score": c.get("netchop_31d_cterm_score", ""),
            "netchop_31d_cleavage_sites": c.get("netchop_31d_cleavage_sites", ""),
            "netchop_processing_status": c.get("netchop_processing_status", "ASSESSED" if c else "UNASSESSED"),
            "mhcflurry_affinity_percentile": str(to_float(pct, 99.0)),
            "mhcflurry_processing_score": str(to_float(proc, 0.0)),
            "mhcflurry_presentation_score": str(to_float(pres, 0.0)),
            "mhcflurry_wt_affinity_percentile": wt_m.get("mhcflurry_affinity_percentile") or evidence_value("mhcflurry_wt_affinity_percentile"),
            "mhcflurry_wt_processing_score": wt_m.get("mhcflurry_processing_score") or evidence_value("mhcflurry_wt_processing_score"),
            "mhcflurry_wt_presentation_score": wt_m.get("mhcflurry_presentation_score") or evidence_value("mhcflurry_wt_presentation_score"),
            "prime_wt_score": evidence_value("prime_wt_score"),
            "prime_wt_rank": evidence_value("prime_wt_rank"),
            "bigmhc_im_wt_score": evidence_value("bigmhc_im_wt_score"),
            "binding_evidence_score": f"{binding:.4f}",
            "presentation_evidence_score": f"{presentation:.4f}",
            "evidence_completeness": f"{complete:.4f}",
            "presentation_evidence_grade": grade(binding, presentation, complete),
        })
    registry = provenance_registry or ProvenanceRegistry()
    if netmhcpan and not registry.has("netmhcpan"):
        registry.register_passthrough("netmhcpan", netmhcpan)
    if mhcflurry and not registry.has("mhcflurry"):
        registry.register_passthrough("mhcflurry", mhcflurry)
    if netmhcstabpan and not registry.has("netmhcstabpan"):
        registry.register_passthrough("netmhcstabpan", netmhcstabpan)
    if netchop and not registry.has("netchop"):
        registry.register_passthrough("netchop", netchop)
    if out:
        summary = registry.tool_summary_fields()
        composite = provenance_derived(
            "presentation_composite",
            out,
            upstream="netmhcpan+mhcflurry+netmhcstabpan+netchop+immunogenicity",
        )
        for row in rows:
            row.update(summary)
            row.update(composite.as_fields())
        write_tsv(out, rows, PRESENTATION_FIELDS)
    return rows
