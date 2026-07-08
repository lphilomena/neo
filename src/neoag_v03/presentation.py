from __future__ import annotations
from pathlib import Path
from .utils import read_tsv, write_tsv, to_float, norm_rank, clamp, safe_id
from .schemas import PRESENTATION_FIELDS
from .evidence_provenance import ProvenanceRecord, ProvenanceRegistry, provenance_derived, attach_provenance

def by_key(rows):
    return {r.get("peptide_hla_key") or safe_id(f"{r.get('peptide','')}_{r.get('hla_allele','')}"): r for r in rows}

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
    provenance_registry: ProvenanceRegistry | None = None,
):
    peptides = read_tsv(raw_peptides)
    net = by_key(read_tsv(netmhcpan)) if netmhcpan else {}
    mhc = by_key(read_tsv(mhcflurry)) if mhcflurry else {}
    stab = by_key(read_tsv(netmhcstabpan)) if netmhcstabpan else {}
    w = profile.get("presentation_weights", {})
    w_ba = float(w.get("netmhcpan_ba", 0.25))
    w_el = float(w.get("netmhcpan_el", 0.35))
    w_mhcf = float(w.get("mhcflurry_presentation", 0.30))
    w_proc = float(w.get("mhcflurry_processing", 0.10))
    rows = []
    for p in peptides:
        key = safe_id(f"{p.get('peptide','')}_{p.get('hla_allele','')}")
        n = net.get(key, {})
        m = mhc.get(key, {})
        s = stab.get(key, {})
        ba = n.get("netmhcpan_ba_rank", "")
        el = n.get("netmhcpan_el_rank", "")
        pct = m.get("mhcflurry_affinity_percentile", "")
        proc = m.get("mhcflurry_processing_score", "")
        pres = m.get("mhcflurry_presentation_score", "")
        ba_s = norm_rank(ba) if ba != "" else None
        el_s = norm_rank(el) if el != "" else None
        pct_s = norm_rank(pct) if pct != "" else None
        proc_s = clamp(to_float(proc, -1)) if proc != "" else None
        pres_s = clamp(to_float(pres, -1)) if pres != "" else None
        binding_parts = [x for x in [ba_s, pct_s] if x is not None]
        binding = max(binding_parts) if binding_parts else norm_rank(p.get("binding_rank", 99))
        num = den = 0.0
        for val, wt in [(ba_s,w_ba),(el_s,w_el),(pres_s,w_mhcf),(proc_s,w_proc)]:
            if val is not None:
                num += val * wt; den += wt
        if den:
            presentation = num / den
            complete = min(1.0, den / (w_ba+w_el+w_mhcf+w_proc))
        else:
            presentation = clamp(to_float(p.get("presentation_score"), 0.0))
            complete = 0.25 if p.get("presentation_score") else 0.0
        rows.append({
            "peptide_id": p.get("peptide_id",""),
            "event_id": p.get("event_id",""),
            "sample_id": p.get("sample_id",""),
            "peptide": p.get("peptide",""),
            "hla_allele": p.get("hla_allele",""),
            "mhc_class": p.get("mhc_class",""),
            "netmhcpan_ba_rank": str(to_float(ba, 99.0)),
            "netmhcpan_el_rank": str(to_float(el, 99.0)),
            "netmhcstabpan_score": str(to_float(s.get("netmhcstabpan_score"), 0.0)) if s else "",
            "netmhcstabpan_rank": str(to_float(s.get("netmhcstabpan_rank"), 99.0)) if s else "",
            "mhcflurry_affinity_percentile": str(to_float(pct, 99.0)),
            "mhcflurry_processing_score": str(to_float(proc, 0.0)),
            "mhcflurry_presentation_score": str(to_float(pres, 0.0)),
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
    if out:
        summary = registry.tool_summary_fields()
        composite = provenance_derived(
            "presentation_composite",
            out,
            upstream="netmhcpan+mhcflurry+netmhcstabpan+immunogenicity",
        )
        for row in rows:
            row.update(summary)
            row.update(composite.as_fields())
        write_tsv(out, rows, PRESENTATION_FIELDS)
    return rows
