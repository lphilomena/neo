from __future__ import annotations

from pathlib import Path
import html
from typing import Mapping

from . import __version__
from .utils import read_tsv


def esc(x):
    return html.escape(str(x if x is not None else ""))


def _read_optional(path):
    return read_tsv(path) if path and Path(path).exists() else []


def _badge(text: str) -> str:
    t = str(text or "")
    cls = "UNASSESSED"
    if any(x in t.upper() for x in ["PASS", "INTACT", "HIGH", "A", "B"]):
        cls = "PASS"
    if any(x in t.upper() for x in ["CAUTION", "REVIEW", "MEDIUM", "LOW", "C"]):
        cls = "CAUTION"
    if any(x in t.upper() for x in ["DEFECT", "REJECT", "FAIL", "LOST", "GLOBAL", "D"]):
        cls = "FAIL"
    if any(x in t.upper() for x in ["UNASSESSED", "INSUFFICIENT", "MISSING", "INCONCLUSIVE"]):
        cls = "UNASSESSED"
    return f"<span class='badge {cls}'>{esc(text)}</span>"


def _table(rows, headers, max_rows=None):
    rows = rows[:max_rows] if max_rows else rows
    out = ["<table><tr>" + "".join(f"<th>{esc(h)}</th>" for h in headers) + "</tr>"]
    for r in rows:
        out.append("<tr>" + "".join(f"<td>{esc(r.get(h, ''))}</td>" for h in headers) + "</tr>")
    out.append("</table>")
    return "\n".join(out)


def _map_by(rows, key):
    return {r.get(key, ""): r for r in rows if r.get(key)}


def make_report_v041(
    path,
    profile: Mapping,
    events,
    peptides,
    appm_summary=None,
    validation_rows=None,
    *,
    appm_gene_status=None,
    appm_module_scores=None,
    appm_submodule_scores=None,
    appm_conflicts=None,
    appm_peptide_modifiers=None,
    immune_escape_summary=None,
    peptide_escape_flags=None,
    peptide_safety=None,
    ccf=None,
):
    appm_summary = appm_summary or {}
    validation_rows = validation_rows or []
    appm_gene_status = appm_gene_status or []
    appm_module_scores = appm_module_scores or []
    appm_submodule_scores = appm_submodule_scores or []
    appm_conflicts = appm_conflicts or []
    appm_peptide_modifiers = appm_peptide_modifiers or []
    immune_escape_summary = immune_escape_summary or []
    peptide_escape_flags = peptide_escape_flags or []
    peptide_safety = peptide_safety or []
    ccf = ccf or []
    mod_by_pep = _map_by(appm_peptide_modifiers, "peptide_id")
    esc_by_pep = _map_by(peptide_escape_flags, "peptide_id")
    safe_by_pep = _map_by(peptide_safety, "peptide_id")
    ccf_by_event = _map_by(ccf, "event_id")

    p = Path(path); p.parent.mkdir(parents=True, exist_ok=True)
    css = """
<style>
body{font-family:Arial,sans-serif;margin:32px;color:#222;line-height:1.35}h1,h2,h3{color:#17324d}.section{margin-top:26px}
table{border-collapse:collapse;width:100%;margin:12px 0 24px}th,td{border:1px solid #ddd;padding:7px;font-size:12px;vertical-align:top}th{background:#f3f6f9}.badge{padding:3px 7px;border-radius:8px;font-size:12px;display:inline-block}.PASS{background:#d6f5d6}.CAUTION{background:#fff1b8}.FAIL{background:#ffd6d6}.UNASSESSED{background:#eee;color:#555}.card{border:1px solid #ddd;border-radius:10px;padding:14px;margin:12px 0;box-shadow:0 1px 4px #eee}.small{color:#555;font-size:13px}.grid{display:grid;grid-template-columns:repeat(3,1fr);gap:12px}.metric{border:1px solid #ddd;border-radius:8px;padding:10px;background:#fafafa}.warn{background:#fff7e6;border-left:4px solid #e6a700;padding:10px;margin:12px 0}.mono{font-family:Menlo,Consolas,monospace;font-size:12px}
</style>
"""
    out = [f"<!doctype html><html><head><meta charset='utf-8'><title>NeoAg v{esc(__version__)} Evidence Report</title>", css, "</head><body>"]
    out.append(f"<h1>NeoAg Event Pipeline v{esc(__version__)} Evidence Report</h1><p><b>Profile:</b> {esc(profile.get('_profile_name'))}</p>")
    out.append("<div class='warn'><b>Boundary:</b> Computational triage only. APPM defects and immune escape flags are mechanism evidence, not clinical resistance diagnosis. Missing evidence is not evidence of intact presentation.</div>")

    out.append("<div class='section'><h2>Sample-level APPM Evidence Card</h2>")
    out.append("<div class='grid'>")
    for label, score_key, status_key in [
        ("MHC-I", "mhc_i_integrity_score", "mhc_i_integrity_status"),
        ("MHC-II", "mhc_ii_integrity_score", "mhc_ii_integrity_status"),
        ("IFNG response", "ifng_response_score", "ifng_response_status"),
    ]:
        out.append(f"<div class='metric'><b>{esc(label)}</b><br>Score: <b>{esc(appm_summary.get(score_key,''))}</b><br>Status: {_badge(appm_summary.get(status_key,''))}</div>")
    out.append("</div>")
    out.append(f"<p><b>APPM call confidence:</b> {_badge(appm_summary.get('appm_call_confidence',''))} score={esc(appm_summary.get('appm_call_confidence_score',''))}<br><b>Reason:</b> <span class='mono'>{esc(appm_summary.get('confidence_reason',''))}</span></p>")
    out.append(f"<p><b>Evidence completeness:</b> {esc(appm_summary.get('appm_evidence_completeness',''))} ({esc(appm_summary.get('appm_evidence_completeness_score',''))}); <b>Conflicts:</b> {esc(appm_summary.get('evidence_conflict_flag',''))}</p>")
    out.append("</div>")

    if appm_submodule_scores:
        out.append("<div class='section'><h2>MHC/APPM Submodule Scores</h2>")
        headers = ["parent_module", "submodule", "score", "status", "defect_severity", "appm_call_confidence", "driver_defects", "action_hint", "confidence_reason"]
        out.append(_table(appm_submodule_scores, headers))
        out.append("</div>")

    if appm_gene_status:
        out.append("<div class='section'><h2>Top APPM Driver Defects</h2>")
        important = [r for r in appm_gene_status if r.get("functional_status") in {"defective", "caution"} or r.get("biallelic_status") == "BIALLELIC_LOSS"]
        if not important:
            important = [r for r in appm_gene_status if r.get("gene") in {"B2M","TAP1","TAP2","NLRC5","CIITA","JAK1","JAK2","HLA-A","HLA-B","HLA-C"}]
        headers = ["gene", "pathway", "biallelic_status", "functional_status", "copy_number_status", "loh_status", "expression_status", "gene_integrity_status", "reason"]
        out.append(_table(important, headers, max_rows=20))
        out.append("</div>")

    if immune_escape_summary:
        ies = immune_escape_summary[0]
        out.append("<div class='section'><h2>Immune Escape Summary</h2>")
        out.append(f"<p><b>Overall risk:</b> {_badge(ies.get('overall_immune_escape_risk',''))}; <b>Mechanisms:</b> {esc(ies.get('mechanism_summary',''))}; <b>Context:</b> {esc(ies.get('therapy_context',''))}</p>")
        out.append(f"<p><b>Burden:</b> {esc(ies.get('escape_burden_summary',''))}</p>")
        out.append("</div>")

    out.append("<div class='section'><h2>Top Events</h2>")
    out.append(_table(events, ["event_id", "event_name", "event_type", "gene", "event_score", "ccf_status", "safety_status"], max_rows=20))
    out.append("</div>")

    out.append("<div class='section'><h2>Top Peptides</h2>")
    headers = ["peptide_id", "peptide", "hla_allele", "gene", "presentation_evidence_grade", "safety_status", "escape_status", "appm_action", "final_priority", "efficacy_score"]
    pep_rows=[]
    for ppt in peptides[:40]:
        pid=ppt.get("peptide_id", "")
        e=esc_by_pep.get(pid, {})
        a=mod_by_pep.get(pid, {})
        s=safe_by_pep.get(pid, {})
        row={**ppt, "escape_status": e.get("escape_status", ppt.get("escape_status", "")), "appm_action": a.get("appm_action", ppt.get("appm_action", "")), "safety_status": s.get("safety_status", ppt.get("safety_status", ""))}
        pep_rows.append(row)
    out.append(_table(pep_rows, headers, max_rows=40))
    out.append("</div>")

    out.append("<div class='section'><h2>Peptide Mechanism Cards</h2>")
    for ppt in peptides[:20]:
        pid=ppt.get("peptide_id", "")
        e=esc_by_pep.get(pid, {})
        a=mod_by_pep.get(pid, {})
        s=safe_by_pep.get(pid, {})
        c=ccf_by_event.get(ppt.get("event_id", ""), {})
        out.append("<div class='card'>")
        out.append(f"<h3>{esc(ppt.get('peptide'))} — {esc(ppt.get('hla_allele'))}</h3>")
        out.append(f"<p><b>Event:</b> {esc(ppt.get('event_type'))} / {esc(ppt.get('gene'))} / <span class='mono'>{esc(ppt.get('event_id'))}</span></p>")
        out.append(f"<p><b>Decision:</b> {_badge(ppt.get('final_priority',''))}; score={esc(ppt.get('efficacy_score'))}; cap={esc(ppt.get('priority_cap',''))}</p>")
        out.append(f"<p><b>APPM:</b> action={esc(a.get('appm_action', ppt.get('appm_action','')))}; multiplier={esc(a.get('appm_multiplier', ppt.get('appm_multiplier','')))}; confidence={esc(a.get('appm_call_confidence',''))}; reason=<span class='mono'>{esc(a.get('appm_multiplier_reason', ppt.get('appm_reason','')))}</span></p>")
        out.append(f"<p><b>Immune escape:</b> {_badge(e.get('escape_status',''))}; multiplier={esc(e.get('escape_multiplier',''))}; reason=<span class='mono'>{esc(e.get('escape_reason',''))}</span></p>")
        out.append(f"<p><b>Safety:</b> {_badge(s.get('safety_status', ppt.get('safety_status','')))}; tier={esc(s.get('safety_tier',''))}; reason=<span class='mono'>{esc(s.get('safety_reason', ppt.get('safety_reason','')))}</span></p>")
        out.append(f"<p><b>CCF:</b> status={esc(c.get('ccf_status', ppt.get('ccf_status','')))}; best={esc(c.get('ccf_best', ppt.get('ccf_estimate','')))}; confidence={esc(c.get('ccf_confidence',''))}; method={esc(c.get('ccf_method',''))}</p>")
        out.append("</div>")
    out.append("</div>")

    out.append("</body></html>")
    p.write_text("\n".join(out), encoding="utf-8")
