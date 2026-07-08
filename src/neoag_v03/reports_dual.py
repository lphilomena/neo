"""Dual-audience HTML reports: patient communication vs research/technical."""

from __future__ import annotations

import html
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

from .utils import read_tsv

REPORT_CSS = """
<style>
body{font-family:Arial,sans-serif;margin:32px;color:#222;line-height:1.45;max-width:1100px}
h1,h2,h3{color:#17324d}.section{margin-top:28px}
table{border-collapse:collapse;width:100%;margin:12px 0 24px}th,td{border:1px solid #ddd;padding:7px;font-size:12px;vertical-align:top}th{background:#f3f6f9}
.badge{padding:3px 7px;border-radius:8px;font-size:12px;display:inline-block}.PASS{background:#d6f5d6}.CAUTION{background:#fff1b8}.FAIL{background:#ffd6d6}.UNASSESSED{background:#eee;color:#555}
.card{border:1px solid #ddd;border-radius:10px;padding:14px;margin:12px 0;box-shadow:0 1px 4px #eee}
.small{color:#555;font-size:13px}.mono{font-family:Menlo,Consolas,monospace;font-size:11px;word-break:break-all}
.warn{background:#fff7e6;border-left:4px solid #e6a700;padding:12px;margin:14px 0}
.info{background:#f0f7ff;border-left:4px solid #3b82f6;padding:12px;margin:14px 0}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px}
.metric{border:1px solid #ddd;border-radius:8px;padding:10px;background:#fafafa}
ul.compact{margin:8px 0 8px 20px}
.patient h1{font-size:1.6rem}.patient .lead{font-size:1.05rem;color:#333}
</style>
"""


def esc(x: Any) -> str:
    return html.escape(str(x if x is not None else ""))


def _read_optional(path: str | Path | None) -> list[dict[str, str]]:
    if not path:
        return []
    p = Path(path)
    return read_tsv(p) if p.is_file() else []


def _read_json_optional(path: str | Path | None) -> dict[str, Any]:
    if not path:
        return {}
    p = Path(path)
    if not p.is_file():
        return {}
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}


def _badge(text: str) -> str:
    t = str(text or "")
    cls = "UNASSESSED"
    if any(x in t.upper() for x in ["PASS", "INTACT", "HIGH", "A", "B"]):
        cls = "PASS"
    if any(x in t.upper() for x in ["CAUTION", "REVIEW", "MEDIUM", "LOW", "C"]):
        cls = "CAUTION"
    if any(x in t.upper() for x in ["DEFECT", "REJECT", "FAIL", "LOST", "GLOBAL", "D"]):
        cls = "FAIL"
    return f"<span class='badge {cls}'>{esc(text)}</span>"


def _table(rows: list[Mapping[str, Any]], headers: list[str], *, max_rows: int | None = None) -> str:
    view = rows[:max_rows] if max_rows else rows
    out = ["<table><tr>" + "".join(f"<th>{esc(h)}</th>" for h in headers) + "</tr>"]
    for row in view:
        out.append("<tr>" + "".join(f"<td>{esc(row.get(h, ''))}</td>" for h in headers) + "</tr>")
    out.append("</table>")
    return "\n".join(out)


def _map_by(rows: list[Mapping[str, Any]], key: str) -> dict[str, dict[str, str]]:
    return {str(r.get(key, "")): dict(r) for r in rows if r.get(key)}


PRIORITY_PATIENT = {
    "A": "优先推荐进一步验证",
    "B": "值得考虑验证",
    "B_CAUTION": "可考虑，但需关注安全性",
    "C": "证据有限，谨慎推进",
    "C_CAUTION": "证据有限且需安全性复核",
    "D": "当前不建议推进",
}

CONSEQUENCE_PATIENT = {
    "missense": "氨基酸改变（点突变）",
    "frameshift": "移码变异产生的新肽段",
    "splice_junction": "剪接异常产生的新肽段",
    "exon_deletion_junction": "外显子缺失/剪接交界肽段",
    "insertion": "插入/缺失交界肽段",
    "fusion": "基因融合交界肽段",
    "other": "其他变异相关肽段",
}

FIELD_GLOSSARY = {
    "efficacy_score": "综合免疫学评分（0–1），整合表达、结合、呈递、安全性等维度。",
    "final_priority": "最终优先级分层：A/B 为优先候选，C 为需更多证据，D 为不建议推进。",
    "presentation_evidence_grade": "HLA 结合与加工呈递证据等级（A 最优）。",
    "appm_multiplier": "抗原加工呈递通路（APPM）完整性对候选的折减系数。",
    "ccf_multiplier": "肿瘤克隆性（CCF）对候选的折减系数。",
    "safety_status": "正常组织表达、自身肽相似性等安全性初筛结果。",
    "escape_status": "免疫逃逸机制（如 HLA 丢失）对候选的影响评估。",
    "validation_mode": "建议的实验验证设计类型（短肽对、长肽、minigene 等）。",
    "recommended_assay": "推荐的体外验证实验类型。",
}


@dataclass
class ReportBundle:
    profile: Mapping[str, Any]
    events: list[dict[str, str]]
    peptides: list[dict[str, str]]
    appm_summary: Mapping[str, Any] = field(default_factory=dict)
    validation_rows: list[dict[str, str]] = field(default_factory=list)
    sample_id: str = ""
    entry_mode: str = ""
    provenance: dict[str, Any] = field(default_factory=dict)
    peptide_safety: list[dict[str, str]] = field(default_factory=list)
    peptide_escape_flags: list[dict[str, str]] = field(default_factory=list)
    immune_escape_summary: list[dict[str, str]] = field(default_factory=list)
    ccf: list[dict[str, str]] = field(default_factory=list)
    appm_gene_status: list[dict[str, str]] = field(default_factory=list)
    appm_peptide_modifiers: list[dict[str, str]] = field(default_factory=list)
    appm_module_scores: list[dict[str, str]] = field(default_factory=list)
    appm_submodule_scores: list[dict[str, str]] = field(default_factory=list)
    appm_conflicts: list[dict[str, str]] = field(default_factory=list)


def load_report_bundle(
    *,
    profile: Mapping[str, Any],
    events: list[dict[str, str]],
    peptides: list[dict[str, str]],
    appm_summary: Mapping[str, Any] | None = None,
    validation_rows: list[dict[str, str]] | None = None,
    outdir: str | Path | None = None,
    provenance: Mapping[str, Any] | None = None,
    sample_id: str = "",
    entry_mode: str = "",
) -> ReportBundle:
    root = Path(outdir) if outdir else None
    prov = dict(provenance or {})
    if root and not prov:
        prov = _read_json_optional(root / "provenance.v03.json")

    def p(*parts: str) -> Path | None:
        return root / Path(*parts) if root else None

    return ReportBundle(
        profile=profile,
        events=events,
        peptides=peptides,
        appm_summary=appm_summary or {},
        validation_rows=validation_rows or [],
        sample_id=sample_id or str(prov.get("sample_id") or (peptides[0].get("sample_id") if peptides else "")),
        entry_mode=entry_mode or str(prov.get("entry_mode") or ""),
        provenance=prov,
        peptide_safety=_read_optional(p("safety", "peptide_safety.tsv") if root else None),
        peptide_escape_flags=_read_optional(p("immune_escape", "peptide_escape_flags.tsv") if root else None),
        immune_escape_summary=_read_optional(p("immune_escape", "immune_escape_summary.tsv") if root else None),
        ccf=_read_optional(p("clonality", "ccf_2.tsv") if root else None) or _read_optional(p("clonality", "ccf_lite.tsv") if root else None),
        appm_gene_status=_read_optional(p("appm", "appm_gene_status.tsv") if root else None),
        appm_peptide_modifiers=_read_optional(p("appm", "appm_peptide_modifiers.tsv") if root else None),
        appm_module_scores=_read_optional(p("appm", "appm_module_scores.tsv") if root else None),
        appm_submodule_scores=_read_optional(p("appm", "appm_submodule_scores.tsv") if root else None),
        appm_conflicts=_read_optional(p("appm", "appm_conflicts.tsv") if root else None),
    )


def _patient_consequence_label(peptide: Mapping[str, Any]) -> str:
    pc = str(peptide.get("peptide_consequence") or "").lower()
    if pc in CONSEQUENCE_PATIENT:
        return CONSEQUENCE_PATIENT[pc]
    et = str(peptide.get("event_type") or "")
    if et == "Fusion":
        return CONSEQUENCE_PATIENT["fusion"]
    return CONSEQUENCE_PATIENT.get(pc, "肿瘤特异性肽段候选")


def _patient_priority_label(priority: str) -> str:
    return PRIORITY_PATIENT.get(str(priority or "").strip(), "需进一步评估")


def _appm_patient_summary(appm_summary: Mapping[str, Any]) -> str:
    i_status = str(appm_summary.get("mhc_i_integrity_status") or "")
    if "DEFECT" in i_status.upper():
        return "样本 MHC-I 抗原呈递通路可能存在缺陷，部分候选肽段的实际呈递能力可能低于计算预测。"
    if "PARTIAL" in i_status.upper() or "CAUTION" in i_status.upper():
        return "样本 MHC-I 抗原呈递通路存在不确定因素，建议结合实验验证解读候选。"
    return "样本 MHC-I 抗原呈递通路完整性评估未见明确缺陷信号（仍依赖输入证据完整度）。"


def _top_hla_alleles(peptides: list[Mapping[str, Any]], limit: int = 6) -> str:
    alleles: list[str] = []
    for p in peptides:
        hla = str(p.get("hla_allele") or "").strip()
        if hla and hla not in alleles:
            alleles.append(hla)
        if len(alleles) >= limit:
            break
    return "、".join(alleles) if alleles else "未提供"


def _val_by_peptide(validation_rows: list[Mapping[str, Any]]) -> dict[str, dict[str, str]]:
    return _map_by(validation_rows, "peptide_id")


def make_patient_report(path: str | Path, bundle: ReportBundle) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    val_map = _val_by_peptide(bundle.validation_rows)
    advancable = [ppt for ppt in bundle.peptides if str(ppt.get("final_priority", "")).upper() not in {"D", ""}]
    top = advancable[:10] if advancable else bundle.peptides[:10]
    genes = []
    for ppt in top:
        g = str(ppt.get("gene") or "")
        if g and g not in genes:
            genes.append(g)

    out = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        f"<title>新抗原分析报告（患者沟通版）— {esc(bundle.sample_id)}</title>",
        REPORT_CSS,
        "</head><body class='patient'>",
        "<h1>新抗原计算分析报告</h1>",
        "<p class='lead'>本报告用于帮助理解肿瘤新抗原候选的<strong>计算筛选结果</strong>，"
        "便于医患沟通与后续实验规划；<strong>不能替代临床诊断或治疗决策</strong>。</p>",
    ]

    out.append("<div class='section'><h2>1. 分析目的</h2>")
    out.append(
        "<p>通过肿瘤 DNA/RNA 与 HLA 分型信息，识别可能由肿瘤变异产生、并被患者 HLA 分子呈递的肽段候选，"
        "用于疫苗设计、T 细胞治疗靶点探索或免疫监测研究的初步筛选。</p>"
    )
    out.append("</div>")

    out.append("<div class='section'><h2>2. 样本与分型概况</h2>")
    out.append("<ul class='compact'>")
    out.append(f"<li><b>样本编号：</b>{esc(bundle.sample_id or '未注明')}</li>")
    out.append(f"<li><b>分析场景：</b>{esc(bundle.entry_mode or '肿瘤新抗原筛选')}</li>")
    out.append(f"<li><b>患者 HLA 分型（来自分析输入）：</b>{esc(_top_hla_alleles(bundle.peptides))}</li>")
    out.append(f"<li><b>评分方案：</b>{esc(bundle.profile.get('_profile_name', 'default'))}</li>")
    out.append("</ul></div>")

    out.append("<div class='section'><h2>3. 关键发现（摘要）</h2><ul class='compact'>")
    out.append(f"<li>共评估 <b>{len(bundle.peptides)}</b> 条肽段候选，其中 <b>{len(advancable)}</b> 条未被判为“不建议推进”。</li>")
    if genes:
        out.append(f"<li>优先关注的基因包括：<b>{esc('、'.join(genes[:8]))}</b>。</li>")
    out.append(f"<li>{esc(_appm_patient_summary(bundle.appm_summary))}</li>")
    if bundle.immune_escape_summary:
        ies = bundle.immune_escape_summary[0]
        risk = str(ies.get("overall_immune_escape_risk") or "")
        if risk and risk.upper() not in {"LOW", "PASS", ""}:
            out.append(f"<li>免疫逃逸风险提示：<b>{esc(risk)}</b>（机制层面证据，非临床耐药结论）。</li>")
    out.append("</ul></div>")

    out.append("<div class='section'><h2>4. 优先候选肽段（Top 10）</h2>")
    out.append("<p class='small'>下表为计算排序靠前的候选，不代表已证实可诱导抗肿瘤免疫反应。</p>")
    headers = ["rank", "gene", "variant_type", "hla", "priority", "meaning", "next_step"]
    rows = []
    for i, ppt in enumerate(top, 1):
        pid = str(ppt.get("peptide_id") or "")
        val = val_map.get(pid, {})
        mode = str(val.get("validation_mode") or "")
        next_step = str(val.get("validation_strategy") or val.get("recommended_assay") or ppt.get("recommended_use") or "")
        if mode in {"frameshift_long", "splice_junction_long", "fusion_junction_long", "insertion_long"}:
            meaning = "可能由肿瘤特异性序列产生，建议用长肽/minigene 验证真实加工呈递"
        else:
            meaning = "突变肽 vs 正常肽对照验证"
        rows.append({
            "rank": str(i),
            "gene": ppt.get("gene", ""),
            "variant_type": _patient_consequence_label(ppt),
            "hla": ppt.get("hla_allele", ""),
            "priority": _patient_priority_label(str(ppt.get("final_priority") or "")),
            "meaning": meaning,
            "next_step": next_step,
        })
    out.append(_table(rows, headers))
    out.append("</div>")

    out.append("<div class='section'><h2>5. 建议的下一步工作</h2><ul class='compact'>")
    out.append("<li><b>实验验证：</b>对优先候选开展体外免疫学实验（如 ELISpot、四聚体或长肽/minigene 刺激），确认 T 细胞识别。</li>")
    out.append("<li><b>移码/剪接/融合候选：</b>不应仅依赖 8–11 氨基酸短肽；需验证更接近体内加工的长肽或 minigene。</li>")
    out.append("<li><b>安全性复核：</b>对“需谨慎”候选关注正常组织表达与自身肽相似性。</li>")
    out.append("<li><b>临床决策：</b>任何治疗或疫苗方案须由临床团队结合病理、分期与指南独立判断。</li>")
    out.append("</ul></div>")

    out.append("<div class='warn'><h2>6. 重要说明与风险边界</h2><ul class='compact'>")
    out.append("<li>本分析为<strong>计算机辅助筛选</strong>，预测结合亲和力不等于体内呈递，更不等于临床疗效。</li>")
    out.append("<li>APPM、CCF、安全性与免疫逃逸评估依赖输入数据完整度；缺失数据不等于“无风险”。</li>")
    out.append("<li>本报告<strong>不包含</strong>原始测序质控、文件路径或生信命令细节；技术细节见科研技术版报告。</li>")
    out.append("<li>不得将本报告直接用于患者诊断、预后判断或个体化治疗处方。</li>")
    out.append("</ul></div>")
    out.append("</body></html>")
    p.write_text("\n".join(out), encoding="utf-8")


def _profile_threshold_section(profile: Mapping[str, Any]) -> str:
    sections = []
    for key in ("gates", "safety", "ccf_lite", "v03_score_weights", "l3_weights", "appm_penalty"):
        block = profile.get(key)
        if isinstance(block, Mapping) and block:
            rows = [{"parameter": k, "value": v} for k, v in block.items()]
            sections.append(f"<h3>{esc(key)}</h3>{_table(rows, ['parameter', 'value'])}")
    return "".join(sections)


def _provenance_section(provenance: Mapping[str, Any]) -> str:
    if not provenance:
        return "<p class='small'>No provenance metadata supplied.</p>"
    out = ["<ul class='compact'>"]
    out.append(f"<li><b>sample_id:</b> <span class='mono'>{esc(provenance.get('sample_id'))}</span></li>")
    out.append(f"<li><b>entry_mode:</b> {esc(provenance.get('entry_mode'))}</li>")
    out.append(f"<li><b>profile:</b> {esc(provenance.get('profile'))}</li>")
    out.append(f"<li><b>created_at:</b> {esc(provenance.get('created_at'))}</li>")
    tools = provenance.get("tools") or {}
    if tools:
        out.append("</ul><h3>Tool provenance</h3><table><tr><th>tool</th><th>status</th><th>version</th><th>file</th><th>mode</th></tr>")
        for name, rec in tools.items():
            if not isinstance(rec, Mapping):
                continue
            out.append(
                "<tr>"
                f"<td>{esc(name)}</td><td>{esc(rec.get('status'))}</td><td>{esc(rec.get('version'))}</td>"
                f"<td class='mono'>{esc(rec.get('file'))}</td><td>{esc(rec.get('mode'))}</td>"
                "</tr>"
            )
        out.append("</table>")
    else:
        out.append("</ul>")
    return "\n".join(out)


def make_technical_report(path: str | Path, bundle: ReportBundle) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    mod_by_pep = _map_by(bundle.appm_peptide_modifiers, "peptide_id")
    esc_by_pep = _map_by(bundle.peptide_escape_flags, "peptide_id")
    safe_by_pep = _map_by(bundle.peptide_safety, "peptide_id")
    ccf_by_event = _map_by(bundle.ccf, "event_id")
    val_map = _val_by_peptide(bundle.validation_rows)

    out = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        f"<title>NeoAg Technical Evidence Report — {esc(bundle.sample_id)}</title>",
        REPORT_CSS,
        "</head><body>",
        "<h1>NeoAg Pipeline — Research / Technical Report</h1>",
        f"<p><b>Profile:</b> {esc(bundle.profile.get('_profile_name'))} &nbsp; "
        f"<b>Sample:</b> <span class='mono'>{esc(bundle.sample_id)}</span> &nbsp; "
        f"<b>Entry mode:</b> {esc(bundle.entry_mode)}</p>",
        "<div class='warn'><b>Boundary:</b> Computational triage only. "
        "Mechanism flags (APPM, escape, safety) are evidence layers, not clinical diagnoses.</div>",
    ]

    out.append("<div class='section'><h2>Pipeline &amp; Provenance</h2>")
    out.append(_provenance_section(bundle.provenance))
    out.append("</div>")

    out.append("<div class='section'><h2>Profile Thresholds &amp; Weights</h2>")
    out.append(_profile_threshold_section(bundle.profile))
    out.append("</div>")

    out.append("<div class='section'><h2>APPM Summary</h2>")
    if bundle.appm_summary:
        rows = [{"field": k, "value": v} for k, v in bundle.appm_summary.items()]
        out.append(_table(rows, ["field", "value"]))
    if bundle.appm_submodule_scores:
        out.append("<h3>APPM Submodule Scores</h3>")
        out.append(_table(bundle.appm_submodule_scores, [
            "parent_module", "submodule", "score", "status", "defect_severity",
            "appm_call_confidence", "driver_defects", "action_hint", "confidence_reason",
        ]))
    if bundle.appm_gene_status:
        out.append("<h3>APPM Gene Status</h3>")
        out.append(_table(bundle.appm_gene_status, [
            "gene", "pathway", "biallelic_status", "functional_status", "copy_number_status",
            "loh_status", "expression_status", "gene_integrity_status", "reason",
        ], max_rows=50))
    if bundle.appm_conflicts:
        out.append("<h3>APPM Conflicts</h3>")
        out.append(_table(bundle.appm_conflicts, list(bundle.appm_conflicts[0].keys())))
    out.append("</div>")

    if bundle.immune_escape_summary:
        out.append("<div class='section'><h2>Immune Escape Summary</h2>")
        out.append(_table(bundle.immune_escape_summary, list(bundle.immune_escape_summary[0].keys())))
        out.append("</div>")

    out.append("<div class='section'><h2>Ranked Events (full)</h2>")
    out.append(_table(bundle.events, [
        "event_id", "event_name", "event_type", "mutation_source", "peptide_consequence", "gene",
        "event_score", "ccf_estimate", "ccf_status", "clonality_multiplier", "safety_status", "safety_reason",
    ]))
    out.append("</div>")

    out.append("<div class='section'><h2>Ranked Peptides (full)</h2>")
    pep_headers = [
        "peptide_id", "event_id", "gene", "peptide", "wildtype_peptide", "peptide_consequence",
        "hla_allele", "mhc_class", "presentation_evidence_grade", "binding_evidence_score",
        "presentation_evidence_score", "netmhcpan_ba_rank", "netmhcpan_el_rank",
        "appm_multiplier", "ccf_multiplier", "safety_status", "escape_status",
        "efficacy_score", "final_priority", "recommended_use",
    ]
    out.append(_table(bundle.peptides, pep_headers))
    out.append("</div>")

    if bundle.validation_rows:
        out.append("<div class='section'><h2>Validation Plan (full)</h2>")
        headers = list(bundle.validation_rows[0].keys())
        out.append(_table(bundle.validation_rows, headers))
        out.append("</div>")

    if bundle.peptide_safety:
        out.append("<div class='section'><h2>Peptide Safety Evidence</h2>")
        out.append(_table(bundle.peptide_safety, list(bundle.peptide_safety[0].keys()), max_rows=100))
        out.append("</div>")

    if bundle.peptide_escape_flags:
        out.append("<div class='section'><h2>Peptide Escape Flags</h2>")
        out.append(_table(bundle.peptide_escape_flags, list(bundle.peptide_escape_flags[0].keys()), max_rows=100))
        out.append("</div>")

    if bundle.ccf:
        out.append("<div class='section'><h2>CCF / Clonality</h2>")
        out.append(_table(bundle.ccf, list(bundle.ccf[0].keys()), max_rows=100))
        out.append("</div>")

    out.append("<div class='section'><h2>Peptide Mechanism Cards</h2>")
    for ppt in bundle.peptides[:25]:
        pid = str(ppt.get("peptide_id") or "")
        e = esc_by_pep.get(pid, {})
        a = mod_by_pep.get(pid, {})
        s = safe_by_pep.get(pid, {})
        c = ccf_by_event.get(str(ppt.get("event_id") or ""), {})
        v = val_map.get(pid, {})
        out.append("<div class='card'>")
        out.append(f"<h3>{esc(ppt.get('peptide'))} — {esc(ppt.get('hla_allele'))}</h3>")
        out.append(
            f"<p><b>IDs:</b> peptide_id=<span class='mono'>{esc(pid)}</span>; "
            f"event_id=<span class='mono'>{esc(ppt.get('event_id'))}</span></p>"
        )
        out.append(
            f"<p><b>Layers:</b> mutation_source={esc(ppt.get('mutation_source'))}; "
            f"peptide_consequence={esc(ppt.get('peptide_consequence'))}; source_tool={esc(ppt.get('source_tool'))}</p>"
        )
        out.append(
            f"<p><b>Presentation:</b> grade={esc(ppt.get('presentation_evidence_grade'))}; "
            f"BA={esc(ppt.get('netmhcpan_ba_rank'))}; EL={esc(ppt.get('netmhcpan_el_rank'))}; "
            f"MHCflurry={esc(ppt.get('mhcflurry_presentation_score'))}</p>"
        )
        out.append(
            f"<p><b>APPM:</b> multiplier={esc(a.get('appm_multiplier', ppt.get('appm_multiplier')))}; "
            f"reason=<span class='mono'>{esc(a.get('appm_multiplier_reason', ''))}</span></p>"
        )
        out.append(
            f"<p><b>Escape:</b> {_badge(e.get('escape_status', ppt.get('escape_status')))}; "
            f"multiplier={esc(e.get('escape_multiplier', ppt.get('escape_multiplier')))}; "
            f"reason=<span class='mono'>{esc(e.get('escape_reason', ''))}</span></p>"
        )
        out.append(
            f"<p><b>Safety:</b> {_badge(s.get('safety_status', ppt.get('safety_status')))}; "
            f"tier={esc(s.get('safety_tier', ''))}; reason=<span class='mono'>{esc(s.get('safety_reason', ''))}</span></p>"
        )
        out.append(
            f"<p><b>CCF:</b> status={esc(c.get('ccf_status', ppt.get('ccf_status')))}; "
            f"estimate={esc(c.get('ccf_estimate', ppt.get('ccf_estimate')))}; "
            f"multiplier={esc(c.get('clonality_multiplier', ppt.get('ccf_multiplier')))}</p>"
        )
        if v:
            out.append(
                f"<p><b>Validation design:</b> mode={esc(v.get('validation_mode'))}; "
                f"assay={esc(v.get('recommended_assay'))}; minigene=<span class='mono'>{esc(v.get('minigene'))}</span></p>"
            )
        out.append(
            f"<p><b>Decision:</b> {_badge(ppt.get('final_priority'))}; "
            f"efficacy_score={esc(ppt.get('efficacy_score'))}; {esc(ppt.get('recommended_use'))}</p>"
        )
        out.append("</div>")
    out.append("</div>")

    out.append("<div class='section'><h2>Field Glossary</h2><table><tr><th>Field</th><th>Description</th></tr>")
    for field_name, desc in FIELD_GLOSSARY.items():
        out.append(f"<tr><td class='mono'>{esc(field_name)}</td><td>{esc(desc)}</td></tr>")
    out.append("</table></div>")
    out.append("</body></html>")
    p.write_text("\n".join(out), encoding="utf-8")


def make_dual_reports(
    reports_dir: str | Path,
    bundle: ReportBundle,
    *,
    patient_name: str = "evidence_report.patient.html",
    technical_name: str = "evidence_report.technical.html",
    legacy_name: str = "evidence_report.v03.html",
) -> dict[str, str]:
    reports_dir = Path(reports_dir)
    patient_path = reports_dir / patient_name
    technical_path = reports_dir / technical_name
    legacy_path = reports_dir / legacy_name
    make_patient_report(patient_path, bundle)
    make_technical_report(technical_path, bundle)
    legacy_path.write_text(technical_path.read_text(encoding="utf-8"), encoding="utf-8")
    return {
        "evidence_report_patient": str(patient_path),
        "evidence_report_technical": str(technical_path),
        "evidence_report": str(legacy_path),
    }
