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
    "cross_platform_status": "WES/WGS 共同检出、低水平支持、检出能力不足或样本特异性等跨平台证据状态。",
    "cross_platform_multiplier": "跨平台 DNA 证据对排序分数的保守调整系数。",
    "source_chain_confidence_tier": "候选来源链置信度 C1–C4：评估事件、转录本/ORF、肽段与HLA结果能否形成可追溯链路；不同于R1–R4实验推荐等级。",
    "source_chain_orthogonal_status": "独立或跨模态确认状态；同一BAM上的多个caller仅算计算一致性，不等同正交确认。",
    "source_chain_requirement_statuses": "SNV、InDel、Fusion或Splice赛道专属要求的逐项状态，区分NOT_APPLICABLE、UNASSESSED、LOW_POWER、NEGATIVE与CONFLICT。",
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
    wes_qc: list[dict[str, str]] = field(default_factory=list)
    wes_wgs_coding_summary: list[dict[str, str]] = field(default_factory=list)
    targeted_pileup_summary: list[dict[str, str]] = field(default_factory=list)


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
        prov = _read_json_optional(root / "provenance.json")

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
        wes_qc=_read_optional(p("qc", "wes", "wes_qc.tsv") if root else None),
        wes_wgs_coding_summary=_read_optional(
            p("qc", "wes_wgs_coding_comparison", "wes_wgs_coding_summary.tsv") if root else None
        ),
        targeted_pileup_summary=(
            _read_optional(
                p("qc", "wes_wgs_coding_comparison", "targeted_pileup", "protein_altering", "discordant_targeted_pileup_summary.tsv")
                if root else None
            )
            or _read_optional(
                p("qc", "wes_wgs_coding_comparison", "targeted_pileup", "discordant_targeted_pileup_summary.tsv")
                if root else None
            )
        ),
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


def _cross_platform_counts(events: list[dict[str, str]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for event in events:
        status = str(event.get("cross_platform_status") or "")
        if status and status not in {"NOT_APPLICABLE", "UNASSESSED_NOT_IN_COMPARISON"}:
            counts[status] = counts.get(status, 0) + 1
    return counts


def _patient_platform_label(status: str) -> str:
    return {
        "CROSS_PLATFORM_PASS_CONCORDANT": "WES/WGS 共同检出",
        "ALT_PRESENT_BELOW_PASS_OR_CALLER_DIFFERENCE": "另一检测可见低水平支持",
        "COVERED_NO_ALT_SAMPLE_OR_ASSAY_DIFFERENCE": "仅当前样本/时间点明确",
        "OTHER_COVERED_BUT_LIMITED_POWER_AT_SOURCE_VAF": "另一检测能力不足，不能判阴性",
        "OTHER_LOW_OR_NO_COVERAGE": "另一检测覆盖不足",
        "SOURCE_INDEL_NOT_REPRODUCED_REASSEMBLY_REQUIRED": "复杂变异需局部重组复核",
        "SOURCE_PASS_NOT_REPRODUCED_BY_PILEUP": "源检测需重新复核",
        "SOURCE_WEAK_EXACT_PILEUP_SUPPORT": "源检测支持较弱",
        "NORMAL_SUPPORT_REVIEW": "正常血液也有支持，暂不推进",
        "NOT_APPLICABLE": "非 DNA 点突变事件",
        "UNASSESSED_NOT_IN_COMPARISON": "尚未完成跨平台评估",
    }.get(str(status or ""), "尚未完成跨平台评估")


def _patient_rna_label(peptide: Mapping[str, Any]) -> str:
    status = str(peptide.get("rna_support_status") or "")
    if status in {"RNA_ALT_SUPPORTED", "RNA_JUNCTION_SUPPORTED"}:
        return "RNA 已支持"
    if status == "RNA_ALT_NOT_DETECTED":
        return "RNA 未检出突变支持"
    if status in {"UNASSESSED", "", "RNA_ONLY_UNRESOLVED"}:
        return "RNA 证据未完整评估"
    return status.replace("_", " ")


def _patient_event_change(event: Mapping[str, Any]) -> str:
    change = str(event.get("combined_protein_change") or event.get("event_name") or event.get("consequence") or "")
    if ":p." in change:
        change = "p." + change.split(":p.", 1)[1]
    return change.replace("%3D", "=") or "蛋白改变待确认"


def _patient_fusion_interpretation(event: Mapping[str, Any]) -> str:
    gene = str(event.get("gene") or "")
    if gene == "EWSR1::WT1":
        return "DSRCT 的特征性驱动融合；仍需确认断点、阅读框及融合肽的真实加工呈递"
    if gene.startswith("HLA-"):
        return "HLA 高多态区域事件，优先排查比对或转录本拼接影响"
    if str(event.get("safety_status") or "") == "CAUTION":
        return "RNA junction 有支持，但正常组织背景或安全性证据仍需复核"
    return "候选融合事件；需用独立方法确认断点和阅读框"


def _metric_value(rows: list[dict[str, str]], metric: str, default: str = "未评估") -> str:
    for row in rows:
        if row.get("metric") == metric:
            return str(row.get("value") or default)
    return default


def _patient_pileup_category(category: str) -> str:
    return {
        "ALT_PRESENT_BELOW_PASS_OR_CALLER_DIFFERENCE": "另一平台存在低水平 ALT 或调用规则差异",
        "COVERED_NO_ALT_SAMPLE_OR_ASSAY_DIFFERENCE": "覆盖充分但未见 ALT，考虑样本/时间点差异",
        "NORMAL_SUPPORT_REVIEW": "正常血液存在支持，需排除胚系或技术因素",
        "OTHER_COVERED_BUT_LIMITED_POWER_AT_SOURCE_VAF": "有覆盖但按源 VAF 统计检出能力不足",
        "OTHER_LOW_OR_NO_COVERAGE": "另一平台低覆盖或无覆盖",
        "SOURCE_INDEL_NOT_REPRODUCED_REASSEMBLY_REQUIRED": "源平台 InDel 未由简单 pileup 复现，需局部重组",
        "SOURCE_PASS_NOT_REPRODUCED_BY_PILEUP": "源平台 PASS 位点未由 pileup 复现",
        "SOURCE_WEAK_EXACT_PILEUP_SUPPORT": "源平台精确 ALT 支持较弱",
        "WEAK_OR_ABSENT_REVIEW": "支持较弱或缺失，需人工复核",
    }.get(str(category or ""), str(category or "未分类"))


def _patient_pileup_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [
        {
            "source": row.get("source", ""),
            "category": _patient_pileup_category(row.get("category", "")),
            "count": row.get("count", ""),
        }
        for row in rows
    ]


def _patient_cancer_context(row: Mapping[str, Any]) -> str:
    status = str(row.get("cancer_driver_context") or "")
    types = str(row.get("cancer_gene_types") or "")
    if status == "DRIVER_CONTEXT":
        return f"癌症基因背景：{types}" if types else "具有癌症基因背景"
    if status == "LISTED_NO_DRIVER_CLASS":
        return "名单收录，但未归类为癌基因/抑癌基因"
    if status == "NOT_LISTED":
        return "未在当前癌症基因表中"
    return "未评估"


def _make_patient_report_legacy(path: str | Path, bundle: ReportBundle) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    val_map = _val_by_peptide(bundle.validation_rows)
    advancable = [ppt for ppt in bundle.peptides if str(ppt.get("final_priority", "")).upper() not in {"D", ""}]
    ranked_pool = advancable if advancable else bundle.peptides
    top = []
    seen_events = set()
    for candidate in ranked_pool:
        event_id = str(candidate.get("event_id") or candidate.get("peptide_id") or "")
        if event_id in seen_events:
            continue
        seen_events.add(event_id)
        top.append(candidate)
        if len(top) >= 10:
            break
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

    out.append("<div class='section'><h2>2. 样本与 HLA 背景</h2>")
    out.append("<ul class='compact'>")
    out.append(f"<li><b>样本编号：</b>{esc(bundle.sample_id or '未注明')}</li>")
    out.append(f"<li><b>分析场景：</b>{esc(bundle.entry_mode or '肿瘤新抗原筛选')}</li>")
    out.append("<li><b>数据层：</b>肿瘤 WES、肿瘤 WGS、配对血液正常样本、肿瘤 RNA/融合、HLA 分型，以及纯度、CNV 和 CCF 证据。</li>")
    out.append(f"<li><b>当前候选使用的 HLA-I 背景：</b>{esc(_top_hla_alleles(bundle.peptides))}</li>")
    out.append("<li><b>HLA LOH 背景：</b>当前检出的丢失信号位于 HLA-II（DQA1/DQB1）；未见 HLA-A/B/C 丢失影响当前 HLA-I 候选，受 HLA LOH 影响的候选肽数为 0。</li>")
    out.append(f"<li><b>评分方案：</b>{esc(bundle.profile.get('_profile_name', 'default'))}</li>")
    out.append("<li><b>临床样本关系：</b>精确取材日期、部位、治疗前后关系仍须以临床样本清单核实，本报告不据测序文件名推断。</li>")
    out.append("</ul><p class='small'>HLA 分型用于判断候选肽可能由哪些 HLA 分子呈递；它本身不证明肿瘤细胞已经加工并展示该肽段。</p></div>")

    out.append("<div class='section'><h2>3. 肿瘤是否具备抗原呈递条件</h2>")
    presentation_rows = [
        {"dimension": "MHC-I 抗原呈递", "conclusion": "部分保留，需谨慎解释", "meaning": "与当前 HLA-I 候选最直接相关；未显示通路完全缺失，但部分环节仍需实验确认"},
        {"dimension": "MHC-II 抗原呈递", "conclusion": "存在不确定因素", "meaning": "DQA1/DQB1 丢失及部分低表达信号提示 HLA-II 背景需要复核"},
        {"dimension": "IFNG 应答", "conclusion": "未见明确缺陷", "meaning": "现有输入支持该应答通路总体保留，但不等同于临床免疫治疗敏感"},
        {"dimension": "HLA-I 等位基因丢失", "conclusion": "未见影响当前候选", "meaning": "未发现 HLA-A/B/C 丢失影响当前纳入排序的 HLA-I 候选"},
        {"dimension": "HLA-II 等位基因丢失", "conclusion": "需要单独复核", "meaning": "当前信号位于 DQA1/DQB1，只进入 HLA-II 背景，不应降低当前 HLA-I 候选"},
        {"dimension": "证据完整性", "conclusion": "目前仍不完整", "meaning": "部分通路证据缺失；缺失不能解释为正常，也不能据此诊断免疫逃逸"},
    ]
    out.append(_table(presentation_rows, ["dimension", "conclusion", "meaning"]))
    out.append("<p><b>综合判断：</b>肿瘤具备部分 HLA-I 抗原呈递条件，并非“完全不能呈递”；但 TAP 等环节存在谨慎信号，且 APPM 输入完整度较低，因此应表述为<b>部分保留、仍需实验确认</b>，不能据此判断免疫治疗敏感或耐药。</p>")
    out.append("<p class='small'>本患者沟通版不展示模型分值；详细评分、状态字段和计算依据保留在科研技术版报告中。</p>")
    out.append("</div>")

    out.append("<div class='section'><h2>4. 关键发现（摘要）</h2><ul class='compact'>")
    priority_counts = {}
    for peptide in bundle.peptides:
        label = str(peptide.get("final_priority") or "未分级")
        priority_counts[label] = priority_counts.get(label, 0) + 1
    out.append(
        f"<li>共评估 <b>{len(bundle.events)}</b> 个候选事件及 <b>{len(bundle.peptides)}</b> 个肽段-HLA 组合；"
        f"其中 <b>{len(advancable)}</b> 个组合未被判为“不建议推进”。</li>"
    )
    out.append(
        "<li><b>当前分层：</b>"
        + "；".join(f"{esc(level)} 级 {count} 个" for level, count in sorted(priority_counts.items()))
        + "。分层数量较多不代表存在同等数量的独立治疗靶点。</li>"
    )
    if genes:
        out.append(f"<li>优先关注的基因包括：<b>{esc('、'.join(genes[:8]))}</b>。</li>")
    out.append(f"<li>{esc(_appm_patient_summary(bundle.appm_summary))}</li>")
    if bundle.immune_escape_summary:
        ies = bundle.immune_escape_summary[0]
        risk = str(ies.get("overall_immune_escape_risk") or "")
        if risk and risk.upper() not in {"LOW", "PASS", ""}:
            out.append(f"<li>免疫逃逸风险提示：<b>{esc(risk)}</b>（机制层面证据，非临床耐药结论）。</li>")
    platform_counts = _cross_platform_counts(bundle.events)
    if platform_counts:
        concordant = platform_counts.get("CROSS_PLATFORM_PASS_CONCORDANT", 0)
        review = sum(platform_counts.values()) - concordant
        out.append(
            f"<li>WES/WGS DNA 交叉复核：<b>{concordant}</b> 个事件由两平台共同检出；"
            f"<b>{review}</b> 个事件存在低水平、检出能力或样本时间点差异，已在排序中标记复核。</li>"
        )
    out.append("</ul></div>")

    dna_events = [event for event in bundle.events if str(event.get("mutation_source") or "") in {"SNV", "InDel"}]
    featured_genes = {"TP53", "KRAS"}
    dna_events.sort(key=lambda event: (
        str(event.get("gene") or "") in featured_genes,
        str(event.get("haplotype_status") or "") == "PHASED_CIS_COMBINED",
        str(event.get("rna_support_status") or "") == "RNA_ALT_SUPPORTED",
        str(event.get("cross_platform_status") or "") == "CROSS_PLATFORM_PASS_CONCORDANT",
        float(event.get("event_score") or 0),
    ), reverse=True)
    out.append("<div class='section'><h2>5. 主要 DNA 突变及 RNA/跨平台证据</h2>")
    out.append("<p class='small'>以下为结合跨平台、RNA 与事件评分选出的主要 DNA 变异；它们不是临床用药清单，也不代表均为肿瘤驱动事件。</p>")
    mutation_rows = []
    for event in dna_events[:12]:
        mutation_rows.append({
            "gene": event.get("gene", ""),
            "cancer_context": _patient_cancer_context(event),
            "protein_change": _patient_event_change(event),
            "type": event.get("mutation_source", ""),
            "rna_evidence": _patient_rna_label(event),
            "wes_wgs_evidence": _patient_platform_label(str(event.get("cross_platform_status") or "")),
            "interpretation": (
                "两相邻变异已完成同一单倍型重构" if event.get("haplotype_status") == "PHASED_CIS_COMBINED"
                else "需结合病理、克隆性和实验验证判断作用"
            ),
        })
    out.append(_table(mutation_rows, ["gene", "cancer_context", "protein_change", "type", "rna_evidence", "wes_wgs_evidence", "interpretation"]))
    out.append("<p><b>重点提示：</b>TP53 等共同检出且有 RNA 支持的变异可信度相对更高；KRAS 等仅在一个肿瘤文库明确检出的变异应按样本/时间点特异结果解释；TBR1 相邻变异按已重构的同一单倍型解读。</p>")
    out.append("</div>")

    fusion_events = [event for event in bundle.events if str(event.get("event_type") or "") == "Fusion"]
    fusion_events.sort(key=lambda event: (str(event.get("gene") or "") == "EWSR1::WT1", float(event.get("event_score") or 0)), reverse=True)
    out.append("<div class='section'><h2>6. 融合基因事件及 DSRCT 背景解释</h2>")
    out.append("<p><b>EWSR1::WT1</b> 是 DSRCT 的标志性融合，本样本存在 RNA junction 支持；但由融合产生的新抗原仍需独立验证阅读框、异常 junction 和实际 HLA 呈递。</p>")
    fusion_rows = []
    for event in fusion_events[:10]:
        fusion_rows.append({
            "fusion": event.get("gene", ""),
            "cancer_context": _patient_cancer_context(event),
            "junction_reads": event.get("rna_junction_reads", ""),
            "rna_status": _patient_rna_label(event),
            "expression": event.get("event_expression", ""),
            "safety": event.get("safety_status", ""),
            "interpretation": _patient_fusion_interpretation(event),
        })
    out.append(_table(fusion_rows, ["fusion", "cancer_context", "junction_reads", "rna_status", "expression", "safety", "interpretation"]))
    out.append("<p class='small'>肝脏高表达基因、HLA 区域或仅有少量 junction reads 的融合可能包含 read-through、正常背景或比对伪影，其证据等级不能与 EWSR1::WT1 等驱动融合等同。</p>")
    out.append("</div>")

    event_by_gene: dict[str, dict[str, str]] = {}
    for event in bundle.events:
        gene = str(event.get("gene") or "")
        if gene and gene not in event_by_gene:
            event_by_gene[gene] = event
    discussable = []
    discussion_notes = {
        "EWSR1::WT1": "DSRCT 标志性驱动融合，疾病相关性最明确；新抗原价值仍取决于 junction 阅读框、加工呈递和功能实验",
        "TP53": "WES/WGS 共同检出并有 RNA 支持，事件真实性较强；仍需验证突变肽相对 WT 的特异性",
        "TBR1": "两个相邻变异已按同一单倍型重构；应只保留少量 combined-mutant 肽，避免重叠窗口重复占位",
        "KRAS": "经典热点但仅在一个肿瘤文库明确；需先确认取材/时间点和独立 DNA/RNA 支持",
        "TOMM34": "跨平台与 RNA 证据较完整，可作为非驱动但可验证的突变肽事件讨论",
        "ACLY": "跨平台与 RNA 证据较完整，可纳入突变肽对照验证候选组",
    }
    for gene in ("EWSR1::WT1", "TP53", "TBR1", "TOMM34", "ACLY", "KRAS"):
        event = event_by_gene.get(gene)
        if not event:
            continue
        discussable.append({
            "event": gene,
            "change": _patient_event_change(event),
            "rna": _patient_rna_label(event),
            "dna": _patient_platform_label(str(event.get("cross_platform_status") or "")),
            "why_discuss": discussion_notes[gene],
        })
    out.append("<div class='section'><h2>7. 目前最值得讨论的候选事件</h2>")
    out.append("<p>这里的“值得讨论”综合考虑疾病生物学、DNA/RNA 可复现性、单倍型、HLA 呈递预测和实验可行性，不等同于自动分数最高。</p>")
    out.append(_table(discussable, ["event", "change", "rna", "dna", "why_discuss"]))
    out.append("</div>")

    manual_statuses = {
        "COVERED_NO_ALT_SAMPLE_OR_ASSAY_DIFFERENCE": "可能反映不同取材/时间点的真实异质性",
        "SOURCE_INDEL_NOT_REPRODUCED_REASSEMBLY_REQUIRED": "复杂 InDel 可能无法由简单 pileup 复现",
        "SOURCE_PASS_NOT_REPRODUCED_BY_PILEUP": "源 caller 给出 PASS，但原始 reads 仍需人工核验",
        "SOURCE_WEAK_EXACT_PILEUP_SUPPORT": "存在少量精确 ALT reads，但证据偏弱",
        "NORMAL_SUPPORT_REVIEW": "正常样本也有支持，需先排除胚系或技术伪影",
    }
    manual_rows = []
    for event in bundle.events:
        status = str(event.get("cross_platform_status") or "")
        gene = str(event.get("gene") or "")
        if status not in manual_statuses:
            continue
        manual_rows.append({
            "event": gene,
            "change": _patient_event_change(event),
            "retain_reason": manual_statuses[status],
            "do_not_auto_advance": _patient_platform_label(status),
            "required_review": "IGV/局部重组、独立测序或匹配时间点复核",
        })
        if len(manual_rows) >= 12:
            break
    ews = event_by_gene.get("EWSR1::WT1")
    if ews:
        manual_rows.insert(0, {
            "event": "EWSR1::WT1",
            "change": _patient_event_change(ews),
            "retain_reason": "DSRCT 标志性驱动融合，生物学优先级高于一般自动分数",
            "do_not_auto_advance": "融合肽呈递、安全性和 WT/正常背景尚未完成",
            "required_review": "断点/阅读框、junction RNA、长肽/minigene 和功能验证",
        })
    out.append("<div class='section'><h2>8. 需要人工保留、但不应按自动分数直接推进的事件</h2>")
    out.append(_table(manual_rows, ["event", "change", "retain_reason", "do_not_auto_advance", "required_review"]))
    out.append("<p class='small'>人工保留表示不应因单一分数或暂时缺证而删除；它也不表示可以绕过证据补充直接进入治疗设计。</p></div>")

    tier_rows = [
        {"tier": "研究层 1A", "scope": "WES/WGS 共同检出、RNA 支持、呈递证据较好的 SNV/InDel", "examples": "TP53、TOMM34、ACLY 及其他满足条件事件", "action": "优先做 MT/WT 成对肽验证"},
        {"tier": "研究层 1B", "scope": "疾病驱动融合及异常 junction", "examples": "EWSR1::WT1", "action": "单独成组，优先长肽/minigene，不以短肽分数替代加工验证"},
        {"tier": "研究层 2", "scope": "样本/时间点特异或另一平台低水平支持", "examples": "KRAS 等", "action": "先确认目标取材中的 DNA/RNA 存在，再决定是否进入免疫学实验"},
        {"tier": "人工复核层", "scope": "复杂 InDel、源检测未复现、弱支持", "examples": "按 targeted pileup 标记的事件", "action": "IGV、局部组装或独立测序后重新评分"},
        {"tier": "暂缓层", "scope": "正常样本支持、明显安全性风险或总体为 D", "examples": "AXDND1 等正常支持事件", "action": "不进入首批实验；先排除胚系、正常组织表达和交叉反应"},
    ]
    out.append("<div class='section'><h2>9. 最值得关注的候选分层</h2>")
    out.append(_table(tier_rows, ["tier", "scope", "examples", "action"]))
    out.append("</div>")

    out.append("<div class='section'><h2>10. 优先候选肽段（按事件去重 Top 10）</h2>")
    out.append("<p class='small'>下表为计算排序靠前的候选，不代表已证实可诱导抗肿瘤免疫反应。</p>")
    headers = ["rank", "gene", "cancer_context", "variant_type", "hla", "source_chain", "priority", "rna_evidence", "dna_evidence", "meaning", "next_step"]
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
            "cancer_context": _patient_cancer_context(ppt),
            "variant_type": _patient_consequence_label(ppt),
            "hla": ppt.get("hla_allele", ""),
            "source_chain": ppt.get("source_chain_confidence_tier", "未评估"),
            "priority": _patient_priority_label(str(ppt.get("final_priority") or "")),
            "rna_evidence": _patient_rna_label(ppt),
            "dna_evidence": _patient_platform_label(str(ppt.get("cross_platform_status") or "")),
            "meaning": meaning,
            "next_step": next_step,
        })
    out.append(_table(rows, headers))
    out.append("</div>")

    out.append("<div class='section'><h2>11. WES 与 WGS 蛋白改变型 SNV/InDel 对比</h2>")
    if bundle.wes_wgs_coding_summary:
        coding_rows = [
            {"metric": "WES 蛋白改变型 SNV/InDel", "value": _metric_value(bundle.wes_wgs_coding_summary, "protein_altering_wes")},
            {"metric": "WGS 蛋白改变型 SNV/InDel", "value": _metric_value(bundle.wes_wgs_coding_summary, "protein_altering_wgs")},
            {"metric": "两平台共同检出", "value": _metric_value(bundle.wes_wgs_coding_summary, "protein_altering_common")},
            {"metric": "WES-only", "value": _metric_value(bundle.wes_wgs_coding_summary, "protein_altering_wes_only")},
            {"metric": "WGS-only", "value": _metric_value(bundle.wes_wgs_coding_summary, "protein_altering_wgs_only")},
            {"metric": "共同位点 VAF 相关性", "value": _metric_value(bundle.wes_wgs_coding_summary, "common_af_pearson")},
        ]
        out.append(_table(coding_rows, ["metric", "value"]))
    else:
        out.append("<p>WES/WGS coding 区域汇总尚未加载。</p>")
    if bundle.targeted_pileup_summary:
        out.append("<h3>差异位点回查结果</h3>")
        out.append(_table(_patient_pileup_rows(bundle.targeted_pileup_summary), ["source", "category", "count"]))
    out.append("<p class='small'>本节只统计 VEP 明确标记为 missense、frameshift、inframe insertion/deletion、start/stop lost 或 stop gained 等会改变蛋白序列的 SNV/InDel；不包括 synonymous、非编码、仅 splice-region 但没有明确蛋白改变的记录，也不包括 fusion/SV。</p>")
    out.append("<p>低重合度不能简单解释为某个平台错误。回查显示主要来源包括低 VAF 下检出能力不足、另一平台存在但未通过 PASS、不同肿瘤文库/时间点差异，以及长 InDel 需要局部重组复核。只有在源 BAM 可复现、另一 BAM 覆盖和统计检出能力充分且无 ALT 时，才视为较强的样本差异证据。</p>")
    out.append("</div>")

    out.append("<div class='section'><h2>12. 为什么当前没有直接进入高优先级的候选</h2><ul class='compact'>")
    if not any(str(p.get("final_priority") or "") in {"A", "B", "B_CAUTION"} for p in bundle.peptides):
        out.append("<li>当前没有 A/B 级候选。主要原因不是“完全没有候选”，而是正常组织安全性、RNA 支持、突变特异性或样本间一致性仍需补证。</li>")
    out.append("<li>C_CAUTION 表示候选具有一定计算证据，但在进入首批实验或治疗设计前必须完成针对性复核。</li>")
    out.append("<li>D 级表示目前不建议推进，常见原因包括安全性风险、正常样本支持、呈递证据不足或关键证据缺失。</li>")
    out.append("<li>同一事件可产生多个长度和多个 HLA 限制性肽段，因此肽段组合数远高于独立变异事件数。</li>")
    out.append("</ul></div>")

    out.append("<div class='section'><h2>13. 当前证据缺口与样本差异</h2><ul class='compact'>")
    if platform_counts:
        out.append(f"<li>WES/WGS 共同检出的事件：<b>{platform_counts.get('CROSS_PLATFORM_PASS_CONCORDANT', 0)}</b> 个。</li>")
        out.append(
            f"<li>另一检测可见低水平 ALT 支持：<b>{platform_counts.get('ALT_PRESENT_BELOW_PASS_OR_CALLER_DIFFERENCE', 0)}</b> 个；"
            f"因检出能力不足不能判阴性：<b>{platform_counts.get('OTHER_COVERED_BUT_LIMITED_POWER_AT_SOURCE_VAF', 0)}</b> 个。</li>"
        )
        out.append(
            f"<li>覆盖充分但呈现样本/时间点差异：<b>{platform_counts.get('COVERED_NO_ALT_SAMPLE_OR_ASSAY_DIFFERENCE', 0)}</b> 个。"
            "这提示不同取材、肿瘤异质性或低纯度影响，不能把一个时间点的结果概括为所有肿瘤组织。</li>"
        )
        source_review = sum(platform_counts.get(key, 0) for key in (
            "SOURCE_INDEL_NOT_REPRODUCED_REASSEMBLY_REQUIRED",
            "SOURCE_PASS_NOT_REPRODUCED_BY_PILEUP",
            "SOURCE_WEAK_EXACT_PILEUP_SUPPORT",
        ))
        out.append(f"<li>源检测仍需局部重组、IGV 或重复测序复核：<b>{source_review}</b> 个事件。</li>")
    out.append("<li>RNA 表达量只能说明基因被表达，不能替代突变 RNA reads 或融合/剪接 junction 的直接支持。</li>")
    out.append("<li>精确取材日期、部位及治疗前后关系需与临床样本记录核对；本报告仅描述已测序文库。</li>")
    out.append("</ul></div>")

    out.append("<div class='section'><h2>14. 建议的下一步验证顺序</h2><ol class='compact'>")
    out.append("<li><b>先核对样本：</b>确认 WES、WGS、RNA 的取材部位、日期、治疗前后关系及肿瘤含量，避免把时间点差异误作技术失败。</li>")
    out.append("<li><b>确认事件真实性：</b>对样本特异、复杂 InDel 和弱支持事件做 IGV、局部组装、靶向深测或独立 PCR；TBR1 保留 read-backed phasing 结论。</li>")
    out.append("<li><b>确认突变转录：</b>SNV/InDel 检查 RNA alt reads/RNA VAF；融合和剪接检查 junction reads、阅读框及异常转录本。</li>")
    out.append("<li><b>确认突变特异性：</b>比较 MT 与 WT 的 HLA 结合、呈递和免疫原性；WT 相当或更强者不进入首批。</li>")
    out.append("<li><b>确认加工呈递：</b>短肽候选做 MT/WT 成对验证；移码、剪接和融合优先长肽或 minigene，并在条件允许时做免疫肽组学。</li>")
    out.append("<li><b>确认免疫功能与安全性：</b>再开展 ELISpot、四聚体/多聚体、细胞毒实验，并补正常组织、HSPC、自身肽和脱靶复核。</li>")
    out.append("</ol></div>")

    experiment_rows = [
        {"group": "A：SNV/InDel 短肽组", "contents": "每个事件 1–2 个最佳 peptide-HLA 组合", "controls": "对应 WT 肽、无关肽、阳性刺激", "purpose": "验证突变特异性和 T 细胞识别"},
        {"group": "B：移码/剪接长肽组", "contents": "覆盖新生尾部或异常 junction 的长肽/minigene", "controls": "正常转录本/WT 构建", "purpose": "验证真实加工而非仅短肽结合"},
        {"group": "C：融合专项组", "contents": "EWSR1::WT1 单独成组，其他融合分开", "controls": "断点阴性、WT 两端序列、无关融合", "purpose": "验证断点、阅读框、呈递和 DSRCT 特异背景"},
        {"group": "D：人工保留复核组", "contents": "样本特异、复杂 InDel、弱支持事件", "controls": "另一平台/另一时间点、正常样本", "purpose": "先解决事件真实性，不与已确认候选混合解释"},
    ]
    out.append("<div class='section'><h2>15. 建议的实验候选组织方式</h2>")
    out.append(_table(experiment_rows, ["group", "contents", "controls", "purpose"]))
    out.append("<p class='small'>每组应预先定义入组条件、排除条件和重复数；同一事件的高度重叠肽应去冗余，避免一个事件因窗口数量多而在实验中被过度代表。</p></div>")

    out.append("<div class='section'><h2>16. 面向患者的核心结论</h2>")
    out.append("<p>本次分析发现了若干值得继续研究的候选，尤其包括与 DSRCT 密切相关的 <b>EWSR1::WT1</b> 融合，以及部分在 DNA、RNA 和 HLA 预测层面得到支持的突变事件。肿瘤的 HLA-I 抗原呈递能力看起来是<b>部分保留</b>的，因此继续做新抗原实验验证具有研究依据。</p>")
    out.append("<p>但目前没有候选达到“仅凭计算结果即可用于治疗”的证据标准。部分事件在 WES/WGS、不同取材或正常样本之间存在差异，正常组织安全性和真实加工呈递也尚未完整验证。当前最重要的下一步是按分层进行事件确认、RNA/断点验证、MT-WT 对照和 T 细胞功能实验，而不是直接按自动排名选择治疗方案。</p>")
    out.append("</div>")

    out.append("<div class='warn'><h2>17. 数据来源与解释边界</h2><ul class='compact'>")
    out.append("<li><b>数据来源：</b>肿瘤 WES、肿瘤 WGS、配对血液正常样本、肿瘤 RNA/融合结果、HLA 分型、纯度/CNV/CCF，以及呈递、免疫原性、APPM/逃逸与正常组织安全性参考。</li>")
    out.append("<li><b>跨平台边界：</b>WES 与 WGS 可能来自不同肿瘤文库或时间点；本报告展示的蛋白改变型 SNV/InDel 差异同时受捕获范围、深度、低纯度、异质性、caller、VEP 转录本选择与局部组装影响。完整 coding/splice PASS 全集保留在技术 QC 中。</li>")
    out.append("<li><b>融合边界：</b>检测到驱动融合不等于其 junction 肽一定被加工、呈递或被 T 细胞识别。</li>")
    out.append("<li>本分析为<strong>计算机辅助筛选</strong>，预测结合亲和力不等于体内呈递，更不等于临床疗效。</li>")
    out.append("<li>APPM、CCF、安全性与免疫逃逸评估依赖输入数据完整度；缺失数据不等于“无风险”。</li>")
    out.append("<li>本报告<strong>不包含</strong>原始测序质控、文件路径或生信命令细节；技术细节见科研技术版报告。</li>")
    out.append("<li>不得将本报告直接用于患者诊断、预后判断或个体化治疗处方。</li>")
    out.append("</ul></div>")
    out.append("</body></html>")
    p.write_text("\n".join(out), encoding="utf-8")


R_GRADE_PATIENT = {
    "R1": ("第一批实验优先", "关键证据较完整，可优先进入研究性实验验证。"),
    "R2": ("值得推进", "总体证据较好，但仍有一项或少量谨慎因素需要补充。"),
    "R3": ("优先补证据", "先补 RNA、事件真实性、安全性或呈递证据，再决定是否进入免疫学实验。"),
    "R4": ("当前暂不推进", "存在硬失败、明确风险、证据明显不足或呈递一致弱等原因。"),
}


def _patient_grade(row: Mapping[str, Any]) -> str:
    for key in ("pipeline_r_grade", "evidence_grade", "r_grade", "final_priority"):
        value = str(row.get(key) or "").strip().upper()
        if value:
            return value
    return "UNASSESSED"


def _patient_track(row: Mapping[str, Any]) -> str:
    text = " ".join(str(row.get(key) or "") for key in (
        "event_type", "mutation_source", "peptide_consequence", "source_type", "event_kind"
    )).lower()
    if "fusion" in text:
        return "Fusion"
    if "splice" in text or "junction" in text or "exon" in text:
        return "Splice"
    if "sv" in text or "structural" in text:
        return "DNA SV"
    if any(token in text for token in ("indel", "frameshift", "insertion", "deletion")):
        return "InDel"
    if any(token in text for token in ("snv", "missense", "substitution")):
        return "SNV"
    return "Other"


def _patient_representatives(rows: list[dict[str, str]], limit: int, track: str | None = None) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []
    seen: set[str] = set()
    for row in rows:
        if track and _patient_track(row) != track:
            continue
        event_id = str(row.get("event_id") or row.get("event_name") or row.get("peptide_id") or "")
        key = event_id or f"{row.get('gene', '')}|{row.get('peptide', '')}"
        if key in seen:
            continue
        seen.add(key)
        selected.append(row)
        if len(selected) >= limit:
            break
    return selected


def _patient_evidence_summary(row: Mapping[str, Any]) -> str:
    evidence: list[str] = []
    if str(row.get("rna_support_status") or "") in {"RNA_ALT_SUPPORTED", "RNA_JUNCTION_SUPPORTED", "RNA_CONFIRMED"}:
        evidence.append("RNA支持")
    if str(row.get("cross_platform_status") or "") == "CROSS_PLATFORM_PASS_CONCORDANT":
        evidence.append("DNA跨平台一致")
    presentation = str(row.get("presentation_consensus_state") or row.get("presentation_evidence_grade") or "")
    if presentation:
        evidence.append(f"呈递证据 {presentation}")
    specificity = str(row.get("mutant_specificity_status") or row.get("mutant_specificity_state") or "")
    if specificity:
        evidence.append(f"MT/WT {specificity}")
    source_chain = str(row.get("source_chain_confidence_tier") or "")
    if source_chain:
        evidence.append(f"来源链 {source_chain}")
    return "；".join(evidence) if evidence else "现有证据尚未完整归一化"


def _patient_limitation(row: Mapping[str, Any]) -> str:
    reasons = str(row.get("hard_failure_codes") or row.get("priority_cap_reason_codes") or row.get("reason_codes") or "")
    if reasons:
        return reasons.replace("|", "；")
    missing: list[str] = []
    if str(row.get("rna_support_status") or "") in {"", "UNASSESSED", "RNA_ONLY_UNRESOLVED"}:
        missing.append("RNA证据未评估")
    if str(row.get("safety_status") or "") in {"", "UNASSESSED", "SAFETY_PARTIAL"}:
        missing.append("安全性证据不完整")
    if str(row.get("hla_loh_status") or row.get("escape_status") or "") in {"", "UNASSESSED"}:
        missing.append("HLA LOH未完整评估")
    return "；".join(missing) if missing else "未见明确限制；仍需实验验证"


def _patient_validation(row: Mapping[str, Any], val_map: Mapping[str, Mapping[str, str]]) -> str:
    val = val_map.get(str(row.get("peptide_id") or ""), {})
    explicit = str(val.get("validation_strategy") or val.get("recommended_assay") or row.get("recommended_validation") or row.get("recommended_use") or "")
    if explicit:
        return explicit
    track = _patient_track(row)
    if track == "SNV":
        return "MT/WT成对短肽与ELISpot/多聚体"
    if track == "InDel":
        return "新生尾部长肽或minigene"
    if track == "Fusion":
        return "RT-PCR/Sanger确认断点，再做融合junction长肽或minigene"
    if track == "Splice":
        return "targeted RNA确认junction，再做异常junction长肽或minigene"
    return "先确认事件真实性，再设计功能实验"


def _patient_hla_rows(peptides: list[dict[str, str]]) -> list[dict[str, str]]:
    loci: dict[str, list[str]] = {"A": [], "B": [], "C": []}
    for row in peptides:
        allele = str(row.get("hla_allele") or "").replace("HLA-", "")
        locus = allele.split("*", 1)[0]
        if locus in loci and allele and allele not in loci[locus]:
            loci[locus].append(allele)
    return [{"位点": key, "推荐等位基因": " / ".join(values[:2]) if values else "未评估", "用途": "限制性HLA-I呈递背景"} for key, values in loci.items()]


def _patient_appm_rows(bundle: ReportBundle) -> list[dict[str, str]]:
    summary = bundle.appm_summary
    dimensions = [
        ("MHC-I核心", "mhc_i_integrity_status"),
        ("MHC-II背景", "mhc_ii_integrity_status"),
        ("IFNG/JAK-STAT", "ifng_response_status"),
        ("APPM证据完整度", "appm_evidence_completeness"),
    ]
    rows = []
    for label, key in dimensions:
        status = str(summary.get(key) or "UNASSESSED")
        rows.append({"维度": label, "结果": status, "通俗解释": "未评估不等于正常或阴性" if status == "UNASSESSED" else "用于判断抗原加工呈递条件，不能单独预测临床疗效"})
    lost = sorted({str(row.get("hla_allele") or row.get("allele") or "") for row in bundle.peptide_escape_flags if str(row.get("restricting_hla_lost") or "").lower() in {"yes", "true", "1"}})
    rows.append({"维度": "限制性HLA-I LOH", "结果": "LOST: " + ", ".join(lost) if lost else "未见影响当前候选或未评估", "通俗解释": "仅HLA-A/B/C丢失可直接影响相应HLA-I候选"})
    return rows


def _patient_tool_rows(provenance: Mapping[str, Any]) -> list[dict[str, str]]:
    tools = provenance.get("tools") or {}
    rows: list[dict[str, str]] = []
    if isinstance(tools, Mapping):
        for name, record in tools.items():
            if not isinstance(record, Mapping):
                continue
            rows.append({"流程/工具": str(name), "版本": str(record.get("version") or "未记录"), "状态": str(record.get("status") or "未评估"), "作用": str(record.get("purpose") or record.get("mode") or "证据生成")})
    return rows


def make_patient_report(path: str | Path, bundle: ReportBundle) -> None:
    """Write the template-aligned, sample-agnostic patient HTML report."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    val_map = _val_by_peptide(bundle.validation_rows)
    ranked = list(bundle.peptides)
    top = _patient_representatives(ranked, 10)
    grade_counts: dict[str, int] = {}
    for row in ranked:
        grade = _patient_grade(row)
        grade_counts[grade] = grade_counts.get(grade, 0) + 1
    track_counts: dict[str, int] = {}
    event_seen: set[str] = set()
    for row in bundle.events or ranked:
        event_id = str(row.get("event_id") or row.get("event_name") or row.get("peptide_id") or "")
        if event_id and event_id in event_seen:
            continue
        if event_id:
            event_seen.add(event_id)
        track = _patient_track(row)
        track_counts[track] = track_counts.get(track, 0) + 1

    out = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        f"<title>肿瘤新抗原筛选报告—{esc(bundle.sample_id)}</title>", REPORT_CSS,
        "</head><body class='patient'><h1>肿瘤新抗原筛选报告</h1><p class='small'>患者沟通版</p>",
        "<div class='info'><b>阅读提示：</b>本报告先给出结论与候选分层，再说明样本、HLA、变异、肽段和验证建议。"
        "R1–R4是研究性证据等级，不是疗效等级。</div>",
        "<div class='warn'><b>重要说明：</b>本报告为研究性计算筛选，不能替代临床诊断或治疗决策。预测候选不等于体内真实呈递、不等于T细胞能够识别，也不等于确定治疗方案；所有候选均需进一步实验和临床专业判断。</div>",
    ]
    if str(bundle.provenance.get("report_role") or "").lower() == "pipeline_snapshot":
        out.append(
            "<div class='warn'><b>报告状态：</b>这是 Pipeline 运行阶段的结果快照，不是最终患者报告。"
            "请在 Evidence-consensus 完成后运行 open-neo-review，并以其 reports/patient_report.html 为准。</div>"
        )

    out.append("<div class='section'><h2>1. 报告摘要</h2>")
    summary = [
        {"项目": "样本编号", "结果": bundle.sample_id or "未注明", "说明": "以运行清单为准"},
        {"项目": "分析入口", "结果": bundle.entry_mode or "未注明", "说明": "可能为VCF、融合、剪接或联合入口"},
        {"项目": "独立事件", "结果": str(len(event_seen) or len(bundle.events)), "说明": "按event_id去重"},
        {"项目": "肽段-HLA组合", "结果": str(len(ranked)), "说明": "同一事件可产生多个组合"},
        {"项目": "评分配置", "结果": str(bundle.profile.get("_profile_name") or "未记录"), "说明": "研究性规则版本"},
    ]
    out.append(_table(summary, ["项目", "结果", "说明"]))
    grade_rows = [{"等级": grade, "数量": count, "含义": R_GRADE_PATIENT.get(grade, ("旧版/其他分层", "见技术报告"))[0], "下一步": R_GRADE_PATIENT.get(grade, ("复核", "按证据缺口处理"))[1]} for grade, count in sorted(grade_counts.items())]
    out.append("<h3>关键结论与R1–R4分层</h3>" + _table(grade_rows, ["等级", "数量", "含义", "下一步"]))
    out.append("<p>候选选择同时考虑事件真实性、RNA支持、HLA呈递、MT/WT突变特异性、克隆性、HLA/APPM和安全性；缺失证据统一视为未评估，不作为阴性结论。</p></div>")

    out.append("<div class='section'><h2>2. 患者样本与测序数据</h2>")
    qc_rows = [
        {"项目": "肿瘤/正常配对", "结果": str(bundle.provenance.get("pairing_status") or "未评估"), "解释": "需由样本清单和指纹核验"},
        {"项目": "肿瘤纯度/倍性", "结果": str(bundle.provenance.get("purity_ploidy") or "未评估"), "解释": "用于CNV、CCF和LOH解释"},
        {"项目": "肿瘤DNA深度", "结果": str(bundle.provenance.get("tumor_dna_depth") or "未评估"), "解释": "低深度会降低变异检出能力"},
        {"项目": "正常DNA深度", "结果": str(bundle.provenance.get("normal_dna_depth") or "未评估"), "解释": "用于排除胚系和正常支持"},
        {"项目": "RNA质量/覆盖", "结果": str(bundle.provenance.get("rna_qc_status") or "未评估"), "解释": "决定RNA未检出的解释强度"},
        {"项目": "参考版本", "结果": str(bundle.provenance.get("genome_build") or "未记录"), "解释": "FASTA、GTF、VEP和坐标必须一致"},
    ]
    out.append(_table(qc_rows, ["项目", "结果", "解释"]))
    out.append("<p class='small'>未提供的QC项目保持“未评估”，不会自动写成正常。</p></div>")

    out.append("<div class='section'><h2>3. HLA分型与抗原呈递条件</h2>")
    out.append("<h3>推荐HLA-I背景</h3>" + _table(_patient_hla_rows(ranked), ["位点", "推荐等位基因", "用途"]))
    out.append("<h3>抗原加工呈递与免疫逃逸</h3>" + _table(_patient_appm_rows(bundle), ["维度", "结果", "通俗解释"]))
    out.append("<p>这些结果说明候选是否具备被加工和呈递的条件，但不能单独判断免疫治疗敏感、耐药或患者获益。</p></div>")

    out.append("<div class='section'><h2>4. 重点变异事件（按类型、按事件去重）</h2>")
    overall = [{"事件类型": track, "事件数": count, "主要复核重点": {"SNV": "DNA深度/VAF、RNA alt、MT/WT", "InDel": "局部重比对、阅读框、NMD和phasing", "Fusion": "精确断点、junction reads、frame和正常read-through", "Splice": "精确junction、PSI/reads、正常isoform和ORF", "DNA SV": "断点与异常转录本"}.get(track, "事件真实性和证据完整性")} for track, count in sorted(track_counts.items())]
    out.append(_table(overall, ["事件类型", "事件数", "主要复核重点"]))
    for track in ("SNV", "InDel", "Fusion", "Splice", "DNA SV"):
        representatives = _patient_representatives(ranked, 5, track)
        if not representatives:
            continue
        rows = []
        for rank, row in enumerate(representatives, 1):
            rows.append({"排名": rank, "基因/事件": row.get("gene") or row.get("event_name") or row.get("event_id", ""), "改变": _patient_event_change(row), "来源链": row.get("source_chain_confidence_tier", "未评估"), "R等级": _patient_grade(row), "关键证据": _patient_evidence_summary(row), "主要限制": _patient_limitation(row)})
        out.append(f"<h3>{esc(track)} Top 5</h3>" + _table(rows, ["排名", "基因/事件", "改变", "来源链", "R等级", "关键证据", "主要限制"]))
    out.append("<p class='small'>不同事件赛道的证据结构不同，Top5用于赛道内审阅，不应仅凭序号直接跨赛道比较。</p></div>")

    out.append("<div class='section'><h2>5. 候选肽段Top 10（跨赛道、按事件去重）</h2>")
    peptide_rows = []
    for rank, row in enumerate(top, 1):
        peptide_rows.append({"排名": rank, "基因": row.get("gene", ""), "类型": _patient_track(row), "肽段-HLA": f"{row.get('peptide', '')} / {row.get('hla_allele', '')}", "等级": _patient_grade(row), "关键证据": _patient_evidence_summary(row), "主要限制": _patient_limitation(row), "建议实验": _patient_validation(row, val_map), "审阅状态": row.get("review_status", "待实验复核")})
    out.append(_table(peptide_rows, ["排名", "基因", "类型", "肽段-HLA", "等级", "关键证据", "主要限制", "建议实验", "审阅状态"]))
    out.append("<p class='small'>排序顺序为：R1–R4证据等级 → 同赛道Pareto排序 → 确定性tie-break → 事件去重。旧版weighted baseline仅作为并行比较，不覆盖证据共识结论。</p></div>")

    out.append("<div class='section'><h2>6. Top候选解读与实验建议</h2>")
    interpretation_rows = []
    for row in top[:5]:
        interpretation_rows.append({"候选": f"{row.get('gene', '')} | {row.get('peptide', '')} | {row.get('hla_allele', '')}", "为什么值得关注": _patient_evidence_summary(row), "当前不确定性": _patient_limitation(row), "建议下一步": _patient_validation(row, val_map)})
    out.append(_table(interpretation_rows, ["候选", "为什么值得关注", "当前不确定性", "建议下一步"]))
    out.append("<ol><li>先确认事件和异常转录本真实性；</li><li>再补RNA alt/VAF或精确junction证据；</li><li>完成MT/WT、正常背景和限制性HLA复核；</li><li>最后开展短肽、长肽、minigene及T细胞功能实验。</li></ol></div>")

    out.append("<div class='section'><h2>7. 分析方法与工具状态</h2>")
    method_rows = [
        {"阶段": "事件输入与质控", "方法": "VCF/Fusion/Splice/SV标准化与来源链检查", "状态": "按当前输入执行"},
        {"阶段": "肽段构建", "方法": "SNV MT/WT；InDel novel tail；Fusion/Splice精确junction ORF", "状态": "无可追溯ORF时不形成高等级证据"},
        {"阶段": "呈递预测", "方法": "核心结合/呈递、稳定性和免疫原性样模型分组汇总", "状态": "工具缺失记为UNASSESSED"},
        {"阶段": "综合排序", "方法": "证据状态、hard fail、priority cap、R1–R4、Pareto和事件去重", "状态": "研究性规则"},
    ]
    out.append(_table(method_rows, ["阶段", "方法", "状态"]))
    tool_rows = _patient_tool_rows(bundle.provenance)
    if tool_rows:
        out.append("<h3>运行工具</h3>" + _table(tool_rows, ["流程/工具", "版本", "状态", "作用"]))
    out.append("</div>")

    out.append("<div class='section'><h2>8. 局限性与总体结论</h2><ul>")
    out.append("<li>DNA/RNA覆盖不足时的未检出属于低检出能力或未评估，不是阴性。</li><li>计算呈递不等于体内真实呈递，体内呈递也不等于T细胞能够识别。</li><li>正常表达、HSPC、正常蛋白组、正常ligandome和正常junction不完整时，安全性必须标记为PARTIAL。</li><li>患者版不输出用药建议、疗效承诺或已确认新抗原结论。</li></ul>")
    out.append("<p><b>总体结论：</b>本次结果提供了可追溯的研究候选和分层验证顺序。应优先围绕事件真实性、RNA支持、MT/WT特异性、HLA/APPM和正常背景补证，再决定首批实验集合。</p></div>")

    grade_rows = [{"等级": grade, "定义": value[0], "核心要求/处理": value[1]} for grade, value in R_GRADE_PATIENT.items()]
    out.append("<div class='section'><h2>附录A：R1–R4证据分层</h2>" + _table(grade_rows, ["等级", "定义", "核心要求/处理"]) + "</div>")
    glossary = [
        {"术语": "HLA", "通俗解释": "细胞表面的‘展示架’，把肽段展示给T细胞。"},
        {"术语": "APPM", "通俗解释": "抗原从蛋白被切割、运输到HLA展示的一整套加工呈递机制。"},
        {"术语": "LOH", "通俗解释": "某个HLA等位基因在肿瘤中丢失，可能使受它限制的候选无法呈递。"},
        {"术语": "MT/WT", "通俗解释": "突变肽与正常肽成对比较，用于判断候选是否真正具有突变特异性。"},
        {"术语": "CCF", "通俗解释": "估计携带该事件的肿瘤细胞比例；低可信结果不能当作精确比例。"},
        {"术语": "RNA alt reads / RNA VAF", "通俗解释": "直接支持突变转录本的RNA reads及其比例，比仅有gene TPM更接近突变表达证据。"},
        {"术语": "junction reads", "通俗解释": "跨越融合或异常剪接连接点的reads，必须精确对应同一junction。"},
        {"术语": "来源链 C1–C4", "通俗解释": "评价事件到转录本、ORF和肽段是否可追溯；与最终R1–R4推荐等级不同。"},
        {"术语": "Pareto排序", "通俗解释": "在多个证据维度间保留没有被全面压倒的候选，避免单一总分掩盖重要短板。"},
    ]
    out.append("<div class='section'><h2>附录B：术语说明</h2>" + _table(glossary, ["术语", "通俗解释"]) + "</div>")
    trace_rows = [
        {"文件": "运行与来源清单", "用途": "记录输入、参考、工具版本和运行身份"},
        {"文件": "全部工具证据表", "用途": "所有工具证据的统一可追溯表"},
        {"文件": "事件级候选排序表", "用途": "按事件去重的研究候选排序"},
        {"文件": "肽段-HLA候选排序表", "用途": "肽段-HLA级研究候选排序"},
        {"文件": "实验验证计划", "用途": "短肽、长肽、minigene和targeted RNA建议"},
        {"文件": "证据冲突清单", "用途": "来源字段冲突和复核原因"},
    ]
    out.append("<div class='section'><h2>附录C：附件与可追溯文件</h2>" + _table(trace_rows, ["文件", "用途"]) + "<p class='small'>患者版仅列逻辑文件名，不展示服务器绝对路径；实际位置和校验和见技术报告与run manifest。</p></div>")
    out.append("</body></html>")
    p.write_text("\n".join(out), encoding="utf-8")


def _profile_threshold_section(profile: Mapping[str, Any]) -> str:
    sections = []
    for key in ("gates", "safety", "ccf_lite", "score_weights", "l3_weights", "appm_penalty"):
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

    parallel_rankings = bundle.provenance.get("parallel_rankings", {})
    if isinstance(parallel_rankings, Mapping) and parallel_rankings.get("evidence_consensus"):
        out.append("<div class='section'><h2>Experimental parallel evidence-consensus ranking</h2>")
        out.append(
            "<div class='warn'><b>Research-only parallel analysis:</b> "
            "This section does not replace the current primary weighted ranking. "
            "The evidence-consensus ranking has not been experimentally calibrated "
            "and is provided only for algorithm comparison and candidate review.</div>"
        )
        metadata_rows = [
            {"field": "primary patient-report ranking", "value": parallel_rankings.get("legacy_weighted", "ranked_peptides.tsv")},
            {"field": "parallel peptide ranking", "value": parallel_rankings.get("evidence_consensus", "")},
            {"field": "parallel event ranking", "value": parallel_rankings.get("event_consensus", "")},
            {"field": "comparison", "value": parallel_rankings.get("comparison", "")},
            {"field": "rules", "value": parallel_rankings.get("rules_name", "")},
            {"field": "rules version", "value": parallel_rankings.get("rules_version", "")},
            {"field": "rules status", "value": parallel_rankings.get("rules_status", "PROVISIONAL_RESEARCH_ONLY")},
        ]
        out.append(_table(metadata_rows, ["field", "value"]))
        out.append("</div>")

    if bundle.wes_qc:
        out.append("<div class='section'><h2>Independent WES QC</h2>")
        wes_headers = [
            "sample_id", "qc_status", "total_reads", "primary_mapping_rate_pct",
            "properly_paired_rate_pct", "duplicate_rate_pct", "target_definition",
            "mean_target_coverage", "pct_target_bases_20x", "pct_target_bases_30x",
            "on_target_rate_pct", "capture_rate_status", "formal_capture_rate_pct",
        ]
        out.append(_table(bundle.wes_qc, wes_headers))
        if any(row.get("capture_rate_status") != "ASSESSED" for row in bundle.wes_qc):
            out.append(
                "<div class='warn'><b>WES capture QC is partial:</b> coverage and "
                "on-target values use a GENCODE CDS proxy. The assay-specific capture "
                "BED is required before reporting a formal capture rate.</div>"
            )
        out.append("</div>")

    platform_counts = _cross_platform_counts(bundle.events)
    if platform_counts:
        out.append("<div class='section'><h2>WES/WGS Cross-platform DNA Evidence</h2>")
        rows = [{"cross_platform_status": status, "event_count": count} for status, count in sorted(platform_counts.items())]
        out.append(_table(rows, ["cross_platform_status", "event_count"]))
        out.append(
            "<div class='warn'><b>Interpretation:</b> Power-limited absence is not treated as a negative result. "
            "Source-unreproduced InDels require local assembly/IGV review. Sample-specific calls describe the "
            "sequenced specimen and must not be generalized to every tumor time point.</div>"
        )
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
        "source_chain_track", "source_chain_confidence_tier", "source_chain_confidence_label",
        "source_chain_orthogonal_status", "source_chain_orthogonal_sources",
        "source_chain_hard_failure", "source_chain_hard_failure_codes",
        "source_chain_reason_codes", "source_chain_missing_requirements",
        "source_chain_low_power_requirements", "source_chain_negative_requirements",
        "source_chain_conflict_requirements",
        "evidence_grade", "pipeline_r_grade", "review_status", "experiment_priority",
        "cancer_gene_list_status", "cancer_gene_symbols", "cancer_gene_types", "cancer_driver_context",
        "oncokb_annotated", "cosmic_cgc_flag", "cancer_gene_source_count", "cancer_gene_sources", "cancer_gene_context",
        "event_score", "raw_ccf", "ccf_estimate", "ccf_status", "ccf_confidence",
        "ccf_method", "ccf_warning", "clonality_multiplier", "safety_status", "safety_reason",
        "safety_evidence_completeness", "safety_missing_layers", "normal_expression_status",
        "normal_hspc_status", "reference_proteome_status", "normal_ligandome_status",
        "phase_group_id", "haplotype_status", "phase_support_reads",
        "phase_total_informative_reads", "phase_confidence", "component_event_ids",
        "combined_protein_change", "comparison_status", "cross_platform_status",
        "cross_platform_confidence", "cross_platform_multiplier", "cross_platform_priority_cap",
        "cross_platform_review_required", "wes_tumor_depth", "wes_tumor_alt_count",
        "wes_tumor_alt_vaf", "wgs_tumor_depth", "wgs_tumor_alt_count",
        "wgs_tumor_alt_vaf", "normal_depth", "normal_alt_count", "normal_alt_vaf",
    ]))
    out.append("</div>")

    out.append("<div class='section'><h2>Ranked Peptides (full)</h2>")
    pep_headers = [
        "peptide_id", "event_id", "gene", "cancer_gene_types", "cancer_driver_context", "cancer_gene_context",
        "source_chain_track", "source_chain_confidence_tier", "source_chain_confidence_label",
        "source_chain_orthogonal_status", "source_chain_orthogonal_sources",
        "source_chain_requirement_statuses", "source_chain_reason_codes",
        "source_chain_hard_failure", "source_chain_hard_failure_codes",
        "peptide", "wildtype_peptide", "peptide_consequence",
        "hla_allele", "mhc_class", "presentation_evidence_grade", "binding_evidence_score",
        "presentation_evidence_score", "netmhcpan_ba_rank", "netmhcpan_el_rank",
        "netmhcpan_wt_rank_el", "agretopicity_el", "mt_wt_el_rank_difference",
        "mhcflurry_mt_wt_presentation_difference", "prime_mt_wt_score_difference",
        "bigmhc_mt_wt_score_difference", "mutation_positions_in_peptide",
        "mutation_anchor_only", "mutation_tcr_facing", "mutant_specificity_status",
        "mutant_specificity_gate_status", "mutant_specificity_reason", "mutant_specificity_multiplier",
        "phase_group_id", "haplotype_status", "phase_support_reads",
        "phase_total_informative_reads", "phase_confidence", "component_event_ids",
        "combined_protein_change", "redundancy_group",
        "comparison_status", "cross_platform_status", "cross_platform_confidence",
        "cross_platform_multiplier", "cross_platform_review_required",
        "appm_multiplier", "ccf_multiplier", "safety_status", "safety_evidence_completeness",
        "safety_missing_layers", "normal_expression_status", "normal_hspc_status",
        "reference_proteome_status", "normal_ligandome_status", "anchor_assessment_status",
        "escape_status",
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
        if ppt.get("source_chain_confidence_tier"):
            out.append(
                f"<p><b>Source-chain confidence:</b> {_badge(ppt.get('source_chain_confidence_tier'))}; "
                f"track={esc(ppt.get('source_chain_track'))}; "
                f"orthogonal={esc(ppt.get('source_chain_orthogonal_status'))}; "
                f"sources={esc(ppt.get('source_chain_orthogonal_sources'))}; "
                f"missing=<span class='mono'>{esc(ppt.get('source_chain_missing_requirements'))}</span>; "
                f"low-power=<span class='mono'>{esc(ppt.get('source_chain_low_power_requirements'))}</span>; "
                f"negative=<span class='mono'>{esc(ppt.get('source_chain_negative_requirements'))}</span>; "
                f"conflict=<span class='mono'>{esc(ppt.get('source_chain_conflict_requirements'))}</span></p>"
            )
        final_grade = ppt.get("evidence_grade") or ppt.get("pipeline_r_grade") or ppt.get("priority")
        if final_grade:
            out.append(
                f"<p><b>Final recommendation tier:</b> {_badge(final_grade)}; "
                f"review={esc(ppt.get('review_status'))}; "
                f"experiment_priority={esc(ppt.get('experiment_priority'))}. "
                "Source-chain confidence and final recommendation answer different questions: "
                "event authenticity alone does not establish neoantigen priority.</p>"
            )
        out.append(
            f"<p><b>Presentation:</b> grade={esc(ppt.get('presentation_evidence_grade'))}; "
            f"BA={esc(ppt.get('netmhcpan_ba_rank'))}; EL={esc(ppt.get('netmhcpan_el_rank'))}; "
            f"MHCflurry={esc(ppt.get('mhcflurry_presentation_score'))}</p>"
        )
        out.append(
            f"<p><b>Mutant specificity:</b> {_badge(ppt.get('mutant_specificity_gate_status'))}; "
            f"status={esc(ppt.get('mutant_specificity_status'))}; "
            f"MT_EL={esc(ppt.get('netmhcpan_mt_rank_el', ppt.get('netmhcpan_el_rank')))}; "
            f"WT_EL={esc(ppt.get('netmhcpan_wt_rank_el'))}; "
            f"agretopicity={esc(ppt.get('agretopicity_el'))}; "
            f"positions={esc(ppt.get('mutation_positions_in_peptide'))}; "
            f"anchor_only={esc(ppt.get('mutation_anchor_only'))}; "
            f"TCR_facing={esc(ppt.get('mutation_tcr_facing'))}; "
            f"reason=<span class='mono'>{esc(ppt.get('mutant_specificity_reason'))}</span></p>"
        )
        out.append(
            f"<p><b>Haplotype:</b> {_badge(ppt.get('haplotype_status'))}; "
            f"phase_group={esc(ppt.get('phase_group_id'))}; "
            f"support={esc(ppt.get('phase_support_reads'))}/{esc(ppt.get('phase_total_informative_reads'))}; "
            f"confidence={esc(ppt.get('phase_confidence'))}; "
            f"components=<span class='mono'>{esc(ppt.get('component_event_ids'))}</span>; "
            f"protein=<span class='mono'>{esc(ppt.get('combined_protein_change'))}</span></p>"
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
            f"tier={esc(s.get('safety_tier', ''))}; completeness={esc(s.get('safety_evidence_completeness', ppt.get('safety_evidence_completeness')))}; "
            f"missing=<span class='mono'>{esc(s.get('safety_missing_layers', ppt.get('safety_missing_layers')))}</span>; "
            f"reason=<span class='mono'>{esc(s.get('safety_reason', ppt.get('safety_reason')))}</span></p>"
        )
        out.append(
            f"<p><b>CCF:</b> status={esc(c.get('ccf_status', ppt.get('ccf_status')))}; "
            f"raw={esc(c.get('raw_ccf', ppt.get('raw_ccf')))}; "
            f"estimate={esc(c.get('ccf_estimate', ppt.get('ccf_estimate')))}; "
            f"confidence={esc(c.get('ccf_confidence', ppt.get('ccf_confidence')))}; "
            f"method={esc(c.get('ccf_method', ppt.get('ccf_method')))}; "
            f"multiplier={esc(c.get('clonality_multiplier', ppt.get('ccf_multiplier')))}; "
            f"warning=<span class='mono'>{esc(c.get('ccf_warning', ppt.get('ccf_warning')))}</span></p>"
        )
        out.append(
            f"<p><b>WES/WGS evidence:</b> {_badge(ppt.get('cross_platform_status'))}; "
            f"confidence={esc(ppt.get('cross_platform_confidence'))}; "
            f"multiplier={esc(ppt.get('cross_platform_multiplier'))}; "
            f"WES={esc(ppt.get('wes_tumor_alt_count'))}/{esc(ppt.get('wes_tumor_depth'))}; "
            f"WGS={esc(ppt.get('wgs_tumor_alt_count'))}/{esc(ppt.get('wgs_tumor_depth'))}; "
            f"normal={esc(ppt.get('normal_alt_count'))}/{esc(ppt.get('normal_depth'))}</p>"
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
    legacy_name: str = "evidence_report.html",
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
