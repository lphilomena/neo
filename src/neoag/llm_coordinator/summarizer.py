from __future__ import annotations

import json
from typing import Any

from .guardrails import computational_boundary_note, validate_final_text
from .model_provider import BaseModelProvider
from .prompts import COORDINATOR_SYSTEM_PROMPT, SUMMARY_PROMPT_TEMPLATE
from .schemas import CoordinatorPlan, IntentResult, SkillCallResult
from .result_collector import collect_result_summaries


def _is_project_overview_request(message: str) -> bool:
    m = message.lower()
    return any(x in m for x in [
        "概述", "介绍项目", "项目功能", "功能介绍", "项目是做什么", "项目做什么",
        "project overview", "what does this project do",
    ])


def project_overview_summary() -> str:
    return """# 项目功能概述

本项目是一个面向肿瘤新抗原分析的综合流程，核心目标是把体细胞变异、HLA 分型、抗原呈递预测、免疫逃逸证据和候选排序整合到一套可复用的分析框架中。

## 主要功能
- 从 SNV/InDel 等变异输入生成候选突变肽，并支持 sliding-window / run-full 等运行模式。
- 调用 NetMHCpan、MHCflurry、PRIME、BigMHC 等工具评估肽段与 HLA 的结合和呈递可能性。
- 整合 HLA typing、HLA LOH、APPM、CCF/克隆性、RNA 表达或 junction 证据等多维信息。
- 对候选 peptide/event 进行综合评分和排序，输出 ranked_peptides、ranked_events 和验证优先级表。
- 生成技术报告和患者沟通报告，便于结果复核、解释和后续实验设计。
- 提供 Nextflow、脚本化运行、Docker/Apptainer 工具容器、安装验收和 LLM agent/web 页面等部署与使用方式。

## 典型输入
- 体细胞 VCF 或事件表
- HLA 分型文件
- 可选 RNA 表达、fusion/splice、HLA LOH、purity/CNV/CCF 等辅助证据

## 典型输出
- ranked_peptides.tsv
- ranked_events.tsv
- validation_plan.tsv
- evidence_report.html
- 各工具证据表、APPM/immune escape/CCF 相关中间结果

## 定位
该项目用于候选新抗原的计算筛选、证据整合和实验验证优先级排序。计算结果不能直接等同于已验证新抗原或临床治疗方案，需要结合实验验证和专业审阅。
"""


def _skill_markdown_report(results: list[SkillCallResult], output_name: str) -> str | None:
    for r in results:
        path = r.outputs.get(output_name)
        if path:
            try:
                return open(path, "r", encoding="utf-8", errors="replace").read()
            except Exception:
                return None
    return None


def _inspection_report(results: list[SkillCallResult]) -> str | None:
    return _skill_markdown_report(results, "result_inspection.md")


def deterministic_summary(message: str, intent: IntentResult, plan: CoordinatorPlan, results: list[SkillCallResult]) -> str:
    if intent.intent in {"general_explanation", "project_overview"} and _is_project_overview_request(message):
        return project_overview_summary()
    if intent.intent == "inspect_results":
        report = _inspection_report(results)
        if report:
            return report + "\n## 边界\n" + computational_boundary_note() + "\n"
    if intent.intent == "purity_cnv_review":
        report = _skill_markdown_report(results, "purity_cnv_review.md")
        if report:
            return report + "\n## 边界\n" + computational_boundary_note() + "\n"
    if intent.intent == "hla_typing_compare":
        report = _skill_markdown_report(results, "hla_typing_compare.md")
        if report:
            return report + "\n## 边界\n" + computational_boundary_note() + "\n"
    if intent.intent == "fusion_rna_run":
        report = _skill_markdown_report(results, "fusion_rna_review.md")
        if report:
            return report + "\n## 边界\n" + computational_boundary_note() + "\n"
    if intent.intent == "rna_fastq_to_tpm":
        report = _skill_markdown_report(results, "rna_expression_summary.md")
        if report:
            return report + "\n## 边界\n" + computational_boundary_note() + "\n"
    lines = ["# LLM-assisted Coordinator Summary", "", f"意图：**{intent.intent}**（来源：{intent.source}，置信度：{intent.confidence:.2f}）", "", "## 执行步骤"]
    for r in results:
        lines.append(f"- {r.skill_name}: {r.status}")
    outputs = []
    for r in results:
        for name, path in r.outputs.items():
            outputs.append((name, path))
    if outputs:
        lines += ["", "## 输出文件"]
        for name, path in outputs:
            lines.append(f"- {name}: `{path}`")
    if plan.approval_required:
        lines += ["", "## 需要人工确认", plan.approval_reason or "该计划包含高风险或写入型操作。"]
    lines += ["", "## 边界", computational_boundary_note()]
    return "\n".join(lines) + "\n"


def summarize_with_llm(message: str, intent: IntentResult, plan: CoordinatorPlan, results: list[SkillCallResult], provider: BaseModelProvider, use_llm: bool = True) -> tuple[str, list[str]]:
    if intent.intent in {"general_explanation", "project_overview"} and _is_project_overview_request(message):
        return validate_final_text(project_overview_summary())
    if intent.intent == "inspect_results":
        report = _inspection_report(results)
        if report:
            return validate_final_text(report + "\n## 边界\n" + computational_boundary_note() + "\n")
    if intent.intent == "purity_cnv_review":
        report = _skill_markdown_report(results, "purity_cnv_review.md")
        if report:
            return validate_final_text(report + "\n## 边界\n" + computational_boundary_note() + "\n")
    if intent.intent == "hla_typing_compare":
        report = _skill_markdown_report(results, "hla_typing_compare.md")
        if report:
            return validate_final_text(report + "\n## 边界\n" + computational_boundary_note() + "\n")
    if intent.intent == "fusion_rna_run":
        report = _skill_markdown_report(results, "fusion_rna_review.md")
        if report:
            return validate_final_text(report + "\n## 边界\n" + computational_boundary_note() + "\n")
    if intent.intent == "rna_fastq_to_tpm":
        report = _skill_markdown_report(results, "rna_expression_summary.md")
        if report:
            return validate_final_text(report + "\n## 边界\n" + computational_boundary_note() + "\n")
    if not use_llm:
        text = deterministic_summary(message, intent, plan, results)
        return validate_final_text(text)
    payload = json.dumps(collect_result_summaries(results), ensure_ascii=False, indent=2)[:8000]
    prompt = SUMMARY_PROMPT_TEMPLATE.format(
        message=message,
        intent=json.dumps(intent.to_dict(), ensure_ascii=False),
        plan=json.dumps(plan.to_dict(), ensure_ascii=False, indent=2),
        skill_summaries=payload,
        boundary_note=computational_boundary_note(),
    )
    try:
        resp = provider.complete([
            {"role": "system", "content": COORDINATOR_SYSTEM_PROMPT},
            {"role": "user", "content": prompt},
        ], temperature=0.1)
        text = resp.text or deterministic_summary(message, intent, plan, results)
    except Exception:
        text = deterministic_summary(message, intent, plan, results)
    return validate_final_text(text)
