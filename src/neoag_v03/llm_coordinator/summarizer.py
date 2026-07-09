from __future__ import annotations

import json
from typing import Any

from .guardrails import computational_boundary_note, validate_final_text
from .model_provider import BaseModelProvider
from .prompts import COORDINATOR_SYSTEM_PROMPT, SUMMARY_PROMPT_TEMPLATE
from .schemas import CoordinatorPlan, IntentResult, SkillCallResult
from .result_collector import collect_result_summaries


def deterministic_summary(message: str, intent: IntentResult, plan: CoordinatorPlan, results: list[SkillCallResult]) -> str:
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
