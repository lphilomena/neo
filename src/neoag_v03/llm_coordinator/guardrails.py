from __future__ import annotations

from typing import Iterable

from .schemas import CoordinatorPlan, PlanStep, SAFE_SKILLS, WRITE_OR_EXECUTION_SKILLS
from neoag_v03.agents.guardrails import boundary_note, sanitize_patient_language

HIGH_RISK_TERMS = [
    "hpc",
    "slurm",
    "sbatch",
    "提交",
    "覆盖",
    "删除",
    "rm -",
    "安装",
    "download",
    "下载",
    "正式报告",
    "发布",
]

CLINICAL_RISK_TERMS = [
    "一定有效",
    "保证有效",
    "确定治疗方案",
    "临床耐药",
    "已确认新抗原",
    "一定获益",
]


def message_needs_approval(message: str) -> tuple[bool, str]:
    m = message.lower()
    for t in HIGH_RISK_TERMS:
        if t.lower() in m:
            return True, f"User request contains high-impact operation term: {t}"
    return False, ""


def step_needs_approval(step: PlanStep) -> tuple[bool, str]:
    if step.approval_required:
        return True, f"Step {step.step_id} marks approval_required"
    if step.skill in WRITE_OR_EXECUTION_SKILLS and step.mode == "execute":
        return True, f"Skill {step.skill} may run workflows or write analysis outputs"
    return False, ""


def apply_plan_guardrails(message: str, plan: CoordinatorPlan) -> CoordinatorPlan:
    reasons: list[str] = []
    ok, reason = message_needs_approval(message)
    if ok:
        reasons.append(reason)
    for step in plan.steps:
        ok, reason = step_needs_approval(step)
        if ok:
            step.approval_required = True
            reasons.append(reason)
    if reasons:
        plan.approval_required = True
        plan.approval_reason = "; ".join(sorted(set(reasons)))
    return plan


def is_skill_safe_to_execute(skill: str) -> bool:
    return skill in SAFE_SKILLS


def validate_final_text(text: str) -> tuple[str, list[str]]:
    warnings: list[str] = []
    out = text
    for term in CLINICAL_RISK_TERMS:
        if term in out:
            warnings.append(f"clinical_boundary_rewrite:{term}")
    out = sanitize_patient_language(out)
    if "候选新抗原" in out and "已验证" in out and "不等同" not in out:
        warnings.append("candidate_vs_confirmed_boundary_missing")
        out += "\n\n边界提示：候选新抗原不等同于已验证新抗原，也不等同于确定治疗方案。"
    return out, warnings


def computational_boundary_note() -> str:
    return boundary_note()
