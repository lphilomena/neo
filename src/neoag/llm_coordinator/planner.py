from __future__ import annotations

import json
import time
from typing import Any

from neoag.agents.intent_schema import INTENT_TO_SKILLS
from neoag.agents.skill_router import classify_intent
from .guardrails import apply_plan_guardrails
from .json_utils import extract_json_object
from .model_provider import BaseModelProvider, parse_intent_response
from .prompts import COORDINATOR_SYSTEM_PROMPT, INTENT_PROMPT_TEMPLATE, PLANNER_PROMPT_TEMPLATE
from .schema_validation import SchemaValidationError, validate_plan_object, validate_plan_payload
from .schemas import CoordinatorPlan, InputState, IntentResult, PlanStep, validate_intent


def rule_intent(message: str) -> IntentResult:
    intent = classify_intent(message)
    requires_approval = any(k in message.lower() for k in ["hpc", "slurm", "sbatch", "安装", "删除", "覆盖"])
    return IntentResult(
        intent=intent,
        confidence=0.74,
        user_goal=message[:240],
        requires_execution=intent not in {"general_explanation", "unknown"},
        requires_human_approval=requires_approval,
        missing_information_questions=[],
        source="rule",
    )


def classify_with_llm(message: str, provider: BaseModelProvider, file_kinds: list[str], min_confidence: float = 0.80) -> IntentResult:
    if getattr(provider, "name", "") == "rule":
        return rule_intent(message)
    prompt = INTENT_PROMPT_TEMPLATE.format(message=message, file_kinds=", ".join(file_kinds) or "none")
    try:
        resp = provider.complete([
            {"role": "system", "content": COORDINATOR_SYSTEM_PROMPT},
            {"role": "user", "content": prompt},
        ], temperature=0.0, response_format="json_object")
        intent = parse_intent_response(resp.text, source=resp.provider)
        if intent and intent.confidence >= min_confidence:
            return intent
    except SchemaValidationError:
        raise
    except Exception:
        pass
    return rule_intent(message)


def _default_steps_for_intent(intent: str) -> list[PlanStep]:
    skills = INTENT_TO_SKILLS.get(intent, INTENT_TO_SKILLS.get("input_check", []))
    if not skills and intent in {"project_overview", "general_explanation"}:
        skills = ["project-overview"]
    return [
        PlanStep(
            step_id=f"s{i+1}",
            skill=s,
            mode="execute",
            reason=f"Run registered skill for intent {intent}" if s != "project-overview" else "Generate a deterministic project overview without external tools",
        )
        for i, s in enumerate(skills)
    ]


def _coerce_plan(obj: dict[str, Any], fallback_intent: str, provider_name: str) -> CoordinatorPlan | None:
    try:
        steps: list[PlanStep] = []
        for i, raw in enumerate(obj.get("steps") or []):
            if not isinstance(raw, dict):
                continue
            skill = str(raw.get("skill", "")).strip()
            if not skill:
                continue
            steps.append(
                PlanStep(
                    step_id=str(raw.get("step_id") or f"s{i+1}"),
                    skill=skill,
                    mode=str(raw.get("mode") or "execute") if str(raw.get("mode") or "execute") in {"plan", "execute", "dry_run"} else "execute",  # type: ignore[arg-type]
                    reason=str(raw.get("reason", "")),
                    required_inputs=[str(x) for x in (raw.get("required_inputs") or [])],
                    expected_outputs=[str(x) for x in (raw.get("expected_outputs") or [])],
                    approval_required=bool(raw.get("approval_required", False)),
                )
            )
        intent = validate_intent(str(obj.get("intent") or fallback_intent))
        return CoordinatorPlan(
            plan_id=str(obj.get("plan_id") or f"plan_{int(time.time())}"),
            goal=str(obj.get("goal") or intent),
            intent=intent,
            steps=steps or _default_steps_for_intent(intent),
            approval_required=bool(obj.get("approval_required", False)),
            approval_reason=str(obj.get("approval_reason", "")),
            missing_inputs=obj.get("missing_inputs") if isinstance(obj.get("missing_inputs"), list) else [],
            questions_to_user=[str(q) for q in (obj.get("questions_to_user") or [])],
            model_provider=provider_name,
        )
    except Exception:
        return None


def plan_with_llm(message: str, intent: IntentResult, input_state: InputState, skills_registry: dict[str, Any], provider: BaseModelProvider, min_valid_steps: int = 1) -> CoordinatorPlan:
    if getattr(provider, "name", "") == "rule":
        plan = CoordinatorPlan(
            plan_id=f"plan_{int(time.time())}",
            goal=intent.user_goal or intent.intent,
            intent=intent.intent,
            steps=_default_steps_for_intent(intent.intent),
            approval_required=intent.requires_human_approval,
            approval_reason="intent requires approval" if intent.requires_human_approval else "",
            missing_inputs=input_state.missing_inputs,
            questions_to_user=intent.missing_information_questions,
            model_provider="rule",
        )
        plan = apply_plan_guardrails(message, plan)
        validate_plan_object(plan, skills_registry)
        return plan

    reg = json.dumps(skills_registry, ensure_ascii=False, indent=2)[:24000]
    istate = json.dumps(input_state.to_dict(), ensure_ascii=False, indent=2)[:6000]
    prompt = PLANNER_PROMPT_TEMPLATE.format(skills_registry=reg, input_state=istate, intent=json.dumps(intent.to_dict(), ensure_ascii=False))
    try:
        resp = provider.complete([
            {"role": "system", "content": COORDINATOR_SYSTEM_PROMPT},
            {"role": "user", "content": prompt},
        ], temperature=0.0, response_format="json_object")
        obj = extract_json_object(resp.text)
        if obj:
            validate_plan_payload(obj, skills_registry)
            plan = _coerce_plan(obj, intent.intent, resp.provider)
            if plan and len(plan.steps) >= min_valid_steps:
                plan = apply_plan_guardrails(message, plan)
                validate_plan_object(plan, skills_registry)
                return plan
    except SchemaValidationError:
        raise
    except Exception:
        pass
    plan = CoordinatorPlan(
        plan_id=f"plan_{int(time.time())}",
        goal=intent.user_goal or intent.intent,
        intent=intent.intent,
        steps=_default_steps_for_intent(intent.intent),
        approval_required=intent.requires_human_approval,
        approval_reason="intent requires approval" if intent.requires_human_approval else "",
        missing_inputs=input_state.missing_inputs,
        questions_to_user=intent.missing_information_questions,
        model_provider="rule_fallback",
    )
    plan = apply_plan_guardrails(message, plan)
    validate_plan_object(plan, skills_registry)
    return plan
