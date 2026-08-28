from __future__ import annotations

from dataclasses import asdict, is_dataclass
from typing import Any

from .schemas import SUPPORTED_INTENTS, CoordinatorPlan, ParsedTask, SkillCallResult

ALLOWED_ACTIONS = {"plan", "run", "assist", "inspect", "status", "explain"}
ALLOWED_RISKS = {"low", "medium", "high"}
ALLOWED_STEP_MODES = {"plan", "execute", "dry_run"}
ALLOWED_STEP_STATUSES = {"PLANNED", "PASS", "FAIL", "SKIPPED", "APPROVAL_REQUIRED"}


class SchemaValidationError(ValueError):
    """Raised when an LLM/tool JSON object is parseable but violates schema."""

    def __init__(self, schema_name: str, issues: list[str], payload: Any | None = None) -> None:
        self.schema_name = schema_name
        self.issues = issues
        self.payload = payload
        super().__init__(f"{schema_name} schema validation failed: " + "; ".join(issues))

    def to_dict(self) -> dict[str, Any]:
        payload = self.payload
        if is_dataclass(payload):
            payload = asdict(payload)
        return {"schema": self.schema_name, "issues": self.issues, "payload": payload}


def _is_str(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _is_list_of_str(value: Any) -> bool:
    return isinstance(value, list) and all(isinstance(x, str) for x in value)


def _is_list_of_dict(value: Any) -> bool:
    return isinstance(value, list) and all(isinstance(x, dict) for x in value)


def _skill_names(registry: dict[str, Any] | None) -> set[str]:
    return set((registry or {}).keys())


def validate_intent_payload(obj: dict[str, Any]) -> None:
    issues: list[str] = []
    if not isinstance(obj, dict):
        raise SchemaValidationError("intent", ["payload must be an object"], obj)
    if obj.get("intent") not in SUPPORTED_INTENTS:
        issues.append("intent must be one of the supported intents")
    try:
        conf = float(obj.get("confidence"))
        if conf < 0 or conf > 1:
            issues.append("confidence must be between 0 and 1")
    except Exception:
        issues.append("confidence must be numeric")
    if "user_goal" in obj and not isinstance(obj.get("user_goal"), str):
        issues.append("user_goal must be a string")
    for key in ["requires_execution", "requires_human_approval"]:
        if key in obj and not isinstance(obj.get(key), bool):
            issues.append(f"{key} must be boolean")
    if "missing_information_questions" in obj and not _is_list_of_str(obj.get("missing_information_questions")):
        issues.append("missing_information_questions must be a list of strings")
    if issues:
        raise SchemaValidationError("intent", issues, obj)


def validate_task_payload(obj: dict[str, Any]) -> None:
    issues: list[str] = []
    if not isinstance(obj, dict):
        raise SchemaValidationError("task", ["payload must be an object"], obj)
    if obj.get("intent") not in SUPPORTED_INTENTS:
        issues.append("intent must be one of the supported intents")
    if "action" in obj and obj.get("action") not in ALLOWED_ACTIONS:
        issues.append(f"action must be one of {sorted(ALLOWED_ACTIONS)}")
    if "risk" in obj and str(obj.get("risk")).lower() not in ALLOWED_RISKS:
        issues.append(f"risk must be one of {sorted(ALLOWED_RISKS)}")
    if "inputs" in obj and not _is_list_of_str(obj.get("inputs")):
        issues.append("inputs must be a list of strings")
    if "parameters" in obj and not isinstance(obj.get("parameters"), dict):
        issues.append("parameters must be an object")
    if "missing_fields" in obj and not _is_list_of_str(obj.get("missing_fields")):
        issues.append("missing_fields must be a list of strings")
    if "needs_status_tracking" in obj and not isinstance(obj.get("needs_status_tracking"), bool):
        issues.append("needs_status_tracking must be boolean")
    if issues:
        raise SchemaValidationError("task", issues, obj)


def validate_parsed_task(task: ParsedTask) -> None:
    validate_task_payload(task.to_dict())


def validate_plan_payload(obj: dict[str, Any], registry: dict[str, Any] | None = None) -> None:
    issues: list[str] = []
    if not isinstance(obj, dict):
        raise SchemaValidationError("plan", ["payload must be an object"], obj)
    if obj.get("intent") not in SUPPORTED_INTENTS:
        issues.append("intent must be one of the supported intents")
    if "goal" in obj and not isinstance(obj.get("goal"), str):
        issues.append("goal must be a string")
    if "approval_required" in obj and not isinstance(obj.get("approval_required"), bool):
        issues.append("approval_required must be boolean")
    if "missing_inputs" in obj and not _is_list_of_dict(obj.get("missing_inputs")):
        issues.append("missing_inputs must be a list of objects")
    if "questions_to_user" in obj and not _is_list_of_str(obj.get("questions_to_user")):
        issues.append("questions_to_user must be a list of strings")
    steps = obj.get("steps")
    if not isinstance(steps, list):
        issues.append("steps must be a list")
    elif not steps and not bool(obj.get("approval_required", False)):
        issues.append("steps must be a non-empty list unless approval is required")
    else:
        known_skills = _skill_names(registry)
        for i, step in enumerate(steps):
            if not isinstance(step, dict):
                issues.append(f"steps[{i}] must be an object")
                continue
            skill = step.get("skill")
            if not _is_str(skill):
                issues.append(f"steps[{i}].skill is required")
            elif known_skills and skill not in known_skills:
                issues.append(f"steps[{i}].skill is not registered: {skill}")
            mode = step.get("mode", "execute")
            if mode not in ALLOWED_STEP_MODES:
                issues.append(f"steps[{i}].mode must be one of {sorted(ALLOWED_STEP_MODES)}")
            for key in ["required_inputs", "expected_outputs"]:
                if key in step and not _is_list_of_str(step.get(key)):
                    issues.append(f"steps[{i}].{key} must be a list of strings")
            if "approval_required" in step and not isinstance(step.get("approval_required"), bool):
                issues.append(f"steps[{i}].approval_required must be boolean")
    if issues:
        raise SchemaValidationError("plan", issues, obj)


def validate_plan_object(plan: CoordinatorPlan, registry: dict[str, Any] | None = None) -> None:
    validate_plan_payload(plan.to_dict(), registry)


def validate_skill_result_payload(obj: dict[str, Any]) -> None:
    issues: list[str] = []
    if not isinstance(obj, dict):
        raise SchemaValidationError("skill_result", ["payload must be an object"], obj)
    if not _is_str(obj.get("skill_name")):
        issues.append("skill_name is required")
    if obj.get("status") not in ALLOWED_STEP_STATUSES:
        issues.append(f"status must be one of {sorted(ALLOWED_STEP_STATUSES)}")
    if "cmd" in obj and not _is_list_of_str(obj.get("cmd")):
        issues.append("cmd must be a list of strings")
    if "returncode" in obj and obj.get("returncode") is not None and not isinstance(obj.get("returncode"), int):
        issues.append("returncode must be integer or null")
    if "outputs" in obj and not isinstance(obj.get("outputs"), dict):
        issues.append("outputs must be an object")
    if issues:
        raise SchemaValidationError("skill_result", issues, obj)


def validate_skill_result(result: SkillCallResult) -> None:
    validate_skill_result_payload(result.to_dict())


def validate_skill_results(results: list[SkillCallResult]) -> None:
    for result in results:
        validate_skill_result(result)
