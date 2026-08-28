from __future__ import annotations

import pytest

from neoag.llm_coordinator.schema_validation import (
    SchemaValidationError,
    validate_plan_payload,
    validate_skill_result,
    validate_task_payload,
)
from neoag.llm_coordinator.schemas import SkillCallResult


def test_plan_schema_rejects_unregistered_skill() -> None:
    payload = {"intent": "sliding_run", "steps": [{"skill": "not-registered"}]}
    with pytest.raises(SchemaValidationError) as exc:
        validate_plan_payload(payload, {"neoag-sliding-run": {}})
    assert exc.value.schema_name == "plan"
    assert "not-registered" in ";".join(exc.value.issues)


def test_task_schema_rejects_unknown_action() -> None:
    payload = {"intent": "sliding_run", "action": "teleport", "risk": "medium"}
    with pytest.raises(SchemaValidationError) as exc:
        validate_task_payload(payload)
    assert exc.value.schema_name == "task"
    assert "action" in ";".join(exc.value.issues)


def test_skill_result_schema_accepts_valid_result() -> None:
    result = SkillCallResult(skill_name="neoag-input-qc", status="PASS", outputs={"input_status.json": "out/input_status.json"})
    validate_skill_result(result)
