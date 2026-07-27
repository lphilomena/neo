from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any

from neoag.controlled_execution.io_utils import now_iso, write_json

MACRO_SCHEMA_VERSION = "open-neo-macro-skill-v1"
TERMINAL_FAILURE = {"BLOCKED", "FAILED", "UNSAFE", "APPROVAL_REQUIRED", "NEEDS_RANKING"}


def _present(value: Any) -> bool:
    if isinstance(value, (list, dict)):
        return bool(value)
    return value not in {None, ""}


def validate_json_schema(data: dict[str, Any], schema: dict[str, Any]) -> list[str]:
    """Validate the JSON-schema subset used by public macro Skills."""
    errors: list[str] = []
    properties = schema.get("properties") or {}
    def matches_type(value: Any, name: str) -> bool:
        if name == "integer":
            return isinstance(value, int) and not isinstance(value, bool)
        if name == "number":
            return isinstance(value, (int, float)) and not isinstance(value, bool)
        return isinstance(value, {
            "string": str, "boolean": bool, "array": list, "object": dict,
        }[name]) if name in {"string", "boolean", "array", "object"} else False
    for key in schema.get("required") or []:
        if key not in data or not _present(data.get(key)):
            errors.append(f"MISSING_REQUIRED:{key}")
    alternatives = schema.get("anyOf") or []
    if alternatives and not any(
        all(key in data and _present(data.get(key)) for key in option.get("required") or [])
        for option in alternatives
    ):
        choices = ["+".join(option.get("required") or []) for option in alternatives]
        errors.append("MISSING_INPUT_CHOICE:" + "|".join(choices))
    for key, value in data.items():
        if not _present(value) or key not in properties:
            continue
        spec = properties[key]
        allowed_types = spec.get("type")
        if allowed_types:
            names = allowed_types if isinstance(allowed_types, list) else [allowed_types]
            if not any(matches_type(value, name) for name in names):
                errors.append(f"INVALID_TYPE:{key}")
                continue
        if spec.get("enum") and value not in spec["enum"]:
            errors.append(f"INVALID_ENUM:{key}:{value}")
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            if "minimum" in spec and value < spec["minimum"]:
                errors.append(f"BELOW_MINIMUM:{key}:{spec['minimum']}")
            if "maximum" in spec and value > spec["maximum"]:
                errors.append(f"ABOVE_MAXIMUM:{key}:{spec['maximum']}")
        if isinstance(value, list):
            if "minItems" in spec and len(value) < spec["minItems"]:
                errors.append(f"TOO_FEW_ITEMS:{key}:{spec['minItems']}")
            item_spec = spec.get("items") or {}
            item_type = item_spec.get("type")
            if item_type and any(not matches_type(item, item_type) for item in value):
                errors.append(f"INVALID_ITEM_TYPE:{key}")
            if item_spec.get("enum") and any(item not in item_spec["enum"] for item in value):
                errors.append(f"INVALID_ITEM_ENUM:{key}")
    return errors


@dataclass
class MacroInput:
    outdir: str = ""
    case_id: str = ""
    extras: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_mapping(cls, data: dict[str, Any]) -> "MacroInput":
        known = {item.name for item in fields(cls) if item.name != "extras"}
        values = {key: data[key] for key in known if key in data and data[key] is not None}
        values["extras"] = {key: value for key, value in data.items() if key not in known}
        return cls(**values)

    def to_mapping(self) -> dict[str, Any]:
        payload = {
            key: value for key, value in asdict(self).items()
            if key != "extras" and _present(value)
        }
        payload.update(self.extras)
        return payload


@dataclass
class InstallCheckInput(MacroInput):
    project_root: str = "."
    release_tarball: str = ""
    deployment_tier: str = "review"
    mode: str = "plan"
    approved: bool = False


@dataclass
class RunInput(MacroInput):
    sample_manifest: str = ""
    mode: str = "plan"
    approved: bool = False
    result_dir: str = ""
    production_manifest: str = ""
    tumor_rna_fastq: list[str] = field(default_factory=list)


@dataclass
class ReviewInput(MacroInput):
    result_dir: str = ""
    top_n: int = 12
    clinical_context: str = ""
    disease_profile: str = ""
    therapy_context: str = "research"
    reports: list[str] = field(default_factory=lambda: ["patient", "technical", "onepage"])


@dataclass
class MacroStep:
    step_id: str
    name: str
    status: str = "PLANNED"
    detail: str = ""
    inputs: dict[str, Any] = field(default_factory=dict)
    outputs: dict[str, Any] = field(default_factory=dict)
    failure_code: str = ""


@dataclass
class MacroResult:
    skill: str
    case_id: str
    run_id: str
    mode: str
    status: str = "PLANNED"
    started_at: str = field(default_factory=now_iso)
    finished_at: str = ""
    steps: list[MacroStep] = field(default_factory=list)
    outputs: dict[str, str] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    blocking_issues: list[str] = field(default_factory=list)
    missing_evidence: list[str] = field(default_factory=list)
    approval_required: bool = False
    approved: bool = False
    provenance: dict[str, Any] = field(default_factory=dict)

    def finish(self, status: str) -> "MacroResult":
        self.status = status
        self.finished_at = now_iso()
        return self

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["schema_version"] = MACRO_SCHEMA_VERSION
        data["approval"] = {
            "required": self.approval_required,
            "approved": self.approved,
        }
        data.pop("approval_required", None)
        data.pop("approved", None)
        return data

    def write(self, path: str | Path) -> Path:
        payload = self.to_dict()
        result_path = write_json(path, payload)
        from .state import persist_result_state

        persist_result_state(result_path, payload)
        return result_path
