from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from neoag.controlled_execution.io_utils import now_iso, write_json

MACRO_SCHEMA_VERSION = "open-neo-macro-skill-v1"
TERMINAL_FAILURE = {"BLOCKED", "FAILED", "UNSAFE", "APPROVAL_REQUIRED"}


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
        return write_json(path, self.to_dict())
