from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Literal

ExecutionMode = Literal["plan", "execute-safe", "execute-with-approval"]
StepMode = Literal["plan", "execute", "dry_run"]
StepStatus = Literal["PLANNED", "PASS", "FAIL", "SKIPPED", "APPROVAL_REQUIRED"]

SUPPORTED_INTENTS = {
    "input_check",
    "tool_check",
    "demo_smoke",
    "sliding_run",
    "run_scoring",
    "fusion_rna_run",
    "rna_fastq_to_tpm",
    "appm_escape_review",
    "hla_typing_compare",
    "ccf_review",
    "purity_cnv_review",
    "ranking_compare",
    "patient_report_update",
    "experiment_design",
    "general_explanation",
    "project_overview",
    "check_status",
    "inspect_results",
    "data_transfer",
    "debug_error",
    "update_docs",
    "git_release",
    "setup_tool",
    "workflow_run_request",
    "release_qc",
    "unknown",
}

SAFE_SKILLS = {
    "neoag-input-qc",
    "neoag-tool-and-reference-qc",
    "neoag-hla-loh-appm-review",
    "neoag-hla-typing-run-and-compare",
    "neoag-ccf-clonality-review",
    "neoag-purity-cnv-run-and-review",
    "neoag-fusion-rna-run",
    "neoag-rna-fastq-to-tpm",
    "neoag-ranking-compare",
    "neoag-patient-report",
}

WRITE_OR_EXECUTION_SKILLS = {
    "neoag-run-demo-and-smoke",
    "neoag-evidence-scoring",
}

@dataclass
class FileRecord:
    path: str
    name: str
    kind: str = "unknown"
    exists: bool = True
    size_bytes: int | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

@dataclass
class InputState:
    sample_id: str | None = None
    files: list[FileRecord] = field(default_factory=list)
    known_files: dict[str, str] = field(default_factory=dict)
    features: dict[str, bool] = field(default_factory=dict)
    missing_inputs: list[dict[str, str]] = field(default_factory=list)
    recommended_workflow: str | None = None
    input_qc_path: str | None = None

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["files"] = [f.to_dict() for f in self.files]
        return d

@dataclass
class ParsedTask:
    intent: str
    task_type: str = "unknown"
    action: str = "plan"
    workflow: str | None = None
    target: str | None = None
    inputs: list[str] = field(default_factory=list)
    parameters: dict[str, Any] = field(default_factory=dict)
    risk: str = "low"
    missing_fields: list[str] = field(default_factory=list)
    needs_status_tracking: bool = True
    user_visible_summary: str = ""
    source: str = "rule"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

@dataclass
class IntentResult:
    intent: str
    confidence: float
    user_goal: str
    requires_execution: bool = False
    requires_human_approval: bool = False
    missing_information_questions: list[str] = field(default_factory=list)
    source: str = "rule"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

@dataclass
class PlanStep:
    step_id: str
    skill: str
    mode: StepMode = "execute"
    reason: str = ""
    required_inputs: list[str] = field(default_factory=list)
    expected_outputs: list[str] = field(default_factory=list)
    approval_required: bool = False

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

@dataclass
class CoordinatorPlan:
    plan_id: str
    goal: str
    intent: str
    steps: list[PlanStep] = field(default_factory=list)
    approval_required: bool = False
    approval_reason: str = ""
    missing_inputs: list[dict[str, str]] = field(default_factory=list)
    questions_to_user: list[str] = field(default_factory=list)
    model_provider: str = "rule"

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["steps"] = [s.to_dict() for s in self.steps]
        return d

@dataclass
class SkillCallResult:
    skill_name: str
    status: StepStatus
    cmd: list[str] = field(default_factory=list)
    returncode: int | None = None
    stdout: str = ""
    stderr: str = ""
    outputs: dict[str, str] = field(default_factory=dict)
    summary: str = ""
    warnings: list[str] = field(default_factory=list)
    failure_code: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

@dataclass
class CaseState:
    case_id: str
    message: str
    mode: ExecutionMode
    intent: IntentResult
    input_state: InputState
    plan: CoordinatorPlan
    skill_results: list[SkillCallResult] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    outputs: list[dict[str, str]] = field(default_factory=list)
    boundary_note: str = ""
    audit_log: list[dict[str, Any]] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "case_id": self.case_id,
            "message": self.message,
            "mode": self.mode,
            "intent": self.intent.to_dict(),
            "input_state": self.input_state.to_dict(),
            "plan": self.plan.to_dict(),
            "skill_results": [r.to_dict() for r in self.skill_results],
            "warnings": self.warnings,
            "outputs": self.outputs,
            "boundary_note": self.boundary_note,
            "audit_log": self.audit_log,
        }


def validate_intent(intent: str) -> str:
    if intent not in SUPPORTED_INTENTS:
        return "unknown"
    return intent
