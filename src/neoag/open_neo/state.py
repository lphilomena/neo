from __future__ import annotations

import json
import hashlib
import re
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from neoag.controlled_execution.io_utils import append_jsonl, ensure_dir, now_iso, write_json


STEP_STATUSES = {
    "PLANNED", "RUNNING", "PASS", "PASS_WITH_WARNINGS", "PARTIAL",
    "REUSED", "SKIPPED", "UNASSESSED", "LOW_CONFIDENCE", "DRY_RUN",
    "QUEUED", "BLOCKED", "FAILED", "UNSAFE", "APPROVAL_REQUIRED",
    "NEEDS_RANKING",
}
REUSABLE_STEP_STATUSES = {"PASS", "PASS_WITH_WARNINGS", "PARTIAL", "REUSED"}


def safe_identifier(value: str, fallback: str = "CASE") -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "")).strip("_.-")
    return cleaned or fallback


def new_run_id(case_id: str, prefix: str = "openneo") -> str:
    stamp = now_iso().replace(":", "").replace("+", "_").replace("-", "")[:15]
    return f"{prefix}-{safe_identifier(case_id)}-{stamp}-{uuid.uuid4().hex[:8]}"


@dataclass(frozen=True)
class RunLayout:
    root: Path
    manifests: Path
    input_qc: Path
    pipeline: Path
    evidence: Path
    ranking: Path
    review: Path
    reports: Path
    logs: Path

    @classmethod
    def create(cls, root: str | Path) -> "RunLayout":
        base = ensure_dir(root)
        names = {n: ensure_dir(base / n) for n in [
            "manifests", "input_qc", "pipeline", "evidence",
            "ranking", "review", "reports", "logs",
        ]}
        return cls(base, **names)

    @property
    def audit_log(self) -> Path:
        return self.root / "audit_log.jsonl"

    @property
    def case_state(self) -> Path:
        return self.root / "case_state.json"

    @property
    def run_manifest(self) -> Path:
        return self.root / "run_manifest.json"

    @property
    def run_state(self) -> Path:
        return self.root / "run_state.json"

    @property
    def skill_result(self) -> Path:
        return self.root / "skill_result.json"


def audit(layout: RunLayout, event: str, status: str, **metadata: Any) -> None:
    append_jsonl(layout.audit_log, {
        "time": now_iso(), "event": event, "status": status, "metadata": metadata,
    })


def update_case_state(layout: RunLayout, **fields: Any) -> dict[str, Any]:
    current: dict[str, Any] = {}
    if layout.case_state.is_file():
        try:
            current = json.loads(layout.case_state.read_text(encoding="utf-8"))
        except Exception:
            current = {}
    current.update(fields)
    current["updated_at"] = now_iso()
    write_json(layout.case_state, current)
    return current


def _file_signature(path: str | Path) -> dict[str, Any]:
    target = Path(path)
    record: dict[str, Any] = {"path": str(target), "exists": target.is_file()}
    if not target.is_file():
        return record
    record["size_bytes"] = target.stat().st_size
    record["mtime_ns"] = target.stat().st_mtime_ns
    if target.stat().st_size > 50 * 1024 * 1024:
        record["sha256"] = "not_computed_large_file"
        return record
    digest = hashlib.sha256()
    with target.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    record["sha256"] = digest.hexdigest()
    return record


def _path_records(value: Any) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    if isinstance(value, dict):
        for nested in value.values():
            records.extend(_path_records(nested))
    elif isinstance(value, (list, tuple)):
        for nested in value:
            records.extend(_path_records(nested))
    elif isinstance(value, (str, Path)) and str(value):
        records.append(_file_signature(value))
    return records


def load_run_state(path: str | Path) -> dict[str, Any]:
    target = Path(path)
    if not target.is_file():
        return {}
    try:
        payload = json.loads(target.read_text(encoding="utf-8"))
        return payload if isinstance(payload, dict) else {}
    except (OSError, json.JSONDecodeError):
        return {}


def persist_result_state(result_path: str | Path, result: dict[str, Any]) -> Path:
    """Persist resumable macro state next to skill_result.json."""
    target = Path(result_path).resolve().parent / "run_state.json"
    steps = []
    for raw in result.get("steps") or []:
        if not isinstance(raw, dict):
            continue
        status = str(raw.get("status") or "PLANNED")
        if status not in STEP_STATUSES:
            status = "UNASSESSED"
        steps.append({
            "step_id": str(raw.get("step_id") or ""),
            "name": str(raw.get("name") or ""),
            "status": status,
            "failure_code": str(raw.get("failure_code") or ""),
            "input_signatures": _path_records(raw.get("inputs") or {}),
            "output_signatures": _path_records(raw.get("outputs") or {}),
            "updated_at": now_iso(),
        })
    state = {
        "schema_version": "open-neo-run-state-v1",
        "run_id": result.get("run_id", ""),
        "case_id": result.get("case_id", ""),
        "skill": result.get("skill", ""),
        "mode": result.get("mode", ""),
        "status": result.get("status", "PLANNED"),
        "started_at": result.get("started_at", ""),
        "finished_at": result.get("finished_at", ""),
        "steps": steps,
        "updated_at": now_iso(),
    }
    return write_json(target, state)


def resume_step_decision(previous_state: dict[str, Any], step_id: str) -> dict[str, Any]:
    """Return a deterministic reuse decision based on status and file signatures."""
    step = next(
        (row for row in previous_state.get("steps") or [] if str(row.get("step_id")) == str(step_id)),
        None,
    )
    if not step:
        return {"decision": "RUN", "reason": "NO_PREVIOUS_STEP"}
    if step.get("status") not in REUSABLE_STEP_STATUSES:
        return {"decision": "RUN", "reason": f"PREVIOUS_STATUS_{step.get('status', 'UNKNOWN')}"}
    for expected in step.get("input_signatures") or []:
        if not expected.get("exists"):
            continue
        path = Path(str(expected.get("path") or ""))
        if not path.is_file():
            return {"decision": "RUN", "reason": f"INPUT_MISSING:{path}"}
        observed = _file_signature(path)
        expected_hash = expected.get("sha256")
        if expected_hash and expected_hash != "not_computed_large_file" and observed.get("sha256") != expected_hash:
            return {"decision": "RUN", "reason": f"INPUT_HASH_CHANGED:{path}"}
        if expected_hash == "not_computed_large_file" and (
            observed.get("size_bytes") != expected.get("size_bytes")
            or observed.get("mtime_ns") != expected.get("mtime_ns")
        ):
            return {"decision": "RUN", "reason": f"INPUT_METADATA_CHANGED:{path}"}
    outputs = step.get("output_signatures") or []
    if not outputs:
        return {"decision": "RUN", "reason": "NO_OUTPUT_SIGNATURES"}
    for expected in outputs:
        path = Path(str(expected.get("path") or ""))
        if not path.is_file():
            return {"decision": "RUN", "reason": f"OUTPUT_MISSING:{path}"}
        observed = _file_signature(path)
        expected_hash = expected.get("sha256")
        if expected_hash and expected_hash != "not_computed_large_file" and observed.get("sha256") != expected_hash:
            return {"decision": "RUN", "reason": f"OUTPUT_HASH_CHANGED:{path}"}
        if expected_hash == "not_computed_large_file" and (
            observed.get("size_bytes") != expected.get("size_bytes")
            or observed.get("mtime_ns") != expected.get("mtime_ns")
        ):
            return {"decision": "RUN", "reason": f"OUTPUT_METADATA_CHANGED:{path}"}
    return {"decision": "REUSE", "reason": "OUTPUT_SIGNATURES_MATCH"}


def build_resume_plan(previous_state: dict[str, Any]) -> list[dict[str, str]]:
    plan: list[dict[str, str]] = []
    for step in previous_state.get("steps") or []:
        step_id = str(step.get("step_id") or "")
        decision = resume_step_decision(previous_state, step_id)
        plan.append({
            "step_id": step_id,
            "name": str(step.get("name") or ""),
            "previous_status": str(step.get("status") or ""),
            "decision": str(decision["decision"]),
            "reason": str(decision["reason"]),
        })
    return plan
