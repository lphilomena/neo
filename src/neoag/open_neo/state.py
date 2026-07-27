from __future__ import annotations

import json
import re
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from neoag.controlled_execution.io_utils import append_jsonl, ensure_dir, now_iso, write_json


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
