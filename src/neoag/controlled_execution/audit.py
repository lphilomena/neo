from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

from .io_utils import append_jsonl, now_iso, sha256_text


@dataclass
class AuditEvent:
    timestamp: str
    event_type: str
    status: str
    message: str = ""
    case_id: str = ""
    request_id: str = ""
    risk_level: str = "LOW"
    command_preview: str = ""
    approval_status: str = "not_required"
    metadata: dict[str, Any] = field(default_factory=dict)


class AuditLogger:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def log(self, event_type: str, status: str, message: str = "", **kwargs: Any) -> None:
        event = AuditEvent(
            timestamp=now_iso(),
            event_type=event_type,
            status=status,
            message=message,
            **kwargs,
        )
        append_jsonl(self.path, event)


def command_preview(cmd: Sequence[str]) -> str:
    return " ".join(str(x) for x in cmd)


def run_command(
    cmd: Sequence[str],
    *,
    cwd: str | Path | None = None,
    timeout: int = 120,
    env: dict[str, str] | None = None,
    logger: AuditLogger | None = None,
    risk_level: str = "LOW",
    allow_execute: bool = True,
) -> dict[str, Any]:
    preview = command_preview(cmd)
    if logger:
        logger.log("command.start", "DRY_RUN" if not allow_execute else "START", command_preview=preview, risk_level=risk_level)
    if not allow_execute:
        result = {"cmd": preview, "returncode": None, "ok": False, "dry_run": True, "stdout": "", "stderr": ""}
        if logger:
            logger.log("command.finish", "DRY_RUN", command_preview=preview, risk_level=risk_level, metadata=result)
        return result
    try:
        proc = subprocess.run(
            list(cmd),
            cwd=str(cwd) if cwd else None,
            env={**os.environ, **(env or {})},
            text=True,
            capture_output=True,
            timeout=timeout,
        )
        result = {
            "cmd": preview,
            "returncode": proc.returncode,
            "ok": proc.returncode == 0,
            "dry_run": False,
            "stdout": proc.stdout[-10000:],
            "stderr": proc.stderr[-10000:],
        }
    except Exception as exc:
        result = {"cmd": preview, "returncode": 999, "ok": False, "dry_run": False, "stdout": "", "stderr": str(exc)}
    if logger:
        logger.log("command.finish", "PASS" if result.get("ok") else "FAIL", command_preview=preview, risk_level=risk_level, metadata=result)
    return result
