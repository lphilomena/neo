from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any


def append_task_event(project_root: str | Path, outdir: str | Path, event: dict[str, Any]) -> None:
    payload = {"ts": time.time(), **event}
    for path in [Path(outdir) / "task_events.jsonl", Path(project_root) / "results" / "agent_task_index.jsonl"]:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(payload, ensure_ascii=False, sort_keys=True) + "\n")


def latest_task_events(project_root: str | Path, limit: int = 20) -> list[dict[str, Any]]:
    path = Path(project_root) / "results" / "agent_task_index.jsonl"
    if not path.exists():
        return []
    rows = []
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines()[-limit:]:
        try:
            rows.append(json.loads(line))
        except Exception:
            pass
    return rows
