from __future__ import annotations

from pathlib import Path
from typing import Any

from .schemas import SkillCallResult


def read_text_head(path: str | Path, limit: int = 2500) -> str:
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
        return text[:limit]
    except Exception:
        return ""


def collect_result_summaries(results: list[SkillCallResult]) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    for r in results:
        item: dict[str, Any] = {
            "skill": r.skill_name,
            "status": r.status,
            "summary": r.summary,
            "failure_code": r.failure_code,
            "outputs": r.outputs,
        }
        for name, path in r.outputs.items():
            if name.endswith(".md") or name.endswith("report.md"):
                item["report_excerpt"] = read_text_head(path)
                break
        if r.stderr and r.status == "FAIL":
            item["stderr_excerpt"] = r.stderr[-1200:]
        summaries.append(item)
    return summaries


def flatten_output_links(results: list[SkillCallResult]) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for r in results:
        for name, path in r.outputs.items():
            out.append({"skill": r.skill_name, "label": name, "path": path})
    return out
