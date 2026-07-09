from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

from neoag_v03.agents.coordinator import build_skill_command
from .guardrails import is_skill_safe_to_execute
from .schemas import CoordinatorPlan, ExecutionMode, SkillCallResult


def load_registry(project_root: str | Path = ".") -> dict[str, Any]:
    p = Path(project_root) / ".agents" / "config" / "skills_registry.json"
    if p.exists():
        return json.loads(p.read_text(encoding="utf-8"))
    return {}


def expected_outputs_for_skill(registry: dict[str, Any], skill: str, skill_outdir: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in registry.get(skill, {}).get("outputs", []) or []:
        path = skill_outdir / item
        if path.exists():
            out[item] = str(path)
    # Also include common markdown outputs that might not be in older registries.
    for path in skill_outdir.glob("*.md"):
        out.setdefault(path.name, str(path))
    for path in skill_outdir.glob("*.json"):
        out.setdefault(path.name, str(path))
    return out


def _run_subprocess(cmd: list[str]) -> dict[str, Any]:
    proc = subprocess.run(cmd, text=True, capture_output=True)
    return {
        "returncode": proc.returncode,
        "stdout": proc.stdout[-6000:],
        "stderr": proc.stderr[-6000:],
    }


def execute_plan(plan: CoordinatorPlan, files: dict[str, str], result_dir: str | None, outdir: Path, project_root: str, mode: ExecutionMode, allow_high_risk: bool = False) -> list[SkillCallResult]:
    registry = load_registry(project_root)
    results: list[SkillCallResult] = []
    for step in plan.steps:
        skill_out = outdir / step.skill
        if mode == "plan" or step.mode in {"plan", "dry_run"}:
            cmd = build_skill_command(step.skill, files, result_dir, outdir, project_root, execute=False)
            results.append(SkillCallResult(skill_name=step.skill, status="PLANNED", cmd=cmd or [], summary=step.reason))
            continue
        if step.approval_required and not allow_high_risk:
            cmd = build_skill_command(step.skill, files, result_dir, outdir, project_root, execute=False)
            results.append(SkillCallResult(skill_name=step.skill, status="APPROVAL_REQUIRED", cmd=cmd or [], summary="Human approval required before execution", failure_code="APPROVAL_REQUIRED"))
            continue
        if mode == "execute-safe" and not is_skill_safe_to_execute(step.skill):
            cmd = build_skill_command(step.skill, files, result_dir, outdir, project_root, execute=False)
            results.append(SkillCallResult(skill_name=step.skill, status="APPROVAL_REQUIRED", cmd=cmd or [], summary="Skill is not safe for automatic execution in execute-safe mode", failure_code="APPROVAL_REQUIRED"))
            continue
        cmd = build_skill_command(step.skill, files, result_dir, outdir, project_root, execute=True)
        if cmd is None:
            results.append(SkillCallResult(skill_name=step.skill, status="SKIPPED", summary="Required inputs not found for this skill", failure_code="MISSING_REQUIRED_INPUT"))
            continue
        run = _run_subprocess(cmd)
        outputs = expected_outputs_for_skill(registry, step.skill, skill_out)
        results.append(SkillCallResult(
            skill_name=step.skill,
            status="PASS" if run["returncode"] == 0 else "FAIL",
            cmd=cmd,
            returncode=run["returncode"],
            stdout=run["stdout"],
            stderr=run["stderr"],
            outputs=outputs,
            summary=f"Executed {step.skill}" if run["returncode"] == 0 else f"{step.skill} failed",
            failure_code=None if run["returncode"] == 0 else "SKILL_RUNTIME_ERROR",
        ))
    return results
