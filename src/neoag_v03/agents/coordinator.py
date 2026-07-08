from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

from .guardrails import boundary_note, needs_human_approval
from .skill_router import classify_intent, find_named_files, skills_for_intent

SKILL_TO_MODULE = {
    "neoag-input-qc": "neoag_v03.agent_skills.input_qc",
    "neoag-tool-and-reference-qc": "neoag_v03.agent_skills.tool_reference_qc",
    "neoag-run-demo-and-smoke": "neoag_v03.agent_skills.run_demo_smoke",
    "neoag-evidence-scoring": "neoag_v03.agent_skills.evidence_scoring",
    "neoag-hla-loh-appm-review": "neoag_v03.agent_skills.appm_review",
    "neoag-ccf-clonality-review": "neoag_v03.agent_skills.ccf_review",
    "neoag-ranking-compare": "neoag_v03.agent_skills.ranking_compare",
    "neoag-patient-report": "neoag_v03.agent_skills.patient_report",
}


def build_skill_command(skill: str, files: dict[str, str], result_dir: str | None, outdir: Path, project_root: str, execute: bool) -> list[str] | None:
    if skill not in SKILL_TO_MODULE:
        return None
    module = SKILL_TO_MODULE[skill]
    skill_out = outdir / skill
    base = [sys.executable, "-m", module]
    if skill == "neoag-input-qc":
        cmd = base + ["--outdir", str(skill_out)]
        if result_dir:
            cmd += ["--result-dir", result_dir]
        return cmd
    if skill == "neoag-tool-and-reference-qc":
        return base + ["--project-root", project_root, "--outdir", str(skill_out)]
    if skill == "neoag-run-demo-and-smoke":
        cmd = base + ["--project-root", project_root, "--outdir", str(skill_out)]
        if execute:
            cmd += ["--execute"]
        return cmd
    if skill == "neoag-ranking-compare":
        if not files.get("recommendation") or not files.get("netmhcpan42"):
            return None
        return base + ["--recommendation", files["recommendation"], "--netmhcpan42", files["netmhcpan42"], "--outdir", str(skill_out)]
    if skill == "neoag-hla-loh-appm-review":
        cmd = base + ["--outdir", str(skill_out)]
        for flag, key in [("--evidence-report", "evidence_report"), ("--hla-loh", "hla_loh"), ("--appm-gene-status", "appm_gene_status"), ("--appm-submodule-scores", "appm_submodule_scores"), ("--ranked-peptides", "recommendation")]:
            if files.get(key):
                cmd += [flag, files[key]]
        return cmd
    if skill == "neoag-ccf-clonality-review":
        cmd = base + ["--outdir", str(skill_out)]
        if files.get("recommendation"):
            cmd += ["--ranked-peptides", files["recommendation"]]
        elif files.get("ranked_peptides"):
            cmd += ["--ranked-peptides", files["ranked_peptides"]]
        if files.get("purity_table"):
            cmd += ["--purity-table", files["purity_table"]]
        return cmd
    if skill == "neoag-patient-report":
        cmd = base + ["--outdir", str(skill_out)]
        for flag, key in [("--recommendation", "recommendation"), ("--netmhcpan42", "netmhcpan42"), ("--evidence-report", "evidence_report")]:
            if files.get(key):
                cmd += [flag, files[key]]
        rc_report = outdir / "neoag-ranking-compare" / "ranking_compare_report.md"
        appm_report = outdir / "neoag-hla-loh-appm-review" / "appm_escape_review.md"
        ccf_report = outdir / "neoag-ccf-clonality-review" / "ccf_clonality_review.md"
        if rc_report.exists():
            cmd += ["--ranking-compare-report", str(rc_report)]
        if appm_report.exists():
            cmd += ["--appm-review", str(appm_report)]
        if ccf_report.exists():
            cmd += ["--ccf-review", str(ccf_report)]
        return cmd
    return None


def run_skill(cmd: list[str], dry_run: bool) -> dict[str, Any]:
    record: dict[str, Any] = {"cmd": cmd, "dry_run": dry_run, "status": "PLANNED"}
    if dry_run:
        return record
    proc = subprocess.run(cmd, text=True, capture_output=True)
    record.update({"returncode": proc.returncode, "stdout": proc.stdout[-4000:], "stderr": proc.stderr[-4000:], "status": "PASS" if proc.returncode == 0 else "FAIL"})
    return record


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="P1 Coordinator Agent for Project B Skills Pack")
    ap.add_argument("--message", required=True, help="User request in natural language")
    ap.add_argument("--result-dir")
    ap.add_argument("--file", action="append", default=[])
    ap.add_argument("--project-root", default=".")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--execute", action="store_true", help="Execute skills; default is dry-run plan")
    args = ap.parse_args(argv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    intent = classify_intent(args.message)
    skills = skills_for_intent(intent)
    files = find_named_files(args.result_dir, args.file)
    approval = needs_human_approval(args.message, skills)
    dry_run = not args.execute or approval

    calls = []
    for skill in skills:
        cmd = build_skill_command(skill, files, args.result_dir, outdir, args.project_root, execute=args.execute)
        if cmd is None:
            reason = "documented workflow; follow .agents/skills/neoag-sliding-run/SKILL.md" if skill == "neoag-sliding-run" else "required inputs not found"
            calls.append({"skill": skill, "status": "SKIPPED", "reason": reason, "known_files": files})
            continue
        rec = run_skill(cmd, dry_run=dry_run)
        rec["skill"] = skill
        calls.append(rec)

    state = {
        "message": args.message,
        "intent": intent,
        "planned_skills": skills,
        "known_files": files,
        "human_approval_required": approval,
        "dry_run": dry_run,
        "calls": calls,
        "boundary_note": boundary_note(),
    }
    (outdir / "case_state.json").write_text(json.dumps(state, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    lines = ["# Coordinator Agent plan", "", f"Intent: **{intent}**", f"Dry run: **{dry_run}**", "", "## Known files"]
    for k, v in files.items():
        lines.append(f"- {k}: `{v}`")
    lines += ["", "## Skill calls"]
    for c in calls:
        lines.append(f"- {c.get('skill')}: {c.get('status')} — {' '.join(c.get('cmd', [])) if c.get('cmd') else c.get('reason','')}")
    lines += ["", "## Boundary", boundary_note()]
    (outdir / "coordinator_plan.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(state, ensure_ascii=False, indent=2))
    return 0 if all(c.get("status") in {"PASS", "PLANNED", "SKIPPED"} for c in calls) else 1

if __name__ == "__main__":
    raise SystemExit(main())
