from __future__ import annotations

from pathlib import Path
from typing import Any

from .io import ensure_dir, write_json


def run_input_qc(args: dict[str, Any]) -> dict[str, Any]:
    from neoag_v03.agent_skills.input_qc import main as input_qc_main
    outdir = ensure_dir(args["outdir"])
    argv = ["--outdir", str(outdir)]
    if args.get("manifest"):
        argv += ["--manifest", str(args["manifest"])]
    if args.get("result_dir") or args.get("manifest_or_result_dir"):
        p = args.get("result_dir") or args.get("manifest_or_result_dir")
        argv += ["--result-dir", str(p)]
    rc = input_qc_main(argv)
    res = {"status": "PASS" if rc == 0 else "FAIL", "skill": "neoag-input-qc", "summary": "Input QC completed", "outputs": {"input_status": str(outdir / "input_status.json")}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_doctor(args: dict[str, Any]) -> dict[str, Any]:
    from neoag_v03.controlled_execution.doctor import run_doctor
    outdir = ensure_dir(args["outdir"])
    run_doctor(project_root=args.get("project_root") or ".", outdir=str(outdir), tools_manifest=args.get("tools_manifest"), reference_manifest=args.get("reference_manifest"), sample_manifest=args.get("sample_manifest"), profile=args.get("profile", "local"), run_demo=False, run_pytest=False, run_nextflow=False, release_audit=True, allow_execute=False)
    res = {"status": "PASS", "skill": "neoag-doctor", "summary": "Doctor read-only check completed", "outputs": {"doctor_status": str(outdir / "doctor_status.json")}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_tool_reference_qc(args: dict[str, Any]) -> dict[str, Any]:
    from neoag_v03.agent_skills.tool_reference_qc import main as tool_qc_main
    outdir = ensure_dir(args["outdir"])
    argv = ["--project-root", str(args.get("project_root") or "."), "--outdir", str(outdir)]
    if args.get("tools_manifest"):
        argv += ["--tools-manifest", str(args["tools_manifest"])]
    if args.get("reference_manifest"):
        argv += ["--reference-manifest", str(args["reference_manifest"])]
    rc = tool_qc_main(argv)
    res = {"status": "PASS" if rc == 0 else "FAIL", "skill": "neoag-tool-reference-qc", "summary": "Tool/reference QC completed", "outputs": {"report": str(outdir / "tool_smoke_report.md")}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_demo_smoke(args: dict[str, Any]) -> dict[str, Any]:
    from neoag_v03.agent_skills.run_demo_smoke import main as run_demo_main
    outdir = ensure_dir(args["outdir"])
    argv = ["--project-root", str(args.get("project_root") or "."), "--outdir", str(outdir)]
    if args.get("run_pytest"):
        argv.append("--run-pytest")
    if args.get("run_nextflow"):
        argv.append("--run-nextflow")
    rc = run_demo_main(argv)
    res = {"status": "PASS" if rc == 0 else "FAIL", "skill": "neoag-run-demo-and-smoke", "summary": "Demo smoke completed", "outputs": {"smoke_report": str(outdir / "smoke_test_report.md")}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_pipeline_full(args: dict[str, Any]) -> dict[str, Any]:
    from neoag_v03.controlled_execution.pipeline_runner import run_pipeline_full
    outdir = ensure_dir(args["outdir"])
    execute = bool(args.get("execute") or False)
    approved = bool(args.get("approved") or False)
    if execute and not approved:
        res = {"status": "APPROVAL_REQUIRED", "skill": "neoag-pipeline-full", "summary": "pipeline-full execution is high risk; rerun with approved=true", "outputs": {}}
        write_json(outdir / "skill_result.json", res)
        return res
    run_pipeline_full(sample_manifest=args.get("sample_manifest"), tools_manifest=args.get("tools_manifest"), reference_manifest=args.get("reference_manifest"), outdir=str(outdir), profile=args.get("profile", "local"), allow_execute=execute and approved)
    res = {"status": "PASS", "skill": "neoag-pipeline-full", "summary": "Pipeline-full plan/execution completed", "outputs": {"pipeline_plan": str(outdir / "pipeline_plan.md")}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_release_qc(args: dict[str, Any]) -> dict[str, Any]:
    from neoag_v03.controlled_execution.release_audit import scan_release_boundary, write_release_audit
    outdir = ensure_dir(args["outdir"])
    root = args.get("project_root") or args.get("input") or "."
    result = scan_release_boundary(root)
    write_release_audit(result, str(outdir))
    res = {"status": "PASS", "skill": "neoag-release-qc", "summary": "Release boundary audit completed", "outputs": {"release_audit": str(outdir)}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_gateway_submit(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    approved = bool(args.get("approved") or False)
    task = args.get("task") or args.get("input") or ""
    risk = str(args.get("risk_level") or "HIGH")
    status = "READY_TO_SUBMIT" if approved else "APPROVAL_REQUIRED"
    job = {"status": status, "task": task, "gateway_url": args.get("gateway_url", ""), "risk_level": risk, "approved": approved}
    write_json(outdir / "gateway_job.json", job)
    (outdir / "gateway_submission.md").write_text(f"# Gateway submission\n\nStatus: {status}\nTask: {task}\nRisk: {risk}\n", encoding="utf-8")
    res = {"status": status, "skill": "neoag-gateway-submit", "summary": "Gateway submission prepared", "outputs": {"gateway_job": str(outdir / "gateway_job.json")}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_hpc_runner(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    approved = bool(args.get("approved") or False)
    cmd = str(args.get("command") or "neoag-v03 pipeline-full --sample-manifest sample.yaml --profile slurm")
    script = "#!/usr/bin/env bash\nset -euo pipefail\n# Dry-run generated by neoag-hpc-runner. Review before submit.\n" + cmd + "\n"
    (outdir / "hpc_dry_run.sh").write_text(script, encoding="utf-8")
    manifest = {"status": "APPROVAL_REQUIRED" if not approved else "READY_TO_SUBMIT", "command": cmd, "approved": approved, "note": "This skill never submits directly; use Execution Gateway or scheduler after approval."}
    write_json(outdir / "hpc_job_manifest.json", manifest)
    (outdir / "hpc_submission.md").write_text("# HPC runner\n\nGenerated dry-run submission script. Human approval is required before scheduler submission.\n", encoding="utf-8")
    res = {"status": manifest["status"], "skill": "neoag-hpc-runner", "summary": "Generated HPC dry-run script", "outputs": {"script": str(outdir / "hpc_dry_run.sh")}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_remote_deploy(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    project_root = Path(args.get("project_root") or ".")
    tier = str(args.get("tier") or args.get("deployment_tier") or "tier0")
    execute = bool(args.get("execute") or False)
    approved = bool(args.get("approved") or False)
    skill_dir = project_root / ".agents" / "skills" / "neoag-remote-deploy"
    steps = [
        "Read .agents/skills/neoag-remote-deploy/SKILL.md",
        "Run scripts/00_preflight.sh",
        "Run scripts/02_install_core.sh only if core entry points are missing or user approved install",
        "Run scripts/03_check_runtime.sh",
        "Run scripts/04_configure_manifests.py",
        "Run scripts/05_run_smoke_tests.sh for the selected tier",
        "Run scripts/06_run_doctor.sh",
        "Run scripts/07_write_deploy_report.py",
    ]
    plan = outdir / "remote_deploy_plan.md"
    plan.write_text(
        "# Remote deployment plan\n\n"
        f"Project root: `{project_root}`\n\n"
        f"Deployment tier: `{tier}`\n\n"
        "## Steps\n\n"
        + "\n".join(f"- {step}" for step in steps)
        + "\n\n## Safety\n\nThis skill does not install external tools, download large references, submit HPC jobs, delete files, or overwrite results without explicit approval.\n",
        encoding="utf-8",
    )
    if not execute:
        res = {
            "status": "DRY_RUN",
            "skill": "neoag-remote-deploy",
            "summary": "Deployment plan generated; rerun with execute=true approved=true for local low-risk checks.",
            "outputs": {"plan": str(plan)},
        }
        write_json(outdir / "skill_result.json", res)
        return res
    if not approved:
        res = {
            "status": "APPROVAL_REQUIRED",
            "skill": "neoag-remote-deploy",
            "summary": "Deployment execution requires approved=true even though bundled scripts are conservative.",
            "outputs": {"plan": str(plan)},
        }
        write_json(outdir / "skill_result.json", res)
        return res
    import subprocess
    script = skill_dir / "scripts" / "00_preflight.sh"
    if script.exists():
        subprocess.run([str(script), "--project-root", str(project_root), "--outdir", str(outdir)], check=False)
    script = skill_dir / "scripts" / "04_configure_manifests.py"
    if script.exists():
        subprocess.run([str(script), "--project-root", str(project_root), "--outdir", "configs/local"], check=False)
    script = skill_dir / "scripts" / "07_write_deploy_report.py"
    if script.exists():
        subprocess.run([str(script), "--project-root", str(project_root), "--workdir", str(outdir)], check=False)
    res = {
        "status": "PASS",
        "skill": "neoag-remote-deploy",
        "summary": "Low-risk deployment preparation completed. Review Doctor/smoke steps before production use.",
        "outputs": {"plan": str(plan), "deployment_report": str(outdir / "deployment_report.md")},
    }
    write_json(outdir / "skill_result.json", res)
    return res
