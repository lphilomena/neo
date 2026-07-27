from __future__ import annotations

import json
import platform
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from neoag.controlled_execution.doctor import run_doctor
from neoag.controlled_execution.io_utils import markdown_table, sha256_file, write_json, write_tsv

from .contracts import MacroResult, MacroStep
from .state import RunLayout, audit, new_run_id, safe_identifier, update_case_state

TIER_TOOLS = {
    "review": ["python", "neoag", "neoag-skill"],
    "core": ["python", "neoag", "neoag-skill"],
    "prediction": ["python", "neoag", "vep", "netmhcpan", "mhcflurry", "pvacseq"],
    "full": ["python", "neoag", "vep", "netmhcpan", "mhcflurry", "pvacseq", "gatk", "bwa", "star", "optitype", "arriba", "facets"],
}


def _environment_inventory(project_root: Path) -> list[dict[str, Any]]:
    commands = {
        "python": sys.executable,
        "java": "java",
        "docker": "docker",
        "apptainer": "apptainer",
        "nextflow": "nextflow",
        "git": "git",
    }
    rows = []
    for name, command in commands.items():
        path = command if Path(command).is_absolute() and Path(command).exists() else shutil.which(command)
        rows.append({"component": name, "status": "FOUND" if path else "MISSING", "path": str(path or "")})
    disk = shutil.disk_usage(project_root)
    rows.append({"component": "disk_free_bytes", "status": "INFO", "path": str(disk.free)})
    rows.append({"component": "python_version", "status": "INFO", "path": platform.python_version()})
    rows.append({"component": "platform", "status": "INFO", "path": platform.platform()})
    return rows


def _write_local_manifest_templates(layout: RunLayout, project_root: Path) -> dict[str, str]:
    tools = layout.manifests / "tools_manifest.local.yaml"
    refs = layout.manifests / "reference_manifest.local.yaml"
    tools.write_text(
        """manifest_version: 1
tools:
  vep:
    executable: vep
    required: false
  netmhcpan:
    executable: netMHCpan
    license_required: true
    distributable: false
  mhcflurry:
    executable: mhcflurry-predict
    required: false
""",
        encoding="utf-8",
    )
    refs.write_text(
        """manifest_version: 1
genome_build: GRCh38
references:
  reference_fasta:
    path: ''
    required: true
  gencode_gtf:
    path: ''
    required: true
  vep_cache:
    path: ''
    required: true
  normal_proteome:
    path: ''
    required: false
""",
        encoding="utf-8",
    )
    paths = layout.manifests / "paths.env"
    paths.write_text(f"export OPEN_NEO_PROJECT_ROOT={project_root}\n", encoding="utf-8")
    return {"tools_manifest_template": str(tools), "reference_manifest_template": str(refs), "paths_env": str(paths)}


def _tier_status(tier: str, doctor_rows: list[Any]) -> tuple[str, list[str]]:
    by_name = {str(r.name).lower(): str(r.status) for r in doctor_rows}
    import_ok = by_name.get("python_import_neoag") == "OK"
    missing: list[str] = []
    for name in TIER_TOOLS[tier]:
        status = by_name.get(name.lower())
        if name in {"neoag", "neoag-skill"} and import_ok:
            status = "OK"
        if name == "neoag-skill" and (shutil.which("neoag-skill") or shutil.which("neoag")):
            status = "OK"
        if status not in {"OK", "INFO"}:
            missing.append(name)
    if not missing:
        return "READY", []
    if tier in {"review", "core"} and any(x in missing for x in ["python", "neoag"]):
        return "BLOCKED", missing
    return "PARTIAL", missing


def _deployment_command(args: dict[str, Any], project_root: Path, layout: RunLayout, *, execute: bool) -> list[str]:
    deploy_root = Path(str(args.get("deploy_root") or "/opt/neoag"))
    script = project_root / ".agents/skills/neoag-remote-deploy/scripts/16_install_new_machine.sh"
    command = [
        "bash", str(script),
        "--project-root", str(project_root),
        "--tools-root", str(args.get("tools_root") or deploy_root / "env_tool"),
        "--reference-root", str(args.get("reference_root") or deploy_root / "refs"),
        "--licensed-root", str(args.get("licensed_root") or deploy_root / "licensed_tools"),
        "--outdir", str(layout.root / "deployment"),
        "--" + str(args.get("installer_profile") or "minimal"),
    ]
    if args.get("asset_source_host"):
        command += ["--asset-source-host", str(args["asset_source_host"])]
    if bool(args.get("allow_download", False)):
        command.append("--allow-download")
    if bool(args.get("no_sync_assets", False)):
        command.append("--no-sync-assets")
    if execute:
        command.append("--execute")
    return command


def _run_deployment(command: list[str], project_root: Path, log_path: Path) -> tuple[bool, str]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(command, cwd=project_root, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    log_path.write_text(proc.stdout or "", encoding="utf-8")
    return proc.returncode == 0, str(log_path)


def run_install_check(args: dict[str, Any]) -> dict[str, Any]:
    project_root = Path(args.get("project_root") or ".").resolve()
    case_id = safe_identifier(str(args.get("case_id") or "INSTALL"))
    mode = str(args.get("mode") or "verify").lower()
    tier = str(args.get("deployment_tier") or "core").lower()
    if tier not in TIER_TOOLS:
        tier = "core"
    layout = RunLayout.create(args.get("outdir") or f"work/open-neo-install-check/{case_id}")
    result = MacroResult("open-neo-install-check", case_id, new_run_id(case_id, "install"), mode, approval_required=mode in {"repair", "install"}, approved=bool(args.get("approved", False)))
    audit(layout, "install_check.start", "START", mode=mode, tier=tier, project_root=str(project_root))

    result.steps.append(MacroStep("01", "environment-preflight"))
    inventory = _environment_inventory(project_root)
    write_tsv(layout.root / "environment_inventory.tsv", inventory)
    result.steps[-1].status = "PASS"
    result.steps[-1].outputs = {"environment_inventory": str(layout.root / "environment_inventory.tsv")}

    if args.get("release_tarball") and not args.get("sha256"):
        result.blocking_issues.append("CHECKSUM_REQUIRED")
        result.steps.append(MacroStep("02", "release-checksum", "BLOCKED", "A release tarball requires --sha256", failure_code="CHECKSUM_REQUIRED"))
        result.finish("BLOCKED").write(layout.skill_result)
        return result.to_dict()
    if args.get("release_tarball") and args.get("sha256"):
        tarball = Path(str(args["release_tarball"]))
        if not tarball.is_file():
            result.blocking_issues.append("CHECKSUM_FAILED")
            result.steps.append(MacroStep("02", "release-checksum", "FAILED", f"file not found: {tarball}", failure_code="CHECKSUM_FAILED"))
            result.finish("BLOCKED").write(layout.skill_result)
            return result.to_dict()
        observed = sha256_file(tarball)
        if observed != str(args["sha256"]):
            result.blocking_issues.append("CHECKSUM_FAILED")
            result.steps.append(MacroStep("02", "release-checksum", "FAILED", f"observed={observed}", failure_code="CHECKSUM_FAILED"))
            result.finish("BLOCKED").write(layout.skill_result)
            return result.to_dict()

    templates = _write_local_manifest_templates(layout, project_root)
    result.outputs.update(templates)
    result.steps.append(MacroStep("02", "local-manifest-templates", "PASS", outputs=templates))

    if mode in {"repair", "install"} and not result.approved:
        result.blocking_issues.append("APPROVAL_REQUIRED")
        result.steps.append(MacroStep("03", "approval-gate", "APPROVAL_REQUIRED", "Installation or repair requires explicit approval", failure_code="APPROVAL_REQUIRED"))
        result.finish("APPROVAL_REQUIRED").write(layout.skill_result)
        return result.to_dict()

    deploy_command = _deployment_command(args, project_root, layout, execute=mode in {"repair", "install"})
    write_json(layout.root / "deployment_command.json", {"command": deploy_command, "execute": mode in {"repair", "install"}})
    result.outputs["deployment_command"] = str(layout.root / "deployment_command.json")
    if mode in {"repair", "install"}:
        if not Path(deploy_command[1]).is_file():
            result.blocking_issues.append("CORE_INSTALL_FAILED")
            result.steps.append(MacroStep("03", "portable-deployment", "BLOCKED", f"installer not found: {deploy_command[1]}", failure_code="CORE_INSTALL_FAILED"))
            result.finish("BLOCKED").write(layout.skill_result)
            return result.to_dict()
        ok, deploy_log = _run_deployment(deploy_command, project_root, layout.logs / "portable_deployment.log")
        result.steps.append(MacroStep("03", "portable-deployment", "PASS" if ok else "FAILED", outputs={"log": deploy_log}, failure_code="" if ok else "CORE_INSTALL_FAILED"))
        result.outputs["deployment_log"] = deploy_log
        if not ok:
            result.blocking_issues.append("CORE_INSTALL_FAILED")
            result.finish("FAILED").write(layout.skill_result)
            return result.to_dict()

    tools_manifest = args.get("tools_manifest")
    reference_manifest = args.get("reference_manifest")
    result.steps.append(MacroStep("04", "doctor-and-mini-smoke"))
    doctor = run_doctor(
        project_root=project_root,
        outdir=layout.root / "doctor",
        tools_manifest=tools_manifest,
        reference_manifest=reference_manifest,
        sample_manifest=args.get("sample_manifest"),
        profile=str(args.get("profile") or "local"),
        run_demo=bool(args.get("run_demo", tier in {"core", "prediction", "full"})),
        run_pytest=bool(args.get("run_pytest", False)),
        run_nextflow=bool(args.get("run_nextflow", False)),
        mini_smoke=bool(args.get("mini_smoke", tier in {"prediction", "full"})),
        release_audit=bool(args.get("release_audit", True)),
        allow_execute=mode != "plan",
    )
    tier_status, missing = _tier_status(tier, doctor.rows)
    result.steps[-1].status = doctor.status
    result.steps[-1].outputs = doctor.outputs
    result.warnings.extend([f"tier_missing:{x}" for x in missing])
    result.outputs.update({f"doctor_{k}": v for k, v in doctor.outputs.items()})

    final_status = tier_status
    if doctor.status == "UNSAFE":
        final_status = "UNSAFE"
    elif doctor.status == "BLOCKED" and tier_status != "READY":
        final_status = "BLOCKED"
    elif missing:
        final_status = "PARTIAL" if final_status != "BLOCKED" else final_status
    elif doctor.status in {"PARTIAL", "BLOCKED"} and tier_status == "READY":
        result.warnings.append(f"doctor_status={doctor.status}; optional tools/references are incomplete for higher deployment tiers")

    report_rows = [{"deployment_tier": tier, "tier_status": final_status, "doctor_status": doctor.status, "missing_tier_requirements": ",".join(missing), "project_root": str(project_root)}]
    write_tsv(layout.root / "deployment_status.tsv", report_rows)
    md = [
        "# Open-Neo installation and environment check",
        "",
        f"- Deployment tier: **{tier}**",
        f"- Tier status: **{final_status}**",
        f"- Doctor status: **{doctor.status}**",
        f"- Project root: `{project_root}`",
        "",
        "## Missing tier requirements",
        "",
        ", ".join(missing) if missing else "None.",
        "",
        "## Environment inventory",
        "",
        markdown_table(inventory, max_rows=30),
        "",
        "## Boundary",
        "",
        "This check does not bypass licensed-tool terms and does not interpret missing tools as biological negative evidence.",
    ]
    (layout.root / "deployment_report.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    result.outputs.update({"deployment_report": str(layout.root / "deployment_report.md"), "deployment_status": str(layout.root / "deployment_status.tsv")})
    result.provenance = {"python": platform.python_version(), "project_root": str(project_root), "deployment_tier": tier}
    update_case_state(layout, case_id=case_id, current_intent="install_check", deployment_tier=tier, status=final_status, outputs=result.outputs)
    audit(layout, "install_check.finish", final_status, missing=missing)
    result.finish(final_status).write(layout.skill_result)
    return result.to_dict()
