from __future__ import annotations

import json
import platform
from dataclasses import asdict
from pathlib import Path
from typing import Any

from neoag.controlled_execution.doctor import run_doctor
from neoag.controlled_execution.io_utils import sha256_file, write_json, write_tsv
from neoag.production_runner import run_production

from .contracts import MacroResult, MacroStep
from .execution_adapters import (
    discover_result_artifacts,
    ensure_parallel_ranking,
    run_cli,
    write_output_manifest,
    write_run_config,
)
from .routing import (
    RoutingResult,
    build_inventory,
    normalize_manifest,
    route_inputs,
    inspect_manifest,
    write_routing_outputs,
)
from .state import RunLayout, audit, new_run_id, safe_identifier, update_case_state


def _result_from_args(args: dict[str, Any], layout: RunLayout) -> RoutingResult:
    sample_manifest = args.get("sample_manifest")
    if sample_manifest:
        return inspect_manifest(sample_manifest)
    data: dict[str, Any] = {
        "case_id": args.get("case_id") or args.get("sample_id") or "CASE001",
        "sample_id": args.get("sample_id") or args.get("case_id") or "SAMPLE001",
        "genome_build": args.get("genome_build") or "GRCh38",
        "profile": args.get("profile") or "default",
        "inputs": {},
    }
    keys = [
        "somatic_vcf", "fusion_tsv", "splice_junction_tsv", "sv_vcf", "capture_bed",
        "peptide_csv", "raw_events", "raw_peptides", "sv_raw_events", "sv_raw_peptides",
        "hla_file", "hla_alleles", "expression_tsv", "transcript_expression_tsv",
        "rna_evidence_tsv", "purity_tsv", "cnv_tsv", "hla_loh_tsv", "normal_expression",
        "normal_hla_ligands", "reference_proteome", "normal_junctions", "reference_fasta",
        "gencode_gtf", "vep_cache", "production_manifest", "result_dir",
        "comprehensive_evidence", "weighted_baseline",
    ]
    project_root = Path(args.get("project_root") or ".").resolve()
    path_like = {
        "somatic_vcf", "fusion_tsv", "splice_junction_tsv", "sv_vcf", "capture_bed",
        "peptide_csv", "raw_events", "raw_peptides", "sv_raw_events", "sv_raw_peptides",
        "hla_file", "expression_tsv", "transcript_expression_tsv", "rna_evidence_tsv",
        "purity_tsv", "cnv_tsv", "hla_loh_tsv", "normal_expression",
        "normal_hla_ligands", "reference_proteome", "normal_junctions", "reference_fasta",
        "gencode_gtf", "vep_cache", "production_manifest", "result_dir",
        "comprehensive_evidence", "weighted_baseline",
    }
    for key in keys:
        value = args.get(key)
        if value is None or value == "" or value == []:
            continue
        if key in path_like:
            values = value if isinstance(value, list) else [value]
            resolved = []
            for item in values:
                path = Path(str(item))
                resolved.append(str(path if path.is_absolute() else (project_root / path).resolve()))
            value = resolved if isinstance(value, list) else resolved[0]
        data["inputs"][key] = value
    data["execution"] = {
        "profile": args.get("profile") or "default",
        "backend": args.get("backend") or "production-run",
        "reuse_existing": bool(args.get("reuse_existing", True)),
    }
    manifest_path = layout.manifests / "sample_manifest.generated.json"
    write_json(manifest_path, data)
    return inspect_manifest(manifest_path)


def _doctor_preflight(args: dict[str, Any], layout: RunLayout, *, execute: bool) -> dict[str, Any]:
    if not bool(args.get("doctor", True)):
        return {"status": "SKIPPED", "outputs": {}}
    doctor = run_doctor(
        project_root=args.get("project_root") or ".",
        outdir=layout.pipeline / "doctor",
        tools_manifest=args.get("tools_manifest"),
        reference_manifest=args.get("reference_manifest"),
        sample_manifest=args.get("sample_manifest"),
        profile=str(args.get("execution_profile") or "local"),
        run_demo=False,
        run_pytest=False,
        run_nextflow=False,
        mini_smoke=bool(args.get("mini_smoke", False)),
        release_audit=bool(args.get("release_audit", False)),
        allow_execute=execute,
    )
    return {"status": doctor.status, "outputs": doctor.outputs, "summary": doctor.summary}


def _run_sv_only(args: dict[str, Any], routing: RoutingResult, layout: RunLayout) -> dict[str, Any]:
    route = next(r for r in routing.routes if r.route in {"sv_wgs", "sv_wes"})
    sv_vcf = routing.inputs.get("sv_vcf")
    values = sv_vcf if isinstance(sv_vcf, list) else [sv_vcf]
    hla = routing.inputs.get("hla_alleles") or routing.inputs.get("hla_file")
    if isinstance(hla, list):
        hla = ",".join(hla)
    command = [
        "sv-run-full-wes" if route.route == "sv_wes" else "sv-run-full",
        "--sample-id", routing.sample_id,
        "--profile", str(routing.inputs.get("profile") or ("sv_wes_phase1_5" if route.route == "sv_wes" else "sv_wgs_phase1")),
        "--outdir", str(layout.pipeline / "result"),
        "--reference-fasta", str(routing.inputs["reference_fasta"]),
        "--gencode-gtf", str(routing.inputs["gencode_gtf"]),
        "--hla", str(hla),
        "--sv-vcf", *[str(v) for v in values if v],
    ]
    if route.route == "sv_wes":
        command += ["--capture-bed", str(routing.inputs["capture_bed"])]
    optional = {
        "expression_tsv": "--expression",
        "splice_junction_tsv": "--rna-junctions",
        "normal_expression": "--normal-expression",
        "normal_hla_ligands": "--normal-hla-ligands",
        "hla_loh_tsv": "--hla-loh",
        "purity_tsv": "--purity",
        "cnv_tsv": "--cnv",
        "reference_proteome": "--reference-proteome",
        "normal_junctions": "--normal-junctions",
    }
    for key, flag in optional.items():
        if routing.inputs.get(key):
            command += [flag, str(routing.inputs[key])]
    if bool(args.get("stub", False)):
        command += ["--binding-stub", "--immunogenicity-stub"]
    return run_cli(command, cwd=args.get("project_root") or ".", log_path=layout.logs / "sv_run_full.log", timeout=int(args.get("timeout", 7200)))


def _run_standard(args: dict[str, Any], routing: RoutingResult, layout: RunLayout) -> dict[str, Any]:
    config = write_run_config(
        layout.manifests / "run.open_neo.generated.toml",
        {**routing.inputs, "enabled_tools": args.get("enabled_tools") or ["netmhcpan", "mhcflurry"], "immunogenicity_stub": args.get("immunogenicity_stub", args.get("stub", False))},
        [asdict(r) for r in routing.routes],
        outdir=layout.pipeline / "result",
        stub=bool(args.get("stub", False)),
    )
    return run_cli(
        ["run-full", "--config", str(config), "--outdir", str(layout.pipeline / "result")],
        cwd=args.get("project_root") or ".",
        log_path=layout.logs / "run_full.log",
        timeout=int(args.get("timeout", 7200)),
    )


def _result_root(layout: RunLayout, production_result: Any | None = None) -> Path:
    if production_result is not None and getattr(production_result, "final_outdir", ""):
        return Path(production_result.final_outdir)
    return layout.pipeline / "result"


def run_open_neo(args: dict[str, Any]) -> dict[str, Any]:
    mode = str(args.get("mode") or ("ranking-only" if args.get("comprehensive_evidence") else "plan")).lower()
    if mode not in {"plan", "dry-run", "execute", "resume", "ranking-only"}:
        mode = "plan"
    approved = bool(args.get("approved", False))
    provisional_out = Path(args.get("outdir") or "work/open-neo-run")
    provisional_out.mkdir(parents=True, exist_ok=True)
    temp_layout = RunLayout.create(provisional_out)
    routing = _result_from_args(args, temp_layout)
    case_id = safe_identifier(routing.case_id)
    layout = temp_layout
    result = MacroResult("open-neo-run", case_id, new_run_id(case_id, "run"), mode, approval_required=mode in {"execute", "resume"}, approved=approved)
    audit(layout, "open_neo_run.start", "START", mode=mode, routes=[r.route for r in routing.routes])

    route_outputs = write_routing_outputs(routing, layout.input_qc)
    result.outputs.update(route_outputs)
    route_status = "PASS" if routing.status == "PASS" else ("PARTIAL" if routing.status == "PARTIAL" else "BLOCKED")
    result.steps.append(MacroStep("01", "input-detection-and-routing", route_status, detail="; ".join(r.reason for r in routing.routes), outputs=route_outputs))
    result.warnings.extend(routing.warnings)
    result.missing_evidence.extend([x["field"] for x in routing.missing])
    if routing.status == "BLOCKED":
        if any(x.get("field") == "hla_alleles_or_hla_file" for x in routing.missing):
            result.blocking_issues.append("HLA_MISSING")
        else:
            result.blocking_issues.append("AMBIGUOUS_INPUT" if routing.ambiguous else "ROUTE_FAILED")
        result.finish("BLOCKED").write(layout.skill_result)
        return result.to_dict()

    if mode in {"execute", "resume"} and not approved:
        result.steps.append(MacroStep("02", "approval-gate", "APPROVAL_REQUIRED", "Execution requires --approved", failure_code="APPROVAL_REQUIRED"))
        result.blocking_issues.append("APPROVAL_REQUIRED")
        result.finish("APPROVAL_REQUIRED").write(layout.skill_result)
        return result.to_dict()

    execute = mode in {"execute", "resume"}
    if mode == "resume" and not (routing.inputs.get("production_manifest") or routing.inputs.get("result_dir")):
        result.steps.append(MacroStep("02", "resume-input-check", "BLOCKED", "Resume requires a production manifest or an existing result directory", failure_code="RESUME_REQUIRES_PRODUCTION_MANIFEST_OR_RESULT_DIR"))
        result.blocking_issues.append("RESUME_REQUIRES_PRODUCTION_MANIFEST_OR_RESULT_DIR")
        result.finish("BLOCKED").write(layout.skill_result)
        return result.to_dict()
    preflight = _doctor_preflight(args, layout, execute=execute)
    result.steps.append(MacroStep("02", "doctor-preflight", preflight["status"], outputs=preflight.get("outputs", {})))
    result.outputs.update({f"doctor_{k}": v for k, v in preflight.get("outputs", {}).items()})
    if preflight["status"] in {"BLOCKED", "UNSAFE"} and execute and not bool(args.get("allow_partial", False)):
        result.blocking_issues.append("DOCTOR_BLOCKED")
        result.finish(preflight["status"]).write(layout.skill_result)
        return result.to_dict()

    # Plan/dry-run always produces a concrete generated run configuration where possible.
    route_dicts = [asdict(r) for r in routing.routes]
    plan_config = None
    if not routing.inputs.get("production_manifest") and not routing.inputs.get("result_dir") and not routing.inputs.get("comprehensive_evidence"):
        plan_config = write_run_config(layout.manifests / "run.open_neo.generated.toml", routing.inputs, route_dicts, outdir=layout.pipeline / "result", stub=bool(args.get("stub", False)))
        result.outputs["generated_run_config"] = str(plan_config)
    result.steps.append(MacroStep("03", "pipeline-plan", "PASS", inputs={"routes": route_dicts}, outputs={"generated_run_config": str(plan_config or "")}))

    if mode in {"plan", "dry-run"}:
        status = "DRY_RUN" if mode == "dry-run" else "PASS"
        plan_md = ["# Open-Neo run plan", "", f"- Case: {case_id}", f"- Mode: {mode}", f"- Input status: {routing.status}", "", "## Routes"]
        plan_md.extend(f"- {r.route} -> {r.skill}: {r.reason} [{r.status}]" for r in routing.routes)
        plan_md += ["", "## Missing inputs/evidence", ", ".join(result.missing_evidence) if result.missing_evidence else "None.", "", "Weighted baseline and evidence-consensus rankings are generated together after execution."]
        (layout.root / "run_plan.md").write_text("\n".join(plan_md) + "\n", encoding="utf-8")
        result.outputs["run_plan"] = str(layout.root / "run_plan.md")
        update_case_state(layout, case_id=case_id, current_intent="run", mode=mode, routes=route_dicts, status=status)
        result.finish(status).write(layout.skill_result)
        return result.to_dict()

    production_result = None
    command_result: dict[str, Any] | None = None
    result.steps.append(MacroStep("04", "pipeline-execution"))
    try:
        if routing.inputs.get("production_manifest"):
            production_result = run_production(
                routing.inputs["production_manifest"],
                outdir=layout.pipeline / "production",
                project_root=args.get("project_root") or ".",
                execute=True,
                force=bool(args.get("force", False)),
                skip_ranking=False,
            )
            success = production_result.status in {"PASS", "LOW_CONFIDENCE", "PARTIAL"}
            result.steps[-1].status = production_result.status
            result.steps[-1].outputs = {"production_outdir": production_result.outdir, "final_outdir": getattr(production_result, "final_outdir", "")}
        elif routing.inputs.get("result_dir") or routing.inputs.get("comprehensive_evidence"):
            success = True
            result.steps[-1].status = "REUSED"
            result.steps[-1].detail = "Existing result/evidence inputs reused"
        elif len(routing.routes) == 1 and routing.routes[0].route in {"sv_wgs", "sv_wes"}:
            command_result = _run_sv_only(args, routing, layout)
            success = bool(command_result["ok"])
            result.steps[-1].status = "PASS" if success else "FAILED"
            result.steps[-1].outputs = {"log": command_result.get("log", "")}
        elif any(r.route in {"sv_wgs", "sv_wes"} for r in routing.routes):
            success = False
            result.steps[-1].status = "BLOCKED"
            result.steps[-1].failure_code = "SV_REQUIRES_PRODUCTION_MANIFEST_OR_PREBUILT_RAW"
            result.steps[-1].detail = "Mixed SV plus other entries require a production manifest or prebuilt sv_raw_events/sv_raw_peptides"
        else:
            command_result = _run_standard(args, routing, layout)
            success = bool(command_result["ok"])
            result.steps[-1].status = "PASS" if success else "FAILED"
            result.steps[-1].outputs = {"log": command_result.get("log", "")}
    except Exception as exc:
        success = False
        result.steps[-1].status = "FAILED"
        result.steps[-1].detail = str(exc)
        result.steps[-1].failure_code = "PIPELINE_STAGE_FAILED"

    if not success:
        result.blocking_issues.append(result.steps[-1].failure_code or "PIPELINE_STAGE_FAILED")
        result.finish("BLOCKED" if result.steps[-1].status == "BLOCKED" else "FAILED").write(layout.skill_result)
        return result.to_dict()

    if routing.inputs.get("result_dir"):
        result_root = Path(routing.inputs["result_dir"])
    elif routing.inputs.get("comprehensive_evidence"):
        result_root = Path(routing.inputs["weighted_baseline"]).parent
    else:
        result_root = _result_root(layout, production_result)

    result.steps.append(MacroStep("05", "dual-ranking-and-artifact-validation"))
    direct_ranking_only = bool(routing.inputs.get("comprehensive_evidence") and not routing.inputs.get("result_dir"))
    ranking_target = layout.root if direct_ranking_only else None
    ranking = ensure_parallel_ranking(
        project_root=args.get("project_root") or ".",
        result_dir=result_root,
        outdir=ranking_target,
        comprehensive_evidence=routing.inputs.get("comprehensive_evidence"),
        weighted_baseline=routing.inputs.get("weighted_baseline"),
        rules=args.get("rules"),
        provenance=args.get("provenance"),
    )
    if direct_ranking_only:
        result_root = layout.root
    artifacts = discover_result_artifacts(result_root)
    required = ["weighted_baseline", "consensus_peptides", "consensus_events"]
    absent = [x for x in required if not artifacts.get(x)]
    if ranking["status"] == "FAILED" or absent:
        result.steps[-1].status = "FAILED"
        result.steps[-1].detail = "missing artifacts: " + ",".join(absent)
        result.steps[-1].failure_code = "CONSENSUS_RANKING_FAILED"
        result.blocking_issues.append("CONSENSUS_RANKING_FAILED")
        result.finish("FAILED").write(layout.skill_result)
        return result.to_dict()
    result.steps[-1].status = "PASS" if ranking["status"] == "PASS" else "REUSED"
    result.steps[-1].outputs = artifacts
    result.outputs.update(artifacts)
    output_manifest = write_output_manifest(result_root, layout.root / "output_manifest.json")
    result.outputs["output_manifest"] = str(output_manifest)

    macro_run_manifest = {
        "schema_version": "open-neo-run-manifest-v1",
        "run_id": result.run_id,
        "case_id": case_id,
        "sample_id": routing.sample_id,
        "genome_build": routing.genome_build,
        "mode": mode,
        "profile": routing.inputs.get("profile"),
        "routes": route_dicts,
        "result_root": str(result_root.resolve()),
        "artifacts": artifacts,
        "input_hashes": [{"key": row.get("input_key"), "path": row.get("path"), "sha256": sha256_file(row["path"]) if row.get("is_file") and Path(row["path"]).stat().st_size < 50 * 1024 * 1024 else "not_computed"} for row in routing.inventory],
        "status": "PASS_WITH_WARNINGS" if result.warnings or routing.missing else "PASS",
    }
    write_json(layout.run_manifest, macro_run_manifest)
    result.outputs["run_manifest"] = str(layout.run_manifest)
    final_status = "PASS_WITH_WARNINGS" if result.warnings or routing.missing or (production_result and production_result.status == "LOW_CONFIDENCE") else "PASS"
    result.provenance = {"python": platform.python_version(), "project_root": str(Path(args.get("project_root") or ".").resolve()), "result_root": str(result_root.resolve())}
    update_case_state(layout, case_id=case_id, current_intent="run", status=final_status, result_root=str(result_root.resolve()), artifacts=artifacts)
    audit(layout, "open_neo_run.finish", final_status, result_root=str(result_root), artifacts=artifacts)
    result.finish(final_status)
    payload = result.to_dict()
    if mode == "ranking-only":
        payload["production_command"] = "neoag evidence-rank"
        payload["algorithm_owner"] = "src/neoag/evidence_consensus.py"
    write_json(layout.skill_result, payload)
    return payload
