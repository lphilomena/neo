from __future__ import annotations

import json
import platform
from dataclasses import asdict
from pathlib import Path
from typing import Any

from neoag.controlled_execution.doctor import run_doctor
from neoag.controlled_execution.io_utils import sha256_file, write_json, write_tsv
from neoag.controlled_execution.pipeline_runner import run_pipeline_full
from neoag.production_runner import run_production
from neoag.schemas import EVENT_FIELDS, PEPTIDE_FIELDS
from neoag.skill_taxonomy.entry_skills import run_sv_wes, run_sv_wgs
from neoag.utils import read_tsv as read_rows

from .contracts import MacroResult, MacroStep
from .errors import FailureCode
from .execution_adapters import (
    discover_result_artifacts,
    ensure_parallel_ranking,
    run_cli,
    submit_gateway_run,
    write_named_output_manifest,
    write_run_config,
)
from .capability_planner import build_automatic_production_plan, is_automatic_production_candidate
from .rna_preprocessing import prepare_rna_evidence
from .rna_fusion_splice_profile import (
    generate_rna_fusion_splice_manifest,
    is_rna_fastq_profile_candidate,
)
from .routing import (
    RoutingResult,
    build_inventory,
    normalize_manifest,
    route_inputs,
    inspect_manifest,
    write_routing_outputs,
)
from .tool_consensus import build_tool_consensus, enrich_all_tool_results
from .state import (
    RunLayout,
    audit,
    build_resume_plan,
    load_run_state,
    new_run_id,
    persist_result_state,
    safe_identifier,
    update_case_state,
)


def _result_from_args(args: dict[str, Any], layout: RunLayout) -> RoutingResult:
    sample_manifest = args.get("sample_manifest")
    data: dict[str, Any] = {
        "case_id": args.get("case_id") or args.get("sample_id") or "CASE001",
        "sample_id": args.get("sample_id") or args.get("case_id") or "SAMPLE001",
        "genome_build": args.get("genome_build") or "GRCh38",
        "profile": args.get("profile") or "default",
        "inputs": {},
    }
    keys = [
        "tumor_dna_bam", "normal_dna_bam", "tumor_rna_bam", "tumor_dna_fastq", "normal_dna_fastq", "tumor_rna_fastq",
        "tumor_sample_id", "normal_sample_id", "evidence_consensus_rules",
        "somatic_vcf", "fusion_tsv", "splice_junction_tsv", "sv_vcf", "capture_bed",
        "peptide_csv", "raw_events", "raw_peptides", "sv_raw_events", "sv_raw_peptides",
        "hla_file", "hla_alleles", "expression_tsv", "transcript_expression_tsv",
        "rna_evidence_tsv", "purity_tsv", "cnv_tsv", "hla_loh_tsv", "normal_expression",
        "normal_hla_ligands", "reference_proteome", "normal_junctions", "reference_fasta",
        "gencode_gtf", "vep_cache", "production_manifest", "result_dir",
        "salmon_index", "tx2gene", "rsem_reference", "star_index", "ctat_genome_lib",
        "easyfuse_ref", "normal_readthrough", "snaf_workflow", "snaf_db", "snaf_python", "altanalyze_image", "splicemutr_workflow",
        "case_root", "asset_root", "predictor_deps", "netmhcpan_home", "netmhcstabpan_home", "sequenza", "purple",
        "rna_fastq1", "rna_fastq2", "rna_bam", "rna_vaf", "easyfuse_star_index", "star_index_build_dir",
        "star_executable", "samtools_executable", "fusion_caller_root", "prime_evidence", "bigmhc_evidence", "deepimmuno_evidence", "python",
        "rna_threads", "star_sjdb_overhang",
        "comprehensive_evidence", "weighted_baseline",
        "input_dir",
    ]
    project_root = Path(args.get("project_root") or ".").resolve()
    path_like = {
        "tumor_dna_bam", "normal_dna_bam", "tumor_rna_bam", "tumor_dna_fastq", "normal_dna_fastq", "tumor_rna_fastq",
        "evidence_consensus_rules",
        "somatic_vcf", "fusion_tsv", "splice_junction_tsv", "sv_vcf", "capture_bed",
        "peptide_csv", "raw_events", "raw_peptides", "sv_raw_events", "sv_raw_peptides",
        "hla_file", "expression_tsv", "transcript_expression_tsv", "rna_evidence_tsv",
        "purity_tsv", "cnv_tsv", "hla_loh_tsv", "normal_expression",
        "normal_hla_ligands", "reference_proteome", "normal_junctions", "reference_fasta",
        "gencode_gtf", "vep_cache", "production_manifest", "result_dir",
        "salmon_index", "tx2gene", "rsem_reference", "star_index", "ctat_genome_lib",
        "easyfuse_ref", "normal_readthrough", "snaf_workflow", "snaf_db", "snaf_python", "splicemutr_workflow",
        "case_root", "asset_root", "predictor_deps", "netmhcpan_home", "netmhcstabpan_home", "sequenza", "purple",
        "rna_fastq1", "rna_fastq2", "rna_bam", "rna_vaf", "easyfuse_star_index", "star_index_build_dir",
        "star_executable", "samtools_executable", "fusion_caller_root", "prime_evidence", "bigmhc_evidence", "deepimmuno_evidence", "python",
        "comprehensive_evidence", "weighted_baseline",
        "input_dir",
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
    for key in ("tumor_sample_id", "normal_sample_id", "assay_type"):
        if args.get(key):
            data[key] = args[key]
    if args.get("tool_results"):
        data["tool_results"] = args["tool_results"]
    data["execution"] = {
        "profile": args.get("profile") or "default",
        "backend": args.get("backend") or "production-run",
        "reuse_existing": bool(args.get("reuse_existing", True)),
    }
    manifest_path = layout.manifests / "sample_manifest.generated.json"
    write_json(manifest_path, data)
    if sample_manifest:
        overrides = normalize_manifest(data, base_dir=project_root)
        inspected = inspect_manifest(sample_manifest, overrides=overrides, input_dir=args.get("input_dir"), output_dir=layout.root)
    else:
        inspected = inspect_manifest(manifest_path, input_dir=args.get("input_dir"), output_dir=layout.root)
    effective_path = layout.manifests / "sample_manifest.effective.json"
    write_json(effective_path, {
        "schema_version": "open-neo-sample-v1",
        "case_id": inspected.case_id,
        "sample_id": inspected.sample_id,
        "genome_build": inspected.genome_build,
        "tumor_type": inspected.inputs.get("tumor_type", ""),
        "inputs": {key: value for key, value in inspected.inputs.items() if key not in {"case_id", "sample_id", "genome_build", "profile", "backend", "reuse_existing", "execution_mode", "tumor_type", "tumor_sample_id", "normal_sample_id", "assay_type", "tool_results"}},
        "tumor_sample_id": inspected.inputs.get("tumor_sample_id", ""),
        "normal_sample_id": inspected.inputs.get("normal_sample_id", ""),
        "assay_type": inspected.inputs.get("assay_type", ""),
        "tool_results": inspected.inputs.get("tool_results", {}),
        "execution": {"mode": args.get("mode") or inspected.inputs.get("execution_mode") or "plan", "profile": inspected.inputs.get("profile", "default"), "backend": inspected.inputs.get("backend", "production-run"), "reuse_existing": inspected.inputs.get("reuse_existing", True)},
    })
    inspected.manifest_path = str(effective_path)
    return inspected


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



def _path_arg(value: Any) -> str:
    if isinstance(value, list):
        return ",".join(str(v) for v in value if v)
    return str(value or "")


def _run_production_case_wrapper(args: dict[str, Any], routing: RoutingResult, layout: RunLayout) -> dict[str, Any]:
    project_root = Path(args.get("project_root") or ".").resolve()
    wrapper = project_root / "scripts" / "run_production_case.sh"
    command = [
        "bash", str(wrapper),
        "--project-root", str(project_root),
        "--sample-id", routing.sample_id,
        "--case-root", str(routing.inputs["case_root"]),
        "--outdir", str(layout.pipeline / "production_case"),
        "--somatic-vcf", str(routing.inputs["somatic_vcf"]),
        "--rna-threads", str(routing.inputs.get("rna_threads") or args.get("rna_threads") or 16),
    ]
    optional_flags = {
        "profile": "--profile",
        "evidence_consensus_rules": "--evidence-consensus-rules",
        "event_top_n": "--event-top-n",
        "candidate_top_n": "--candidate-top-n",
        "asset_root": "--asset-root",
        "reference_fasta": "--reference-fasta",
        "gencode_gtf": "--gencode-gtf",
        "sequenza": "--sequenza",
        "purple": "--purple",
        "expression_tsv": "--expression",
        "transcript_expression_tsv": "--transcript-expression",
        "rna_bam": "--rna-bam",
        "rna_vaf": "--rna-vaf",
        "star_index": "--star-index",
        "easyfuse_star_index": "--easyfuse-star-index",
        "star_index_build_dir": "--star-index-build-dir",
        "star_sjdb_overhang": "--star-sjdb-overhang",
        "star_executable": "--star-executable",
        "samtools_executable": "--samtools-executable",
        "normal_readthrough": "--normal-readthrough",
        "prime_evidence": "--prime-evidence",
        "bigmhc_evidence": "--bigmhc-evidence",
        "deepimmuno_evidence": "--deepimmuno-evidence",
        "predictor_deps": "--pred-deps",
        "netmhcpan_home": "--netmhcpan-home",
        "netmhcstabpan_home": "--netmhcstabpan-home",
        "python": "--python",
    }
    for key, flag in optional_flags.items():
        value = routing.inputs.get(key, args.get(key))
        if value is None or value == "" or value == []:
            continue
        if key == "profile" and str(value) == "default":
            continue
        command += [flag, _path_arg(value)]
    if routing.inputs.get("rna_fastq1") or routing.inputs.get("rna_fastq2"):
        command += ["--rna-fastq1", _path_arg(routing.inputs.get("rna_fastq1")), "--rna-fastq2", _path_arg(routing.inputs.get("rna_fastq2"))]
    for caller_root in routing.inputs.get("fusion_caller_root") or []:
        if caller_root:
            command += ["--fusion-caller-root", str(caller_root)]
    return run_cli(command, cwd=project_root, log_path=layout.logs / "production_case_wrapper.log", timeout=int(args.get("timeout", 7200)))


def _run_standard(args: dict[str, Any], routing: RoutingResult, layout: RunLayout) -> dict[str, Any]:
    effective_inputs = dict(routing.inputs)
    sv_routes = [route for route in routing.routes if route.route in {"sv_wgs", "sv_wes"}]
    if sv_routes:
        handler = run_sv_wes if sv_routes[0].route == "sv_wes" else run_sv_wgs
        sv_values = effective_inputs.get("sv_vcf") if isinstance(effective_inputs.get("sv_vcf"), list) else [effective_inputs.get("sv_vcf")]
        merged_events: list[dict[str, str]] = []
        merged_peptides: list[dict[str, str]] = []
        for index, sv_vcf in enumerate(sv_values, 1):
            sv_dir = layout.pipeline / "entry_adapters" / sv_routes[0].route / f"input_{index:03d}"
            sv_result = handler({"outdir": str(sv_dir), "sv_vcf": sv_vcf, "capture_bed": effective_inputs.get("capture_bed"), "sample_id": routing.sample_id})
            if sv_result.get("status") != "PASS":
                return {"ok": False, "returncode": 2, "stderr": str(sv_result), "stdout": "", "log": ""}
            merged_events.extend(read_rows(sv_result["outputs"]["raw_events"]))
            merged_peptides.extend(read_rows(sv_result["outputs"]["raw_peptides"]))
        merged_sv = layout.pipeline / "entry_adapters" / sv_routes[0].route / "merged"
        merged_sv.mkdir(parents=True, exist_ok=True)
        event_fields = list(dict.fromkeys([*EVENT_FIELDS, *(key for row in merged_events for key in row)]))
        peptide_fields = list(dict.fromkeys([*PEPTIDE_FIELDS, *(key for row in merged_peptides for key in row)]))
        write_tsv(merged_sv / "raw_events.tsv", merged_events, event_fields)
        write_tsv(merged_sv / "raw_peptides.tsv", merged_peptides, peptide_fields)
        effective_inputs["sv_raw_events"] = str(merged_sv / "raw_events.tsv")
        effective_inputs["sv_raw_peptides"] = str(merged_sv / "raw_peptides.tsv")
    config = write_run_config(
        layout.manifests / "run.open_neo.generated.toml",
        {**effective_inputs, "enabled_tools": args.get("enabled_tools") or ["netmhcpan", "mhcflurry"], "immunogenicity_stub": args.get("immunogenicity_stub", args.get("stub", False))},
        [asdict(r) for r in routing.routes],
        outdir=layout.pipeline / "result",
        stub=bool(args.get("stub", False)),
    )
    return run_cli(
        ["run-full", "--config", str(config), "--outdir", str(layout.pipeline / "result"), "--reports", "technical"],
        cwd=args.get("project_root") or ".",
        log_path=layout.logs / "run_full.log",
        timeout=int(args.get("timeout", 7200)),
    )


def _result_root(layout: RunLayout, production_result: Any | None = None) -> Path:
    if production_result is not None and getattr(production_result, "final_outdir", ""):
        return Path(production_result.final_outdir)
    return layout.pipeline / "result"


def _register_production_tool_outputs(inputs: dict[str, Any], production_result: Any | None) -> None:
    if production_result is None:
        return
    tool_results = inputs.setdefault("tool_results", {})
    stage_domains = {
        "sample_identity_bam_matcher": ("sample_identity", "bam-matcher"),
        "hla_optitype": ("hla_typing", "optitype"),
        "hla_hla_la": ("hla_typing", "hla-la"),
        "hla_spechla": ("hla_typing", "spechla"),
        "hla_loh_lohhla": ("hla_loh", "lohhla"),
        "hla_loh_multi_tool": ("hla_loh", "multi-tool"),
        "purity_facets": ("purity_cnv", "facets"),
        "purity_sequenza": ("purity_cnv", "sequenza"),
        "purity_purple": ("purity_cnv", "purple"),
        "purity_ascat": ("purity_cnv", "ascat"),
        "rna_expression": ("rna_expression", "primary"),
        "rsem_expression_crosscheck": ("rna_expression", "rsem"),
        "easyfuse_discovery": ("fusion", "easyfuse"),
        "star_fusion_discovery": ("fusion", "star-fusion"),
        "fusioncatcher_discovery": ("fusion", "fusioncatcher"),
        "arriba_discovery": ("fusion", "arriba"),
        "junction_extraction": ("splice_rna", "regtools"),
        "snaf_discovery": ("splice_neoantigen", "snaf"),
        "splicemutr_discovery": ("splice_neoantigen", "splicemutr"),
    }
    for stage in getattr(production_result, "stages", []):
        if stage.status not in {"PASS", "REUSED"} or stage.name not in stage_domains:
            continue
        domain, tool = stage_domains[stage.name]
        values = [str(value) for value in stage.outputs.values() if isinstance(value, str) and Path(value).exists()]
        if values:
            tool_results.setdefault(domain, {})[tool] = values[0]
        if stage.name == "hla_loh_multi_tool":
            for tool, key in (("lohhla", "lohhla_hla_loh"), ("spechla", "spechla_hla_loh")):
                path = str(stage.outputs.get(key) or "")
                if path and Path(path).is_file():
                    tool_results.setdefault("hla_loh", {})[tool] = path
            consensus = str(stage.outputs.get("hla_loh_consensus") or "")
            if consensus and Path(consensus).is_file():
                tool_results.setdefault("hla_loh", {})["provided_consensus"] = consensus


def run_open_neo(args: dict[str, Any]) -> dict[str, Any]:
    provisional_out = Path(args.get("outdir") or "work/open-neo-run")
    provisional_out.mkdir(parents=True, exist_ok=True)
    temp_layout = RunLayout.create(provisional_out)
    routing = _result_from_args(args, temp_layout)
    mode = str(args.get("mode") or routing.inputs.get("execution_mode") or ("ranking-only" if args.get("comprehensive_evidence") else "plan")).lower()
    if mode not in {"plan", "dry-run", "execute", "resume", "ranking-only"}:
        mode = "plan"
    approved = bool(args.get("approved", False))
    case_id = safe_identifier(routing.case_id)
    layout = temp_layout
    previous_state = load_run_state(layout.run_state)
    result = MacroResult("open-neo-run", case_id, new_run_id(case_id, "run"), mode, approval_required=mode in {"execute", "resume"}, approved=approved)
    if mode == "resume":
        resume_plan = build_resume_plan(previous_state)
        write_tsv(layout.manifests / "resume_plan.tsv", resume_plan)
        result.outputs["resume_plan"] = str(layout.manifests / "resume_plan.tsv")
        result.provenance["previous_run_id"] = str(previous_state.get("run_id") or "")
        result.provenance["resume_reuse_steps"] = sum(1 for row in resume_plan if row["decision"] == "REUSE")
    audit(layout, "open_neo_run.start", "START", mode=mode, routes=[r.route for r in routing.routes])

    rna_modes = [
        label for label, present in (
            ("RNA_FASTQ", bool(routing.inputs.get("tumor_rna_fastq") or routing.inputs.get("rna_fastq1"))),
            ("RNA_BAM", bool(routing.inputs.get("tumor_rna_bam") or routing.inputs.get("rna_bam"))),
            ("RNA_VAF", bool(routing.inputs.get("rna_evidence_tsv") or routing.inputs.get("rna_vaf"))),
        ) if present
    ]
    if len(rna_modes) > 1:
        result.steps.append(MacroStep(
            "00", "rna-evidence-input-mode", "BLOCKED",
            "Use only one RNA allele-evidence input mode: FASTQ pair, RNA BAM, or existing RNA VAF; observed=" + ",".join(rna_modes),
            failure_code=FailureCode.AMBIGUOUS_INPUT.value,
        ))
        result.blocking_issues.append(FailureCode.AMBIGUOUS_INPUT.value)
        result.finish("BLOCKED").write(layout.skill_result)
        return result.to_dict()

    route_outputs = write_routing_outputs(routing, layout.input_qc)
    result.outputs.update(route_outputs)
    automatic_candidate = is_automatic_production_candidate(routing.inputs)
    route_status = "PASS" if routing.status == "PASS" else ("PARTIAL" if routing.status == "PARTIAL" or automatic_candidate else "BLOCKED")
    result.steps.append(MacroStep("01", "input-detection-and-routing", route_status, detail="; ".join(r.reason for r in routing.routes), outputs=route_outputs))
    result.warnings.extend(routing.warnings)
    result.missing_evidence.extend([x["field"] for x in routing.missing])
    if routing.status == "BLOCKED" and not automatic_candidate:
        if any(x.get("field") == "hla_alleles_or_hla_file" for x in routing.missing):
            result.blocking_issues.append(FailureCode.HLA_MISSING.value)
        else:
            result.blocking_issues.append(FailureCode.AMBIGUOUS_INPUT.value if routing.ambiguous else FailureCode.ROUTE_FAILED.value)
        result.finish("BLOCKED").write(layout.skill_result)
        return result.to_dict()
    if routing.status == "BLOCKED" and automatic_candidate:
        result.warnings.append(
            "Input QC found prerequisites that may be generated by the capability-aware production DAG; final blocking status is deferred to the capability planner"
        )

    if mode in {"execute", "resume"} and not approved:
        result.steps.append(MacroStep("02", "approval-gate", "APPROVAL_REQUIRED", "Execution requires --approved", failure_code=FailureCode.APPROVAL_REQUIRED.value))
        result.blocking_issues.append(FailureCode.APPROVAL_REQUIRED.value)
        result.finish("APPROVAL_REQUIRED").write(layout.skill_result)
        return result.to_dict()

    production_case_candidate = any(r.route == "production_case_wrapper" for r in routing.routes)
    if mode in {"execute", "resume"} and not bool(args.get("_gateway_dispatched", False)) and not production_case_candidate:
        gateway_url = str(args.get("gateway_url") or "")
        if not gateway_url:
            result.steps.append(MacroStep("02", "gateway-execution-boundary", "BLOCKED", "Direct heavy execution is disabled; submit through NeoAg Gateway with --gateway-url", failure_code=FailureCode.GATEWAY_REQUIRED.value))
            result.blocking_issues.append(FailureCode.GATEWAY_REQUIRED.value)
            result.finish("BLOCKED").write(layout.skill_result)
            return result.to_dict()
        gateway_payload = {key: value for key, value in args.items() if not key.startswith("_gateway") and key not in {"gateway_url", "gateway_wait"} and value is not None and value != ""}
        gateway_payload["mode"] = mode
        gateway_payload["approved"] = True
        submitted = submit_gateway_run(gateway_url, gateway_payload, wait=bool(args.get("gateway_wait", False)), timeout=int(args.get("timeout", 7200)))
        result.steps.append(MacroStep("02", "gateway-submit", str(submitted.get("status") or "FAILED"), outputs={"job_id": str(submitted.get("job_id") or "")}, detail=str(submitted.get("message") or "")))
        result.outputs["gateway_job_id"] = str(submitted.get("job_id") or "")
        result.outputs["gateway_response"] = str(layout.root / "gateway_response.json")
        write_json(layout.root / "gateway_response.json", submitted)
        final = "QUEUED" if submitted.get("status") in {"QUEUED", "RUNNING"} else str(submitted.get("status") or "FAILED")
        result.finish(final).write(layout.skill_result)
        return result.to_dict()

    execute = mode in {"execute", "resume"}
    if mode == "resume" and not (
        routing.inputs.get("production_manifest")
        or routing.inputs.get("result_dir")
        or automatic_candidate
        or is_rna_fastq_profile_candidate(routing.inputs)
    ):
        result.steps.append(MacroStep("02", "resume-input-check", "BLOCKED", "Resume requires a production manifest or an existing result directory", failure_code=FailureCode.RESUME_INPUT_REQUIRED.value))
        result.blocking_issues.append(FailureCode.RESUME_INPUT_REQUIRED.value)
        result.finish("BLOCKED").write(layout.skill_result)
        return result.to_dict()

    automatic_plan = None
    if is_automatic_production_candidate(routing.inputs):
        automatic_manifest = layout.manifests / "capability_aware.production.toml"
        automatic_plan = build_automatic_production_plan(
            routing.inputs,
            automatic_manifest,
            project_root=args.get("project_root") or ".",
            outdir=layout.pipeline / "automatic_production",
            tools_manifest=args.get("tools_manifest"),
            reference_manifest=args.get("reference_manifest"),
            policy=str(args.get("automatic_tool_policy") or "all-available"),
        )
        routing.inputs["production_manifest"] = automatic_plan.manifest
        result.outputs.update(automatic_plan.outputs)
        result.outputs["production_results_manifest"] = automatic_plan.manifest
        generated_domains = {
            decision.domain for decision in automatic_plan.decisions
            if decision.status == "SELECTED"
        }
        if "hla_typing" in generated_domains:
            result.missing_evidence = [
                item for item in result.missing_evidence
                if item not in {"hla_file", "hla_alleles_or_hla_file"}
            ]
        if "somatic_variant" in generated_domains:
            result.missing_evidence = [item for item in result.missing_evidence if item != "somatic_vcf"]
        if routing.inputs.get("tumor_rna_fastq"):
            result.outputs["generated_production_manifest"] = automatic_plan.manifest
            result.outputs["rna_fusion_splice_requirements"] = automatic_plan.outputs["capability_decisions"]
        result.steps.append(MacroStep(
            "02a", "capability-aware-production-plan", automatic_plan.status,
            detail="missing_required=" + ",".join(automatic_plan.missing_required),
            outputs=automatic_plan.outputs,
        ))
        result.missing_evidence.extend(automatic_plan.missing_required)
        if execute and automatic_plan.status == "BLOCKED":
            result.blocking_issues.append(FailureCode.PRODUCTION_MANIFEST_REQUIRED.value)
            result.finish("BLOCKED").write(layout.skill_result)
            return result.to_dict()

    auto_rna_profile: dict[str, Any] | None = None
    if automatic_plan is None and is_rna_fastq_profile_candidate(routing.inputs):
        profile_manifest = layout.manifests / "rna_fusion_splice.production.toml"
        requirements_tsv = layout.manifests / "rna_fusion_splice.requirements.tsv"
        auto_rna_profile = generate_rna_fusion_splice_manifest(
            routing.inputs,
            profile_manifest,
            project_root=args.get("project_root") or ".",
            outdir=layout.pipeline / "rna_fusion_splice_run",
        )
        routing.inputs["production_manifest"] = str(profile_manifest)
        write_tsv(requirements_tsv, auto_rna_profile["requirements"])
        result.outputs.update({
            "generated_production_manifest": str(profile_manifest),
            "rna_fusion_splice_requirements": str(requirements_tsv),
        })
        profile_status = "PASS" if auto_rna_profile["ready_for_execute"] else "PARTIAL"
        result.steps.append(MacroStep(
            "02a", "rna-fastq-fusion-splice-profile", profile_status,
            detail="missing_required=" + ",".join(auto_rna_profile["missing_required"]),
            outputs={
                "production_manifest": str(profile_manifest),
                "requirements": str(requirements_tsv),
            },
        ))
        if execute and not auto_rna_profile["ready_for_execute"]:
            result.blocking_issues.extend(
                f"{FailureCode.RNA_FUSION_SPLICE_MISSING.value}:{field}" for field in auto_rna_profile["missing_required"]
            )
            result.finish("BLOCKED").write(layout.skill_result)
            return result.to_dict()

    preflight = _doctor_preflight(args, layout, execute=execute)
    result.steps.append(MacroStep("02", "doctor-preflight", preflight["status"], outputs=preflight.get("outputs", {})))
    result.outputs.update({f"doctor_{k}": v for k, v in preflight.get("outputs", {}).items()})
    if preflight["status"] in {"BLOCKED", "UNSAFE"} and execute and not bool(args.get("allow_partial", False)):
        result.blocking_issues.append(FailureCode.DOCTOR_BLOCKED.value)
        result.finish(preflight["status"]).write(layout.skill_result)
        return result.to_dict()

    rna_preprocessing = (
        {"status": "PLANNED_IN_PRODUCTION_DAG", "outputs": {}}
        if automatic_plan is not None
        else prepare_rna_evidence(
            routing.inputs,
            project_root=args.get("project_root") or ".",
            outdir=layout.pipeline / "rna_preprocessing",
            execute=execute and auto_rna_profile is None,
            method=str(args.get("rna_quant_method") or "auto"),
            timeout=int(args.get("timeout", 7200)),
        )
    )
    rna_outputs = rna_preprocessing.get("outputs") or {}
    if automatic_plan is not None:
        rna_status_dir = layout.pipeline / "rna_preprocessing"
        write_tsv(rna_status_dir / "rna_preprocessing_status.tsv", [{
            "stage": "automatic_production_dag", "status": "PLANNED",
            "command_preview": "see capability_aware.production.toml", "output": "",
            "message": "RNA stages are owned by the unified production DAG",
        }])
        write_json(rna_status_dir / "rna_preprocessing_summary.json", rna_preprocessing)
    if rna_outputs.get("gene_tpm"):
        routing.inputs["expression_tsv"] = rna_outputs["gene_tpm"]
    if rna_outputs.get("transcript_tpm"):
        routing.inputs["transcript_expression_tsv"] = rna_outputs["transcript_tpm"]
    if rna_outputs.get("rna_alt_vaf"):
        routing.inputs["rna_evidence_tsv"] = rna_outputs["rna_alt_vaf"]
    rna_files = {
        "rna_preprocessing_status": str(layout.pipeline / "rna_preprocessing" / "rna_preprocessing_status.tsv"),
        "rna_preprocessing_summary": str(layout.pipeline / "rna_preprocessing" / "rna_preprocessing_summary.json"),
        **{f"rna_{key}": value for key, value in rna_outputs.items()},
    }
    result.outputs.update(rna_files)
    result.steps.append(MacroStep("02b", "rna-expression-and-allele-evidence", str(rna_preprocessing.get("status") or "UNASSESSED"), outputs=rna_files))

    # Plan/dry-run always produces a concrete generated run configuration where possible.
    route_dicts = [asdict(r) for r in routing.routes]
    plan_config = None
    if not routing.inputs.get("production_manifest") and not routing.inputs.get("result_dir") and not routing.inputs.get("comprehensive_evidence"):
        plan_config = write_run_config(layout.manifests / "run.open_neo.generated.toml", routing.inputs, route_dicts, outdir=layout.pipeline / "result", stub=bool(args.get("stub", False)))
        result.outputs["generated_run_config"] = str(plan_config)
    result.steps.append(MacroStep("03", "pipeline-plan", "PASS", inputs={"routes": route_dicts}, outputs={"generated_run_config": str(plan_config or "")}))

    pipeline_plan = run_pipeline_full(
        sample_manifest=routing.manifest_path,
        tools_manifest=args.get("tools_manifest"),
        reference_manifest=args.get("reference_manifest"),
        project_root=args.get("project_root") or ".",
        outdir=layout.pipeline / "pipeline_full_plan",
        profile=str(args.get("execution_profile") or "local"),
        dry_run=True,
        allow_partial=True,
        run_doctor_first=False,
    )
    plan_outputs = {
        "pipeline_plan": str(Path(pipeline_plan.output_dir) / "pipeline_plan.md"),
        "pipeline_status": str(Path(pipeline_plan.output_dir) / "pipeline_status.tsv"),
        "pipeline_run_manifest": str(Path(pipeline_plan.output_dir) / "run_manifest.json"),
    }
    result.outputs.update(plan_outputs)
    result.steps.append(MacroStep("03b", "pipeline-full-dry-run", pipeline_plan.status, outputs=plan_outputs))
    initial_consensus = build_tool_consensus(routing.inputs, layout.pipeline / "tool_consensus")
    result.outputs.update({f"consensus_{key.removesuffix('.tsv')}": value for key, value in initial_consensus.items()})

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
        if any(r.route == "production_case_wrapper" for r in routing.routes):
            command_result = _run_production_case_wrapper(args, routing, layout)
            success = bool(command_result["ok"])
            result.steps[-1].status = "PASS" if success else "FAILED"
            result.steps[-1].outputs = {"log": command_result.get("log", ""), "production_outdir": str(layout.pipeline / "production_case")}
            if success:
                routing.inputs["result_dir"] = str(layout.pipeline / "production_case")
        elif routing.inputs.get("production_manifest"):
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
            _register_production_tool_outputs(routing.inputs, production_result)
        elif routing.inputs.get("result_dir") or routing.inputs.get("comprehensive_evidence"):
            success = True
            result.steps[-1].status = "REUSED"
            result.steps[-1].detail = "Existing result/evidence inputs reused"
        elif all(r.route == "production_inputs" for r in routing.routes):
            success = False
            result.steps[-1].status = "BLOCKED"
            result.steps[-1].failure_code = FailureCode.PRODUCTION_MANIFEST_REQUIRED.value
            result.steps[-1].detail = "BAM/FASTQ execution requires an explicit production_manifest with approved stage commands"
        elif len(routing.routes) == 1 and routing.routes[0].route in {"sv_wgs", "sv_wes"}:
            command_result = _run_sv_only(args, routing, layout)
            success = bool(command_result["ok"])
            result.steps[-1].status = "PASS" if success else "FAILED"
            result.steps[-1].outputs = {"log": command_result.get("log", "")}
        elif any(r.route in {"sv_wgs", "sv_wes"} for r in routing.routes):
            command_result = _run_standard(args, routing, layout)
            success = bool(command_result["ok"])
            result.steps[-1].status = "PASS" if success else "FAILED"
            result.steps[-1].detail = "Raw SV was normalized through the SV entry adapter and merged with the other entry routes"
            result.steps[-1].outputs = {"log": command_result.get("log", "")}
        else:
            command_result = _run_standard(args, routing, layout)
            success = bool(command_result["ok"])
            result.steps[-1].status = "PASS" if success else "FAILED"
            result.steps[-1].outputs = {"log": command_result.get("log", "")}
    except Exception as exc:
        success = False
        result.steps[-1].status = "FAILED"
        result.steps[-1].detail = str(exc)
        result.steps[-1].failure_code = FailureCode.PIPELINE_STAGE_FAILED.value

    if not success:
        result.blocking_issues.append(result.steps[-1].failure_code or FailureCode.PIPELINE_STAGE_FAILED.value)
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
        raw_events=routing.inputs.get("raw_events"),
        raw_peptides=routing.inputs.get("raw_peptides"),
        expression_evidence=routing.inputs.get("expression_tsv"),
        rna_junction_evidence=routing.inputs.get("rna_evidence_tsv"),
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
        result.steps[-1].failure_code = FailureCode.CONSENSUS_RANKING_FAILED.value
        result.blocking_issues.append(FailureCode.CONSENSUS_RANKING_FAILED.value)
        result.finish("FAILED").write(layout.skill_result)
        return result.to_dict()
    result.steps[-1].status = "PASS" if ranking["status"] == "PASS" else "REUSED"
    result.steps[-1].outputs = artifacts
    result.outputs.update(artifacts)
    evidence_for_consensus = artifacts.get("consensus_peptides") or artifacts.get("comprehensive_evidence")
    consensus_outputs = build_tool_consensus(routing.inputs, layout.pipeline / "tool_consensus", evidence_path=evidence_for_consensus)
    result.outputs.update({f"consensus_{key.removesuffix('.tsv')}": value for key, value in consensus_outputs.items()})
    if artifacts.get("all_tool_results"):
        enrich_all_tool_results(artifacts["all_tool_results"], consensus_outputs["tool_consensus_summary.tsv"])
    output_manifest = write_named_output_manifest(result.outputs, layout.root / "output_manifest.json")
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
    persist_result_state(layout.skill_result, payload)
    return payload
