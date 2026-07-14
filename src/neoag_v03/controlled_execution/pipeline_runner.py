from __future__ import annotations

import argparse
import os
import subprocess
import sys
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .audit import AuditLogger, run_command
from .doctor import run_doctor
from .io_utils import ensure_dir, markdown_table, now_iso, sha256_file, write_json, write_tsv
from .manifests import infer_sample_id, load_manifest, manifest_paths, normalize_tools_manifest, tool_container_digests, validate_manifests
from .version import ENHANCEMENT_VERSION


@dataclass
class PipelineStep:
    step_id: str
    name: str
    status: str = "PENDING"
    mode: str = "dry_run"
    reason: str = ""
    command: str = ""
    outputs: dict[str, str] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    failure_reason: str = ""
    requires_execution: bool = False


@dataclass
class PipelineRun:
    run_id: str
    sample_id: str
    started_at: str
    status: str
    profile: str
    dry_run: bool
    steps: list[PipelineStep]
    output_dir: str
    warnings: list[str] = field(default_factory=list)
    finished_at: str = ""


def _manifest_file_hashes(data: dict[str, Any]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for key, val in manifest_paths(data):
        p = Path(os.path.expandvars(os.path.expanduser(val)))
        if p.is_file():
            try:
                rows.append({"key": key, "path": str(p), "sha256": sha256_file(p), "size_bytes": str(p.stat().st_size)})
            except Exception as exc:
                rows.append({"key": key, "path": str(p), "sha256": "", "size_bytes": "", "error": str(exc)})
    return rows


def _git_sha(project_root: str | Path) -> str:
    try:
        proc = subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(project_root), text=True, capture_output=True, timeout=10)
        return proc.stdout.strip() if proc.returncode == 0 else ""
    except Exception:
        return ""


def _tool_versions(tools_data: dict[str, Any]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for name, spec in normalize_tools_manifest(tools_data).items():
        exe = str(spec.get("executable") or spec.get("path") or "")
        mode = str(spec.get("mode") or "")
        rows.append({"tool": name, "mode": mode, "executable": exe, "version": "UNASSESSED", "version_command": f"{exe} --version" if exe else ""})
    return rows


def _find_existing(sample_data: dict[str, Any], *keys: str) -> str | None:
    def walk(obj: Any, trail: list[str]) -> str | None:
        if isinstance(obj, dict):
            for k, v in obj.items():
                low = k.lower()
                if low in keys and isinstance(v, str) and v:
                    return v
                r = walk(v, trail + [k])
                if r:
                    return r
        elif isinstance(obj, list):
            for v in obj:
                if isinstance(v, str) and v:
                    # caller decides when list-like key is wanted
                    continue
                r = walk(v, trail)
                if r:
                    return r
        return None
    return walk(sample_data, [])


def _exists(path: str | None) -> bool:
    if not path:
        return False
    return Path(os.path.expandvars(os.path.expanduser(path))).exists()


def plan_pipeline(sample_data: dict[str, Any], *, profile: str = "local") -> list[PipelineStep]:
    steps: list[PipelineStep] = []
    def add(step_id: str, name: str, reason: str, mode: str = "planned", requires_execution: bool = False):
        steps.append(PipelineStep(step_id=step_id, name=name, status="PLANNED", mode=mode, reason=reason, requires_execution=requires_execution))

    add("01", "input-qc", "Validate sample manifest and available evidence")
    add("02", "doctor", "Check package, tools, references and release boundary")
    hla = _find_existing(sample_data, "hla_typing", "hla", "hla_tsv", "hla_file")
    hla_loh = _find_existing(sample_data, "hla_loh", "hla_loh_tsv", "spechla", "lohhla")
    raw_events = _find_existing(sample_data, "raw_events")
    raw_peptides = _find_existing(sample_data, "raw_peptides")
    presentation = _find_existing(sample_data, "presentation", "presentation_evidence")
    ranked_recommendation = _find_existing(sample_data, "ranked_peptides_recommendation", "recommendation", "ranked_peptides")
    ranked_netmhcpan42 = _find_existing(sample_data, "ranked_peptides_netmhcpan42", "netmhcpan42")
    evidence_report = _find_existing(sample_data, "evidence_report", "report", "html_report")
    somatic_vcf = _find_existing(sample_data, "somatic_vcf", "vcf", "snv_vcf")
    sv_vcf = _find_existing(sample_data, "sv_vcf", "manta_vcf", "gridss_vcf")
    fusion = _find_existing(sample_data, "fusion_tsv", "fusions", "fusion")
    expression = _find_existing(sample_data, "expression_tsv", "expression", "gene_expression", "tpm")
    purity = _find_existing(sample_data, "purity_tsv", "purity")
    cnv = _find_existing(sample_data, "cnv_segments", "cnv_tsv", "cnv")

    has_ranked_declared = bool(ranked_recommendation or ranked_netmhcpan42 or evidence_report)
    has_ranked_results = _exists(ranked_recommendation) or _exists(ranked_netmhcpan42) or _exists(evidence_report)
    ranked_mode = "result_review" if has_ranked_results else ("optional_unassessed" if has_ranked_declared else "planned")

    add("03", "hla-typing-hla-consensus", "HLA typing / HLA consensus from supplied HLA typing or typing tools", "skip_if_existing" if _exists(hla) else "requires_execution", requires_execution=not _exists(hla))
    add("04", "hla-loh-consensus", "HLA LOH consensus from SpecHLA/LOHHLA/FACETS/PURPLE evidence", "skip_if_existing" if _exists(hla_loh) else "optional_unassessed")
    if _exists(raw_events) and _exists(raw_peptides):
        add("05", "raw-intermediates-snv-dna-sv-fusion-splice", "Use supplied raw_events/raw_peptides", "skip_compute")
    elif _exists(somatic_vcf):
        add("05", "raw-intermediates-snv-dna-sv-fusion-splice", "Build raw_events/raw_peptides from somatic VCF and VEP annotation", "requires_execution", requires_execution=True)
    elif _exists(sv_vcf):
        add("05", "raw-intermediates-snv-dna-sv-fusion-splice", "Build raw_events/raw_peptides from DNA SV VCF", "requires_execution", requires_execution=True)
    elif _exists(fusion):
        add("05", "raw-intermediates-snv-dna-sv-fusion-splice", "Build raw_events/raw_peptides from fusion/splice evidence", "requires_execution", requires_execution=True)
    elif has_ranked_declared:
        add("05", "raw-intermediates-snv-dna-sv-fusion-splice", "Raw intermediates unavailable, but ranked/report outputs are declared; continue in result-review mode and mark missing evidence as UNASSESSED", "unassessed_existing_results")
    else:
        add("05", "raw-intermediates-snv-dna-sv-fusion-splice", "No raw_events/raw_peptides, variant source, ranked peptides or report found", "blocked")
    if _exists(presentation) or _exists(ranked_netmhcpan42) or _exists(ranked_recommendation):
        presentation_mode = "skip_if_existing"
        presentation_requires_execution = False
    elif presentation or ranked_netmhcpan42 or ranked_recommendation:
        presentation_mode = "optional_unassessed"
        presentation_requires_execution = False
    else:
        presentation_mode = "requires_execution"
        presentation_requires_execution = True
    add("06", "presentation-prediction", "Use existing presentation/ranked evidence or run HLA binding/presentation tools", presentation_mode, requires_execution=presentation_requires_execution)
    add("07", "expression-rna-evidence", "Integrate gene TPM, RNA alt reads and junction evidence", "skip_if_existing" if _exists(expression) else "optional_unassessed")
    add("08", "purity-cnv-ccf", "Build CCF/clonality from purity and CNV", "requires_execution" if _exists(purity) or _exists(cnv) else "optional_unassessed", requires_execution=bool(_exists(purity) or _exists(cnv)))
    add("09", "appm-immune-escape", "Build APPM and immune escape peptide flags", ranked_mode)
    add("10", "safety-gate", "Build normal-proteome/ligandome/junction safety evidence", ranked_mode)
    add("11", "recommendation-ranking", "Score candidates, or compare supplied recommendation and NetMHCpan ranked outputs", ranked_mode)
    add("12", "validation-plan", "Generate short peptide / long peptide / minigene validation plan", ranked_mode)
    add("13", "reports", "Generate patient and technical reports, or review supplied evidence report", ranked_mode)
    add("14", "audit-provenance", "Write output manifest, provenance and pipeline summary")
    return steps


def _write_plan_markdown(run: PipelineRun, outdir: Path) -> None:
    rows = [s.__dict__ for s in run.steps]
    md = [
        "# NeoAg pipeline-full plan",
        "",
        f"Run ID: `{run.run_id}`",
        f"Sample ID: `{run.sample_id}`",
        f"Status: **{run.status}**",
        f"Dry run: `{run.dry_run}`",
        "",
        "## Steps",
        markdown_table(rows, ["step_id", "name", "status", "mode", "requires_execution", "reason", "failure_reason"], max_rows=100),
        "",
        "## Warnings",
    ]
    if run.warnings:
        md.extend(f"- {w}" for w in run.warnings)
    else:
        md.append("No warnings.")
    (outdir / "pipeline_plan.md").write_text("\n".join(md) + "\n", encoding="utf-8")


def run_pipeline_full(
    *,
    sample_manifest: str | Path,
    outdir: str | Path,
    tools_manifest: str | Path | None = None,
    reference_manifest: str | Path | None = None,
    project_root: str | Path = ".",
    profile: str = "local",
    dry_run: bool = True,
    allow_partial: bool = True,
    run_doctor_first: bool = True,
    run_demo_for_fixture: bool = False,
) -> PipelineRun:
    od = ensure_dir(outdir)
    logger = AuditLogger(od / "audit_log.jsonl")
    sample_res = load_manifest(sample_manifest)
    reference_res = load_manifest(reference_manifest)
    tools_res = load_manifest(tools_manifest)
    sample_data = sample_res.data
    reference_data = reference_res.data
    tools_data = tools_res.data
    sample_id = infer_sample_id(sample_data)
    run_id = f"neoag-{sample_id}-{uuid.uuid4().hex[:10]}"
    logger.log("pipeline.start", "START", case_id=sample_id, request_id=run_id, metadata={"dry_run": dry_run})
    steps = plan_pipeline(sample_data, profile=profile)
    warnings: list[str] = []
    if not sample_res.exists or sample_res.error:
        warnings.append(f"sample_manifest_problem: {sample_res.error or 'not found'}")
    file_hashes = _manifest_file_hashes(sample_data)
    reference_hashes = _manifest_file_hashes(reference_data)
    tool_versions = _tool_versions(tools_data)
    container_digests = tool_container_digests(tools_data)
    manifest_validation = validate_manifests(sample_data, reference_data, tools_data)
    write_tsv(od / "input_file_hashes.tsv", file_hashes, ["key", "path", "sha256", "size_bytes", "error"])
    write_tsv(od / "reference_hashes.tsv", reference_hashes, ["key", "path", "sha256", "size_bytes", "error"])
    write_tsv(od / "tool_versions.tsv", tool_versions, ["tool", "mode", "executable", "version", "version_command"])
    write_tsv(od / "container_digests.tsv", container_digests, ["tool", "mode", "image", "container_digest", "executable"])
    write_tsv(od / "manifest_validation.tsv", manifest_validation, ["manifest", "field", "status", "message"])

    if any(r.get("status") == "FAIL" for r in manifest_validation):
        warnings.append("manifest_validation_has_failures")

    if run_doctor_first:
        doc_out = od / "doctor"
        doc = run_doctor(
            project_root=project_root,
            outdir=doc_out,
            tools_manifest=tools_manifest,
            reference_manifest=reference_manifest,
            sample_manifest=sample_manifest,
            profile=profile,
            run_demo=False,
            run_pytest=False,
            run_nextflow=False,
            release_audit=True,
            allow_execute=not dry_run,
        )
        steps[1].status = doc.status
        steps[1].outputs = doc.outputs
        if doc.status in {"BLOCKED", "UNSAFE"} and not allow_partial:
            for s in steps[2:]:
                s.status = "SKIPPED"
                s.failure_reason = f"doctor_status={doc.status}"
            status = "BLOCKED"
            run = PipelineRun(run_id, sample_id, now_iso(), status, profile, dry_run, steps, str(od), warnings, now_iso())
            _finalize_run(run, od, sample_manifest=sample_manifest, reference_manifest=reference_manifest, tools_manifest=tools_manifest, project_root=project_root, input_hashes=file_hashes, reference_hashes=reference_hashes, tool_versions=tool_versions, container_digests=container_digests, manifest_validation=manifest_validation)
            return run
        if doc.status in {"PARTIAL", "BLOCKED", "UNSAFE"}:
            warnings.append(f"doctor_status={doc.status}; pipeline continues only in partial/unassessed mode")

    # Phase 2 POC: deterministic safe execution. Heavy bioinformatics steps are planned
    # but not executed unless demo fixture is requested. Existing result review paths are
    # surfaced to downstream Skills/Agents.
    for step in steps:
        if step.step_id in {"01", "02"}:
            if step.status == "PLANNED":
                step.status = "PASS"
            continue
        if step.mode == "blocked":
            step.status = "BLOCKED"
            step.failure_reason = step.reason
            if not allow_partial:
                break
        elif step.mode == "requires_execution":
            step.status = "DRY_RUN" if dry_run else "SKIPPED"
            step.failure_reason = "requires external tool/workflow execution; use Gateway/HPC runner after Doctor is READY"
        elif step.mode in {"optional_unassessed", "requires_input", "unassessed_existing_results"}:
            step.status = "UNASSESSED"
        else:
            step.status = "PASS"

    if run_demo_for_fixture:
        demo_out = od / "demo_fixture"
        res = run_command([sys.executable, "-m", "neoag_v03.cli", "run-demo", "--sample-id", sample_id or "PIPELINE", "--outdir", str(demo_out)], cwd=project_root, timeout=300, logger=logger, risk_level="LOW", allow_execute=not dry_run)
        status = "DRY_RUN" if res.get("dry_run") else ("PASS" if res.get("ok") else "FAIL")
        steps.append(PipelineStep(step_id="15", name="optional-run-demo-fixture", status=status, mode="execute" if not dry_run else "dry_run", reason="Optional fixture smoke run", command=res.get("cmd", ""), outputs={"demo_outdir": str(demo_out)}, failure_reason="" if res.get("ok") else res.get("stderr", "")))

    if any(s.status == "BLOCKED" for s in steps):
        status = "BLOCKED"
    elif any(s.status in {"FAIL"} for s in steps):
        status = "FAIL"
    elif any(s.status in {"UNASSESSED", "DRY_RUN", "SKIPPED"} for s in steps):
        status = "PARTIAL" if not dry_run else "DRY_RUN"
    else:
        status = "PASS"
    run = PipelineRun(run_id, sample_id, now_iso(), status, profile, dry_run, steps, str(od), warnings, now_iso())
    _finalize_run(run, od, sample_manifest=sample_manifest, reference_manifest=reference_manifest, tools_manifest=tools_manifest, project_root=project_root, input_hashes=file_hashes, reference_hashes=reference_hashes, tool_versions=tool_versions, container_digests=container_digests, manifest_validation=manifest_validation)
    logger.log("pipeline.finish", status, case_id=sample_id, request_id=run_id, metadata={"warnings": warnings})
    return run


def _finalize_run(run: PipelineRun, outdir: Path, *, sample_manifest: str | Path, reference_manifest: str | Path | None, tools_manifest: str | Path | None, project_root: str | Path, input_hashes: list[dict[str, str]], reference_hashes: list[dict[str, str]], tool_versions: list[dict[str, str]], container_digests: list[dict[str, str]], manifest_validation: list[dict[str, str]]) -> None:
    step_rows = [s.__dict__ for s in run.steps]
    write_tsv(outdir / "pipeline_status.tsv", step_rows, ["step_id", "name", "status", "mode", "requires_execution", "reason", "command", "failure_reason"])
    output_manifest = {"run_id": run.run_id, "status": run.status, "outdir": str(outdir), "files": sorted(str(p.relative_to(outdir)) for p in outdir.rglob("*") if p.is_file())}
    write_json(outdir / "output_manifest.json", output_manifest)
    run_manifest = {
        "run_id": run.run_id,
        "sample_id": run.sample_id,
        "pipeline_version": ENHANCEMENT_VERSION,
        "git_sha": _git_sha(project_root),
        "profile": run.profile,
        "dry_run": run.dry_run,
        "sample_manifest": str(sample_manifest),
        "reference_manifest": str(reference_manifest or ""),
        "tools_manifest": str(tools_manifest or ""),
        "input_hashes": input_hashes,
        "reference_hashes": reference_hashes,
        "tool_versions": tool_versions,
        "container_digests": container_digests,
        "manifest_validation": manifest_validation,
        "start_time": run.started_at,
        "end_time": run.finished_at,
        "status": run.status,
        "warnings": run.warnings,
        "steps": step_rows,
        "output_manifest": str(outdir / "output_manifest.json"),
    }
    write_json(outdir / "run_manifest.json", run_manifest)
    _write_plan_markdown(run, outdir)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="NeoAg pipeline-full: controlled-execution manifest-driven full pipeline runner")
    ap.add_argument("--sample-manifest", required=True)
    ap.add_argument("--tools-manifest")
    ap.add_argument("--reference-manifest")
    ap.add_argument("--project-root", default=".")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--profile", default="local")
    ap.add_argument("--execute", action="store_true", help="Allow safe local execution. Heavy external steps are still planned unless implemented by Gateway/HPC runner.")
    ap.add_argument("--strict", action="store_true", help="Stop when Doctor is BLOCKED/UNSAFE instead of continuing partial mode")
    ap.add_argument("--run-demo-fixture", action="store_true")
    args = ap.parse_args(argv)
    run = run_pipeline_full(
        sample_manifest=args.sample_manifest,
        tools_manifest=args.tools_manifest,
        reference_manifest=args.reference_manifest,
        project_root=args.project_root,
        outdir=args.outdir,
        profile=args.profile,
        dry_run=not args.execute,
        allow_partial=not args.strict,
        run_demo_for_fixture=args.run_demo_fixture,
    )
    print(f"NeoAg pipeline-full status: {run.status}")
    print(f"  run_id: {run.run_id}")
    print(f"  output_dir: {run.output_dir}")
    print(f"  run_manifest: {Path(run.output_dir) / 'run_manifest.json'}")
    return 0 if run.status in {"PASS", "PARTIAL", "DRY_RUN"} else 2


if __name__ == "__main__":
    raise SystemExit(main())
