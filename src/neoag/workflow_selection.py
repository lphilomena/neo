from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from .controlled_execution.manifests import load_manifest, normalize_tools_manifest


@dataclass(frozen=True)
class PlannedStage:
    order: int
    name: str
    status: str
    required: bool
    risk_level: str
    runner: str
    reason: str


@dataclass(frozen=True)
class WorkflowSelection:
    sample_id: str
    workflow: str
    status: str
    entry_mode: str
    stages: list[PlannedStage]
    detected_inputs: dict[str, dict[str, Any]]
    missing_required: list[str]
    recommended_command: str
    execution_allowed: bool = False


def _mapping(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _pick(*values: Any) -> str:
    for value in values:
        if isinstance(value, str) and value.strip():
            return value.strip()
        if isinstance(value, list):
            for item in value:
                if isinstance(item, str) and item.strip():
                    return item.strip()
    return ""


def _resolve(path: str, manifest_dir: Path) -> Path | None:
    if not path:
        return None
    candidate = Path(path).expanduser()
    return candidate if candidate.is_absolute() else (manifest_dir / candidate).resolve()


def _input_state(path: str, manifest_dir: Path, *, show_paths: bool) -> dict[str, Any]:
    resolved = _resolve(path, manifest_dir)
    return {
        "provided": bool(path),
        "exists": bool(resolved and resolved.exists()),
        "path": str(resolved) if show_paths and resolved else (resolved.name if resolved else ""),
    }


def _tool_declared(tools: dict[str, dict[str, Any]], *names: str) -> bool:
    normalized = {name.lower().replace("-", "_") for name in tools}
    return any(name.lower().replace("-", "_") in normalized for name in names)


def select_workflow(
    sample_manifest: str | Path,
    *,
    tools_manifest: str | Path | None = None,
    reference_manifest: str | Path | None = None,
    outdir: str | Path = "work/workflow_selection",
    show_paths: bool = False,
) -> WorkflowSelection:
    sample_result = load_manifest(sample_manifest)
    if not sample_result.exists or sample_result.error:
        raise ValueError(f"Cannot load sample manifest: {sample_result.error or 'file_not_found'}")

    tools_result = load_manifest(tools_manifest)
    references_result = load_manifest(reference_manifest)
    sample = sample_result.data
    inputs = _mapping(sample.get("inputs"))
    tumor = _mapping(sample.get("tumor"))
    normal = _mapping(sample.get("normal"))
    manifest_dir = Path(sample_result.path).resolve().parent

    paths = {
        "ranked_peptides": _pick(
            inputs.get("ranked_peptides_evidence_consensus"),
            inputs.get("ranked_peptides_recommendation"),
            inputs.get("ranked_peptides"),
        ),
        "raw_events": _pick(inputs.get("raw_events")),
        "raw_peptides": _pick(inputs.get("raw_peptides")),
        "somatic_vcf": _pick(inputs.get("somatic_vcf"), inputs.get("vcf")),
        "hla_typing": _pick(inputs.get("hla_typing"), inputs.get("hla_file")),
        "tumor_dna_bam": _pick(tumor.get("dna_bam"), tumor.get("wgs_bam"), tumor.get("wes_bam")),
        "normal_dna_bam": _pick(normal.get("dna_bam"), normal.get("wgs_bam"), normal.get("wes_bam")),
        "tumor_rna_bam": _pick(tumor.get("rna_bam"), tumor.get("wts_bam")),
        "rna_fastq": _pick(tumor.get("rna_fastq1"), tumor.get("wts_fastq1"), inputs.get("rna_fastq1")),
        "expression": _pick(inputs.get("expression_tsv"), inputs.get("gene_tpm")),
        "fusion": _pick(inputs.get("fusion_tsv")),
        "splice": _pick(inputs.get("splice_tsv"), inputs.get("junction_tsv")),
        "purity": _pick(inputs.get("purity_tsv")),
        "hla_loh": _pick(inputs.get("hla_loh")),
    }
    states = {name: _input_state(path, manifest_dir, show_paths=show_paths) for name, path in paths.items()}
    has = {name: state["provided"] and state["exists"] for name, state in states.items()}

    tools = normalize_tools_manifest(tools_result.data) if tools_result.exists and not tools_result.error else {}
    paired_bam = has["tumor_dna_bam"] and has["normal_dna_bam"]
    hla_source = has["hla_typing"] or has["normal_dna_bam"] or has["rna_fastq"]

    if has["ranked_peptides"]:
        workflow, entry_mode = "result_review", "ranked_results"
    elif has["raw_events"] and has["raw_peptides"]:
        workflow, entry_mode = "intermediates_to_ranking", "intermediates"
    elif has["somatic_vcf"] and has["hla_typing"]:
        workflow, entry_mode = "vcf_to_ranking", "somatic_vcf"
    elif has["somatic_vcf"] and hla_source:
        workflow, entry_mode = "vcf_with_hla_typing", "somatic_vcf"
    elif paired_bam:
        workflow, entry_mode = "paired_bam_production", "paired_bam"
    else:
        workflow, entry_mode = "insufficient_inputs", "unknown"

    stages: list[PlannedStage] = []

    def add(name: str, status: str, required: bool, risk: str, runner: str, reason: str) -> None:
        stages.append(PlannedStage(len(stages) + 1, name, status, required, risk, runner, reason))

    add("input_qc", "READY", True, "LOW", "neoag-input-qc", "sample manifest loaded")
    doctor_ready = tools_result.exists and references_result.exists and not tools_result.error and not references_result.error
    add(
        "doctor",
        "READY" if doctor_ready else "NEEDS_MANIFESTS",
        True,
        "LOW",
        "neoag-doctor",
        "tools and reference manifests supplied" if doctor_ready else "supply tools and reference manifests before execution",
    )

    if workflow == "result_review":
        add("result_review", "READY", True, "LOW", "neoag-result-inspector", "ranked peptide results already exist")
        add("ranking_compare", "OPTIONAL", False, "LOW", "neoag-ranking-compare", "compare weighted and evidence-consensus rankings when both exist")
    else:
        if has["hla_typing"]:
            add("hla_typing", "REUSE", True, "LOW", "neoag-hla-typing-run-and-compare", "existing HLA typing input")
        elif hla_source:
            add("hla_typing", "PLANNED", True, "MEDIUM", "neoag-hla-typing-run-and-compare", "derive consensus HLA before presentation prediction")
        else:
            add("hla_typing", "BLOCKED", True, "MEDIUM", "neoag-hla-typing-run-and-compare", "no HLA result or suitable normal/RNA input")

        if has["raw_events"] and has["raw_peptides"]:
            add("candidate_ingestion", "REUSE", True, "LOW", "neoag-peptide-csv", "canonical event and peptide tables exist")
        elif has["somatic_vcf"]:
            add("candidate_ingestion", "PLANNED", True, "MEDIUM", "neoag-vcf", "parse and annotate supplied somatic VCF")
        elif paired_bam:
            add("variant_calling", "PLANNED", True, "HIGH", "neoag-pipeline-full", "paired BAM variant calling requires approval")
        else:
            add("candidate_ingestion", "BLOCKED", True, "MEDIUM", "neoag-vcf", "no canonical intermediates, somatic VCF, or paired BAM")

        if has["purity"]:
            add("purity_cnv_consensus", "REUSE", False, "LOW", "neoag-purity-cnv-run-and-review", "existing purity evidence")
        elif paired_bam and any(_tool_declared(tools, name) for name in ("facets", "purple", "sequenza", "ascat")):
            add("purity_cnv_consensus", "PLANNED", False, "HIGH", "neoag-purity-cnv-run-and-review", "run declared purity tools; preserve per-tool failures and UNASSESSED values")
        else:
            add("purity_cnv_consensus", "UNASSESSED", False, "MEDIUM", "neoag-purity-cnv-run-and-review", "no reusable result or runnable paired-BAM tool declaration")

        if has["hla_loh"]:
            add("hla_loh_consensus", "REUSE", False, "LOW", "neoag-hla-loh-appm-review", "existing HLA LOH evidence")
        elif paired_bam and (_tool_declared(tools, "lohhla") or _tool_declared(tools, "spechla")):
            add("hla_loh_consensus", "PLANNED", False, "HIGH", "neoag-hla-loh-appm-review", "run declared HLA LOH tools after HLA/purity review")
        else:
            add("hla_loh_consensus", "UNASSESSED", False, "MEDIUM", "neoag-hla-loh-appm-review", "HLA LOH evidence unavailable")

        rna_available = has["tumor_rna_bam"] or has["rna_fastq"]
        add(
            "rna_expression",
            "REUSE" if has["expression"] else ("PLANNED" if rna_available else "UNASSESSED"),
            False,
            "MEDIUM",
            "neoag-rna-fastq-to-tpm",
            "existing expression evidence" if has["expression"] else ("tumor RNA input available" if rna_available else "no tumor RNA input"),
        )
        if has["fusion"] or has["splice"]:
            add("fusion_splice", "REUSE", False, "LOW", "neoag-fusion-rna-run", "existing fusion/splice evidence")
        elif rna_available:
            add("fusion_splice", "PLANNED", False, "HIGH", "neoag-fusion-rna-run", "RNA discovery execution requires approval")
        else:
            add("fusion_splice", "UNASSESSED", False, "MEDIUM", "neoag-fusion-rna-run", "no tumor RNA input")
        add("unified_ranking", "PLANNED", True, "MEDIUM", "open-neo-run", "merge candidate sources and preserve missing evidence as UNASSESSED")

    missing_required = [stage.name for stage in stages if stage.required and stage.status in {"BLOCKED", "NEEDS_MANIFESTS"}]
    if any(stage.status == "BLOCKED" for stage in stages if stage.required):
        status = "BLOCKED"
    elif missing_required or any(stage.status == "UNASSESSED" for stage in stages):
        status = "PARTIAL"
    else:
        status = "READY_TO_PLAN"

    command = (
        f"neoag-production-run --manifest configs/workflows/production_workflow.private.toml --outdir {Path(outdir)}"
        if workflow != "result_review"
        else f"neoag-skill run neoag-result-inspector --outdir {Path(outdir) / 'review'}"
    )
    return WorkflowSelection(
        sample_id=str(sample.get("sample_id") or sample.get("case_id") or "UNKNOWN_SAMPLE"),
        workflow=workflow,
        status=status,
        entry_mode=entry_mode,
        stages=stages,
        detected_inputs=states,
        missing_required=missing_required,
        recommended_command=command,
    )


def write_selection(selection: WorkflowSelection, outdir: str | Path) -> dict[str, str]:
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    json_path = output / "workflow_selection.json"
    md_path = output / "workflow_selection.md"
    payload = asdict(selection)
    json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    lines = [
        "# NeoAg workflow selection",
        "",
        f"- Sample: `{selection.sample_id}`",
        f"- Selected workflow: `{selection.workflow}`",
        f"- Entry mode: `{selection.entry_mode}`",
        f"- Readiness: `{selection.status}`",
        "- Execution: disabled; this command only writes a plan",
        "",
        "| Order | Stage | Status | Required | Risk | Runner | Reason |",
        "|---:|---|---|---|---|---|---|",
    ]
    for stage in selection.stages:
        lines.append(
            f"| {stage.order} | {stage.name} | {stage.status} | {str(stage.required).lower()} | "
            f"{stage.risk_level} | {stage.runner} | {stage.reason} |"
        )
    lines.extend(["", "## Recommended dry-run entry", "", f"`{selection.recommended_command}`", ""])
    if selection.missing_required:
        lines.extend(["## Blocking or required preparation", "", *[f"- {item}" for item in selection.missing_required], ""])
    md_path.write_text("\n".join(lines), encoding="utf-8")
    return {"json": str(json_path), "markdown": str(md_path)}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Select a NeoAg workflow from manifests without executing tools")
    parser.add_argument("--sample-manifest", required=True)
    parser.add_argument("--tools-manifest")
    parser.add_argument("--reference-manifest")
    parser.add_argument("--outdir", default="work/workflow_selection")
    parser.add_argument("--show-paths", action="store_true", help="Include resolved paths in local plan outputs")
    args = parser.parse_args(argv)
    selection = select_workflow(
        args.sample_manifest,
        tools_manifest=args.tools_manifest,
        reference_manifest=args.reference_manifest,
        outdir=args.outdir,
        show_paths=args.show_paths,
    )
    outputs = write_selection(selection, args.outdir)
    print(json.dumps({"status": selection.status, "workflow": selection.workflow, "outputs": outputs}, ensure_ascii=False))


if __name__ == "__main__":
    main()
