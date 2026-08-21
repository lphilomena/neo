from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
import tomllib
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from .adapters.pvactools_parser import parse_pvactools_outputs
from .provenance import (
    CONFLICT_FIELDS,
    PROVENANCE_FIELDS,
    merge_rows_preserving_provenance,
)
from .schemas import EVENT_FIELDS, PEPTIDE_FIELDS
from .utils import read_tsv, write_tsv


@dataclass
class StageResult:
    name: str
    status: str
    required: bool
    source: str = ""
    command: str = ""
    log: str = ""
    outputs: dict[str, Any] = field(default_factory=dict)
    message: str = ""


@dataclass
class ProductionResult:
    sample_id: str
    status: str
    outdir: str
    dry_run: bool
    stages: list[StageResult]
    source_status: str = "UNASSESSED"
    detected_sources: list[str] = field(default_factory=list)
    missing_sources: list[str] = field(default_factory=list)
    generated_config: str = ""
    final_outdir: str = ""
    provenance_outputs: dict[str, str] = field(default_factory=dict)


def _variant_candidate_for_mtwt(row: dict[str, str]) -> bool:
    tokens = " ".join(
        str(row.get(key) or "")
        for key in ("event_type", "mutation_source", "source_chain_track", "peptide_consequence")
    ).upper()
    if any(token in tokens for token in ("FUSION", "SPLICE", "SV")):
        return False
    return any(token in tokens for token in ("SNV", "INDEL", "MISSENSE", "FRAMESHIFT", "VARIANT"))


def _assess_mtwt_completeness(ranked_peptides: Path) -> tuple[bool, dict[str, Any]]:
    summary: dict[str, Any] = {
        "ranked_peptides": str(ranked_peptides),
        "variant_rows": 0,
        "with_wildtype_peptide": 0,
        "with_wt_binding": 0,
        "with_mutant_specificity": 0,
        "status": "NOT_ASSESSED",
    }
    if not ranked_peptides.is_file():
        summary["status"] = "MISSING_RANKED_PEPTIDES"
        return False, summary
    wt_binding_fields = (
        "netmhcpan_wt_rank_el",
        "netmhcpan_wt_affinity",
        "mhcflurry_wt_affinity_percentile",
        "mhcflurry_wt_affinity",
        "mt_wt_el_rank_difference",
        "mt_wt_fold_change",
    )
    specificity_fields = ("mutant_specificity_status", "mutant_specificity_state")
    with ranked_peptides.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="	")
        for row in reader:
            if not _variant_candidate_for_mtwt(row):
                continue
            summary["variant_rows"] += 1
            if str(row.get("wildtype_peptide") or "").strip():
                summary["with_wildtype_peptide"] += 1
            if any(str(row.get(field) or "").strip() for field in wt_binding_fields):
                summary["with_wt_binding"] += 1
            if any(str(row.get(field) or "").strip().upper() not in {"", "UNASSESSED", "NA"} for field in specificity_fields):
                summary["with_mutant_specificity"] += 1
    if int(summary["variant_rows"]) == 0:
        summary["status"] = "NO_VARIANT_CANDIDATES"
        return True, summary
    ok = (
        int(summary["with_wildtype_peptide"]) > 0
        and int(summary["with_wt_binding"]) > 0
        and int(summary["with_mutant_specificity"]) > 0
    )
    summary["status"] = "PASS" if ok else "FAILED"
    if not ok:
        summary["message"] = (
            "SNV/InDel candidates require MT/WT assessment by default. "
            "Ensure VEP Wildtype/Frameshift plugins are configured via NEOAG_VEP_PLUGINS "
            "or refs.vep_plugins, then rerun candidate upstream and final ranking."
        )
    return ok, summary


def load_production_manifest(path: str | Path) -> dict[str, Any]:
    return tomllib.loads(Path(path).read_text(encoding="utf-8"))


def _expand(value: str, context: dict[str, str]) -> str:
    expanded = os.path.expanduser(os.path.expandvars(value))
    for key, replacement in context.items():
        expanded = expanded.replace("${" + key + "}", replacement).replace("$" + key, replacement)
    import re

    protected_values: dict[str, str] = {}

    def protect_shell_param(match: re.Match[str]) -> str:
        token = f"__OPEN_NEO_SHELL_PARAM_{len(protected_values)}__"
        protected_values[token] = match.group(0)
        return token

    protected = re.sub(r"\$\{[^{}]*\}", protect_shell_param, expanded)
    try:
        rendered = protected.format_map(context)
    except KeyError as exc:
        raise ValueError(f"Unknown production manifest placeholder: {exc.args[0]}") from exc
    for token, original in protected_values.items():
        rendered = rendered.replace(token, original)
    return rendered


def _expand_value(value: Any, context: dict[str, str]) -> Any:
    if isinstance(value, str):
        return _expand(value, context)
    if isinstance(value, list):
        return [_expand_value(item, context) for item in value]
    if isinstance(value, dict):
        return {str(key): _expand_value(item, context) for key, item in value.items()}
    return value


def _flatten_paths(outputs: dict[str, Any]) -> list[Path]:
    paths: list[Path] = []
    for value in outputs.values():
        values = value if isinstance(value, list) else [value]
        for item in values:
            if isinstance(item, str) and item:
                paths.append(Path(item))
    return paths


def _outputs_ready(outputs: dict[str, Any]) -> bool:
    paths = _flatten_paths(outputs)
    return bool(paths) and all(path.exists() for path in paths)


def _tool_group_audit(requirements, stages, results, output):
    """Audit cross-tool production requirements against sample outputs."""
    by_name = {row.name: row for row in results}
    rows = []
    passed = True
    for domain, raw_rule in requirements.items():
        rule = raw_rule if isinstance(raw_rule, dict) else {}
        tools = [str(value) for value in rule.get("tools", [])]
        minimum = int(rule.get("min_successful", len(tools)))
        successful, declared = [], []
        for tool in tools:
            key = tool.lower().replace("-", "_")
            matches = []
            for name, spec in stages.items():
                result = by_name.get(name)
                if not result:
                    continue
                configured = str(spec.get("tool") or "").lower().replace("-", "_")
                output_keys = {str(value).lower().replace("-", "_") for value in result.outputs}
                if configured == key or key in name.lower().replace("-", "_") or any(key in value for value in output_keys):
                    matches.append(result)
            if matches:
                declared.append(tool)
            for result in matches:
                if result.status not in {"PASS", "REUSED"}:
                    continue
                matching = {name: value for name, value in result.outputs.items() if key in str(name).lower().replace("-", "_")}
                if _outputs_ready(matching or result.outputs):
                    successful.append(tool)
                    break
        reasons = []
        undeclared = [tool for tool in tools if tool not in declared]
        if bool(rule.get("require_all_declared", True)) and undeclared:
            reasons.append("undeclared=" + ",".join(undeclared))
        if len(successful) < minimum:
            reasons.append(f"successful={len(successful)}<{minimum}")
        status = "FAIL" if reasons else "PASS"
        passed = passed and status == "PASS"
        rows.append({"domain": str(domain), "status": status, "required_tools": ",".join(tools), "declared_tools": ",".join(declared), "successful_tools": ",".join(successful), "min_successful": str(minimum), "reason": "; ".join(reasons) or "cross-tool requirement satisfied"})
    write_tsv(output, rows, ["domain", "status", "required_tools", "declared_tools", "successful_tools", "min_successful", "reason"])
    return passed, rows


def _ordered_stages(stages: dict[str, dict[str, Any]]) -> list[str]:
    ordered: list[str] = []
    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(name: str) -> None:
        if name in visited:
            return
        if name in visiting:
            raise ValueError(f"Production stage dependency cycle at: {name}")
        if name not in stages:
            raise ValueError(f"Unknown production stage dependency: {name}")
        visiting.add(name)
        for dependency in stages[name].get("depends_on") or []:
            visit(str(dependency))
        visiting.remove(name)
        visited.add(name)
        ordered.append(name)

    for stage_name in stages:
        visit(stage_name)
    return ordered


def _run_stage(
    name: str,
    spec: dict[str, Any],
    *,
    context: dict[str, str],
    logs_dir: Path,
    execute: bool,
    force: bool,
) -> StageResult:
    required = bool(spec.get("required", False))
    source = str(spec.get("source") or "")
    outputs = _expand_value(spec.get("outputs") or {}, context)
    command = _expand(str(spec.get("command") or ""), context).strip()
    log_path = logs_dir / f"{name}.log"

    if _outputs_ready(outputs) and not force:
        return StageResult(name, "REUSED", required, source, command, str(log_path), outputs)
    if not execute:
        status = "PLANNED" if command else ("BLOCKED" if required else "LOW_CONFIDENCE")
        message = "command planned" if command else "outputs missing and no command configured"
        return StageResult(name, status, required, source, command, str(log_path), outputs, message)
    if not command:
        status = "FAILED" if required else "LOW_CONFIDENCE"
        return StageResult(
            name,
            status,
            required,
            source,
            command,
            str(log_path),
            outputs,
            "outputs missing and no command configured",
        )

    logs_dir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        ["bash", "-lc", command],
        cwd=context["project_root"],
        text=True,
        capture_output=True,
    )
    log_path.write_text(
        proc.stdout + ("\n--- STDERR ---\n" if proc.stderr else "") + proc.stderr,
        encoding="utf-8",
    )
    if proc.returncode != 0:
        status = "FAILED" if required else "LOW_CONFIDENCE"
        return StageResult(
            name,
            status,
            required,
            source,
            command,
            str(log_path),
            outputs,
            f"command returned {proc.returncode}",
        )
    if not _outputs_ready(outputs):
        status = "FAILED" if required else "LOW_CONFIDENCE"
        return StageResult(
            name,
            status,
            required,
            source,
            command,
            str(log_path),
            outputs,
            "command completed but declared outputs are missing",
        )
    return StageResult(name, "PASS", required, source, command, str(log_path), outputs)


def _deduplicate(
    rows: list[dict[str, str]],
    fields: list[str],
    key_fields: tuple[str, ...],
) -> list[dict[str, str]]:
    """Compatibility wrapper around the provenance-preserving v0.4.4 merge."""

    merged, _, _ = merge_rows_preserving_provenance(
        rows,
        fields,
        key_fields,
        entity_type="production_entity",
    )
    return merged


def _annotate_stage_provenance(
    rows: list[dict[str, str]],
    *,
    stage: StageResult,
    source_file: str,
    entity_id_field: str,
) -> list[dict[str, str]]:
    annotated: list[dict[str, str]] = []
    for row_number, raw in enumerate(rows, 1):
        row = dict(raw)
        row["_provenance_stage_name"] = stage.name
        row["_provenance_stage_source"] = stage.source
        row["_provenance_source_file"] = source_file
        if not row.get("source_tools") and stage.source:
            row["source_tools"] = stage.source
        if not row.get("source_tool") and entity_id_field == "peptide_id" and stage.source:
            row["source_tool"] = stage.source
        if not row.get("source_file"):
            row["source_file"] = source_file
        if not row.get("source_row_number"):
            row["source_row_number"] = str(row_number)
        if not row.get("source_record_id"):
            row["source_record_id"] = str(row.get(entity_id_field) or f"{stage.name}:{row_number}")
        if not row.get("source_records"):
            row["source_records"] = row["source_record_id"]
        if not row.get("provenance_record_count"):
            row["provenance_record_count"] = "1"
        if not row.get("evidence_conflict_status"):
            row["evidence_conflict_status"] = "NONE"
        annotated.append(row)
    return annotated


def _candidate_rows(
    stage: StageResult,
    *,
    sample_id: str,
    profile: str,
    normalized_dir: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    outputs = stage.outputs
    raw_events = outputs.get("raw_events")
    raw_peptides = outputs.get("raw_peptides")
    if raw_events and raw_peptides and Path(raw_events).is_file() and Path(raw_peptides).is_file():
        events = _annotate_stage_provenance(
            read_tsv(raw_events),
            stage=stage,
            source_file=str(raw_events),
            entity_id_field="event_id",
        )
        peptides = _annotate_stage_provenance(
            read_tsv(raw_peptides),
            stage=stage,
            source_file=str(raw_peptides),
            entity_id_field="peptide_id",
        )
        return events, peptides

    pvac = outputs.get("pvac_files") or outputs.get("pvac_file") or outputs.get("pvac")
    if not pvac:
        return [], []
    paths = pvac if isinstance(pvac, list) else [pvac]
    paths = [str(path) for path in paths if Path(path).exists()]
    if not paths:
        return [], []
    stage_dir = normalized_dir / stage.name
    events_path = stage_dir / "raw_events.tsv"
    peptides_path = stage_dir / "raw_peptides.tsv"
    events, peptides = parse_pvactools_outputs(
        paths,
        sample_id,
        profile,
        events_path,
        peptides_path,
    )
    return (
        _annotate_stage_provenance(
            events,
            stage=stage,
            source_file=";".join(paths),
            entity_id_field="event_id",
        ),
        _annotate_stage_provenance(
            peptides,
            stage=stage,
            source_file=";".join(paths),
            entity_id_field="peptide_id",
        ),
    )


def _read_hla_alleles(run_cfg: dict[str, Any], stage_results: list[StageResult]) -> list[str]:
    configured = [str(value) for value in (run_cfg.get("hla_alleles") or []) if str(value).strip()]
    if configured:
        return configured
    hla_file = str(run_cfg.get("hla_file") or "")
    if not hla_file:
        for stage in stage_results:
            candidate = stage.outputs.get("hla_file")
            if isinstance(candidate, str) and candidate:
                hla_file = candidate
                break
    if not hla_file or not Path(hla_file).is_file():
        return []
    import re

    values = re.split(r"[\s,;]+", Path(hla_file).read_text(encoding="utf-8", errors="replace"))
    result: list[str] = []
    for value in values:
        allele = value.strip().strip('"').strip("'")
        if allele.upper().startswith("HLA-") and "*" in allele and allele not in result:
            result.append(allele)
    return result


def _write_final_config(
    path: Path,
    *,
    sample_id: str,
    profile: str,
    tools_stub: bool,
    immunogenicity_stub: bool,
    enabled_predictors: list[str],
    required_predictors: list[str],
    hla_alleles: list[str],
    raw_events: Path,
    raw_peptides: Path,
    evidence: dict[str, Any],
    reuse_immunogenicity_sources: list[str] | None = None,
) -> None:
    def toml_value(value: Any) -> str:
        if isinstance(value, bool):
            return str(value).lower()
        if isinstance(value, list):
            return "[" + ", ".join(json.dumps(str(item)) for item in value) + "]"
        return json.dumps(str(value))

    lines = [
        "[sample]",
        f"id = {toml_value(sample_id)}",
        f"profile = {toml_value(profile)}",
        "",
        "[tools]",
        f"stub = {toml_value(tools_stub)}",
        f"enabled = {toml_value(enabled_predictors)}",
        f"immunogenicity_stub = {toml_value(immunogenicity_stub)}",
        f"reuse_immunogenicity_sources = {toml_value(reuse_immunogenicity_sources or [])}",
        "",
        "[inputs]",
        'entry_mode = "intermediates"',
        f"raw_events = {toml_value(str(raw_events))}",
        f"raw_peptides = {toml_value(str(raw_peptides))}",
        f"hla_alleles = {toml_value(hla_alleles)}",
        f"required_presentation_predictors = {toml_value(required_predictors)}",
        "extract_appm_from_vcf = false",
    ]
    allowed_evidence = {
        "vep_appm",
        "expression",
        "transcript_expression",
        "rna_vaf",
        "rna_junction_tsv",
        "cross_site_rna_evidence",
        "hla_loh",
        "purity",
        "cnv",
        "normal_expression",
        "normal_hla_ligands",
        "reference_proteome",
        "normal_junctions",
        "evidence_consensus_rules",
        "netmhcpan",
        "mhcflurry",
        "netmhcstabpan",
        "netchop",
    }
    for key, value in evidence.items():
        if key in allowed_evidence and value and Path(str(value)).exists():
            lines.append(f"{key} = {toml_value(value)}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


_REUSABLE_IMMUNOGENICITY_OUTPUTS = {
    "prime_evidence": ("prime", "prime_evidence.tsv"),
    "bigmhc_evidence": ("bigmhc_im", "bigmhc_im_evidence.tsv"),
    "deepimmuno_evidence": ("deepimmuno", "deepimmuno_evidence.tsv"),
}


def _materialize_reusable_immunogenicity_outputs(
    final_outdir: Path, evidence: dict[str, Any]
) -> list[str]:
    """Materialize explicitly supplied normalized predictor evidence for run-full."""
    presentation = final_outdir / "presentation"
    reused: list[str] = []
    for evidence_key, (tool_key, filename) in _REUSABLE_IMMUNOGENICITY_OUTPUTS.items():
        raw_source = str(evidence.get(evidence_key) or "").strip()
        if not raw_source:
            continue
        source = Path(raw_source)
        if not source.is_file() or source.stat().st_size == 0:
            raise FileNotFoundError(f"Reusable {tool_key} evidence is missing or empty: {source}")
        presentation.mkdir(parents=True, exist_ok=True)
        destination = presentation / filename
        if source.resolve() != destination.resolve():
            shutil.copy2(source, destination)
        reused.append(tool_key)
    return reused


def _write_result(result: ProductionResult, outdir: Path) -> None:
    (outdir / "production_run_summary.json").write_text(
        json.dumps(asdict(result), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_tsv(
        outdir / "production_stage_status.tsv",
        [asdict(stage) for stage in result.stages],
        ["name", "status", "required", "source", "command", "log", "outputs", "message"],
    )


def run_production(
    manifest_path: str | Path,
    *,
    outdir: str | Path | None = None,
    project_root: str | Path = ".",
    execute: bool = False,
    force: bool = False,
    skip_ranking: bool = False,
    reports_only: bool = False,
) -> ProductionResult:
    manifest = load_production_manifest(manifest_path)
    run_cfg = manifest.get("run") or {}
    root = Path(project_root).resolve()
    sample_id = str(run_cfg.get("sample_id") or "SAMPLE001")
    profile = str(run_cfg.get("profile") or "default")
    run_outdir = Path(outdir or run_cfg.get("outdir") or f"results/{sample_id}_production").resolve()
    run_outdir.mkdir(parents=True, exist_ok=True)
    context = {
        "project_root": str(root),
        "outdir": str(run_outdir),
        "sample_id": sample_id,
        "manifest_dir": str(Path(manifest_path).resolve().parent),
        "OPEN_NEO_REFERENCE_ROOT": os.environ.get("OPEN_NEO_REFERENCE_ROOT", ""),
        "OPEN_NEO_TOOLS_ROOT": os.environ.get("OPEN_NEO_TOOLS_ROOT", ""),
        "OPEN_NEO_LICENSED_ROOT": os.environ.get("OPEN_NEO_LICENSED_ROOT", ""),
    }
    expanded_run = _expand_value(run_cfg, context)
    expanded_evidence = _expand_value(manifest.get("evidence") or {}, context)
    if reports_only:
        from .report_from_final import materialize_hla_loh_layout, write_reports_from_final

        final_outdir = run_outdir / "final"
        materialize_hla_loh_layout(
            final_outdir,
            hla_loh=str(expanded_evidence.get("hla_loh") or ""),
            manifest=manifest,
        )
        outputs = write_reports_from_final(final_outdir)
        stage = StageResult("rebuild_reports", "PASS", True, outputs=outputs)
        result = ProductionResult(
            sample_id=sample_id,
            status="PASS",
            outdir=str(run_outdir),
            dry_run=False,
            stages=[stage],
            generated_config=str(run_outdir / "run.production.generated.toml"),
            final_outdir=str(final_outdir),
        )
        _write_result(result, run_outdir)
        return result
    stage_specs = manifest.get("stages") or {}
    logs_dir = run_outdir / "logs"
    stage_results: list[StageResult] = []
    by_name: dict[str, StageResult] = {}

    for name in _ordered_stages(stage_specs):
        spec = stage_specs[name]
        blocked_dependencies = [
            dep for dep in (spec.get("depends_on") or [])
            if by_name[str(dep)].status in {"FAILED", "BLOCKED"}
        ]
        if blocked_dependencies:
            required = bool(spec.get("required", False))
            status = "BLOCKED" if required else "LOW_CONFIDENCE"
            result = StageResult(
                name,
                status,
                required,
                str(spec.get("source") or ""),
                outputs=_expand_value(spec.get("outputs") or {}, context),
                message="blocked dependencies: " + ", ".join(str(dep) for dep in blocked_dependencies),
            )
        else:
            result = _run_stage(
                name,
                spec,
                context=context,
                logs_dir=logs_dir,
                execute=execute,
                force=force,
            )
        stage_results.append(result)
        by_name[name] = result

    required_failures = [stage for stage in stage_results if stage.required and stage.status in {"FAILED", "BLOCKED"}]
    if required_failures or not execute:
        status = "BLOCKED" if required_failures else "DRY_RUN"
        result = ProductionResult(sample_id, status, str(run_outdir), not execute, stage_results)
        _write_result(result, run_outdir)
        return result

    requirements = expanded_run.get("required_tool_groups") or {}
    if requirements:
        gate_path = run_outdir / "production_release_gate.tsv"
        gate_ok, gate_rows = _tool_group_audit(requirements, stage_specs, stage_results, gate_path)
        gate = StageResult("production_release_gate", "PASS" if gate_ok else "FAILED", True, outputs={"release_gate": str(gate_path)}, message="" if gate_ok else "; ".join(f"{row['domain']}:{row['reason']}" for row in gate_rows if row["status"] == "FAIL"))
        stage_results.append(gate)
        if not gate_ok:
            result = ProductionResult(sample_id, "BLOCKED", str(run_outdir), False, stage_results)
            _write_result(result, run_outdir)
            return result

    normalized_dir = run_outdir / "normalized_sources"
    all_events: list[dict[str, str]] = []
    all_peptides: list[dict[str, str]] = []
    detected_sources: list[str] = []
    for stage in stage_results:
        if not stage.source or stage.status not in {"PASS", "REUSED"}:
            continue
        events, peptides = _candidate_rows(
            stage,
            sample_id=sample_id,
            profile=profile,
            normalized_dir=normalized_dir,
        )
        if events:
            all_events.extend(events)
        if peptides:
            all_peptides.extend(peptides)
        if events or peptides:
            detected_sources.append(stage.source)

    merged_dir = run_outdir / "merged"
    merged_events = merged_dir / "raw_events.tsv"
    merged_peptides = merged_dir / "raw_peptides.tsv"
    merged_event_rows, event_provenance, event_conflicts = merge_rows_preserving_provenance(
        all_events,
        EVENT_FIELDS,
        ("event_id",),
        entity_type="production_event",
    )
    # Biological peptide identity is event + sequence + restricting HLA. Caller-
    # specific peptide IDs are provenance, not a reason to duplicate the entity.
    merged_peptide_rows, peptide_provenance, peptide_conflicts = merge_rows_preserving_provenance(
        all_peptides,
        PEPTIDE_FIELDS,
        ("event_id", "peptide", "hla_allele"),
        entity_type="production_peptide",
    )
    write_tsv(merged_events, merged_event_rows, EVENT_FIELDS)
    write_tsv(merged_peptides, merged_peptide_rows, PEPTIDE_FIELDS)
    event_provenance_path = merged_dir / "event_provenance.tsv"
    peptide_provenance_path = merged_dir / "peptide_provenance.tsv"
    conflicts_path = merged_dir / "evidence_conflicts.tsv"
    write_tsv(event_provenance_path, event_provenance, PROVENANCE_FIELDS)
    write_tsv(peptide_provenance_path, peptide_provenance, PROVENANCE_FIELDS)
    write_tsv(conflicts_path, event_conflicts + peptide_conflicts, CONFLICT_FIELDS)
    secondary_events = str(expanded_evidence.get("secondary_rna_events") or "")
    if secondary_events:
        from .cross_site_rna import build_cross_site_rna_evidence

        cross_site_path = run_outdir / "evidence" / "cross_site_rna_evidence.tsv"
        cross_site_path.parent.mkdir(parents=True, exist_ok=True)
        cross_site_summary = build_cross_site_rna_evidence(
            merged_events,
            secondary_events,
            cross_site_path,
            secondary_sample_id=str(expanded_evidence.get("secondary_sample_id") or "SECONDARY_RNA"),
            identity_status=str(expanded_evidence.get("secondary_identity_status") or "UNASSESSED"),
        )
        expanded_evidence["cross_site_rna_evidence"] = str(cross_site_path)
        cross_site_summary_path = cross_site_path.with_suffix(".summary.json")
        cross_site_summary_path.write_text(json.dumps(cross_site_summary, indent=2) + "\n", encoding="utf-8")
        stage_results.append(StageResult(
            "cross_site_rna_evidence", "PASS", False, source="SECONDARY_SITE_RNA",
            outputs={"cross_site_rna_evidence": str(cross_site_path), "summary": str(cross_site_summary_path)},
            message="secondary-site evidence is annotation-only unless sample identity is confirmed",
        ))
    provenance_outputs = {
        "event_provenance": str(event_provenance_path),
        "peptide_provenance": str(peptide_provenance_path),
        "evidence_conflicts": str(conflicts_path),
    }
    if not merged_peptide_rows:
        stage_results.append(StageResult(
            "candidate_peptide_gate",
            "FAILED",
            True,
            outputs={"raw_events": str(merged_events), "raw_peptides": str(merged_peptides)},
            message=(
                "events were retained for technical review, but no peptide-generating ORF or explicit "
                "junction peptide was available; MHC prediction and patient ranking were not started"
            ),
        ))
        result = ProductionResult(
            sample_id,
            "BLOCKED",
            str(run_outdir),
            False,
            stage_results,
            detected_sources=detected_sources,
            provenance_outputs=provenance_outputs,
        )
        _write_result(result, run_outdir)
        return result

    expected_sources = [str(source) for source in (expanded_run.get("expected_peptide_sources") or [])]
    detected_folded = {source.casefold() for source in detected_sources}
    missing_sources = [source for source in expected_sources if source.casefold() not in detected_folded]
    source_status = "LOW_CONFIDENCE" if missing_sources else "COMPLETE"
    coverage_path = run_outdir / "peptide_source_coverage.tsv"
    write_tsv(
        coverage_path,
        [{
            "status": source_status,
            "expected_sources": ",".join(expected_sources),
            "detected_sources": ",".join(detected_sources),
            "missing_sources": ",".join(missing_sources),
        }],
        ["status", "expected_sources", "detected_sources", "missing_sources"],
    )

    hla_alleles = _read_hla_alleles(expanded_run, stage_results)
    if not hla_alleles:
        result = ProductionResult(
            sample_id,
            "BLOCKED",
            str(run_outdir),
            False,
            stage_results,
            source_status,
            detected_sources,
            missing_sources,
            provenance_outputs=provenance_outputs,
        )
        result.stages.append(StageResult("unified_ranking", "FAILED", True, message="no HLA alleles available"))
        _write_result(result, run_outdir)
        return result

    def normalize_predictor_name(tool: str) -> str:
        return "bigmhc_im" if tool == "bigmhc" else tool

    enabled_predictors = [
        normalize_predictor_name(str(tool))
        for tool in (expanded_run.get("presentation_predictors") or ["netmhcpan", "mhcflurry"])
    ]
    required_predictors = [
        normalize_predictor_name(str(tool))
        for tool in (expanded_run.get("required_presentation_predictors") or enabled_predictors)
    ]
    if "netchop" in enabled_predictors:
        netchop_path = Path(str(expanded_evidence.get("netchop") or run_outdir / "processing/netchop_evidence.tsv"))
        if netchop_path.is_file() and not force:
            netchop_stage = StageResult("netchop_processing", "REUSED", "netchop" in required_predictors, outputs={"netchop": str(netchop_path)})
        else:
            binary = str(expanded_run.get("netchop_executable") or os.environ.get("NEOAG_NETCHOP_BIN") or "netChop")
            home = str(expanded_run.get("netchop_home") or os.environ.get("NETCHOP_HOME") or "")
            command = [sys.executable, str(root / "scripts/run_netchop_evidence.py"), "--raw-peptides", str(merged_peptides), "--output", str(netchop_path), "--binary", binary]
            if home:
                command += ["--home", home]
            proc = subprocess.run(command, cwd=root, text=True, capture_output=True)
            log_path = logs_dir / "netchop_processing.log"
            logs_dir.mkdir(parents=True, exist_ok=True)
            log_path.write_text(proc.stdout + ("\n--- STDERR ---\n" if proc.stderr else "") + proc.stderr, encoding="utf-8")
            ok = proc.returncode == 0 and netchop_path.is_file()
            netchop_stage = StageResult("netchop_processing", "PASS" if ok else "FAILED", "netchop" in required_predictors, command=shlex.join(command), log=str(log_path), outputs={"netchop": str(netchop_path)}, message="" if ok else f"NetChop returned {proc.returncode}")
        stage_results.append(netchop_stage)
        if netchop_stage.required and netchop_stage.status == "FAILED":
            result = ProductionResult(sample_id, "BLOCKED", str(run_outdir), False, stage_results)
            _write_result(result, run_outdir)
            return result
        expanded_evidence["netchop"] = str(netchop_path)
    tools_stub = bool(expanded_run.get("tools_stub", False))
    immunogenicity_stub = bool(expanded_run.get("immunogenicity_stub", False))
    final_outdir = run_outdir / "final"
    try:
        reused_immunogenicity_sources = _materialize_reusable_immunogenicity_outputs(
            final_outdir, expanded_evidence
        )
    except (FileNotFoundError, OSError) as exc:
        stage_results.append(StageResult(
            "immunogenicity_evidence_reuse", "FAILED", True, message=str(exc)
        ))
        result = ProductionResult(sample_id, "BLOCKED", str(run_outdir), False, stage_results)
        _write_result(result, run_outdir)
        return result
    if reused_immunogenicity_sources:
        stage_results.append(StageResult(
            "immunogenicity_evidence_reuse",
            "PASS",
            False,
            outputs={
                tool: str(final_outdir / "presentation" / filename)
                for _, (tool, filename) in _REUSABLE_IMMUNOGENICITY_OUTPUTS.items()
                if tool in reused_immunogenicity_sources
            },
            message="reused normalized evidence: " + ", ".join(reused_immunogenicity_sources),
        ))
    config_path = run_outdir / "run.production.generated.toml"
    _write_final_config(
        config_path,
        sample_id=sample_id,
        profile=profile,
        tools_stub=tools_stub,
        immunogenicity_stub=immunogenicity_stub,
        enabled_predictors=enabled_predictors,
        required_predictors=required_predictors,
        hla_alleles=hla_alleles,
        raw_events=merged_events,
        raw_peptides=merged_peptides,
        evidence=expanded_evidence,
        reuse_immunogenicity_sources=reused_immunogenicity_sources,
    )
    if skip_ranking:
        final_stage = StageResult("unified_ranking", "SKIPPED", True, outputs={"config": str(config_path)})
        final_status = "PARTIAL"
    else:
        py_prefix = Path(sys.executable).resolve().parent.parent
        py_lib = py_prefix / "lib"
        lib_export = (
            f"export LD_LIBRARY_PATH={shlex.quote(str(py_lib))}:"
            '"${LD_LIBRARY_PATH:-}"; '
            if py_lib.is_dir()
            else ""
        )
        command = (
            f"{lib_export}"
            f"source {shlex.quote(str(root / 'conf/tools.env.sh'))}; "
            f"{lib_export}"
            f"{shlex.quote(sys.executable)} -m neoag.cli run-full "
            f"--config {shlex.quote(str(config_path))} --outdir {shlex.quote(str(final_outdir))} "
            f"--reports {shlex.quote(str(expanded_run.get('reports') or 'patient,technical'))}"
        )
        log_path = logs_dir / "unified_ranking.log"
        proc = subprocess.run(["bash", "-lc", command], cwd=root, text=True, capture_output=True)
        logs_dir.mkdir(parents=True, exist_ok=True)
        log_path.write_text(
            proc.stdout + ("\n--- STDERR ---\n" if proc.stderr else "") + proc.stderr,
            encoding="utf-8",
        )
        expected_ranked = final_outdir / "scoring/ranked_peptides.evidence_consensus.tsv"
        final_status = "PASS" if proc.returncode == 0 and expected_ranked.is_file() else "FAILED"
        final_stage = StageResult(
            "unified_ranking",
            final_status,
            True,
            command=command,
            log=str(log_path),
            outputs={"ranked_peptides": str(expected_ranked), "ranked_peptides_weighted_baseline": str(final_outdir / "scoring/ranked_peptides.tsv"), "config": str(config_path)},
            message="" if final_status == "PASS" else f"run-full returned {proc.returncode}",
        )
        if final_status == "PASS":
            mtwt_ok, mtwt_summary = _assess_mtwt_completeness(expected_ranked)
            mtwt_stage = StageResult(
                "mtwt_assessment_gate",
                "PASS" if mtwt_ok else "FAILED",
                True,
                outputs=mtwt_summary,
                message=str(mtwt_summary.get("message") or ""),
            )
            stage_results.append(mtwt_stage)
            if not mtwt_ok:
                final_status = "FAILED"
                final_stage.message = (
                    final_stage.message + "; MT/WT assessment missing for SNV/InDel candidates"
                ).strip("; ")
            else:
                from .report_from_final import materialize_hla_loh_layout, write_reports_from_final

                materialize_hla_loh_layout(
                    final_outdir,
                    hla_loh=str(expanded_evidence.get("hla_loh") or ""),
                    manifest=manifest,
                )
                try:
                    write_reports_from_final(final_outdir)
                    final_stage.message = (final_stage.message + "; reports rebuilt from consensus tables").strip("; ")
                except Exception as exc:  # pragma: no cover - keep production PASS if ranking succeeded
                    final_stage.message = f"{final_stage.message}; report rebuild warning: {exc}".strip("; ")
    stage_results.append(final_stage)
    status = "PASS" if final_status == "PASS" and source_status == "COMPLETE" else (
        "LOW_CONFIDENCE" if final_status == "PASS" else final_status
    )
    result = ProductionResult(
        sample_id=sample_id,
        status=status,
        outdir=str(run_outdir),
        dry_run=False,
        stages=stage_results,
        source_status=source_status,
        detected_sources=detected_sources,
        missing_sources=missing_sources,
        generated_config=str(config_path),
        final_outdir=str(final_outdir),
        provenance_outputs=provenance_outputs,
    )
    _write_result(result, run_outdir)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Execute the manifest-driven NeoAg production workflow")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--outdir")
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--execute", action="store_true", help="Execute stages; default is dry-run")
    parser.add_argument("--force", action="store_true", help="Rerun stages even when declared outputs exist")
    parser.add_argument("--skip-ranking", action="store_true")
    parser.add_argument("--reports-only", action="store_true", help="Rebuild patient/technical reports from existing final/ outputs")
    args = parser.parse_args(argv)
    result = run_production(
        args.manifest,
        outdir=args.outdir,
        project_root=args.project_root,
        execute=args.execute,
        force=args.force,
        skip_ranking=args.skip_ranking,
        reports_only=args.reports_only,
    )
    print(json.dumps(asdict(result), indent=2, ensure_ascii=False))
    return 0 if result.status in {"PASS", "LOW_CONFIDENCE", "PARTIAL", "DRY_RUN"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
