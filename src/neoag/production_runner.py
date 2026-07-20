from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
import tomllib
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from .adapters.pvactools_parser import parse_pvactools_outputs
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


def load_production_manifest(path: str | Path) -> dict[str, Any]:
    return tomllib.loads(Path(path).read_text(encoding="utf-8"))


def _expand(value: str, context: dict[str, str]) -> str:
    expanded = os.path.expanduser(os.path.expandvars(value))
    try:
        return expanded.format_map(context)
    except KeyError as exc:
        raise ValueError(f"Unknown production manifest placeholder: {exc.args[0]}") from exc


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


def _deduplicate(rows: list[dict[str, str]], fields: list[str], key_fields: tuple[str, ...]) -> list[dict[str, str]]:
    seen: set[tuple[str, ...]] = set()
    result: list[dict[str, str]] = []
    for row in rows:
        key = tuple(str(row.get(field) or "") for field in key_fields)
        if key in seen:
            continue
        seen.add(key)
        result.append({field: str(row.get(field) or "") for field in fields})
    return result


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
        return read_tsv(raw_events), read_tsv(raw_peptides)

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
    return events, peptides


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
        "rna_junction_tsv",
        "hla_loh",
        "purity",
        "cnv",
        "normal_expression",
        "normal_hla_ligands",
        "netmhcpan",
        "mhcflurry",
    }
    for key, value in evidence.items():
        if key in allowed_evidence and value and Path(str(value)).exists():
            lines.append(f"{key} = {toml_value(value)}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


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
    }
    expanded_run = _expand_value(run_cfg, context)
    expanded_evidence = _expand_value(manifest.get("evidence") or {}, context)
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
        if peptides:
            all_events.extend(events)
            all_peptides.extend(peptides)
            detected_sources.append(stage.source)

    merged_dir = run_outdir / "merged"
    merged_events = merged_dir / "raw_events.tsv"
    merged_peptides = merged_dir / "raw_peptides.tsv"
    write_tsv(merged_events, _deduplicate(all_events, EVENT_FIELDS, ("event_id",)), EVENT_FIELDS)
    write_tsv(
        merged_peptides,
        _deduplicate(all_peptides, PEPTIDE_FIELDS, ("peptide_id", "event_id", "peptide", "hla_allele")),
        PEPTIDE_FIELDS,
    )

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
        )
        result.stages.append(StageResult("unified_ranking", "FAILED", True, message="no HLA alleles available"))
        _write_result(result, run_outdir)
        return result

    enabled_predictors = [str(tool) for tool in (expanded_run.get("presentation_predictors") or ["netmhcpan", "mhcflurry"])]
    required_predictors = [str(tool) for tool in (expanded_run.get("required_presentation_predictors") or enabled_predictors)]
    tools_stub = bool(expanded_run.get("tools_stub", False))
    immunogenicity_stub = bool(expanded_run.get("immunogenicity_stub", False))
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
    )
    final_outdir = run_outdir / "final"
    if skip_ranking:
        final_stage = StageResult("unified_ranking", "SKIPPED", True, outputs={"config": str(config_path)})
        final_status = "PARTIAL"
    else:
        command = (
            f"source {shlex.quote(str(root / 'conf/tools.env.sh'))}; "
            f"{shlex.quote(sys.executable)} -m neoag.cli run-full "
            f"--config {shlex.quote(str(config_path))} --outdir {shlex.quote(str(final_outdir))}"
        )
        log_path = logs_dir / "unified_ranking.log"
        proc = subprocess.run(["bash", "-lc", command], cwd=root, text=True, capture_output=True)
        logs_dir.mkdir(parents=True, exist_ok=True)
        log_path.write_text(
            proc.stdout + ("\n--- STDERR ---\n" if proc.stderr else "") + proc.stderr,
            encoding="utf-8",
        )
        expected_ranked = final_outdir / "scoring/ranked_peptides.tsv"
        final_status = "PASS" if proc.returncode == 0 and expected_ranked.is_file() else "FAILED"
        final_stage = StageResult(
            "unified_ranking",
            final_status,
            True,
            command=command,
            log=str(log_path),
            outputs={"ranked_peptides": str(expected_ranked), "config": str(config_path)},
            message="" if final_status == "PASS" else f"run-full returned {proc.returncode}",
        )
    stage_results.append(final_stage)
    status = "PASS" if final_status == "PASS" and source_status == "COMPLETE" else (
        "LOW_CONFIDENCE" if final_status == "PASS" else final_status
    )
    result = ProductionResult(
        sample_id,
        status,
        str(run_outdir),
        False,
        stage_results,
        source_status,
        detected_sources,
        missing_sources,
        str(config_path),
        str(final_outdir),
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
    args = parser.parse_args(argv)
    result = run_production(
        args.manifest,
        outdir=args.outdir,
        project_root=args.project_root,
        execute=args.execute,
        force=args.force,
        skip_ranking=args.skip_ranking,
    )
    print(json.dumps(asdict(result), indent=2, ensure_ascii=False))
    return 0 if result.status in {"PASS", "LOW_CONFIDENCE", "PARTIAL", "DRY_RUN"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
