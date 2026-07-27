from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable

from neoag.controlled_execution.io_utils import sha256_file, write_json


ARTIFACT_NAMES: dict[str, tuple[str, ...]] = {
    "weighted_baseline": ("ranked_peptides.weighted_baseline.tsv", "ranked_peptides.tsv", "ranked_peptides.v03.tsv"),
    "consensus_peptides": ("ranked_peptides.evidence_consensus.tsv",),
    "consensus_events": ("ranked_events.evidence_consensus.tsv",),
    "comprehensive_evidence": ("comprehensive_peptide_evidence.tsv", "all_tool_results.tsv"),
    "all_tool_results": ("all_tool_results.tsv",),
    "validation_plan": ("validation_plan.tsv", "validation_plan.v03.tsv"),
    "patient_html": ("evidence_report.patient.html",),
    "technical_html": ("evidence_report.technical.html", "evidence_report.html"),
    "provenance": ("provenance.json", "provenance.v03.json"),
    "run_manifest": ("run_manifest.json", "production_run_summary.json"),
    "comparison_md": ("ranking_compare_weighted_vs_consensus.md",),
    "comparison_tsv": ("ranking_compare_weighted_vs_consensus.tsv", "weighted_vs_consensus_comparison.tsv"),
}


def run_cli(args: list[str], *, cwd: str | Path, log_path: str | Path | None = None, env: dict[str, str] | None = None, timeout: int = 3600) -> dict[str, Any]:
    command = [sys.executable, "-m", "neoag.cli", *args]
    merged_env = os.environ.copy()
    if env:
        merged_env.update(env)
    proc = subprocess.run(command, cwd=str(cwd), text=True, capture_output=True, env=merged_env, timeout=timeout)
    if log_path:
        p = Path(log_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(proc.stdout + ("\n--- STDERR ---\n" if proc.stderr else "") + proc.stderr, encoding="utf-8")
    return {
        "ok": proc.returncode == 0,
        "returncode": proc.returncode,
        "command": command,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "log": str(log_path or ""),
    }


def find_artifact(result_dir: str | Path, names: Iterable[str]) -> Path | None:
    root = Path(result_dir)
    if not root.exists():
        return None
    preferred = [root / "final" / "scoring", root / "scoring", root / "final", root]
    for directory in preferred:
        for name in names:
            candidate = directory / name
            if candidate.is_file():
                return candidate
    matches: list[Path] = []
    for name in names:
        matches.extend(root.rglob(name))
    if not matches:
        return None
    def score(path: Path) -> tuple[int, int, str]:
        text = str(path)
        preferred_score = 0
        if "/final/scoring/" in text:
            preferred_score = 4
        elif "/scoring/" in text:
            preferred_score = 3
        elif "/final/" in text:
            preferred_score = 2
        return (-preferred_score, len(path.parts), text)
    return sorted(set(matches), key=score)[0]


def discover_result_artifacts(result_dir: str | Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for key, names in ARTIFACT_NAMES.items():
        path = find_artifact(result_dir, names)
        if path:
            out[key] = str(path.resolve())
    return out


def ensure_parallel_ranking(
    *,
    project_root: str | Path,
    result_dir: str | Path,
    outdir: str | Path | None = None,
    comprehensive_evidence: str | Path | None = None,
    weighted_baseline: str | Path | None = None,
    rules: str | Path | None = None,
    provenance: str | Path | None = None,
) -> dict[str, Any]:
    artifacts = discover_result_artifacts(outdir if outdir is not None else result_dir)
    if artifacts.get("consensus_peptides") and artifacts.get("consensus_events"):
        return {"status": "REUSED", "artifacts": artifacts}
    source_artifacts = discover_result_artifacts(result_dir)
    comprehensive = Path(comprehensive_evidence or source_artifacts.get("comprehensive_evidence") or "")
    weighted = Path(weighted_baseline or source_artifacts.get("weighted_baseline") or "")
    if not comprehensive.is_file() or not weighted.is_file():
        return {
            "status": "FAILED",
            "failure_code": "CONSENSUS_RANKING_FAILED",
            "message": "Comprehensive evidence and weighted baseline are required",
            "artifacts": artifacts,
        }
    target = Path(outdir or weighted.parent)
    target.mkdir(parents=True, exist_ok=True)
    argv = [
        "evidence-rank",
        "--comprehensive-evidence", str(comprehensive),
        "--weighted-baseline", str(weighted),
        "--outdir", str(target),
    ]
    if rules:
        argv += ["--rules", str(rules)]
    prov = provenance or artifacts.get("provenance")
    if prov:
        argv += ["--provenance", str(prov)]
    result = run_cli(argv, cwd=project_root, log_path=target / "evidence_rank.log")
    return {
        "status": "PASS" if result["ok"] else "FAILED",
        "command": result,
        "artifacts": discover_result_artifacts(target if outdir is not None else result_dir),
    }


def _toml_string(value: str) -> str:
    return json.dumps(str(value), ensure_ascii=False)


def _toml_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    if isinstance(value, list):
        return "[" + ", ".join(_toml_value(v) for v in value) + "]"
    return _toml_string(str(value))


def write_run_config(path: str | Path, inputs: dict[str, Any], routes: list[dict[str, Any]], *, outdir: str | Path, stub: bool = False) -> Path:
    route_names = {str(r.get("route")) for r in routes}
    if route_names == {"peptide_only"}:
        entry_mode = "peptide_only"
    elif route_names <= {"intermediates", "sv_intermediates"}:
        entry_mode = "intermediates" if "intermediates" in route_names else "sv"
    elif route_names == {"splice"}:
        entry_mode = "splice_junction"
    elif route_names == {"fusion"}:
        entry_mode = "fusion"
    elif route_names == {"snv_indel"}:
        entry_mode = "snv_indel"
    else:
        entry_mode = "e2e"

    enabled = inputs.get("enabled_tools") or ["netmhcpan", "mhcflurry"]
    lines = [
        "# Auto-generated by open-neo-run. Review before production execution.",
        "[sample]",
        f"id = {_toml_value(inputs.get('sample_id', 'SAMPLE001'))}",
        f"profile = {_toml_value(inputs.get('profile', 'default'))}",
        f"outdir = {_toml_value(str(Path(outdir).resolve()))}",
        "",
        "[tools]",
        f"stub = {_toml_value(bool(stub))}",
        f"enabled = {_toml_value(enabled)}",
        f"immunogenicity_stub = {_toml_value(bool(inputs.get('immunogenicity_stub', stub)))}",
        "",
        "[inputs]",
        f"entry_mode = {_toml_value(entry_mode)}",
        f"hla_alleles = {_toml_value(inputs.get('hla_alleles') or [])}",
    ]
    mappings = {
        "somatic_vcf": "variants_vcf",
        "fusion_tsv": "fusion_tsv",
        "splice_junction_tsv": "splice_junction_tsv",
        "peptide_csv": "peptide_table",
        "raw_events": "raw_events",
        "raw_peptides": "raw_peptides",
        "sv_raw_events": "sv_raw_events",
        "sv_raw_peptides": "sv_raw_peptides",
        "expression_tsv": "expression",
        "transcript_expression_tsv": "transcript_expression_tsv",
        "rna_evidence_tsv": "rna_vaf_tsv",
        "purity_tsv": "purity",
        "cnv_tsv": "cnv",
        "hla_loh_tsv": "hla_loh",
        "normal_expression": "normal_expression",
        "normal_hla_ligands": "normal_hla_ligands",
        "reference_proteome": "reference_proteome",
        "normal_junctions": "normal_junctions",
        "reference_fasta": "reference_fasta",
        "gencode_gtf": "gencode_gtf",
        "vep_cache": "vep_cache",
    }
    for source, target in mappings.items():
        value = inputs.get(source)
        if value:
            lines.append(f"{target} = {_toml_value(value)}")
    lines += [
        "auto_vep_annotate = true",
        "pass_only = true",
    ]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def write_output_manifest(result_dir: str | Path, path: str | Path) -> Path:
    root = Path(result_dir)
    files = []
    for p in sorted(root.rglob("*")):
        if p.is_file():
            files.append({
                "relative_path": str(p.relative_to(root)),
                "size_bytes": p.stat().st_size,
                "sha256": sha256_file(p) if p.stat().st_size < 50 * 1024 * 1024 else "not_computed_large_file",
            })
    return write_json(path, {"result_dir": str(root.resolve()), "files": files})
