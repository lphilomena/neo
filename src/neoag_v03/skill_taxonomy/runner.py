from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

from .entry_skills import run_fusion, run_peptide_csv, run_splice, run_sv_wes, run_sv_wgs, run_vcf
from .evidence_skills import run_appm_escape, run_ccf, run_expression, run_hla_typing_loh, run_presentation, run_ranking, run_rna_evidence, run_safety
from .governance_skills import run_demo_smoke, run_doctor, run_gateway_submit, run_hpc_runner, run_input_qc, run_pipeline_full, run_release_qc, run_remote_deploy, run_tool_reference_qc
from .registry import SKILLS_BY_NAME, SkillSpec
from .review_skills import run_concept_explainer, run_experiment_design, run_patient_report, run_ranking_compare, run_technical_report
from .io import ensure_dir, write_json

HANDLERS: dict[str, Callable[[dict[str, Any]], dict[str, Any]]] = {
    "vcf": run_vcf,
    "fusion": run_fusion,
    "splice": run_splice,
    "sv_wgs": run_sv_wgs,
    "sv_wes": run_sv_wes,
    "peptide_csv": run_peptide_csv,
    "hla_typing_loh": run_hla_typing_loh,
    "presentation": run_presentation,
    "expression": run_expression,
    "rna_evidence": run_rna_evidence,
    "ccf": run_ccf,
    "appm_escape": run_appm_escape,
    "safety": run_safety,
    "ranking": run_ranking,
    "ranking_compare": run_ranking_compare,
    "experiment_design": run_experiment_design,
    "patient_report": run_patient_report,
    "technical_report": run_technical_report,
    "concept_explainer": run_concept_explainer,
    "input_qc": run_input_qc,
    "doctor": run_doctor,
    "tool_reference_qc": run_tool_reference_qc,
    "run_demo_smoke": run_demo_smoke,
    "pipeline_full": run_pipeline_full,
    "release_qc": run_release_qc,
    "remote_deploy": run_remote_deploy,
    "gateway_submit": run_gateway_submit,
    "hpc_runner": run_hpc_runner,
}


def get_skill(name: str) -> SkillSpec:
    if name not in SKILLS_BY_NAME:
        known = ", ".join(sorted(SKILLS_BY_NAME))
        raise KeyError(f"Unknown skill {name!r}. Known skills: {known}")
    return SKILLS_BY_NAME[name]


def run_skill(name: str, args: dict[str, Any], dry_run: bool = False) -> dict[str, Any]:
    spec = get_skill(name)
    outdir = ensure_dir(args.get("outdir") or f"work/{name}")
    args = dict(args)
    args["outdir"] = str(outdir)
    if dry_run:
        result = {"status": "DRY_RUN", "skill": name, "category": spec.category, "handler": spec.handler, "required_inputs": spec.required_inputs, "outputs": spec.outputs, "risk_level": spec.risk_level, "approval_required": spec.approval_required, "arguments": args}
        write_json(outdir / "skill_result.json", result)
        return result
    handler = HANDLERS.get(spec.handler)
    if handler is None:
        result = {"status": "FAIL", "skill": name, "failure_reason": f"No handler registered for {spec.handler}"}
        write_json(outdir / "skill_result.json", result)
        return result
    return handler(args)


def validate_skill_dirs(root: str | Path = ".") -> dict[str, Any]:
    base = Path(root) / ".agents" / "skills"
    rows = []
    missing = []
    for name, spec in sorted(SKILLS_BY_NAME.items()):
        skill_dir = base / name
        exists = skill_dir.exists()
        has_md = (skill_dir / "SKILL.md").exists()
        rows.append({"skill": name, "category": spec.category, "dir_exists": exists, "skill_md_exists": has_md})
        if not exists or not has_md:
            missing.append(name)
    return {"status": "PASS" if not missing else "INCOMPLETE", "missing": missing, "rows": rows}
