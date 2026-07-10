from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from typing import Any

from .context_builder import build_file_records
from .schemas import FileRecord, InputState
from neoag_v03.agents.skill_router import find_named_files

FEATURE_KEYS = [
    "has_ranked_peptides_recommendation",
    "has_ranked_peptides_netmhcpan42",
    "has_ranked_peptides",
    "has_evidence_report",
    "has_hla_loh",
    "has_purity",
    "has_appm",
    "has_ccf",
    "has_expression",
    "has_hla",
    "has_somatic_vcf",
    "has_sv_vcf",
]


def lightweight_features(named_files: dict[str, str], records: list[FileRecord]) -> dict[str, bool]:
    kinds = {r.kind for r in records}
    return {
        "has_ranked_peptides_recommendation": bool(named_files.get("recommendation")) or "ranked_peptides_recommendation" in kinds,
        "has_ranked_peptides_netmhcpan42": bool(named_files.get("netmhcpan42")) or "ranked_peptides_netmhcpan42" in kinds,
        "has_ranked_peptides": bool(named_files.get("ranked_peptides") or named_files.get("recommendation") or named_files.get("netmhcpan42")),
        "has_evidence_report": bool(named_files.get("evidence_report")) or "evidence_report" in kinds,
        "has_hla_loh": bool(named_files.get("hla_loh")) or "hla_loh" in kinds,
        "has_purity": bool(named_files.get("purity_table")) or "purity_or_cnv" in kinds,
        "has_appm": any(k.startswith("appm") for k in kinds) or bool(named_files.get("appm_gene_status") or named_files.get("appm_submodule_scores")),
        "has_ccf": "ccf" in kinds,
        "has_expression": any("expression" in r.name.lower() or "tpm" in r.name.lower() for r in records),
        "has_hla": any("hla" in r.name.lower() for r in records),
        "has_somatic_vcf": any(r.kind == "vcf" for r in records),
        "has_sv_vcf": any("sv" in r.name.lower() and r.kind == "vcf" for r in records),
    }


def build_input_state(files: list[str] | None = None, result_dir: str | None = None, outdir: str | Path | None = None, execute_input_qc: bool = False) -> InputState:
    records = build_file_records(files, result_dir)
    named = find_named_files(result_dir, files)
    state = InputState(files=records, known_files=named, features=lightweight_features(named, records))
    if execute_input_qc and outdir:
        qc_out = Path(outdir) / "neoag-input-qc"
        cmd = [sys.executable, "-m", "neoag_v03.agent_skills.input_qc", "--outdir", str(qc_out)]
        if result_dir:
            cmd += ["--result-dir", result_dir]
        proc = subprocess.run(cmd, text=True, capture_output=True)
        status_path = qc_out / "input_status.json"
        if status_path.exists():
            data = json.loads(status_path.read_text(encoding="utf-8"))
            if isinstance(data.get("features"), dict):
                state.features.update({str(k): bool(v) for k, v in data["features"].items()})
            if isinstance(data.get("missing_inputs"), list):
                state.missing_inputs = data["missing_inputs"]
            state.recommended_workflow = data.get("recommended_workflow")
            state.input_qc_path = str(status_path)
        else:
            state.missing_inputs.append({"input": "input_qc", "severity": "warning", "reason": f"input-qc failed returncode={proc.returncode}: {proc.stderr[-300:]}"})
    else:
        if state.features.get("has_ranked_peptides_recommendation") and state.features.get("has_ranked_peptides_netmhcpan42"):
            state.recommended_workflow = "ranking_compare"
        elif state.features.get("has_ranked_peptides"):
            state.recommended_workflow = "result_review"
        elif state.features.get("has_sv_vcf"):
            state.recommended_workflow = "dna_sv_workflow"
        elif state.features.get("has_somatic_vcf") and state.features.get("has_hla"):
            state.recommended_workflow = "snv_indel_workflow"
        else:
            state.recommended_workflow = "input_incomplete"
    return state
