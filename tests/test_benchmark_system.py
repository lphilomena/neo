from __future__ import annotations

from neoag.benchmark_system import run_system_benchmark
from neoag.cli import build_parser
from neoag.utils import read_tsv, write_tsv


def test_system_benchmark_writes_synthetic_sensitivity_and_pending_ligandome(tmp_path):
    out = run_system_benchmark(outdir=tmp_path, mode="all")
    assert out["synthetic_perturbation"].is_file()
    assert out["sensitivity_analysis"].is_file()
    assert out["ligandome_ms_validation"].is_file()
    synthetic = read_tsv(out["synthetic_perturbation"])
    by_scenario = {r["scenario"]: r for r in synthetic}
    assert by_scenario["b2m_biallelic_loss"]["final_priority"] == "D"
    assert by_scenario["appm_unassessed"]["observed_behavior"] == "score_preserved"
    assert all(r["claim_scope"] == "synthetic_system_behavior_check" for r in synthetic)
    sensitivity = read_tsv(out["sensitivity_analysis"])
    assert sensitivity
    assert {"appm_multiplier", "escape_multiplier", "ccf_multiplier", "priority_cap"} <= {r["parameter"] for r in sensitivity}
    ligandome = read_tsv(out["ligandome_ms_validation"])
    assert ligandome[0]["benchmark_status"] == "external_required"
    assert ligandome[0]["claim_scope"] == "pending_external_ligandome_ms_validation"


def test_ligandome_ms_validation_computes_metrics_when_labels_are_available(tmp_path):
    ligandome = tmp_path / "ligandome.tsv"
    write_tsv(ligandome, [
        {"peptide": "AAAA", "hla_allele": "HLA-A*02:01", "presented": "1", "presentation_evidence_score": "0.95"},
        {"peptide": "BBBB", "hla_allele": "HLA-A*02:01", "presented": "1", "presentation_evidence_score": "0.80"},
        {"peptide": "CCCC", "hla_allele": "HLA-A*02:01", "presented": "0", "presentation_evidence_score": "0.20"},
        {"peptide": "DDDD", "hla_allele": "HLA-A*02:01", "presented": "0", "presentation_evidence_score": "0.10"},
    ])
    out = run_system_benchmark(outdir=tmp_path / "out", mode="ligandome-ms", ligandome_ms=ligandome)
    rows = read_tsv(out["ligandome_ms_validation"])
    assert rows[0]["benchmark_status"] == "completed"
    assert rows[0]["n"] == "4"
    assert rows[0]["n_presented"] == "2"
    assert rows[0]["auroc"] == "1.0000"


def test_benchmark_system_cli_is_registered():
    parser = build_parser()
    args = parser.parse_args(["benchmark-system", "--outdir", "work/bench", "--mode", "synthetic"])
    assert args.mode == "synthetic"
    assert callable(args.func)
