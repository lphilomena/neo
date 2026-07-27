from __future__ import annotations

import csv
import json
from pathlib import Path

from neoag.open_neo.install_check import run_install_check
from neoag.open_neo.review import run_review
from neoag.open_neo.routing import inspect_manifest
from neoag.open_neo.run import run_open_neo
from neoag.skill_taxonomy.registry import SKILLS_BY_NAME
from neoag.skill_taxonomy.runner import run_skill


def test_public_macro_registry_and_internal_composition():
    for name in ("open-neo-install-check", "open-neo-run", "open-neo-review"):
        spec = SKILLS_BY_NAME[name]
        assert spec.category == "M"
        assert spec.visibility == "public"
        assert spec.entrypoint.startswith("open-neo ")
        assert spec.composes
    assert SKILLS_BY_NAME["neoag-vcf"].category == "A"
    assert SKILLS_BY_NAME["neoag-ranking"].category == "B"
    assert SKILLS_BY_NAME["neoag-experiment-design"].category == "C"
    assert SKILLS_BY_NAME["neoag-doctor"].category == "D"


def test_input_detection_multi_route(tmp_path: Path):
    project = Path.cwd()
    manifest = tmp_path / "sample.yaml"
    manifest.write_text(
        f"""schema_version: open-neo-sample-v1
case_id: CASE_ROUTE
sample_id: SAMPLE_ROUTE
genome_build: GRCh38
inputs:
  fusion_tsv: {project / 'data/fixtures/easyfuse_fusions.pass.tsv'}
  splice_junction_tsv: {project / 'data/fixtures/regtools_splice_junctions.tsv'}
  hla_alleles:
    - HLA-A*02:01
    - HLA-B*07:02
execution:
  profile: default
""",
        encoding="utf-8",
    )
    result = inspect_manifest(manifest)
    assert result.status == "PASS"
    assert {route.route for route in result.routes} == {"fusion", "splice"}
    assert not result.missing


def test_install_check_review_tier_source_checkout(tmp_path: Path):
    result = run_install_check({
        "project_root": str(Path.cwd()),
        "deployment_tier": "review",
        "mode": "plan",
        "release_audit": False,
        "outdir": str(tmp_path / "install"),
    })
    assert result["status"] in {"READY", "PARTIAL"}
    assert Path(result["outputs"]["deployment_report"]).is_file()


def test_install_check_release_requires_checksum(tmp_path: Path):
    tarball = tmp_path / "release.tar"
    tarball.write_bytes(b"release")
    result = run_install_check({
        "project_root": str(Path.cwd()),
        "release_tarball": str(tarball),
        "mode": "plan",
        "outdir": str(tmp_path / "install"),
    })
    assert result["status"] == "BLOCKED"
    assert "CHECKSUM_REQUIRED" in result["blocking_issues"]


def test_open_neo_run_plan_writes_route_and_config(tmp_path: Path):
    result = run_open_neo({
        "mode": "plan",
        "project_root": str(Path.cwd()),
        "sample_id": "PLAN",
        "fusion_tsv": str(Path.cwd() / "data/fixtures/easyfuse_fusions.pass.tsv"),
        "splice_junction_tsv": str(Path.cwd() / "data/fixtures/regtools_splice_junctions.tsv"),
        "hla_alleles": ["HLA-A*02:01", "HLA-B*07:02"],
        "doctor": False,
        "outdir": str(tmp_path / "plan"),
    })
    assert result["status"] == "PASS"
    assert Path(result["outputs"]["route_plan"]).is_file()
    config = Path(result["outputs"]["generated_run_config"])
    assert config.is_file()
    assert 'entry_mode = "e2e"' in config.read_text(encoding="utf-8")


def test_open_neo_run_missing_hla_is_blocked(tmp_path: Path):
    result = run_open_neo({
        "mode": "plan",
        "project_root": str(Path.cwd()),
        "fusion_tsv": str(Path.cwd() / "data/fixtures/easyfuse_fusions.pass.tsv"),
        "doctor": False,
        "outdir": str(tmp_path / "plan"),
    })
    assert result["status"] == "BLOCKED"
    assert "HLA_MISSING" in result["blocking_issues"]


def _write_minimal_evidence(tmp_path: Path) -> tuple[Path, Path]:
    comprehensive = tmp_path / "comprehensive_peptide_evidence.tsv"
    comprehensive.write_text(
        "peptide_id\tevent_id\tevent_type\tpeptide\thla_allele\t"
        "cross_platform_status\trna_alt_reads\trna_vaf\t"
        "netmhcpan_mt_rank_el\tmhcflurry_presentation_score\t"
        "mutant_specificity_status\tccf_estimate\tccf_confidence\t"
        "appm_multiplier\thla_loh_status\trestricting_hla_lost\t"
        "safety_status\tsafety_evidence_completeness\n"
        "P1\tE1\tSNV\tSYFPEITHI\tHLA-A*02:01\t"
        "CROSS_PLATFORM_PASS_CONCORDANT\t8\t0.20\t0.2\t0.8\t"
        "MT_SPECIFIC\t0.9\thigh\t1.0\tRETAINED\tfalse\tPASS\tCOMPLETE\n",
        encoding="utf-8",
    )
    weighted = tmp_path / "ranked_peptides.tsv"
    weighted.write_text("peptide_id\tefficacy_score\tfinal_priority\nP1\t0.8\tB\n", encoding="utf-8")
    return comprehensive, weighted


def test_open_neo_run_ranking_only_keeps_baseline_and_emits_consensus(tmp_path: Path):
    comprehensive, weighted = _write_minimal_evidence(tmp_path)
    outdir = tmp_path / "ranking"
    result = run_open_neo({
        "mode": "ranking-only",
        "project_root": str(Path.cwd()),
        "comprehensive_evidence": str(comprehensive),
        "weighted_baseline": str(weighted),
        "doctor": False,
        "outdir": str(outdir),
    })
    assert result["status"] == "PASS"
    assert result["production_command"] == "neoag evidence-rank"
    assert weighted.read_text(encoding="utf-8").startswith("peptide_id")
    for name in (
        "ranked_peptides.weighted_baseline.tsv",
        "ranked_peptides.evidence_consensus.tsv",
        "ranked_events.evidence_consensus.tsv",
        "all_tool_results.tsv",
    ):
        assert (outdir / name).is_file()


def _write_review_fixture(root: Path) -> None:
    scoring = root / "scoring"
    scoring.mkdir(parents=True)
    (scoring / "ranked_peptides.weighted_baseline.tsv").write_text(
        "peptide_id\tefficacy_score\tfinal_priority\nP1\t0.8\tB\nP2\t0.7\tC\n", encoding="utf-8"
    )
    (scoring / "ranked_peptides.evidence_consensus.tsv").write_text(
        "evidence_rank\tpeptide_id\tevent_id\tevent_type\tgene\tpeptide\thla_allele\tevidence_grade\tpareto_front\n"
        "1\tP1\tE1\tSNV\tGENE1\tSYFPEITHI\tHLA-A*02:01\tR1\t1\n"
        "2\tP2\tE2\tfusion\tGENE2::GENE3\tABCDEFGHI\tHLA-B*07:02\tR3\t1\n",
        encoding="utf-8",
    )
    (scoring / "ranked_events.evidence_consensus.tsv").write_text(
        "event_evidence_rank\tevent_group_id\tevent_id\tgene\tevent_type\tbest_evidence_grade\tbest_pareto_front\tmanual_review_required\trecommended_next_steps\tevent_consensus_trace\t"
        "representative_1_peptide_id\trepresentative_1_peptide\trepresentative_1_hla_allele\n"
        "1\tEVENT:E1\tE1\tGENE1\tSNV\tR1\t1\tno\tMT/WT peptide assay\tcomplete\tP1\tSYFPEITHI\tHLA-A*02:01\n"
        "2\tEVENT:E2\tE2\tGENE2::GENE3\tFusion\tR3\t1\tyes\ttargeted RNA and second caller\tRNA-only\tP2\tABCDEFGHI\tHLA-B*07:02\n",
        encoding="utf-8",
    )


def test_open_neo_review_is_event_level_and_non_mutating(tmp_path: Path):
    result_dir = tmp_path / "result"
    _write_review_fixture(result_dir)
    weighted = result_dir / "scoring/ranked_peptides.weighted_baseline.tsv"
    before = weighted.read_bytes()
    outdir = tmp_path / "review"
    result = run_review({"result_dir": str(result_dir), "top_n": 2, "outdir": str(outdir)})
    assert result["status"] == "PASS_WITH_WARNINGS"
    assert weighted.read_bytes() == before
    review_rows = list(csv.DictReader((outdir / "review/candidate_review.tsv").open(), delimiter="\t"))
    assert [row["gene"] for row in review_rows] == ["GENE1", "GENE2::GENE3"]
    assert (outdir / "review/first_batch_experiment_set.tsv").is_file()
    first_batch = list(csv.DictReader((outdir / "review/first_batch_experiment_set.tsv").open(), delimiter="\t"))
    assert [row["gene"] for row in first_batch] == ["GENE1"]
    completion = list(csv.DictReader((outdir / "review/evidence_completion_queue.tsv").open(), delimiter="\t"))
    assert [row["gene"] for row in completion] == ["GENE2::GENE3"]
    assert (outdir / "reports/patient_report.md").is_file()


def test_open_neo_review_allows_missing_weighted_baseline(tmp_path: Path):
    result_dir = tmp_path / "result"
    _write_review_fixture(result_dir)
    (result_dir / "scoring/ranked_peptides.weighted_baseline.tsv").unlink()
    result = run_review({"result_dir": str(result_dir), "outdir": str(tmp_path / "review")})
    assert result["status"] == "PASS_WITH_WARNINGS"


def test_open_neo_review_blocks_source_output_mutation(tmp_path: Path):
    result_dir = tmp_path / "result"
    _write_review_fixture(result_dir)
    result = run_review({"result_dir": str(result_dir), "outdir": str(result_dir)})
    assert result["status"] == "BLOCKED"
    assert "REPORT_BOUNDARY_VIOLATION" in result["blocking_issues"]


def test_macro_skills_run_through_skill_runner(tmp_path: Path):
    comprehensive, weighted = _write_minimal_evidence(tmp_path)
    result = run_skill("open-neo-run", {
        "comprehensive_evidence": str(comprehensive),
        "weighted_baseline": str(weighted),
        "doctor": False,
        "outdir": str(tmp_path / "skill"),
    })
    assert result["status"] == "PASS"
    assert result["algorithm_owner"] == "src/neoag/evidence_consensus.py"
