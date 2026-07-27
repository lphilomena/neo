from __future__ import annotations

import csv
import json
from pathlib import Path

from neoag.open_neo.install_check import run_install_check
from neoag.open_neo.review import build_review_rows, run_review, select_first_batch
from neoag.open_neo.routing import inspect_manifest
from neoag.open_neo.run import run_open_neo
from neoag.open_neo.rna_preprocessing import prepare_rna_evidence
from neoag.open_neo.tool_consensus import build_tool_consensus
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
    assert Path(result["outputs"]["pipeline_plan"]).is_file()
    assert Path(result["outputs"]["consensus_tool_run_status"]).is_file()
    assert Path(result["outputs"]["consensus_tool_consensus_summary"]).is_file()


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


def test_input_directory_header_detection_and_bam_index_qc(tmp_path: Path):
    data_dir = tmp_path / "inputs"
    data_dir.mkdir()
    fusion = data_dir / "caller.tsv"
    fusion.write_text("gene1\tgene2\tjunction_reads\nA\tB\t5\n", encoding="utf-8")
    hla = data_dir / "typing.txt"
    hla.write_text("HLA-A*02:01 HLA-B*07:02\n", encoding="utf-8")
    manifest = tmp_path / "sample.yaml"
    manifest.write_text(f"case_id: AUTO\nsample_id: AUTO\ninputs:\n  input_dir: {data_dir}\n", encoding="utf-8")
    result = inspect_manifest(manifest, input_dir=data_dir, output_dir=tmp_path / "out")
    assert result.status == "PASS"
    assert result.inputs["fusion_tsv"] == str(fusion)
    assert result.inputs["hla_file"] == str(hla)


def test_bam_without_index_blocks_production_route(tmp_path: Path):
    bam = tmp_path / "tumor.bam"
    bam.write_bytes(b"BAM")
    manifest = tmp_path / "sample.yaml"
    manifest.write_text(f"case_id: BAM\nsample_id: BAM\ninputs:\n  tumor_dna_bam: {bam}\n", encoding="utf-8")
    result = inspect_manifest(manifest, output_dir=tmp_path / "out")
    assert result.status == "BLOCKED"
    assert any(item["field"] == "tumor_dna_bam_index" for item in result.missing)


def test_tool_consensus_emits_domain_outputs(tmp_path: Path):
    optitype = tmp_path / "sample_optitype_result.txt"
    spechla = tmp_path / "sample_spechla_result.txt"
    optitype.write_text("HLA-A*02:01 HLA-B*07:02 HLA-C*07:02\n", encoding="utf-8")
    spechla.write_text("HLA-A*02:01 HLA-B*07:02 HLA-C*07:02\n", encoding="utf-8")
    outputs = build_tool_consensus({
        "sample_id": "S1",
        "tool_results": {"hla_typing": {"optitype": str(optitype), "spechla": str(spechla)}},
    }, tmp_path / "consensus")
    assert Path(outputs["hla_typing_consensus.tsv"]).is_file()
    text = Path(outputs["tool_consensus_summary.tsv"]).read_text(encoding="utf-8")
    assert "hla_typing\tCONSISTENT" in text
    assert Path(outputs["tool_evidence.long.tsv"]).is_file()


def test_rna_preprocessing_plans_gene_transcript_tpm_and_alt_vaf(tmp_path: Path):
    fastq1 = tmp_path / "tumor_R1.fastq.gz"
    fastq2 = tmp_path / "tumor_R2.fastq.gz"
    bam = tmp_path / "tumor.rna.bam"
    vcf = tmp_path / "somatic.vcf"
    tx2gene = tmp_path / "tx2gene.tsv"
    salmon_index = tmp_path / "salmon_index"
    salmon_index.mkdir()
    for path in (fastq1, fastq2, bam):
        path.write_bytes(b"fixture")
    tx2gene.write_text("transcript_id\tgene_id\nTX1\tG1\n", encoding="utf-8")
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t1\t.\tA\tT\t.\tPASS\t.\n", encoding="utf-8")
    result = prepare_rna_evidence(
        {
            "sample_id": "RNA1",
            "tumor_rna_fastq": [str(fastq1), str(fastq2)],
            "tumor_rna_bam": str(bam),
            "somatic_vcf": str(vcf),
            "salmon_index": str(salmon_index),
            "tx2gene": str(tx2gene),
        },
        project_root=Path.cwd(),
        outdir=tmp_path / "rna",
        execute=False,
    )
    stages = {row["stage"]: row for row in result["stages"]}
    assert stages["rna_quantification"]["status"] == "PLANNED"
    assert "run_salmon_fastq_to_tpm.sh" in stages["rna_quantification"]["command_preview"]
    assert stages["rna_alt_vaf"]["status"] == "PLANNED"
    assert "rna_allele_counts_pysam.py" in stages["rna_alt_vaf"]["command_preview"]
    assert Path(tmp_path / "rna/rna_preprocessing_status.tsv").is_file()


def test_direct_execute_requires_gateway(tmp_path: Path):
    comprehensive, weighted = _write_minimal_evidence(tmp_path)
    result = run_open_neo({
        "mode": "execute", "approved": True, "doctor": False,
        "comprehensive_evidence": str(comprehensive), "weighted_baseline": str(weighted),
        "outdir": str(tmp_path / "execute"),
    })
    assert result["status"] == "BLOCKED"
    assert "GATEWAY_REQUIRED" in result["blocking_issues"]


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
    all_tool_header = (outdir / "all_tool_results.tsv").read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "legacy_weighted_rank" in all_tool_header
    assert "evidence_rank" in all_tool_header
    assert (outdir / "pipeline/tool_consensus/presentation_consensus.tsv").is_file()


def _write_review_fixture(root: Path) -> None:
    scoring = root / "scoring"
    scoring.mkdir(parents=True)
    (scoring / "ranked_peptides.weighted_baseline.tsv").write_text(
        "peptide_id\tefficacy_score\tfinal_priority\nP1\t0.8\tB\nP2\t0.7\tC\n", encoding="utf-8"
    )
    (scoring / "ranked_peptides.evidence_consensus.tsv").write_text(
        "evidence_rank\tpeptide_id\tevent_id\tevent_type\tgene\tpeptide\thla_allele\tevidence_grade\tpareto_front\trna_support_state\tsafety_state\tpresentation_consensus_state\tmutant_specificity_state\tclonality_state\tccf_confidence\thla_appm_state\tevidence_completeness_state\thard_failure\n"
        "1\tP1\tE1\tSNV\tGENE1\tSYFPEITHI\tHLA-A*02:01\tR1\t1\tRNA_CONFIRMED\tSAFETY_PASS\tPRESENTATION_CONSISTENT_STRONG\tMT_SPECIFIC\tCLONAL_LIKE\thigh\tHLA_APPM_RETAINED\tCOMPLETE\tno\n"
        "2\tP2\tE2\tfusion\tGENE2::GENE3\tABCDEFGHI\tHLA-B*07:02\tR3\t1\tRNA_UNASSESSED\tSAFETY_PARTIAL\tPRESENTATION_SINGLE_TOOL\tNOVEL_JUNCTION\tunresolved\tlow\tHLA_LOH_UNASSESSED\tPARTIAL\tno\n",
        encoding="utf-8",
    )
    (scoring / "all_tool_results.tsv").write_text((scoring / "ranked_peptides.evidence_consensus.tsv").read_text(encoding="utf-8"), encoding="utf-8")
    (scoring / "validation_plan.tsv").write_text("event_id\trecommended_validation\nE1\tMT/WT pair\nE2\ttargeted RNA\n", encoding="utf-8")
    (scoring / "run_manifest.json").write_text(json.dumps({"schema_version": "1.0", "algorithm": "discrete_state_grade_track_pareto_v2", "status": "PROVISIONAL_RESEARCH_ONLY", "outputs": {}}) + "\n", encoding="utf-8")
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
    source_hashes = {path: path.read_bytes() for path in (result_dir / "scoring").iterdir() if path.is_file()}
    outdir = tmp_path / "review"
    result = run_review({"result_dir": str(result_dir), "top_n": 2, "outdir": str(outdir)})
    assert result["status"] == "PASS_WITH_WARNINGS"
    assert weighted.read_bytes() == source_hashes[weighted]
    assert all(path.read_bytes() == content for path, content in source_hashes.items())
    review_rows = list(csv.DictReader((outdir / "review/candidate_review.tsv").open(), delimiter="\t"))
    assert [row["gene"] for row in review_rows] == ["GENE1", "GENE2::GENE3"]
    assert (outdir / "review/first_batch_experiment_set.tsv").is_file()
    first_batch = list(csv.DictReader((outdir / "review/first_batch_experiment_set.tsv").open(), delimiter="\t"))
    assert [row["gene"] for row in first_batch] == ["GENE1", "GENE2::GENE3"]
    assert first_batch[1]["experiment_priority"] == "TARGETED_RNA_FIRST"
    completion = list(csv.DictReader((outdir / "review/evidence_completion_queue.tsv").open(), delimiter="\t"))
    assert [row["gene"] for row in completion] == ["GENE2::GENE3"]
    assert (outdir / "reports/patient_report.md").is_file()
    assert (outdir / "review/integrity/review_integrity.json").is_file()
    assert (outdir / "review/experiment_design/short_peptide_pool.tsv").is_file()
    assert (outdir / "review/experiment_design/targeted_rna_validation_plan.tsv").is_file()
    assert (outdir / "review/hla_loh_appm_review/appm_escape_review.md").is_file()
    assert (outdir / "review/ccf_clonality_review/ccf_clonality_review.md").is_file()
    assert (outdir / "reports/technical_report.md").is_file()
    assert review_rows[0]["pipeline_r_grade"] == "R1"
    assert review_rows[0]["experiment_priority"] == "EXPERIMENT_PRIORITY_HIGH"


def test_open_neo_review_blocks_incomplete_required_result_set(tmp_path: Path):
    result_dir = tmp_path / "result"
    _write_review_fixture(result_dir)
    (result_dir / "scoring/ranked_peptides.weighted_baseline.tsv").unlink()
    result = run_review({"result_dir": str(result_dir), "outdir": str(tmp_path / "review")})
    assert result["status"] == "BLOCKED"
    assert "REVIEW_INTEGRITY_BLOCKED" in result["blocking_issues"]


def test_open_neo_review_returns_needs_ranking_without_event_consensus(tmp_path: Path):
    result_dir = tmp_path / "result"
    _write_review_fixture(result_dir)
    (result_dir / "scoring/ranked_events.evidence_consensus.tsv").unlink()
    result = run_review({"result_dir": str(result_dir), "outdir": str(tmp_path / "review")})
    assert result["status"] == "NEEDS_RANKING"
    assert "NEEDS_RANKING" in result["blocking_issues"]


def test_open_neo_review_blocks_promoted_hard_failure(tmp_path: Path):
    result_dir = tmp_path / "result"
    _write_review_fixture(result_dir)
    peptide = result_dir / "scoring/ranked_peptides.evidence_consensus.tsv"
    text = peptide.read_text(encoding="utf-8").replace("\tCOMPLETE\tno\n", "\tCOMPLETE\tyes\n", 1)
    peptide.write_text(text, encoding="utf-8")
    (result_dir / "scoring/all_tool_results.tsv").write_text(text, encoding="utf-8")
    result = run_review({"result_dir": str(result_dir), "outdir": str(tmp_path / "review")})
    assert result["status"] == "BLOCKED"


def test_review_priority_uses_evidence_reason_codes_and_phase_deduplication():
    events = [
        {"event_evidence_rank": "1", "event_group_id": "G1", "event_id": "E1", "gene": "FUS", "event_type": "Fusion", "best_evidence_grade": "R3", "representative_1_peptide_id": "P1", "representative_1_peptide": "AAAA", "representative_1_hla_allele": "HLA-A*02:01"},
        {"event_evidence_rank": "2", "event_group_id": "G2", "event_id": "E2", "gene": "SPL", "event_type": "Splice", "best_evidence_grade": "R3", "phase_group_id": "PH1", "representative_1_peptide_id": "P2", "representative_1_peptide": "BBBB", "representative_1_hla_allele": "HLA-B*07:02"},
        {"event_evidence_rank": "3", "event_group_id": "G3", "event_id": "E3", "gene": "DRV", "event_type": "SNV", "best_evidence_grade": "R4", "manual_review_required": "yes", "phase_group_id": "PH1", "representative_1_peptide_id": "P3", "representative_1_peptide": "CCCC", "representative_1_hla_allele": "HLA-C*07:02"},
    ]
    peptides = [
        {"peptide_id": "P1", "rna_support_state": "RNA_UNASSESSED", "safety_state": "SAFETY_PARTIAL"},
        {"peptide_id": "P2", "rna_support_state": "RNA_UNASSESSED", "safety_state": "SAFETY_PARTIAL"},
        {"peptide_id": "P3", "rna_support_state": "RNA_CONFIRMED", "safety_state": "SAFETY_PASS"},
    ]
    all_tool = [{"peptide_id": "P1", "fusion_consensus_status": "SINGLE_CALLER_SUPPORTED"}]
    rows = build_review_rows(events, peptides, all_tool)
    assert rows[0]["experiment_priority"] == "FUSION_CONFIRMATION_FIRST"
    assert rows[0]["review_reason"] == "REVIEW_FUSION_SINGLE_CALLER"
    assert rows[1]["experiment_priority"] == "TARGETED_RNA_FIRST"
    assert rows[2]["experiment_priority"] == "MANUAL_REVIEW_ONLY"
    selected = select_first_batch(rows, 12)
    assert len([row for row in selected if row.get("phase_group_id") == "PH1"]) <= 1


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
