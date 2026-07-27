import json
from pathlib import Path

from neoag.workflow_selection import select_workflow, write_selection


def _write(path: Path, content: str) -> Path:
    path.write_text(content, encoding="utf-8")
    return path


def _support_manifests(tmp_path: Path) -> tuple[Path, Path]:
    tools = _write(
        tmp_path / "tools.yaml",
        """tools:
  facets:
    mode: conda
    executable: runFACETS.R
  lohhla:
    mode: local_or_container
    executable: LOHHLA
""",
    )
    references = _write(
        tmp_path / "references.yaml",
        """genome_build: GRCh38
reference_fasta: reference.fa
gencode_gtf: gencode.gtf
vep_cache: vep_cache
""",
    )
    return tools, references


def test_selects_vcf_to_ranking_without_executing(tmp_path):
    vcf = _write(tmp_path / "sample.vcf", "##fileformat=VCFv4.2\n")
    hla = _write(tmp_path / "hla.txt", "HLA-A*02:01\n")
    sample = _write(
        tmp_path / "sample.yaml",
        f"""sample_id: VCF1
inputs:
  somatic_vcf: {vcf}
  hla_typing: {hla}
""",
    )
    tools, references = _support_manifests(tmp_path)

    result = select_workflow(sample, tools_manifest=tools, reference_manifest=references, outdir=tmp_path / "out")

    assert result.workflow == "vcf_to_ranking"
    assert not result.execution_allowed
    assert any(stage.name == "candidate_ingestion" and stage.status == "PLANNED" for stage in result.stages)
    assert result.detected_inputs["somatic_vcf"]["path"] == "sample.vcf"


def test_paired_bam_plans_consensus_qc_without_fake_results(tmp_path):
    tumor = _write(tmp_path / "tumor.bam", "fixture")
    normal = _write(tmp_path / "normal.bam", "fixture")
    sample = _write(
        tmp_path / "sample.yaml",
        f"""sample_id: BAM1
tumor:
  dna_bam: {tumor}
normal:
  dna_bam: {normal}
inputs: {{}}
""",
    )
    tools, references = _support_manifests(tmp_path)

    result = select_workflow(sample, tools_manifest=tools, reference_manifest=references)
    by_name = {stage.name: stage for stage in result.stages}

    assert result.workflow == "paired_bam_production"
    assert by_name["variant_calling"].risk_level == "HIGH"
    assert by_name["purity_cnv_consensus"].status == "PLANNED"
    assert by_name["hla_loh_consensus"].status == "PLANNED"
    assert all(stage.status != "PASS" for stage in result.stages)


def test_existing_ranked_results_route_to_review(tmp_path):
    ranked = _write(tmp_path / "ranked.tsv", "peptide\thla_allele\nAAAAAAAAL\tHLA-A*02:01\n")
    sample = _write(
        tmp_path / "sample.yaml",
        f"""sample_id: REVIEW1
inputs:
  ranked_peptides_recommendation: {ranked}
""",
    )
    tools, references = _support_manifests(tmp_path)
    result = select_workflow(sample, tools_manifest=tools, reference_manifest=references)

    assert result.workflow == "result_review"
    assert result.entry_mode == "ranked_results"
    assert result.status == "READY_TO_PLAN"


def test_missing_hla_source_is_blocked_and_outputs_are_deterministic(tmp_path):
    vcf = _write(tmp_path / "sample.vcf", "##fileformat=VCFv4.2\n")
    sample = _write(
        tmp_path / "sample.yaml",
        f"""sample_id: BLOCK1
inputs:
  somatic_vcf: {vcf}
""",
    )
    tools, references = _support_manifests(tmp_path)
    result = select_workflow(sample, tools_manifest=tools, reference_manifest=references)
    outputs = write_selection(result, tmp_path / "out")

    assert result.status == "BLOCKED"
    assert "hla_typing" in result.missing_required
    payload = json.loads(Path(outputs["json"]).read_text(encoding="utf-8"))
    assert payload["workflow"] == "insufficient_inputs"
    assert not payload["execution_allowed"]
