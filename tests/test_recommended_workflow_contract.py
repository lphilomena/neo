from pathlib import Path
import csv
import subprocess
import sys

from neoag.agent_skills import hla_typing_compare, purity_cnv_review


ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path: Path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_deploy_skill_installs_bam_matcher_to_shared_prefix():
    script = (ROOT / ".agents/skills/neoag-remote-deploy/scripts/13_install_readme_tools.sh").read_text()
    assert "--bam-matcher" in script
    assert "INSTALL_BAM_MATCHER" in script
    assert "NEOAG_BAM_MATCHER_ENV_PREFIX" in script
    installer = (ROOT / "scripts/install_bam_matcher.sh").read_text()
    assert 'create -y -p "$ENV_PREFIX"' in installer
    assert 'exec "$ENV_PREFIX/bin/python2"' in installer
    assert "Cheetah3==3.2.6.post2" in installer
    assert "new = '\"-i\", \"-X\", \"-u\", \"--min-coverage\"'" in installer
    assert "bam_matcher.identity.hg38.vcf.gz" in installer
    assert "contamination.common.hg38.vcf.gz" in installer


def test_real_input_qc_contract_is_not_manifest_only():
    script = (ROOT / "scripts/run_case_input_qc.py").read_text()
    for value in ("samtools", "quickcheck", "flagstat", "coverage", "bam-matcher", "GetPileupSummaries", "CalculateContamination", "sample_pairing.tsv", "vcf_qc.tsv", "UNASSESSED"):
        assert value in script


def test_hla_compare_outputs_two_field_high_resolution_and_recommendation(tmp_path):
    spechla = tmp_path / "spechla" / "sample.hla.result.txt"
    hlala = tmp_path / "hla-la" / "R1_bestguess_G.txt"
    spechla.parent.mkdir()
    hlala.parent.mkdir()
    text = "HLA-A*02:06:01 HLA-A*30:01:01 HLA-B*13:02:01 HLA-B*48:01:01 HLA-C*06:02:01 HLA-C*08:01:01\n"
    spechla.write_text(text)
    hlala.write_text(text)
    out = tmp_path / "out"
    assert hla_typing_compare.main(["--result-dir", str(tmp_path), "--outdir", str(out)]) == 0
    assert (out / "hla_typing.standardized.tsv").is_file()
    assert (out / "recommended_hla.txt").read_text().splitlines() == [
        "A*02:06", "A*30:01", "B*13:02", "B*48:01", "C*06:02", "C*08:01"
    ]
    rows = read_tsv(out / "hla_typing_consensus.tsv")
    assert rows[0]["consensus_2field"] == "A*02:06 / A*30:01"
    assert rows[0]["consensus_high_resolution"] == "A*02:06:01 / A*30:01:01"


def test_purity_review_writes_recommended_outputs(tmp_path):
    facets = tmp_path / "sample_facets_robust"
    facets.mkdir()
    (facets / "purity.tsv").write_text("sample_id\tpurity\tploidy\nS\t0.40\t2.10\n")
    (facets / "cnv_segments.tsv").write_text("chromosome\tstart\tend\tmajor_cn\tminor_cn\nchr1\t1\t10\t1\t1\n")
    out = tmp_path / "review"
    assert purity_cnv_review.main(["--result-dir", str(facets), "--sample-id", "S", "--outdir", str(out)]) == 0
    for filename in ("purity_cnv_consensus.tsv", "recommended_purity.tsv", "recommended_cnv_segments.tsv", "purity_cnv_review.md"):
        assert (out / filename).is_file()
    segments = read_tsv(out / "recommended_cnv_segments.tsv")
    assert segments[0]["major_cn"] == "1"
    assert segments[0]["minor_cn"] == "1"
    assert segments[0]["loh_status"] == "RETAINED"
    assert (out / "hla_6p21_cnv_tool_evidence.tsv").is_file()
    assert (out / "hla_6p21_cnv_consensus.tsv").is_file()


def test_hla_loh_consensus_preserves_all_four_states(tmp_path):
    lohhla = tmp_path / "lohhla.tsv"
    spechla = tmp_path / "spechla.tsv"
    lohhla.write_text("hla_allele\tloh_status\nHLA-A*02:06\tloh\nHLA-B*13:02\tno\nHLA-C*06:02\tloh\n")
    spechla.write_text("hla_allele\tloh_status\nHLA-A*02:06\tloh\nHLA-B*13:02\tno\nHLA-C*06:02\tno\nHLA-DRB1*01:01\tloh\n")
    out = tmp_path / "out"
    subprocess.run([
        sys.executable, str(ROOT / "scripts/build_hla_loh_consensus.py"), "--sample-id", "S",
        "--lohhla", str(lohhla), "--spechla", str(spechla), "--outdir", str(out),
    ], check=True)
    states = {row["hla_allele"]: row["consensus_status"] for row in read_tsv(out / "hla_loh_consensus.tsv")}
    assert states == {"HLA-A*02:06": "CONSENSUS_LOST", "HLA-B*13:02": "CONSENSUS_RETAINED", "HLA-C*06:02": "DISCORDANT"}
    assert "DRB1" not in (out / "hla_loh_consensus.tsv").read_text()


def test_sample_runners_require_real_completion_outputs():
    contracts = {
        "run_optitype_sample.sh": ("_result.tsv", ".complete"),
        "run_spechla_sample.sh": ("hla.result.txt", ".complete"),
        "run_purple_sample.sh": ("purple.purity.tsv", "purple.cnv.somatic.tsv", ".complete"),
    }
    for filename, required in contracts.items():
        text = (ROOT / "scripts" / filename).read_text()
        for value in required:
            assert value in text
    optitype = (ROOT / "scripts/run_optitype_sample.sh").read_text()
    assert '"$OPTITYPE_BIN" run -i "$r1" -i "$r2"' in optitype
    assert "Reusing existing MHC FASTQs" in optitype


def test_recommended_driver_covers_all_steps_and_dna_only_semantics():
    driver = (ROOT / "scripts/run_recommended_workflow.sh").read_text()
    for value in (
        "run_case_input_qc.py",
        "BAM_MATCHER_REFERENCE",
        "run_spechla_sample.sh",
        "run_hla_la_sample.sh",
        "run_optitype_sample.sh",
        "FACETS_CVAL_PRE=50",
        "FACETS_CVAL_PROC=300",
        "FACETS_MIN_NHET=10",
        "FACETS_TARGET_ROWS=1000000",
        "run_sequenza_sample_by_chrom.sh",
        "run_purple_sample.sh",
        "run_hla_loh_multi_tool.sh",
        "generate_recommended_run_config.py",
        "neoag.cli run-full",
    ):
        assert value in driver
    generator = (ROOT / "scripts/generate_recommended_run_config.py").read_text()
    assert 'profile = "sarcoma"' in generator
    assert "if value:" in generator  # RNA inputs are omitted rather than imputed in DNA-only mode.
