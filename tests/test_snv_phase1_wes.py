from pathlib import Path

from neoag.snv_call.pipeline import run_snv_wes_full, write_snv_run_config

ROOT = Path(__file__).resolve().parents[1]
FX = ROOT / "data" / "fixtures_snv"


def test_write_snv_run_config(tmp_path):
    cfg = write_snv_run_config(
        path=tmp_path / "run.toml",
        sample_id="SNVMINI",
        profile="default",
        tumor_vcf=FX / "mini_somatic.vcf",
        hla_alleles=["HLA-A*02:01", "HLA-B*07:02"],
        tumor_sample_name="SNVMINI_TUMOR",
        normal_sample_name="SNVMINI_NORMAL",
        stub=True,
    )
    text = cfg.read_text(encoding="utf-8")
    assert "tumor_vcf" in text
    assert "SNVMINI_TUMOR" in text
    assert "stub = true" in text


def test_snv_run_full_wes_fixture_stub(tmp_path):
    out = run_snv_wes_full(
        outdir=tmp_path,
        sample_id="SNVMINI",
        profile="default",
        hla_alleles=["HLA-A*02:01", "HLA-B*07:02"],
        tumor_sample_name="SNVMINI_TUMOR",
        normal_sample_name="SNVMINI_NORMAL",
        somatic_vcf=FX / "mini_somatic.vcf",
        upstream_stub=True,
        immunogenicity_stub=True,
    )
    assert Path(out["mutect2_filtered_vcf"]).is_file()
    assert Path(out["ranked_peptides"]).is_file()
    assert Path(out["evidence_report"]).is_file()
    assert Path(out["run_config"]).is_file()
