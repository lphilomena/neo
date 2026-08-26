from pathlib import Path

from neoag.sv.wes_adapter import WESAdapter, build_sv_wes_phase1_5_raw
from neoag.sv.wes_filter import classify_wes_tier
from neoag.utils import read_tsv

ROOT = Path(__file__).resolve().parents[1]
FX = ROOT / "data" / "fixtures_sv"


def test_classify_wes_tier():
    assert classify_wes_tier({"rna_junction_reads": 3}) == "WES_UNINTERPRETABLE"
    assert classify_wes_tier({"rna_junction_reads": 3, "event_confidence_tier": "Tier1", "capture_interpretability": "HIGH"}) == "WES_Tier1"
    assert classify_wes_tier({"rna_junction_reads": 1, "event_confidence_tier": "Tier2", "capture_interpretability": "HIGH"}) == "WES_Tier2"


def test_wes_adapter_builds_raw_tables(tmp_path):
    out = build_sv_wes_phase1_5_raw(
        sample_id="SVMINI",
        sv_vcfs=[FX / "mini_sv.vcf"],
        callers=["GRIDSS2"],
        tumor_sample_name="TUMOR",
        normal_sample_name="NORMAL",
        reference_fasta=FX / "mini_ref.fa",
        gencode_gtf=FX / "mini.gtf",
        hla=FX / "hla.txt",
        outdir=tmp_path,
        expression_tsv=FX / "expression.tsv",
        rna_junction_tsv=FX / "rna_junctions_exact.tsv",
        expressed_products_tsv=FX / "expressed_products.tsv",
        capture_bed=FX / "capture.bed",
        normal_expression_tsv=FX / "normal_expression.tsv",
        normal_hla_ligands_tsv=FX / "normal_hla_ligands.tsv",
    )
    raw_events = read_tsv(out["raw_events"])
    sv_events = read_tsv(out["sv_events_full"])
    assert out["provenance"].endswith("provenance.sv_wes_phase1_5.json")
    assert len(raw_events) == 1
    assert sv_events[0]["event_confidence_tier"] == "WES_Tier1"
    assert sv_events[0]["final_sv_confidence"] == "WES_Tier1"


def test_wes_adapter_class_wrapper(tmp_path):
    adapter = WESAdapter()
    out = adapter.build_raw(
        sample_id="SVMINI",
        sv_vcfs=[FX / "mini_sv.vcf"],
        callers=["GRIDSS2"],
        tumor_sample_name="TUMOR",
        normal_sample_name="NORMAL",
        reference_fasta=FX / "mini_ref.fa",
        gencode_gtf=FX / "mini.gtf",
        hla=FX / "hla.txt",
        outdir=tmp_path,
        rna_junction_tsv=FX / "rna_junctions_exact.tsv",
        expressed_products_tsv=FX / "expressed_products.tsv",
        capture_bed=FX / "capture.bed",
    )
    assert Path(out["raw_peptides"]).is_file()
