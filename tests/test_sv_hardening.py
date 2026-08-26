from pathlib import Path

import pytest

from neoag.sv.exact_evidence import load_exact_junction_evidence
from neoag.sv.identity import canonical_breakpoint_key
from neoag.sv.peptide_builder import build_mhc1_peptides
from neoag.sv.phase1 import build_sv_phase1_raw
from neoag.sv.sv_callset import parse_vcf_records
from neoag.sv.sv_merge import cluster_sv_records
from neoag.utils import read_tsv


ROOT = Path(__file__).resolve().parents[1]
FX = ROOT / "data" / "fixtures_sv"


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_multisample_vcf_never_guesses_tumor_normal():
    with pytest.raises(ValueError, match="requires explicit"):
        parse_vcf_records(FX / "mini_sv.vcf", "GRIDSS2")


def test_nonpass_records_are_excluded_by_default(tmp_path):
    vcf = _write(
        tmp_path / "filtered.vcf",
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\tx\tN\tN]chr2:20]\t.\tLowQual\tSVTYPE=BND\n",
    )
    assert parse_vcf_records(vcf, "GRIDSS2") == []
    assert len(parse_vcf_records(vcf, "GRIDSS2", include_nonpass=True)) == 1


def test_swapped_breakends_have_one_identity_and_cluster(tmp_path):
    vcf = _write(
        tmp_path / "mates.vcf",
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\ta\tN\tN]chr2:20]\t.\tPASS\tSVTYPE=BND\n"
        "chr2\t20\tb\tN\t]chr1:10]N\t.\tPASS\tSVTYPE=BND\n",
    )
    records = parse_vcf_records(vcf, "GRIDSS2")
    assert len(cluster_sv_records(records, distance=0)) == 1
    assert canonical_breakpoint_key("GRCh38", "chr1", 10, "+", "chr2", 20, "-") == canonical_breakpoint_key("grch38", "2", 20, "-", "1", 10, "+")


def test_gene_pair_only_rna_input_is_rejected(tmp_path):
    path = _write(tmp_path / "gene_pair.tsv", "gene1\tgene2\tjunction_reads\nA\tB\t9\n")
    with pytest.raises(ValueError, match="chrom1"):
        load_exact_junction_evidence(path, default_build="GRCh38")


def test_dna_only_event_catalog_generates_no_peptides(tmp_path):
    out = build_sv_phase1_raw(
        sample_id="SVMINI",
        sv_vcfs=[FX / "mini_sv.vcf"],
        callers=["GRIDSS2"],
        tumor_sample_name="TUMOR",
        normal_sample_name="NORMAL",
        reference_fasta=FX / "mini_ref.fa",
        gencode_gtf=FX / "mini.gtf",
        hla=FX / "hla.txt",
        outdir=tmp_path,
    )
    assert read_tsv(out["raw_peptides"]) == []
    sv_events = read_tsv(out["sv_events_full"])
    assert sv_events[0]["reconstruction_status"] == "hypothesis_only"


def test_peptide_ids_are_stable_and_sha_based():
    from neoag.sv.exact_evidence import load_expressed_products

    product = next(iter(load_expressed_products(FX / "expressed_products.tsv", default_build="GRCh38").values()))
    protein = product.to_reconstruction("EVENT1", "S1")
    first = build_mhc1_peptides([protein], ["HLA-A*02:01"])
    second = build_mhc1_peptides([protein], ["HLA-A*02:01"])
    assert [p.peptide_id for p in first] == [p.peptide_id for p in second]
    assert all(len(p.peptide_id.split("_")[-1]) == 16 for p in first)


def test_wes_requires_capture_bed(tmp_path):
    with pytest.raises(ValueError, match="requires capture_bed"):
        build_sv_phase1_raw(
            sample_id="SVMINI",
            sv_vcfs=[FX / "mini_sv.vcf"],
            callers=["GRIDSS2"],
            tumor_sample_name="TUMOR",
            normal_sample_name="NORMAL",
            reference_fasta=FX / "mini_ref.fa",
            gencode_gtf=FX / "mini.gtf",
            hla=FX / "hla.txt",
            outdir=tmp_path,
            wes_mode=True,
        )
