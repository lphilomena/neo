from pathlib import Path

from neoag_v03.sv.bnd_parser import parse_bnd_alt
from neoag_v03.sv.phase1 import build_sv_phase1_raw
from neoag_v03.utils import read_tsv

ROOT = Path(__file__).resolve().parents[1]
FX = ROOT / "data" / "fixtures_sv"


def test_bnd_parser_four_forms():
    examples = ["N]chr2:123]", "N[chr2:123[", "]chr2:123]N", "[chr2:123[N"]
    for alt in examples:
        parsed = parse_bnd_alt(alt)
        assert parsed is not None
        assert parsed.mate_chrom == "chr2"
        assert parsed.mate_pos == 123
        assert parsed.strand1 in {"+", "-"}
        assert parsed.strand2 in {"+", "-"}


def test_sv_phase1_builds_raw_tables(tmp_path):
    out = build_sv_phase1_raw(
        sample_id="SVMINI",
        sv_vcfs=[FX / "mini_sv.vcf"],
        callers=["GRIDSS2"],
        reference_fasta=FX / "mini_ref.fa",
        gencode_gtf=FX / "mini.gtf",
        hla=FX / "hla.txt",
        outdir=tmp_path,
        expression_tsv=FX / "expression.tsv",
        rna_junction_tsv=FX / "rna_junctions.tsv",
        normal_expression_tsv=FX / "normal_expression.tsv",
        normal_hla_ligands_tsv=FX / "normal_hla_ligands.tsv",
    )
    raw_events = read_tsv(out["raw_events"])
    raw_peptides = read_tsv(out["raw_peptides"])
    sv_events = read_tsv(out["sv_events_full"])
    proteins = read_tsv(out["sv_protein_reconstruction"])
    assert len(raw_events) == 1
    assert raw_events[0]["event_type"] == "SV_Fusion"
    assert raw_events[0]["gene"] == "GENE1::GENE2"
    assert len(raw_peptides) > 0
    assert all(p["hla_allele"] in {"HLA-A*02:01", "HLA-B*07:02"} for p in raw_peptides)
    assert sv_events[0]["rna_support_status"] == "RNA_JUNCTION_SUPPORTED"
    assert proteins[0]["reconstruction_method"] == "heuristic_cds_prefix_suffix"
