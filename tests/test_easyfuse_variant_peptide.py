"""Tests for unified EasyFuse catalog builder."""

from pathlib import Path

from neoag.adapters.easyfuse_variant_peptide import (
    EasyFusePeptideConfig,
    build_easyfuse_catalog,
    build_fusion_centered_minigene,
    easyfuse_to_variant_peptide_rows,
    sliding_fusion_neo_peptides,
    write_easyfuse_qc_tables,
)
from neoag.adapters.variant_peptide_adapter import (
    _catalog_rows_to_raw_peptides,
    run_variant_peptide_upstream,
)

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "data/fixtures/easyfuse_fusions.pass.tsv"

CSQ_HEADER = (
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
    "Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|"
    "HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|"
    "Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|MANE_SELECT|"
    'TSL|HGVS_OFFSET|FrameshiftSequence|WildtypeProtein">'
)


def _write_mini_vcf(path: Path) -> None:
    path.write_text(
        "\n".join([
            CSQ_HEADER,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chr1\t100\t.\tA\tT\t.\tPASS\tCSQ=T|missense_variant|MODERATE|GENE1|ENSG0001|"
            "Transcript|ENST00001.1|protein_coding|1/2||ENST00001.1:c.10A>T|"
            "ENSP00001.1:p.Lys4Asn|10|10|4|K/N|||1||HGNC|HGNC:1|YES||MANE_Select|"
            "NM_0001.1|1|||ACDEFGHIKLMNPQRSTVWYACDEFG",
        ]) + "\n",
        encoding="utf-8",
    )


def test_sliding_fusion_neo_peptides_dedupes_windows():
    windows = sliding_fusion_neo_peptides("ABCDEFGHIJK", (8, 9), bp_pos=5)
    peptides = {w["mutant_peptide"] for w in windows}
    assert "ABCDEHGI" not in peptides
    assert len(peptides) == len(windows)
    assert all(8 <= len(p) <= 9 for p in peptides)




def test_build_fusion_centered_minigene():
    mg = build_fusion_centered_minigene("ABCDEFGHIJKLMN", peptide_start_aa=4, peptide_end_aa=11, mini_len=2)
    assert mg == "BC|DEFGHIJK|LM"
    assert mg.split("|")[1] == "DEFGHIJK"


def test_easyfuse_catalog_minigene_centers_each_fusion_peptide():
    result = build_easyfuse_catalog(FIXTURE, "S1", "default", lengths=(8,), peptide_cfg=EasyFusePeptideConfig(junction_only=False, mini_len=3))
    assert result.catalog_rows
    assert all(row["minigene"].split("|")[1] == row["mutant_peptide"] for row in result.catalog_rows)

def test_build_easyfuse_catalog_strict_defaults():
    result = build_easyfuse_catalog(FIXTURE, "S1", "default", lengths=(8, 9))
    assert len(result.events) == 2
    assert result.catalog_rows
    assert all(r["multi_aa_flag"] == "fusion_neo" for r in result.catalog_rows)
    assert all(r["peptide_source"] == "easyfuse" for r in result.catalog_rows)
    assert all(
        r["fusion_generation_method"].startswith("easyfuse_neo_peptide_sliding_window")
        for r in result.catalog_rows
    )
    assert result.summary_qc["junction_only"] == "yes"
    assert result.summary_qc["dedup_per_event"] == "yes"
    assert len(result.filter_qc) == 4
    assert len(result.collapse_qc) == 2


def test_easyfuse_to_variant_peptide_rows_wrapper_has_no_legacy_peptides():
    rows, events, legacy = easyfuse_to_variant_peptide_rows(FIXTURE, "S1", "default", lengths=(8, 9))
    assert events
    assert rows
    assert legacy == []


def test_isoform_strategy_max_junction_reads(tmp_path):
    dup = tmp_path / "dup.tsv"
    dup.write_text(
        "\t".join([
            "BPID", "FTID", "Fusion_Gene", "Breakpoint1", "Breakpoint2", "type", "frame",
            "neo_peptide_sequence", "neo_peptide_sequence_bp", "prediction_prob",
            "prediction_class", "ft_junc_cnt", "ft_anch_cnt",
        ]) + "\n"
        + "\t".join([
            "bp1", "FT_short", "GENE1_GENE2", "chr1:1:+", "chr2:2:+", "trans", "in_frame",
            "SHORTSEQ", "3", "0.9", "positive", "5", "12",
        ]) + "\n"
        + "\t".join([
            "bp1", "FT_long", "GENE1_GENE2", "chr1:1:+", "chr2:2:+", "trans", "in_frame",
            "SHORTSEQEXTRA", "3", "0.7", "positive", "20", "12",
        ]) + "\n",
        encoding="utf-8",
    )
    result = build_easyfuse_catalog(
        dup,
        "S1",
        "default",
        lengths=(8,),
        peptide_cfg=EasyFusePeptideConfig(
            junction_only=False,
            dedup_per_event=True,
            isoform_strategy="max_junction_reads",
        ),
    )
    assert result.collapse_qc[0]["selected_ftid"] == "FT_long"


def test_write_easyfuse_qc_tables(tmp_path):
    result = build_easyfuse_catalog(FIXTURE, "S1", "default", lengths=(8,))
    paths = write_easyfuse_qc_tables(result, tmp_path)
    assert Path(paths["filter_qc"]).is_file()
    assert Path(paths["collapse_qc"]).is_file()
    assert Path(paths["summary_qc"]).is_file()


def test_catalog_rows_to_raw_peptides_links_fusion_events():
    result = build_easyfuse_catalog(FIXTURE, "S1", "default", lengths=(8,))
    event_map = {e["event_id"]: e for e in result.events}
    peptides = _catalog_rows_to_raw_peptides(
        result.catalog_rows[:2],
        event_map,
        sample_id="S1",
        hla_alleles=["HLA-A*02:01"],
        source_tool="EasyFuse",
    )
    assert peptides
    assert peptides[0]["event_type"] == "Fusion"
    assert peptides[0]["source_tool"] == "EasyFuse"
    assert peptides[0]["crosses_junction"] == "yes"


def test_run_variant_peptide_upstream_merges_easyfuse(tmp_path):
    vcf = tmp_path / "mini.vcf"
    _write_mini_vcf(vcf)
    cfg = {
        "inputs": {
            "variant_peptide_extraction": True,
            "variant_peptide_lengths": "8",
            "variant_peptide_mini_len": 8,
            "variant_peptide_filter_normal_proteome": False,
            "easyfuse_pass_csv": str(FIXTURE),
            "tumor_sample_name": "TUMOR",
        }
    }
    parsed = tmp_path / "parsed"
    tools = tmp_path / "tools"
    parsed.mkdir()
    tools.mkdir()

    outs = run_variant_peptide_upstream(
        cfg,
        variants_vcf=vcf,
        parsed_dir=parsed,
        tools_dir=tools,
        sample_id="S1",
        profile_name="default",
        hla_alleles=["HLA-A*02:01"],
    )
    assert int(outs["variant_peptide_rows_vcf"]) > 0
    assert int(outs["variant_peptide_rows_easyfuse"]) > 0
    assert Path(outs["easyfuse_variant_peptides"]).is_file()
    assert Path(outs["filter_qc"]).is_file()
    assert Path(outs["summary_qc"]).is_file()
    catalog = (tmp_path / "tools" / "variant_peptides.tsv").read_text(encoding="utf-8")
    assert "fusion_neo" in catalog
    assert "\teasyfuse\t" in catalog or catalog.endswith("easyfuse\n")
