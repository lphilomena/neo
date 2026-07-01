from pathlib import Path

from neoag_v03.adapters.pvacseq_enrich import (
    enrich_pvacseq_aggregated,
    filter_raw_peptides_normal_proteome,
    parse_pvac_variant_id,
    pvacseq_enrich_enabled,
    refresh_raw_peptides_from_enriched,
    resolve_pvacseq_enrich_options,
    variant_lookup_keys,
)
from neoag_v03.schemas import PEPTIDE_FIELDS
from neoag_v03.utils import read_tsv, write_tsv

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


def test_parse_pvac_variant_id():
    assert parse_pvac_variant_id("chr20-59878018-59878019-T-G") == (
        "chr20", "59878018", "T", "G"
    )
    assert parse_pvac_variant_id("chr1-100-100-A-T") == ("chr1", "100", "A", "T")


def test_variant_lookup_keys():
    keys = variant_lookup_keys("chr1", "100", "A", "T")
    assert "1:100:A>T" in keys
    assert "chr1:100:A>T" in keys


def test_pvacseq_enrich_enabled_defaults(tmp_path):
    vcf = tmp_path / "annotated.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    cfg = {"inputs": {}}
    assert pvacseq_enrich_enabled(cfg, vcf, has_pvacseq_output=True) is True
    cfg_off = {"inputs": {"pvacseq_enrich": False}}
    assert pvacseq_enrich_enabled(cfg_off, vcf, has_pvacseq_output=True) is False


def test_enrich_pvacseq_aggregated_minigene(tmp_path):
    vcf = tmp_path / "mini.vcf"
    _write_mini_vcf(vcf)
    agg = tmp_path / "pvacseq_aggregated.tsv"
    agg.write_text(
        "ID\tIndex\tGene\tBest Peptide\tAllele\n"
        "chr1-100-100-A-T\t1.GENE1.ENST00001.1.missense.4K/N\tGENE1\tCDENGHIK\tHLA-A*02:01\n",
        encoding="utf-8",
    )
    out = tmp_path / "pvacseq_enriched.tsv"
    summary = enrich_pvacseq_aggregated(
        agg,
        vcf,
        out,
        enrich_minigene=True,
        mini_len=2,
        normal_proteome_fasta=None,
    )
    assert summary["variants_matched"] == 1
    text = out.read_text(encoding="utf-8")
    assert "minigene" in text
    assert "multi_aa_flag" in text
    assert "|" in text  # minigene delimiter


def test_enrich_pvacseq_aggregated_hla_columns(tmp_path):
    vcf = tmp_path / "mini.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    agg = tmp_path / "pvacseq_aggregated.tsv"
    agg.write_text(
        "ID\tGene\tBest Peptide\tAllele\tIC50 MT\tIC50 WT\t%ile MT\t%ile WT\n"
        "chr1-100-100-A-T\tGENE1\tCDENGHIK\tHLA-A*02:06\t50.0\t80.0\t1.2\t2.0\n",
        encoding="utf-8",
    )
    out = tmp_path / "pvacseq_enriched.tsv"
    enrich_pvacseq_aggregated(
        agg,
        vcf,
        out,
        enrich_minigene=False,
        hla_alleles=["HLA-A*02:06", "HLA-B*13:02"],
    )
    text = out.read_text(encoding="utf-8")
    assert "sample_hla_alleles" in text
    assert "netmhcpan_mt_ic50" in text


def test_filter_raw_peptides_normal_proteome(tmp_path):
    proteome = tmp_path / "ref.fa"
    proteome.write_text(">p1\nACDEFGHIKLMN\n", encoding="utf-8")
    raw = tmp_path / "raw_peptides.tsv"
    write_tsv(
        raw,
        [
            {
                "peptide_id": "p1",
                "event_id": "e1",
                "sample_id": "S1",
                "peptide": "CDEFGHIK",
                "gene": "G1",
                "hla_allele": "HLA-A*02:01",
            },
            {
                "peptide_id": "p2",
                "event_id": "e2",
                "sample_id": "S1",
                "peptide": "ZZZZZZZZ",
                "gene": "G2",
                "hla_allele": "HLA-A*02:01",
            },
        ],
        PEPTIDE_FIELDS,
    )
    summary = filter_raw_peptides_normal_proteome(
        raw,
        normal_proteome_fasta=proteome,
        annotate_only=False,
    )
    assert summary["peptides_filtered_normal_proteome"] == 1
    kept = raw.read_text(encoding="utf-8")
    assert "ZZZZZZZZ" in kept
    assert "CDEFGHIK" not in kept


def test_refresh_raw_peptides_wildtype_from_enriched(tmp_path):
    enriched = tmp_path / "pvacseq_enriched.tsv"
    enriched.write_text(
        "ID\tIndex\tGene\tBest Peptide\tAllele\tWT Epitope Seq\t"
        "netmhcpan_wt_rank_ba\tnetmhcpan_wt_ic50\n"
        "chr1-100-100-A-T\t1.GENE1.ENST00001.1.missense.4K/N\tGENE1\t"
        "CDENGHIK\tHLA-A*02:06\tCDEFGHIJK\t0.42\t1200.5\n",
        encoding="utf-8",
    )
    raw = tmp_path / "raw_peptides.tsv"
    write_tsv(
        raw,
        [
            {
                "peptide_id": "p1",
                "event_id": "1.GENE1.ENST00001.1.missense.4K/N",
                "sample_id": "S1",
                "peptide": "CDENGHIK",
                "gene": "GENE1",
                "hla_allele": "HLA-A*02:06",
                "wildtype_peptide": "",
                "wildtype_binding_rank": "99",
            }
        ],
        PEPTIDE_FIELDS,
    )
    refresh_raw_peptides_from_enriched(
        enriched,
        sample_id="S1",
        profile_name="default",
        raw_peptides_tsv=raw,
    )
    rows = read_tsv(raw)
    assert rows[0]["wildtype_peptide"] == "CDEFGHIJK"
    assert float(rows[0]["wildtype_binding_rank"]) == 0.42
    assert rows[0]["netmhcpan_wt_ic50"] == "1200.5"


def test_resolve_pvacseq_enrich_options_env(monkeypatch):
    monkeypatch.setenv("NEOAG_NORMAL_PROTEOME_FASTA", "/ref/proteome.fa")
    opts = resolve_pvacseq_enrich_options({"inputs": {"pvacseq_enrich_minigene": False}})
    assert opts["enrich_minigene"] is False
    assert opts["normal_proteome_fasta"] == "/ref/proteome.fa"
    assert opts["filter_normal_proteome"] is True
