from pathlib import Path

from neoag_v03.adapters.variant_peptide_adapter import (
    variant_peptide_extraction_enabled,
    variant_peptide_rows_to_raw_tables,
)
from neoag_v03.tools.upstream import run_upstream

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


def test_variant_peptide_rows_to_raw_tables_expands_hla():
    row = {
        "gene": "GENE1",
        "variant_key": "GENE1|chr1|100|A|T",
        "consequence": "missense_variant",
        "mutant_peptide": "CDENGHIK",
        "wildtype_peptide": "CDEFGHIK",
        "hgvsp": "p.Lys4Asn",
        "chrom": "chr1",
        "pos": "100",
        "ref": "A",
        "alt": "T",
        "transcript_id": "ENST00001.1",
        "multi_aa_flag": "single_aa",
    }
    events, peptides = variant_peptide_rows_to_raw_tables(
        [row],
        sample_id="S1",
        profile_name="default",
        hla_alleles=["HLA-A*02:01", "HLA-B*07:02"],
    )
    assert len(events) == 1
    assert events[0]["source"] == "extract-variant-peptides"
    assert len(peptides) == 2
    assert {p["hla_allele"] for p in peptides} == {"HLA-A*02:01", "HLA-B*07:02"}


def test_variant_peptide_extraction_default_for_snv_indel(tmp_path):
    vcf = tmp_path / "dummy.vcf"
    vcf.write_text("x\n", encoding="utf-8")
    cfg = {"inputs": {"entry_mode": "snv_indel", "variants_vcf": str(vcf)}}
    assert variant_peptide_extraction_enabled(cfg, vcf) is True
    cfg_pvac = {"inputs": {"entry_mode": "snv_indel", "pvac_files": ["a.tsv"]}}
    assert variant_peptide_extraction_enabled(cfg_pvac, vcf) is False


def test_run_upstream_variant_peptide_extraction_stub(tmp_path, monkeypatch):
    vcf = tmp_path / "mini.vcf"
    _write_mini_vcf(vcf)
    cfg_path = tmp_path / "run.toml"
    cfg_path.write_text(
        f"""
[sample]
id = "VPTEST"
profile = "default"

[tools]
stub = true
enabled = ["netmhcpan"]

[inputs]
entry_mode = "snv_indel"
variant_peptide_extraction = true
variants_vcf = "{vcf}"
hla_alleles = ["HLA-A*02:01"]
extract_appm_from_vcf = false
""",
        encoding="utf-8",
    )
    outs = run_upstream(cfg_path, tmp_path / "up")
    assert outs.get("peptide_source") == "extract-variant-peptides"
    assert Path(outs["raw_peptides"]).is_file()
    assert Path(outs["variant_peptides"]).is_file()
    assert Path(outs["netmhcpan"]).is_file()


def test_unique_peptide_hla_pairs_includes_wildtype_peptides(tmp_path):
    from neoag_v03.tools.prep import unique_peptide_hla_pairs
    from neoag_v03.utils import write_tsv

    raw = tmp_path / "raw_peptides.tsv"
    write_tsv(
        raw,
        [{
            "peptide": "MTPEPTIDE",
            "wildtype_peptide": "WTPEPTIDE",
            "hla_allele": "HLA-A*02:06",
        }],
        ["peptide", "wildtype_peptide", "hla_allele"],
    )
    assert unique_peptide_hla_pairs(raw) == [
        ("MTPEPTIDE", "HLA-A*02:06"),
        ("WTPEPTIDE", "HLA-A*02:06"),
    ]
