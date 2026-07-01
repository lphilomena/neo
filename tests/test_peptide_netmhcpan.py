from pathlib import Path

from neoag_v03.adapters.peptide_netmhcpan import (
    annotate_raw_peptides_tsv,
    annotate_variant_peptide_row,
    annotate_variant_peptide_rows,
    format_sample_hla_alleles,
    netmhcpan_columns_from_pvac_row,
    resolve_peptide_source,
)
from neoag_v03.schemas import PEPTIDE_FIELDS
from neoag_v03.utils import read_tsv, write_tsv
from neoag_v03.adapters.pvacseq_enrich import enrich_pvacseq_aggregated


def test_format_sample_hla_alleles():
    assert format_sample_hla_alleles(["HLA-A*02:06", "HLA-B*13:02"]) == "HLA-A*02:06,HLA-B*13:02"


def test_resolve_peptide_source():
    assert resolve_peptide_source({"peptide_source": "easyfuse"}) == "easyfuse"
    assert resolve_peptide_source({"multi_aa_flag": "fusion_neo"}) == "easyfuse"
    assert resolve_peptide_source({"multi_aa_flag": "single_aa"}) == "snv"
    assert resolve_peptide_source({
        "generation_method": "easyfuse_neo_peptide_sliding_window;k=8",
    }) == "easyfuse"


def test_netmhcpan_columns_from_pvac_row():
    mapped = netmhcpan_columns_from_pvac_row({
        "Allele": "HLA-A*02:06",
        "IC50 MT": "85.63",
        "IC50 WT": "96.22",
        "%ile MT": "0.69",
        "%ile WT": "0.77",
    })
    assert mapped["hla_allele"] == "HLA-A*02:06"
    assert mapped["netmhcpan_mt_ic50"] == "85.63"
    assert mapped["netmhcpan_wt_ic50"] == "96.22"
    assert mapped["netmhcpan_mt_rank_ba"] == "0.69"
    assert mapped["netmhcpan_wt_rank_ba"] == "0.77"


def test_annotate_variant_peptide_rows_from_mhcflurry_index(tmp_path):
    csv_path = tmp_path / "mhcflurry.csv"
    csv_path.write_text(
        "peptide,allele,mhcflurry_affinity,mhcflurry_affinity_percentile,"
        "mhcflurry_processing_score,mhcflurry_presentation_score\n"
        "CDENGHIK,HLA-A*02:06,50.0,0.5,0.7,0.8\n"
        "CDEFGHIK,HLA-A*02:06,80.0,1.2,0.6,0.7\n",
        encoding="utf-8",
    )
    rows = [{
        "mutant_peptide": "CDENGHIK",
        "wildtype_peptide": "CDEFGHIK",
        "hla_allele": "HLA-A*02:06",
    }]
    out = annotate_variant_peptide_rows(
        rows,
        ["HLA-A*02:06"],
        mhcflurry_csv=csv_path,
    )
    assert out[0]["mhcflurry_mt_affinity"] == "50.0"
    assert float(out[0]["mhcflurry_mt_affinity_percentile"]) == 0.5
    assert out[0]["mhcflurry_wt_affinity"] == "80.0"
    assert float(out[0]["mhcflurry_wt_affinity_percentile"]) == 1.2


def test_annotate_variant_peptide_rows_adds_stabpan_prime_bigmhc(tmp_path):
    xls = tmp_path / "netmhcpan.xls"
    xls.write_text(
        "Pos\tPeptide\tHLA\tScore_EL\t%Rank_EL\tScore_BA\t%Rank_BA\tBindLevel\n"
        "0\tCDENGHIK\tHLA-A*02:06\t0.1\t0.5\t50.0\t1.2\t\n",
        encoding="utf-8",
    )
    stab = tmp_path / "stab.tsv"
    stab.write_text(
        "Peptide\tHLA\tscore\tpercentile_rank\n"
        "CDENGHIK\tHLA-A*02:06\t1.5\t0.8\n",
        encoding="utf-8",
    )
    prime = tmp_path / "prime.tsv"
    prime.write_text(
        "Peptide\tScore_A0206\t%Rank_A0206\n"
        "CDENGHIK\t0.42\t1.1\n",
        encoding="utf-8",
    )
    bigmhc = tmp_path / "bigmhc.tsv"
    bigmhc.write_text(
        "peptide\thla_allele\tbigmhc_im_score\n"
        "CDENGHIK\tHLA-A*02:06\t0.8123\n",
        encoding="utf-8",
    )
    rows = [{"mutant_peptide": "CDENGHIK", "wildtype_peptide": "CDEFGHIK"}]
    out = annotate_variant_peptide_rows(
        rows,
        ["HLA-A*02:06"],
        netmhcpan_xls=xls,
        netmhcstabpan_tsv=stab,
        prime_tsv=prime,
        bigmhc_im_tsv=bigmhc,
    )
    assert out[0]["netmhcstabpan_score"] == "1.5"
    assert out[0]["prime_score"] == "0.42"
    assert out[0]["bigmhc_im_score"] == "0.8123"
    assert out[0]["iedb_immunogenicity_score"]


def test_annotate_variant_peptide_rows_adds_iedb_from_evidence_tsv(tmp_path):
    iedb = tmp_path / "iedb.tsv"
    iedb.write_text(
        "peptide\thla_allele\tiedb_immunogenicity_score\n"
        "CDENGHIK\tHLA-A*02:06\t0.12345\n",
        encoding="utf-8",
    )
    rows = [{"mutant_peptide": "CDENGHIK", "wildtype_peptide": "CDEFGHIK", "hla_allele": "HLA-A*02:06"}]
    out = annotate_variant_peptide_rows(
        rows,
        ["HLA-A*02:06"],
        iedb_immunogenicity_tsv=iedb,
    )
    assert out[0]["iedb_immunogenicity_score"] == "0.12345"


def test_annotate_variant_peptide_rows_from_index(tmp_path):
    xls = tmp_path / "netmhcpan.xls"
    xls.write_text(
        "Pos\tPeptide\tHLA\tScore_EL\t%Rank_EL\tScore_BA\t%Rank_BA\tBindLevel\n"
        "0\tCDENGHIK\tHLA-A*02:06\t0.1\t0.5\t50.0\t1.2\t\n"
        "0\tCDEFGHIK\tHLA-A*02:06\t0.2\t1.0\t80.0\t2.0\t\n",
        encoding="utf-8",
    )
    rows = [{
        "mutant_peptide": "CDENGHIK",
        "wildtype_peptide": "CDEFGHIK",
        "gene": "GENE1",
    }]
    out = annotate_variant_peptide_rows(
        rows,
        ["HLA-A*02:06", "HLA-A*30:01"],
        netmhcpan_xls=xls,
    )
    assert out[0]["sample_hla_alleles"] == "HLA-A*02:06,HLA-A*30:01"
    assert out[0]["peptide_source"] == "snv"
    assert out[0]["hla_allele"] == "HLA-A*02:06"
    assert out[0]["netmhcpan_mt_ic50"] == "50.0"
    assert out[0]["netmhcpan_wt_ic50"] == "80.0"


def test_enrich_pvacseq_aggregated_adds_hla_mhcflurry_and_netmhcpan(tmp_path):
    vcf = tmp_path / "mini.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
    agg = tmp_path / "pvacseq_aggregated.tsv"
    agg.write_text(
        "ID\tGene\tBest Peptide\tAllele\tIC50 MT\tIC50 WT\t%ile MT\t%ile WT\n"
        "chr1-100-100-A-T\tGENE1\tCDENGHIK\tHLA-A*02:06\t50.0\t80.0\t1.2\t2.0\n",
        encoding="utf-8",
    )
    out = tmp_path / "pvacseq_enriched.tsv"
    mhc = tmp_path / "mhcflurry.csv"
    mhc.write_text(
        "peptide,allele,mhcflurry_affinity,mhcflurry_affinity_percentile,"
        "mhcflurry_processing_score,mhcflurry_presentation_score\n"
        "CDENGHIK,HLA-A*02:06,50.0,1.2,0.7,0.8\n"
        "CDEFGHIK,HLA-A*02:06,80.0,2.0,0.6,0.7\n",
        encoding="utf-8",
    )
    enrich_pvacseq_aggregated(
        agg,
        vcf,
        out,
        enrich_minigene=False,
        hla_alleles=["HLA-A*02:06", "HLA-A*30:01"],
        mhcflurry_csv=mhc,
    )
    text = out.read_text(encoding="utf-8")
    assert "sample_hla_alleles" in text
    assert "HLA-A*02:06,HLA-A*30:01" in text
    assert "netmhcpan_mt_ic50" in text
    assert "netmhcpan_wt_ic50" in text
    assert "mhcflurry_mt_affinity_percentile" in text
    assert "mhcflurry_wt_affinity_percentile" in text


def test_annotate_raw_peptides_tsv(tmp_path):
    xls = tmp_path / "netmhcpan.xls"
    xls.write_text(
        "Pos\tPeptide\tHLA\tScore_EL\t%Rank_EL\tScore_BA\t%Rank_BA\tBindLevel\n"
        "0\tMTPEPTIDE\tHLA-A*02:06\t0.1\t0.5\t100\t0.5\t\n"
        "1\tWTPEPTIDE\tHLA-A*02:06\t0.2\t1.2\t200\t1.2\t\n",
        encoding="utf-8",
    )
    raw = tmp_path / "raw_peptides.tsv"
    write_tsv(
        raw,
        [{
            "peptide_id": "p1",
            "event_id": "e1",
            "sample_id": "S1",
            "peptide": "MTPEPTIDE",
            "wildtype_peptide": "WTPEPTIDE",
            "hla_allele": "HLA-A*02:06",
            "gene": "G1",
        }],
        PEPTIDE_FIELDS,
    )
    mhc = tmp_path / "mhcflurry.csv"
    mhc.write_text(
        "peptide,allele,mhcflurry_affinity,mhcflurry_affinity_percentile,"
        "mhcflurry_processing_score,mhcflurry_presentation_score\n"
        "MTPEPTIDE,HLA-A*02:06,100,0.5,0.7,0.8\n"
        "WTPEPTIDE,HLA-A*02:06,200,1.2,0.6,0.7\n",
        encoding="utf-8",
    )
    annotate_raw_peptides_tsv(raw, xls, mhc)
    rows = read_tsv(raw)
    assert rows[0]["wildtype_peptide"] == "WTPEPTIDE"
    assert float(rows[0]["wildtype_binding_rank"]) == 1.2
    assert float(rows[0]["netmhcpan_mt_rank_ba"]) == 0.5
    assert float(rows[0]["mhcflurry_mt_affinity_percentile"]) == 0.5
    assert float(rows[0]["mhcflurry_wt_affinity_percentile"]) == 1.2


def test_annotate_variant_peptide_row_adds_wt_external_tool_scores():
    row = {
        "mutant_peptide": "MTPEPTIDE",
        "wildtype_peptide": "WTPEPTIDE",
        "hla_allele": "HLA-A*02:06",
    }
    out = annotate_variant_peptide_row(
        row,
        ["HLA-A*02:06"],
        prime_pair_index={
            ("MTPEPTIDE", "HLA-A*02:06"): {"prime_score": "0.9", "prime_rank": "1.0"},
            ("WTPEPTIDE", "HLA-A*02:06"): {"prime_score": "0.2", "prime_rank": "50.0"},
        },
        bigmhc_im_index={
            ("MTPEPTIDE", "HLA-A*02:06"): {"bigmhc_im_score": "0.8"},
            ("WTPEPTIDE", "HLA-A*02:06"): {"bigmhc_im_score": "0.1"},
        },
    )
    assert out["prime_score"] == "0.9"
    assert out["prime_rank"] == "1.0"
    assert out["prime_wt_score"] == "0.2"
    assert out["prime_wt_rank"] == "50.0"
    assert out["bigmhc_im_score"] == "0.8"
    assert out["bigmhc_im_wt_score"] == "0.1"


def test_annotate_variant_peptide_row_sets_wt_tool_scores_na_for_fusion_without_wt():
    row = {
        "mutant_peptide": "FUSIYPEP",
        "wildtype_peptide": "",
        "hla_allele": "HLA-A*02:06",
        "multi_aa_flag": "fusion_neo",
    }
    out = annotate_variant_peptide_row(row, ["HLA-A*02:06"])
    assert out["peptide_source"] == "easyfuse"
    assert out["prime_wt_score"] == "NA"
    assert out["prime_wt_rank"] == "NA"
    assert out["bigmhc_im_wt_score"] == "NA"
