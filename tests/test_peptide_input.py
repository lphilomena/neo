from pathlib import Path

import pytest

from neoag_v03.adapters.peptide_input import (
    convert_peptide_input,
    normalize_hla_allele,
    unique_peptide_hla_records,
    PeptideHlaRecord,
)
from neoag_v03.tools.prep import unique_peptide_hla_pairs
from neoag_v03.cli import main


def test_normalize_hla_formats():
    assert normalize_hla_allele("A*02:01") == "HLA-A*02:01"
    assert normalize_hla_allele("HLA-A0201") == "HLA-A*02:01"
    assert normalize_hla_allele("HLA-A*02:01") == "HLA-A*02:01"


def test_pair_dedup_keeps_same_peptide_different_hla():
    records = [
        PeptideHlaRecord(peptide="SIINFEKL", hla_allele="HLA-A*02:01"),
        PeptideHlaRecord(peptide="SIINFEKL", hla_allele="HLA-B*07:02"),
        PeptideHlaRecord(peptide="SIINFEKL", hla_allele="HLA-A*02:01"),
    ]
    pairs = unique_peptide_hla_records(records)
    assert len(pairs) == 2
    keys = {(p.peptide, p.hla_allele) for p in pairs}
    assert ("SIINFEKL", "HLA-A*02:01") in keys
    assert ("SIINFEKL", "HLA-B*07:02") in keys


def test_convert_flexible_csv(tmp_path):
    inp = tmp_path / "input.csv"
    inp.write_text(
        "seq,hla,gene\n"
        "SIINFEKL,HLA-A*02:01,GENE1\n"
        "SIINFEKL,A*11:01,GENE1\n"
        "SIINFEKL,HLA-A*02:01,GENE1\n",
        encoding="utf-8",
    )
    summary = convert_peptide_input(inp, tmp_path / "out", sample_id="S1")
    assert summary.input_rows == 3
    assert summary.pair_rows == 2
    assert summary.unique_peptides == 1
    assert summary.unique_hla == 2

    pairs = unique_peptide_hla_pairs(summary.raw_peptides_tsv)
    assert ("SIINFEKL", "HLA-A*02:01") in pairs
    assert ("SIINFEKL", "HLA-A*11:01") in pairs


def test_peptide_predict_stub_cli(tmp_path):
    inp = tmp_path / "pairs.tsv"
    inp.write_text(
        "peptide\thla_allele\n"
        "VVVGADGVGK\tHLA-A*11:01\n"
        "AAAAAAA\tHLA-A*02:01\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "predict"
    main(
        [
            "peptide-predict",
            "-i",
            str(inp),
            "-o",
            str(outdir),
            "--stub",
            "--skip-stabpan",
        ]
    )
    assert (outdir / "parsed/raw_peptides.tsv").exists()
    assert (outdir / "peptide_predictions.tsv").exists()
