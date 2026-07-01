from neoag_v03.adapters.deepimmuno import (
    is_valid_peptide,
    normalize_hla_for_deepimmuno,
    normalize_hla_allele,
    parse_deepimmuno,
    predict_pairs,
    predict_pair_stub,
)
from neoag_v03.adapters.peptide_input import convert_peptide_input
from neoag_v03.config import load_profile
from neoag_v03.immunogenicity_composite import component_scores, resolve_immunogenicity_score
from neoag_v03.tools.runner import run_deepimmuno
from neoag_v03.tools.registry import RunContext


def test_normalize_hla_for_deepimmuno():
    assert normalize_hla_for_deepimmuno("HLA-A*02:01") == "HLA-A*0201"
    assert normalize_hla_for_deepimmuno("A*02:01") == "HLA-A*0201"


def test_deepimmuno_valid_peptide_lengths():
    assert is_valid_peptide("YMDGTMSQV")
    assert is_valid_peptide("YMDGTMSQVA")
    assert not is_valid_peptide("SIINFEKL")
    assert not is_valid_peptide("YMDGTMSQVAT")


def test_predict_pairs_stub():
    rows = predict_pairs([("YMDGTMSQV", "HLA-A*02:01"), ("YMDGTMSQV", "HLA-B*07:02")])
    assert len(rows) == 2
    assert rows[0]["deepimmuno_score"] != ""
    assert rows[0]["hla_allele"] == normalize_hla_allele("HLA-A*02:01")


def test_deepimmuno_in_composite_profile():
    profile = load_profile("immunogenicity_extended")
    row = {
        **predict_pairs([("YMDGTMSQV", "HLA-A*02:01")])[0],
        "prime_score": "0.1",
        "prime_rank": "1.0",
        "bigmhc_im_score": "0.2",
    }
    parts = component_scores(row, profile)
    assert "deepimmuno" in parts
    score, meta = resolve_immunogenicity_score({}, row, profile)
    assert score > 0
    assert "deepimmuno" in meta["immunogenicity_source"]


def test_parse_deepimmuno_evidence_tsv(tmp_path):
    path = tmp_path / "deepimmuno.tsv"
    path.write_text(
        "sample_id\tpeptide\thla_allele\tdeepimmuno_score\tsource_file\n"
        "S1\tYMDGTMSQV\tHLA-A*02:01\t0.48081678\t/tools/deepimmuno.tsv\n",
        encoding="utf-8",
    )
    rows = parse_deepimmuno(path, sample_id="S1")
    assert len(rows) == 1
    assert rows[0]["deepimmuno_score"] == "0.48081678"


def test_run_deepimmuno_stub(tmp_path):
    inp = tmp_path / "pairs.tsv"
    inp.write_text("peptide\thla_allele\nYMDGTMSQV\tHLA-A*02:01\n", encoding="utf-8")
    summary = convert_peptide_input(inp, tmp_path / "work", sample_id="S1")
    ctx = RunContext(
        sample_id="S1",
        outdir=tmp_path,
        stub=True,
        raw_peptides=summary.raw_peptides_tsv,
    )
    out = tmp_path / "deepimmuno.tsv"
    run_deepimmuno(ctx, out)
    text = out.read_text(encoding="utf-8")
    assert "deepimmuno_score" in text
    assert "YMDGTMSQV" in text
