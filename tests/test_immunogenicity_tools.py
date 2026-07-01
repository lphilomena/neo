from neoag_v03.adapters.bigmhc_im import parse_bigmhc_im
from neoag_v03.adapters.prime import extract_prime_pair, prime_allele_tag, read_prime_wide_rows


def test_prime_allele_tag():
    assert prime_allele_tag("HLA-A*02:01") == "A0201"
    assert prime_allele_tag("HLA-B*45:01") == "B4501"


def test_extract_prime_pair_from_wide_output():
    row = {
        "Peptide": "YMDGTMSQV",
        "Score_A0201": "0.045642",
        "%Rank_A0201": "1.905",
        "Score_bestAllele": "0.999",
    }
    score, rank = extract_prime_pair(row, "YMDGTMSQV", "HLA-A*02:01")
    assert score == "0.045642"
    assert rank == "1.905"


def test_read_prime_wide_rows_skips_comment_header(tmp_path):
    raw = tmp_path / "prime_out.tsv"
    raw.write_text(
        "####################\n"
        "# Output from PRIME (v2.1)\n"
        "Peptide\t%Rank_bestAllele\tScore_bestAllele\tBestAllele\t%Rank_A0201\tScore_A0201\n"
        "YMDGTMSQV\t0.033\t0.156524\tC0501\t0.092\t0.125796\n",
        encoding="utf-8",
    )
    rows = read_prime_wide_rows(raw)
    assert len(rows) == 1
    score, rank = extract_prime_pair(rows[0], "YMDGTMSQV", "HLA-A*02:01")
    assert score == "0.125796"
    assert rank == "0.092"


def test_parse_bigmhc_im_uses_bigmhc_im_column(tmp_path):
    path = tmp_path / "input.csv.prd"
    path.write_text(
        "pep,mhc,BigMHC_IM\nSIINFEKL,HLA-A*02:01,0.8123\n",
        encoding="utf-8",
    )
    rows = parse_bigmhc_im(path, "HCC1395")
    assert rows
    assert rows[0]["bigmhc_im_score"] == "0.8123"
