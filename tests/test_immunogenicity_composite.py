from neoag.config import load_profile
from neoag.immunogenicity_composite import (
    apply_immunogenicity_evidence,
    combine_component_scores,
    component_scores,
    effective_source_weights,
    resolve_immunogenicity_score,
)


def test_composite_weighted_prime_and_bigmhc():
    profile = load_profile("default")
    row = {
        "prime_score": "0.40000",
        "bigmhc_im_score": "0.80000",
        "hla_allele": "HLA-A*02:01",
    }
    parts = component_scores(row, profile)
    assert "prime" in parts
    assert "bigmhc_im" in parts
    score, meta = resolve_immunogenicity_score({}, row, profile)
    assert score == combine_component_scores(parts, profile, hla_allele="HLA-A*02:01")
    assert meta["immunogenicity_source"] == "bigmhc_im+prime"
    assert parts["prime"] < score < parts["bigmhc_im"]


def test_weighted_composite_ignores_hla_multipliers():
    profile = load_profile("default")
    parts = {"prime": 0.8, "bigmhc_im": 0.2}
    a02 = combine_component_scores(parts, profile, hla_allele="HLA-A*02:01")
    b58 = combine_component_scores(parts, profile, hla_allele="HLA-B*58:01")
    assert a02 == b58
    weights_a02 = effective_source_weights("HLA-A*02:01", parts, profile)
    weights_b58 = effective_source_weights("HLA-B*58:01", parts, profile)
    assert weights_a02 == weights_b58


def test_extended_weighted_composite_ignores_hla_multipliers():
    profile = load_profile("immunogenicity_extended")
    parts = {"prime": 0.5, "bigmhc_im": 0.5, "deepimmuno": 0.9}
    a02 = combine_component_scores(parts, profile, hla_allele="HLA-A*02:01")
    a11 = combine_component_scores(parts, profile, hla_allele="HLA-A*11:01")
    assert a02 == a11
    weights = effective_source_weights("HLA-A*11:01", parts, profile)
    assert weights["deepimmuno"] > 0.0


def test_iedb_in_default_profile_sources():
    profile = load_profile("default")
    assert "iedb" in profile["immunogenicity"]["sources"]
    row = {
        "iedb_immunogenicity_score": "0.30000",
        "hla_allele": "HLA-A*02:06",
    }
    parts = component_scores(row, profile)
    assert "iedb" in parts


def test_stub_evidence_rows_do_not_enter_composite(tmp_path):
    profile = load_profile("default")
    prime = tmp_path / "prime_evidence.tsv"
    prime.write_text(
        "sample_id\tpeptide\thla_allele\tprime_score\tprime_rank\tevidence_status\n"
        "S\tPEPTIDEAA\tHLA-A*02:01\t0.99\t0.01\tstub\n",
        encoding="utf-8",
    )
    rows = [{"peptide": "PEPTIDEAA", "hla_allele": "HLA-A*02:01"}]
    apply_immunogenicity_evidence(rows, {"prime": prime, "bigmhc_im": None, "deepimmuno": None, "iedb": None}, profile)
    assert rows[0]["prime_score"] == ""
    assert rows[0]["prime_rank"] == ""
    assert rows[0]["immunogenicity_composite_score"] == ""
    assert rows[0]["immunogenicity_source"] == "none"


def test_mixed_real_and_stub_evidence_uses_only_real_component(tmp_path):
    profile = load_profile("default")
    prime = tmp_path / "prime_evidence.tsv"
    prime.write_text(
        "sample_id\tpeptide\thla_allele\tprime_score\tprime_rank\tevidence_status\n"
        "S\tPEPTIDEAA\tHLA-A*02:01\t0.99\t0.01\tstub\n",
        encoding="utf-8",
    )
    bigmhc = tmp_path / "bigmhc_im_evidence.tsv"
    bigmhc.write_text(
        "sample_id\tpeptide\thla_allele\tbigmhc_im_score\tevidence_status\n"
        "S\tPEPTIDEAA\tHLA-A*02:01\t0.42\treal\n",
        encoding="utf-8",
    )
    rows = [{"peptide": "PEPTIDEAA", "hla_allele": "HLA-A*02:01"}]
    apply_immunogenicity_evidence(rows, {"prime": prime, "bigmhc_im": bigmhc, "deepimmuno": None, "iedb": None}, profile)
    assert rows[0]["prime_score"] == ""
    assert rows[0]["bigmhc_im_score"] == "0.42"
    assert rows[0]["immunogenicity_source"] == "bigmhc_im"
    assert rows[0]["immunogenicity_composite_score"] == "0.4200"
