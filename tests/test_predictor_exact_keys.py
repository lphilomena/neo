from neoag.adapters.bigmhc_im import bigmhc_by_pair
from neoag.adapters.prime import prime_by_pair
from neoag.presentation import by_key
from neoag.utils import safe_id


def test_presentation_exact_pair_ignores_stale_caller_key_and_normalizes_hla():
    row = {
        "peptide_hla_key": "STALE_CALLER_ID",
        "peptide": "aaaaaaaaa",
        "hla_allele": "A0201",
        "netmhcpan_el_rank": "0.4",
    }
    indexed = by_key([row])
    assert indexed[safe_id("AAAAAAAAA_HLA-A*02:01")] is row
    assert "STALE_CALLER_ID" not in indexed


def test_immunogenicity_pair_indexes_normalize_hla_notation():
    prime = prime_by_pair([{"peptide": "AAAAAAAAA", "hla_allele": "A0201", "prime_score": "0.7"}])
    bigmhc = bigmhc_by_pair([{"peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "bigmhc_im_score": "0.8"}])
    assert prime[("AAAAAAAAA", "HLA-A*02:01")]["prime_score"] == "0.7"
    assert bigmhc[("AAAAAAAAA", "HLA-A*02:01")]["bigmhc_im_score"] == "0.8"
