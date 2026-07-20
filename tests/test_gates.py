from neoag.config import load_profile
from neoag.gates import evaluate_presentation_gate
from neoag.scoring import (
    compute_peptide_efficacy,
    effective_weights,
    is_placeholder_immunogenicity,
)


def test_presentation_gate_pass():
    profile = load_profile("default")
    pres = {
        "presentation_evidence_grade": "A",
        "netmhcpan_el_rank": "0.5",
        "netmhcstabpan_score": "8.0",
    }
    gate = evaluate_presentation_gate({}, {}, pres, profile)
    assert gate["presentation_gate_status"] == "PASS"
    assert gate["presentation_gate_multiplier"] == "1.0000"


def test_presentation_gate_fail_grade():
    profile = load_profile("default")
    pres = {
        "presentation_evidence_grade": "C",
        "netmhcpan_el_rank": "0.5",
        "netmhcstabpan_score": "8.0",
    }
    gate = evaluate_presentation_gate({}, {}, pres, profile)
    assert gate["presentation_gate_status"] == "FAIL"
    assert gate["presentation_gate_multiplier"] == "0.2500"
    assert "grade" in gate["presentation_gate_reason"]


def test_immunogenicity_placeholder_redistribution():
    profile = load_profile("default")
    peptide = {"immunogenicity_score": "0.5"}
    w = effective_weights(profile, peptide)
    assert w["immunogenicity"] == 0.0
    assert w["presentation_evidence"] > profile["score_weights"]["presentation_evidence"]
    assert w["binding_evidence"] > profile["score_weights"]["binding_evidence"]


def test_efficacy_uses_immunogenicity_composite():
    profile = load_profile("default")
    peptide = {"immunogenicity_score": "0.5", "safety_status": "PASS"}
    event = {"event_score": "0.6000", "safety_status": "PASS"}
    pres = {
        "binding_evidence_score": "0.8000",
        "presentation_evidence_score": "0.7000",
        "evidence_completeness": "1.0000",
        "prime_score": "0.8200",
        "prime_rank": "1.5",
        "bigmhc_im_score": "0.7100",
        "presentation_evidence_grade": "A",
        "netmhcpan_el_rank": "0.5",
        "netmhcstabpan_score": "8.0",
    }
    out = compute_peptide_efficacy(peptide, event, pres, profile)
    assert out["immunogenicity_resolved"] == "yes"
    assert out["immunogenicity_source"] == "bigmhc_im+prime"
    assert float(out["efficacy_score"]) > 0.0
    assert is_placeholder_immunogenicity(peptide)
