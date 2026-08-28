from neoag.mutant_specificity import evaluate_mutant_specificity


PROFILE = {
    "mutant_specificity": {
        "near_equal_el_rank_difference": 0.01,
        "positive_agretopicity_ratio": 2.0,
        "positive_el_rank_difference": 0.10,
        "caution_priority_cap": "C_CAUTION",
    }
}


def evaluate(mt, wt, mt_rank, wt_rank):
    return evaluate_mutant_specificity(
        {"peptide": mt, "wildtype_peptide": wt},
        {"netmhcpan_mt_rank_el": str(mt_rank), "netmhcpan_wt_rank_el": str(wt_rank)},
        PROFILE,
    )


def test_wildtype_better_is_capped():
    result = evaluate("PADFVAGSL", "RADFVAGSL", 1.066, 0.055)
    assert result["mutant_specificity_status"] == "WT_BETTER"
    assert result["mutant_specificity_priority_cap"] == "C_CAUTION"


def test_two_very_strong_nearly_identical_ranks_are_similar():
    result = evaluate("FRKPKVTEI", "FRKPKVIEI", 0.002, 0.001)
    assert result["mutant_specificity_status"] == "MT_WT_SIMILAR"
    assert result["mutant_specificity_gate_status"] == "CAUTION"


def test_clear_mutant_advantage_passes():
    result = evaluate("ABCDEFGHI", "ABXDEFGHI", 0.2, 3.0)
    assert result["mutant_specificity_status"] == "MT_SPECIFIC"
    assert result["mutant_specificity_gate_status"] == "PASS"


def test_missing_wildtype_is_unassessed_without_penalty():
    result = evaluate_mutant_specificity({"peptide": "ABCDEFGHI", "wildtype_peptide": ""}, {}, PROFILE)
    assert result["mutant_specificity_status"] == "UNASSESSED"
    assert result["mutant_specificity_multiplier"] == "1.0000"
    assert result["mutant_specificity_priority_cap"] == ""


def test_anchor_and_tcr_facing_positions_are_reported():
    anchor = evaluate("ABCDEFGHI", "AXCDEFGHI", 0.2, 1.0)
    assert anchor["mutation_anchor_only"] == "yes"
    assert anchor["mutation_tcr_facing"] == "no"
    exposed = evaluate("ABCDEFGHI", "ABCXEFGHI", 0.2, 1.0)
    assert exposed["mutation_anchor_only"] == "no"
    assert exposed["mutation_tcr_facing"] == "yes"
    assert anchor["mutation_position_role"] == "PRIMARY_HLA_ANCHOR"
    assert exposed["mutation_position_role"] == "PUTATIVE_TCR_FACING"


def test_small_mt_advantage_with_strong_wt_binding_is_explicit_safety_caution():
    result = evaluate_mutant_specificity(
        {"peptide": "ABCXEFGHI", "wildtype_peptide": "ABCDEFGHI", "mhc_class": "I"},
        {
            "netmhcpan_mt_rank_el": "0.40",
            "netmhcpan_wt_rank_el": "0.70",
            "netmhcpan_wt_rank_ba": "0.90",
            "netmhcpan_wt_ic50": "42",
        },
        PROFILE,
    )
    assert result["mutant_specificity_status"] == "MARGINAL_MT_ADVANTAGE"
    assert result["mutant_specificity_gate_status"] == "CAUTION"
    assert result["wt_self_reactivity_risk_status"] == "WT_STRONG_BINDING_REVIEW"
    assert result["mutation_position_role"] == "PUTATIVE_TCR_FACING"
    assert "not independent evidence" in result["mt_wt_interpretation_caution"]


def test_mhc_ii_position_role_stays_unresolved_without_binding_register():
    result = evaluate_mutant_specificity(
        {"peptide": "ABCDEFGHIJKLMNO", "wildtype_peptide": "ABCXEFGHIJKLMNO", "mhc_class": "II"},
        {"netmhcpan_mt_rank_el": "0.2", "netmhcpan_wt_rank_el": "3.0"},
        PROFILE,
    )
    assert result["mutation_position_role"] == "STRUCTURAL_ROLE_UNCERTAIN"
    assert result["mutation_anchor_only"] == "unknown"


def test_mhc_ii_is_inferred_from_hla_d_allele_when_class_is_missing():
    result = evaluate_mutant_specificity(
        {"peptide": "ABCDEFGHIJKLMNO", "wildtype_peptide": "ABCXEFGHIJKLMNO", "hla_allele": "HLA-DRB1*04:01"},
        {"netmhcpan_mt_rank_el": "0.2", "netmhcpan_wt_rank_el": "3.0"},
        PROFILE,
    )
    assert result["mutation_position_role"] == "STRUCTURAL_ROLE_UNCERTAIN"
