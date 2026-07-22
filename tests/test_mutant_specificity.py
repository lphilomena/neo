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
    result = evaluate("ABCDEFGHI", "ABXDEFGHI", 0.2, 1.0)
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
