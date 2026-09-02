from pathlib import Path

from neoag.cohort_rules import (
    discover_matching_cohort_contract,
    load_cohort_rule_contract,
    validate_cohort_rule_pair,
)


ROOT = Path(__file__).parents[1]


def test_dsrct_contract_locks_the_canonical_profile_and_rules() -> None:
    contract = load_cohort_rule_contract(ROOT / "configs/cohorts/dsrct_v1.toml")
    assert contract["id"] == "DSRCT-OPENNEO-2026.09"
    assert contract["report_contract_version"] == "open-neo-patient-v3"
    assert Path(contract["disease_knowledge_file"]).is_file()
    assert Path(contract["disease_knowledge_file"]).name == "dsrct.json"
    assert validate_cohort_rule_pair(
        contract,
        ranking_profile=ROOT / "profiles/sarcoma_rna_supported_v2_provisional.toml",
        evidence_consensus_rules=ROOT / "configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml",
    ) == []


def test_dsrct_contract_rejects_a_different_rules_file() -> None:
    contract = load_cohort_rule_contract(ROOT / "configs/cohorts/dsrct_v1.toml")
    mismatches = validate_cohort_rule_pair(
        contract,
        ranking_profile=ROOT / "profiles/sarcoma_rna_supported_v2_provisional.toml",
        evidence_consensus_rules=ROOT / "configs/ranking/sarcoma_evidence_consensus_v1.toml",
    )
    assert any("evidence_consensus_rules=" in item for item in mismatches)


def test_contract_can_be_discovered_for_legacy_canonical_pair() -> None:
    contract = discover_matching_cohort_contract(
        ROOT / "profiles/sarcoma_rna_supported_v2_provisional.toml",
        ROOT / "configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml",
    )
    assert contract is not None
    assert contract["id"] == "DSRCT-OPENNEO-2026.09"
