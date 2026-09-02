from pathlib import Path

from neoag.cohort_rules import load_cohort_rule_contract
from neoag.open_neo.review import _read_context
from neoag.report_from_final import enrich_report_provenance
from neoag.reports_dual import load_report_bundle

ROOT = Path(__file__).parents[1]

def test_literal_disease_profile_is_structured_context() -> None:
    assert _read_context("DSRCT", None)["disease"] == "DSRCT"


def test_path_like_missing_disease_profile_is_not_a_diagnosis() -> None:
    assert _read_context("missing/dsrct.yaml", None) == {}


def test_analysis_profile_is_not_promoted_to_clinical_disease(tmp_path: Path) -> None:
    final_dir = tmp_path / "production/final"
    final_dir.mkdir(parents=True)
    provenance = enrich_report_provenance(
        final_dir,
        {},
        manifest={"run": {"profile": "/profiles/sarcoma_rna_supported_v2_provisional.toml"}},
        generated={},
    )
    assert "disease" not in provenance
    assert provenance["analysis_profile"] == "sarcoma_rna_supported_v2_provisional"


def test_dsrct_context_auto_loads_molecular_knowledge() -> None:
    bundle = load_report_bundle(
        profile={"_profile_name": "sarcoma_evidence_consensus_v3_source_chain"},
        events=[],
        peptides=[],
        provenance={"disease": "DSRCT"},
    )
    assert bundle.disease_knowledge["status"] == "LOADED"
    assert bundle.disease_knowledge["disease_id"] == "DSRCT"


def test_contract_selects_knowledge_without_promoting_a_clinical_diagnosis(tmp_path: Path) -> None:
    contract = load_cohort_rule_contract(ROOT / "configs/cohorts/dsrct_v1.toml")
    final_dir = tmp_path / "production/final"
    final_dir.mkdir(parents=True)
    provenance = enrich_report_provenance(
        final_dir,
        {},
        manifest={
            "run": {
                "profile": contract["ranking_profile"],
                "cohort_rule_set": contract["path"],
            },
            "evidence": {"evidence_consensus_rules": contract["evidence_consensus_rules"]},
        },
        generated={},
    )
    assert "disease" not in provenance
    assert provenance["molecular_knowledge_selection"]["source"] == "cohort_rule_contract"
    assert provenance["disease_knowledge_file"] == contract["disease_knowledge_file"]
    bundle = load_report_bundle(profile={}, events=[], peptides=[], provenance=provenance)
    assert bundle.disease_knowledge["status"] == "LOADED"


def test_manifest_clinical_context_is_preserved_separately(tmp_path: Path) -> None:
    final_dir = tmp_path / "production/final"
    final_dir.mkdir(parents=True)
    provenance = enrich_report_provenance(
        final_dir,
        {},
        manifest={"run": {"clinical_context": {"disease": "DSRCT"}}},
        generated={},
    )
    assert provenance["clinical_context"] == {"disease": "DSRCT"}
    assert provenance["cohort_rule_contract"]["comparability_status"] == "UNASSESSED_NO_MATCHING_CONTRACT"
