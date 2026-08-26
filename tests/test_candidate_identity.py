from neoag.candidate_identity import candidate_identity, normalize_hla
from neoag.reports_dual import _patient_representatives


def test_four_level_identity_keeps_events_but_collapses_same_peptide_hla():
    first = {
        "event_id": "FUSION_EWSR1_WT1_BP1",
        "protein_change": "ORF1",
        "peptide": "SSYGQQSEK",
        "hla_allele": "HLA-A*11:02",
        "event_type": "Fusion",
    }
    second = {
        "event_id": "FUSION_EWSR1_WT1_BP2",
        "protein_change": "ORF2",
        "peptide": "ssygqqsek",
        "hla_allele": "A1102",
        "event_type": "Fusion",
    }
    first_ids = candidate_identity(first)
    second_ids = candidate_identity(second)
    assert first_ids["event_identity_id"] != second_ids["event_identity_id"]
    assert first_ids["protein_change_identity_id"] != second_ids["protein_change_identity_id"]
    assert first_ids["peptide_sequence_id"] == second_ids["peptide_sequence_id"]
    assert first_ids["peptide_hla_id"] == second_ids["peptide_hla_id"]
    assert len(_patient_representatives([first, second], 10)) == 1


def test_hla_normalization_is_stable():
    assert normalize_hla("A1102") == "HLA-A*11:02"
    assert normalize_hla("HLA-A*11:02") == "HLA-A*11:02"
