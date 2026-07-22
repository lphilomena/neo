from neoag.config import load_profile
from neoag.gates import evaluate_presentation_gate
from neoag.safety import apply_event_safety, load_normal_expression, safety_multiplier
from neoag.utils import write_tsv


def test_v2_treats_normal_expression_and_hpa_ncpm_as_caution(tmp_path):
    normal = tmp_path / "normal.tsv"
    write_tsv(normal, [{
        "gene": "GENE1",
        "normal_tissue_max_tpm": "50",
        "critical_tissue_max_tpm": "20",
        "critical_tissue_hit": "yes",
        "normal_hspc_tpm": "250",
        "normal_hspc_unit": "HPA_nCPM",
        "normal_expression_status": "ASSESSED",
        "normal_hspc_status": "ASSESSED",
    }])
    profile = load_profile("sarcoma_rna_supported_v2_provisional")
    event = apply_event_safety({"gene": "GENE1"}, profile, load_normal_expression(normal))
    assert event["safety_status"] == "CAUTION"
    assert event["normal_hspc_unit"] == "HPA_nCPM"
    assert "critical_tissue_expression" in event["safety_reason"]
    assert "normal_HSPC_expression" in event["safety_reason"]


def test_v2_caution_multiplier_is_configurable_without_changing_default():
    profile = load_profile("sarcoma_rna_supported_v2_provisional")
    assert safety_multiplier("CAUTION") == 0.45
    assert safety_multiplier("CAUTION", profile) == 0.85


def test_v2_allele_expression_gate_uses_gene_tpm_times_rna_vaf():
    profile = load_profile("sarcoma_rna_supported_v2_provisional")
    result = evaluate_presentation_gate(
        {},
        {
            "event_expression": "10",
            "gene_expression_tpm": "10",
            "rna_alt_reads": "3",
            "rna_vaf": "0.02",
            "peptide_consequence": "missense",
            "mutation_source": "SNV",
        },
        {"presentation_evidence_grade": "A", "netmhcpan_el_rank": "0.2"},
        profile,
    )
    assert result["presentation_gate_status"] == "FAIL"
    assert "allele_expression=0.2000" in result["presentation_gate_reason"]


def test_v2_fusion_caps_are_explicit():
    profile = load_profile("sarcoma_rna_supported_v2_provisional")
    assert profile["fusion"]["rna_only_priority_cap"] == "C_CAUTION"
    assert profile["fusion"]["normal_junction_unassessed_priority_cap"] == "C_CAUTION"
