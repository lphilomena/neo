from pathlib import Path

from neoag_v03.config import load_profile
from neoag_v03.scoring_v03 import resolve_appm_peptide_modifiers_tsv, score_v03
from neoag_v03.utils import write_tsv


def test_resolve_appm_peptide_modifiers_from_summary_sibling(tmp_path):
    appm = tmp_path / "appm"
    appm.mkdir()
    summary = appm / "appm_summary.tsv"
    modifiers = appm / "appm_peptide_modifiers.tsv"
    summary.write_text("sample_id\nS1\n", encoding="utf-8")
    modifiers.write_text("peptide_id\tappm_multiplier\nP1\t0.5000\n", encoding="utf-8")
    resolved = resolve_appm_peptide_modifiers_tsv(appm_summary_tsv=summary)
    assert resolved == str(modifiers)


def test_score_v03_applies_appm_peptide_modifiers(tmp_path):
    events = tmp_path / "events.tsv"
    peptides = tmp_path / "peptides.tsv"
    presentation = tmp_path / "presentation.tsv"
    appm = tmp_path / "appm"
    appm.mkdir()
    write_tsv(events, [{
        "event_id": "E1", "sample_id": "S1", "event_type": "SNV", "gene": "B2M",
        "event_confidence": "0.7", "event_expression": "10", "driver_relevance": "0.5",
        "clonality": "0.5", "persistence": "0.5", "tumor_specificity": "0.7",
        "safety_status": "PASS",
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1", "event_type": "SNV",
        "gene": "B2M", "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "mhc_class": "I",
        "binding_rank": "1.0", "el_rank": "1.0", "presentation_score": "0.8",
        "immunogenicity_score": "0.6", "wildtype_binding_rank": "5.0",
        "self_similarity_score": "0.1", "normal_hla_ligand_overlap": "no",
    }])
    write_tsv(presentation, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1", "peptide": "AAAAAAAAA",
        "hla_allele": "HLA-A*02:01", "mhc_class": "I",
        "binding_evidence_score": "0.8", "presentation_evidence_score": "0.7",
        "presentation_evidence_grade": "A",
    }])
    write_tsv(appm / "appm_summary.tsv", [{
        "mhc_i_integrity_score": "1.0", "mhc_ii_integrity_score": "1.0",
        "mhc_i_integrity_status": "MHC_I_INTACT", "mhc_ii_integrity_status": "MHC_II_INTACT",
    }])
    write_tsv(appm / "appm_peptide_modifiers.tsv", [{
        "peptide_id": "P1", "event_id": "E1", "hla_allele": "HLA-A*02:01", "mhc_class": "I",
        "appm_multiplier": "0.2500", "appm_multiplier_reason": "mhc_i_defective",
        "appm_integrity_status": "MHC_I_DEFECTIVE", "appm_evidence_completeness": "LOW",
        "appm_review_required": "yes", "priority_cap": "C", "appm_action": "CAP",
    }])
    out_events = tmp_path / "ranked_events.tsv"
    out_peptides = tmp_path / "ranked_peptides.tsv"
    _, peps = score_v03(
        events, peptides, presentation, appm / "appm_summary.tsv", None, None, None,
        load_profile("default"), out_events, out_peptides,
        appm_peptide_modifiers_tsv=appm / "appm_peptide_modifiers.tsv",
    )
    assert peps[0]["appm_multiplier"] == "0.2500"
    assert peps[0]["appm_multiplier_reason"] == "mhc_i_defective"
    assert peps[0]["appm_integrity_status"] == "MHC_I_DEFECTIVE"
    assert peps[0]["priority_cap"] == "C"
