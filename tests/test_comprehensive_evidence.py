from pathlib import Path

from neoag.comprehensive_evidence import (
    AUTHORITATIVE_FIELDS,
    EVIDENCE_SOURCE_PRECEDENCE_VERSION,
    build_comprehensive_peptide_evidence,
)
from neoag.utils import read_tsv, write_tsv


def test_comprehensive_evidence_preserves_annotation_and_adds_all_layers(tmp_path):
    annotated = tmp_path / "variant_peptides.annotated.tsv"
    ranked = tmp_path / "ranked.tsv"
    presentation = tmp_path / "presentation.tsv"
    appm = tmp_path / "appm.tsv"
    ccf = tmp_path / "ccf.tsv"
    safety = tmp_path / "safety.tsv"
    escape = tmp_path / "escape.tsv"
    validation = tmp_path / "validation.tsv"
    raw_events = tmp_path / "events.tsv"
    expression = tmp_path / "expression.tsv"
    rna = tmp_path / "rna.tsv"
    output = tmp_path / "comprehensive.tsv"

    write_tsv(annotated, [{
        "peptide_id": "P1",
        "gene": "GENE1",
        "transcript_id": "ENST1",
        "hgvsp": "p.A1V",
        "mutant_peptide": "AAAAAAAA",
        "hla_allele": "HLA-A*02:01",
        "netmhcpan_mt_rank_el": "0.2",
        "mhcflurry_mt_presentation_score": "0.8",
        "prime_score": "0.7",
        "bigmhc_im_score": "0.6",
    }])
    write_tsv(ranked, [{
        "peptide_id": "P1",
        "event_id": "E1",
        "gene": "GENE1",
        "peptide": "AAAAAAAA",
        "hla_allele": "HLA-A*02:01",
        "efficacy_score": "0.91",
        "final_priority": "A",
        "recommended_use": "short peptide",
    }])
    write_tsv(presentation, [{
        "peptide_id": "P1", "presentation_evidence_grade": "A",
        "presentation_evidence_score": "0.94",
    }])
    write_tsv(appm, [{
        "peptide_id": "P1", "appm_multiplier": "0.8", "appm_action": "retain",
        "appm_call_confidence": "high",
    }])
    write_tsv(ccf, [{
        "event_id": "E1", "ccf_best": "0.85", "clonality_status": "clonal",
        "purity": "0.72", "total_cn": "2",
    }])
    write_tsv(safety, [{
        "peptide_id": "P1", "safety_tier": "PASS", "reference_proteome_exact_match": "false",
    }])
    write_tsv(escape, [{
        "peptide_id": "P1", "restricting_hla_lost": "false", "escape_status": "PASS",
    }])
    write_tsv(validation, [{
        "peptide_id": "P1", "validation_mode": "short_peptide", "recommended_assay": "ELISpot",
    }])
    write_tsv(raw_events, [{
        "event_id": "E1", "tumor_vaf": "0.21", "tumor_depth": "100", "rna_vaf": "0.18",
    }])
    write_tsv(expression, [{
        "event_id": "E1", "expression_tpm": "12.5", "expression_source": "Salmon",
    }])
    write_tsv(rna, [{
        "event_id": "E1", "junction_reads": "22", "rna_support_status": "supported",
    }])

    summary = build_comprehensive_peptide_evidence(
        output_tsv=output,
        annotated_peptides=annotated,
        ranked_peptides=ranked,
        raw_events=raw_events,
        presentation_evidence=presentation,
        appm_peptide_modifiers=appm,
        ccf_2=ccf,
        expression_evidence=expression,
        rna_junction_evidence=rna,
        peptide_safety=safety,
        peptide_escape_flags=escape,
        validation_plan=validation,
    )

    row = read_tsv(output)[0]
    assert summary["base_source"] == "annotated_peptides"
    assert row["transcript_id"] == "ENST1"
    assert row["netmhcpan_mt_rank_el"] == "0.2"
    assert row["final_priority"] == "A"
    assert row["presentation_evidence_grade"] == "A"
    assert row["appm_call_confidence"] == "high"
    assert row["ccf_best"] == "0.85"
    assert row["expression_tpm"] == "12.5"
    assert row["junction_reads"] == "22"
    assert row["safety_tier"] == "PASS"
    assert row["escape_status"] == "PASS"
    assert row["validation_mode"] == "short_peptide"
    assert row["comprehensive_evidence_status"] == "COMPLETE"


def test_comprehensive_evidence_falls_back_to_ranked_table(tmp_path):
    ranked = tmp_path / "ranked.tsv"
    output = tmp_path / "comprehensive.tsv"
    write_tsv(ranked, [{
        "peptide_id": "F1", "event_id": "FE1", "event_type": "Fusion",
        "peptide": "FUSIONPEP", "hla_allele": "HLA-A*02:01", "final_priority": "B",
    }])

    summary = build_comprehensive_peptide_evidence(
        output_tsv=output,
        ranked_peptides=ranked,
    )

    row = read_tsv(output)[0]
    assert summary["base_source"] == "ranked_peptides"
    assert row["event_type"] == "Fusion"
    assert row["comprehensive_evidence_status"] == "PARTIAL"


def test_authoritative_sources_override_stale_ranked_copies_and_record_conflicts(tmp_path):
    ranked = tmp_path / "ranked.tsv"
    raw_events = tmp_path / "events.tsv"
    presentation = tmp_path / "presentation.tsv"
    expression = tmp_path / "expression.tsv"
    ccf = tmp_path / "ccf.tsv"
    safety = tmp_path / "safety.tsv"
    event_safety = tmp_path / "event_safety.tsv"
    output = tmp_path / "comprehensive.tsv"
    write_tsv(ranked, [{
        "peptide_id": "P1", "event_id": "E1", "gene": "STALE_GENE",
        "presentation_evidence_score": "0.1", "gene_expression_tpm": "0.1",
        "ccf_estimate": "0.2", "safety_status": "CAUTION",
        "efficacy_score": "0.91", "final_priority": "B",
    }])
    write_tsv(raw_events, [{"event_id": "E1", "gene": "AUTH_GENE", "event_type": "SNV"}])
    write_tsv(presentation, [{"peptide_id": "P1", "presentation_evidence_score": "0.94"}])
    write_tsv(expression, [{"event_id": "E1", "gene_expression_tpm": "12.5"}])
    write_tsv(ccf, [{"event_id": "E1", "ccf_estimate": "0.85"}])
    write_tsv(safety, [{"peptide_id": "P1", "safety_status": "PASS"}])
    write_tsv(event_safety, [{"event_id": "E1", "safety_status": "CAUTION"}])

    summary = build_comprehensive_peptide_evidence(
        output_tsv=output,
        ranked_peptides=ranked,
        raw_events=raw_events,
        presentation_evidence=presentation,
        expression_evidence=expression,
        ccf_2=ccf,
        peptide_safety=safety,
        event_safety=event_safety,
    )
    row = read_tsv(output)[0]
    assert row["gene"] == "AUTH_GENE"
    assert row["presentation_evidence_score"] == "0.94"
    assert row["gene_expression_tpm"] == "12.5"
    assert row["ccf_estimate"] == "0.85"
    assert row["safety_status"] == "PASS"
    assert row["event_safety_status"] == "CAUTION"
    assert row["efficacy_score"] == "0.91"
    assert row["evidence_source_precedence_version"] == EVIDENCE_SOURCE_PRECEDENCE_VERSION
    assert {"gene", "presentation_evidence_score", "gene_expression_tpm", "ccf_estimate", "safety_status"} <= set(row["evidence_conflict_fields"].split(","))
    conflicts = read_tsv(summary["conflicts_tsv"])
    assert summary["conflicts"] >= 5
    assert {item["selected_source"] for item in conflicts} >= {
        "raw_events", "presentation_evidence", "expression_evidence", "ccf_2", "peptide_safety",
    }
    assert "presentation_evidence_score" in AUTHORITATIVE_FIELDS["presentation_evidence"]
