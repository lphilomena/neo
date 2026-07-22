from neoag.config import load_profile
from neoag.peptide_safety_gate import build_peptide_safety_gate
from neoag.safety import load_normal_expression
from neoag.scoring import score
from neoag.utils import write_tsv


def _inputs(tmp_path):
    events = tmp_path / "events.tsv"
    peptides = tmp_path / "peptides.tsv"
    write_tsv(events, [{
        "event_id": "E1", "sample_id": "S1", "event_type": "SNV",
        "mutation_source": "SNV", "gene": "KRAS", "event_confidence": "0.9",
        "event_expression": "10", "driver_relevance": "0.8", "clonality": "0.8",
        "persistence": "0.8", "tumor_specificity": "0.8",
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1",
        "event_type": "SNV", "mutation_source": "SNV", "gene": "KRAS",
        "peptide": "VVVGAGGVG", "wildtype_peptide": "VVVGAGGVG",
        "hla_allele": "HLA-A*02:01", "mhc_class": "I",
        "binding_rank": "0.2", "el_rank": "0.2", "presentation_score": "0.9",
        "immunogenicity_score": "0.7", "wildtype_binding_rank": "0.2",
        "self_similarity_score": "0.1",
    }])
    return events, peptides


def test_missing_safety_resources_are_partial_not_pass(tmp_path):
    events, peptides = _inputs(tmp_path)
    rows, event_rows = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
        out_event_safety=tmp_path / "event_safety.tsv",
        profile=load_profile("default"),
    )
    row = rows[0]
    assert row["safety_status"] == "SAFETY_PARTIAL"
    assert row["normal_tissue_max_tpm"] == ""
    assert row["reference_proteome_exact_match"] == "not_assessed"
    assert "normal_expression" in row["safety_missing_layers"]
    assert "reference_proteome" in row["safety_missing_layers"]
    assert event_rows[0]["event_safety_status"] == "SAFETY_PARTIAL"


def test_zero_expression_is_assessed_when_gene_is_present(tmp_path):
    events, peptides = _inputs(tmp_path)
    normal = tmp_path / "normal.tsv"
    write_tsv(normal, [{
        "gene": "KRAS", "normal_tissue_max_tpm": "0",
        "normal_hspc_tpm": "0", "critical_tissue_hit": "no",
    }])
    rows, _ = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
        normal_expression=normal,
    )
    assert rows[0]["normal_expression_status"] == "ASSESSED"
    assert rows[0]["normal_hspc_status"] == "ASSESSED"
    assert "normal_expression" not in rows[0]["safety_missing_layers"]


def test_event_caution_propagates_to_ranked_peptide(tmp_path):
    events, peptides = _inputs(tmp_path)
    presentation = tmp_path / "presentation.tsv"
    peptide_safety = tmp_path / "peptide_safety.tsv"
    event_safety = tmp_path / "event_safety.tsv"
    write_tsv(presentation, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1",
        "peptide": "VVVGAGGVG", "hla_allele": "HLA-A*02:01", "mhc_class": "I",
        "binding_evidence_score": "0.9", "presentation_evidence_score": "0.9",
        "presentation_evidence_grade": "A",
    }])
    write_tsv(peptide_safety, [{
        "peptide_id": "P1", "event_id": "E1", "safety_status": "PASS",
        "safety_reason": "no_major_signal", "safety_evidence_completeness": "1.0000",
    }])
    write_tsv(event_safety, [{
        "event_id": "E1", "event_safety_status": "CAUTION",
        "event_safety_reason": "critical_tissue_expression",
        "safety_evidence_completeness": "0.8000",
    }])
    _, ranked = score(
        events, peptides, presentation, None, None, None, None,
        load_profile("default"), tmp_path / "ranked_events.tsv", tmp_path / "ranked_peptides.tsv",
        peptide_safety_tsv=peptide_safety, event_safety_tsv=event_safety,
    )
    assert ranked[0]["safety_status"] == "CAUTION"
    assert "critical_tissue_expression" in ranked[0]["safety_reason"]


def test_peptide_local_proteome_match_does_not_fail_whole_event(tmp_path):
    events, peptides = _inputs(tmp_path)
    ref = tmp_path / "proteome.fa"
    ref.write_text(">normal\nXXVVVGAGGVGXX\n", encoding="utf-8")
    rows, event_rows = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
        out_event_safety=tmp_path / "event_safety.tsv",
        profile=load_profile("default"),
        reference_proteome=ref,
    )
    assert rows[0]["safety_status"] == "FAIL"
    assert event_rows[0]["event_safety_status"] == "SAFETY_PARTIAL"


def test_indel_crosses_junction_does_not_require_normal_rna_junctions(tmp_path):
    events, peptides = _inputs(tmp_path)
    event_rows = [{
        "event_id": "E1", "sample_id": "S1", "event_type": "InDel",
        "mutation_source": "InDel", "gene": "GENE1",
    }]
    peptide_rows = [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1",
        "event_type": "InDel", "mutation_source": "InDel", "gene": "GENE1",
        "peptide_consequence": "insertion", "crosses_junction": "yes",
        "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "mhc_class": "I",
    }]
    write_tsv(events, event_rows)
    write_tsv(peptides, peptide_rows)
    rows, _ = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
    )
    assert rows[0]["normal_junction_assessment_status"] == "NOT_APPLICABLE"
    assert "normal_junction" not in rows[0]["safety_missing_layers"]


def test_ctat_normal_fusion_pair_matches_double_dash_catalog(tmp_path):
    events = tmp_path / "events.tsv"
    peptides = tmp_path / "peptides.tsv"
    junctions = tmp_path / "normal_fusions.tsv"
    write_tsv(events, [{
        "event_id": "EF1", "sample_id": "S1", "event_type": "Fusion",
        "mutation_source": "SV", "gene": "GENE1::GENE2",
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "EF1", "sample_id": "S1",
        "event_type": "Fusion", "mutation_source": "SV", "gene": "GENE1::GENE2",
        "peptide_consequence": "fusion", "peptide": "AAAAAAAAA",
        "hla_allele": "HLA-A*02:01", "mhc_class": "I",
    }])
    write_tsv(junctions, [{
        "gene_pair": "GENE1--GENE2", "normal_sample_count": "4",
        "source": "GTEx_recurrent_StarF2019", "tissue": "Liver",
        "junction_class": "normal_recurrent_fusion",
    }])
    rows, _ = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
        profile=load_profile("default"),
        normal_junctions=junctions,
    )
    assert rows[0]["normal_junction_assessment_status"] == "ASSESSED"
    assert rows[0]["normal_junction_seen"] == "yes"
    assert rows[0]["normal_junction_max_reads"] == "4"
    assert rows[0]["safety_status"] == "FAIL"


def test_fusion_catalog_does_not_mark_splice_junction_assessed(tmp_path):
    events = tmp_path / "events.tsv"
    peptides = tmp_path / "peptides.tsv"
    junctions = tmp_path / "normal_fusions.tsv"
    write_tsv(events, [{
        "event_id": "S1", "sample_id": "S1", "event_type": "InDel",
        "mutation_source": "InDel", "gene": "GENE1",
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "S1", "sample_id": "S1",
        "event_type": "InDel", "mutation_source": "InDel", "gene": "GENE1",
        "peptide_consequence": "splice_junction", "peptide": "AAAAAAAAA",
        "hla_allele": "HLA-A*02:01", "mhc_class": "I",
    }])
    write_tsv(junctions, [{
        "gene_pair": "OTHER1::OTHER2", "normal_sample_count": "4",
        "junction_class": "normal_recurrent_fusion",
    }])
    rows, _ = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
        normal_junctions=junctions,
    )
    assert rows[0]["normal_junction_assessment_status"] == "NOT_APPLICABLE"
    assert "normal_junction" not in rows[0]["safety_missing_layers"]


def test_observed_splice_junction_uses_gtex_coordinate_catalog(tmp_path):
    events = tmp_path / "events.tsv"
    peptides = tmp_path / "peptides.tsv"
    junctions = tmp_path / "normal_splice_junctions.tsv"
    write_tsv(events, [{
        "event_id": "SJ1", "sample_id": "S1", "event_type": "splice_junction",
        "mutation_source": "splice", "gene": "GENE1",
        "junction_id": "chr1:100-200:+",
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "SJ1", "sample_id": "S1",
        "event_type": "splice_junction", "mutation_source": "splice", "gene": "GENE1",
        "junction_id": "chr1:100-200:+", "peptide_consequence": "splice_junction",
        "peptide": "AAAAAAAAA", "hla_allele": "HLA-A*02:01", "mhc_class": "I",
    }])
    write_tsv(junctions, [{
        "junction_id": "chr1:100-200:+", "normal_sample_count": "1",
        "source": "GTEx_v11_exon_junctions", "junction_class": "normal_splice_junction",
    }])
    rows, _ = build_peptide_safety_gate(
        raw_events=events,
        raw_peptides=peptides,
        out_peptide_safety=tmp_path / "peptide_safety.tsv",
        normal_junctions=junctions,
    )
    assert rows[0]["normal_junction_assessment_status"] == "ASSESSED"
    assert rows[0]["normal_junction_seen"] == "yes"
    assert rows[0]["normal_junction_source"] == "GTEx_v11_exon_junctions"


def test_normal_expression_merges_gene_symbol_and_ensembl_alias_layers(tmp_path):
    normal = tmp_path / "normal.tsv"
    write_tsv(normal, [
        {
            "gene": "GENE1", "ensembl_gene_id": "ENSG000001",
            "normal_tissue_max_tpm": "12", "normal_expression_status": "ASSESSED",
            "normal_hspc_tpm": "0", "normal_hspc_status": "UNASSESSED",
        },
        {
            "gene": "ENSG000001", "ensembl_gene_id": "ENSG000001",
            "normal_tissue_max_tpm": "0", "normal_expression_status": "UNASSESSED",
            "normal_hspc_tpm": "23", "normal_hspc_status": "ASSESSED",
            "normal_hspc_unit": "HPA_nCPM",
        },
    ])
    loaded = load_normal_expression(normal)
    assert loaded["GENE1"]["normal_expression_status"] == "ASSESSED"
    assert loaded["GENE1"]["normal_hspc_status"] == "ASSESSED"
    assert loaded["GENE1"]["normal_hspc_tpm"] == 23
    assert loaded["ENSG000001"]["normal_tissue_max_tpm"] == 12
