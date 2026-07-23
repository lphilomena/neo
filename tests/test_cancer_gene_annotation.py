from pathlib import Path

from neoag.cancer_gene_annotation import (
    CancerGeneIndex,
    annotate_cancer_gene_context,
    propagate_to_peptide,
)
from neoag.scoring import score, score_event
from neoag.config import load_profile
from neoag.utils import write_tsv


def rows():
    return [
        {"Hugo Symbol": "EWSR1", "Gene Type": "ONCOGENE", "OncoKB Annotated": "Yes", "COSMIC CGC (v99)": "Yes", "# of occurrence within resources (Column K-P)": "5", "Gene Aliases": "EWS"},
        {"Hugo Symbol": "WT1", "Gene Type": "ONCOGENE_AND_TSG", "OncoKB Annotated": "Yes", "COSMIC CGC (v99)": "Yes", "# of occurrence within resources (Column K-P)": "6", "Gene Aliases": "WAGR"},
        {"Hugo Symbol": "ALB", "Gene Type": "NEITHER", "OncoKB Annotated": "Yes", "COSMIC CGC (v99)": "No", "# of occurrence within resources (Column K-P)": "1", "Gene Aliases": "HSA"},
    ]


def test_fusion_and_neither_annotations_propagate_without_score_change():
    index = CancerGeneIndex(rows())
    base = {
        "event_id": "E1", "gene": "EWSR1::WT1", "event_type": "Fusion",
        "event_confidence": "0.8", "event_expression": "10", "driver_relevance": "0.4",
        "clonality": "0.5", "persistence": "0.5", "tumor_specificity": "0.8",
        "safety_status": "PASS", "clonality_multiplier": "1", "cross_platform_multiplier": "1",
    }
    before = score_event(dict(base), load_profile("default"))["event_score"]
    event = annotate_cancer_gene_context(base, index)
    after = score_event(dict(event), load_profile("default"))["event_score"]
    assert before == after
    assert event["cancer_gene_symbols"] == "EWSR1;WT1"
    assert event["cancer_driver_context"] == "DRIVER_CONTEXT"
    peptide = propagate_to_peptide({"peptide_id": "P1"}, event)
    assert peptide["cancer_gene_types"] == "EWSR1:ONCOGENE;WT1:ONCOGENE_AND_TSG"

    alb = annotate_cancer_gene_context({"gene": "ALB"}, index)
    assert alb["cancer_gene_list_status"] == "ANNOTATED"
    assert alb["cancer_driver_context"] == "LISTED_NO_DRIVER_CLASS"


def test_alias_and_missing_statuses():
    index = CancerGeneIndex(rows())
    alias = annotate_cancer_gene_context({"gene": "EWS"}, index)
    assert alias["cancer_gene_symbols"] == "EWSR1"
    assert alias["cancer_gene_match_basis"] == "alias"
    missing = annotate_cancer_gene_context({"gene": "TBR1"}, index)
    assert missing["cancer_gene_list_status"] == "NOT_LISTED"
    unassessed = annotate_cancer_gene_context({"gene": "KRAS"}, None)
    assert unassessed["cancer_gene_list_status"] == "UNASSESSED"


def test_score_integration_writes_event_and_peptide_annotations(tmp_path: Path):
    events = tmp_path / "events.tsv"
    peptides = tmp_path / "peptides.tsv"
    presentation = tmp_path / "presentation.tsv"
    genes = tmp_path / "cancer_genes.tsv"
    write_tsv(events, [{
        "event_id": "E1", "sample_id": "S1", "event_type": "SNV", "mutation_source": "SNV",
        "gene": "ALB", "event_confidence": "0.8", "event_expression": "10",
        "driver_relevance": "0.2", "clonality": "0.5", "persistence": "0.5",
        "tumor_specificity": "0.8", "safety_status": "PASS",
    }])
    write_tsv(peptides, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1", "event_type": "SNV",
        "mutation_source": "SNV", "gene": "ALB", "peptide": "AAAAAAAAA",
        "hla_allele": "HLA-A*02:01", "mhc_class": "I",
    }])
    write_tsv(presentation, [{
        "peptide_id": "P1", "event_id": "E1", "sample_id": "S1", "peptide": "AAAAAAAAA",
        "hla_allele": "HLA-A*02:01", "mhc_class": "I", "binding_evidence_score": "0.8",
        "presentation_evidence_score": "0.8", "presentation_evidence_grade": "A",
    }])
    write_tsv(genes, rows())
    out_events = tmp_path / "ranked_events.tsv"
    out_peptides = tmp_path / "ranked_peptides.tsv"
    scored_events, scored_peptides = score(
        events, peptides, presentation, None, None, None, None, load_profile("default"),
        out_events, out_peptides, cancer_gene_list_tsv=genes,
    )
    assert scored_events[0]["cancer_driver_context"] == "LISTED_NO_DRIVER_CLASS"
    assert scored_peptides[0]["cancer_gene_symbols"] == "ALB"
    assert scored_peptides[0]["cancer_gene_context"].startswith("listed without")
