from neoag.schemas import PEPTIDE_FIELDS
from neoag.splice_prefilter import prefilter_splice_peptides
from neoag.utils import read_tsv, write_tsv


def splice_row(event_id: str, **updates):
    row = {field: "" for field in PEPTIDE_FIELDS}
    row.update({
        "peptide_id": event_id + "-P", "event_id": event_id, "event_type": "Splice",
        "peptide_consequence": "splice_junction", "peptide": "ABCDEFGHI",
        "hla_allele": "HLA-A*02:01", "junction_resolution_status": "RESOLVED",
        "unique_junction_reads": "5", "junction_total_coverage": "20", "splice_psi": "0.2",
        "matched_normal_junction_status": "NOT_DETECTED",
        "normal_cohort_junction_status": "NOT_DETECTED",
        "annotated_normal_isoform_status": "NOVEL", "splice_orf_status": "IN_FRAME",
        "splice_nmd_status": "NOT_PREDICTED", "crosses_junction": "yes",
        "reference_proteome_status": "NOT_DETECTED", "rna_junction_reads": "5",
    })
    row.update(updates)
    return row


def test_splice_prefilter_runs_before_prediction_and_keeps_bounded_review_lane(tmp_path):
    path = tmp_path / "raw_peptides.tsv"
    non_splice = {field: "" for field in PEPTIDE_FIELDS}
    non_splice.update({"peptide_id": "V1-P", "event_id": "V1", "event_type": "SNV", "peptide": "ABCDEFGHI"})
    write_tsv(path, [non_splice, splice_row("PASS"), splice_row("REJECT", unique_junction_reads="1"), splice_row(
        "REVIEW", unique_junction_reads="", junction_total_coverage="", splice_psi="",
        matched_normal_junction_status="", normal_cohort_junction_status="",
    )], PEPTIDE_FIELDS)
    profile = {"gates": {
        "min_splice_unique_reads": 3, "min_splice_total_coverage": 10,
        "min_splice_psi": 0.05, "max_splice_review_prediction_events": 1,
        "max_splice_peptide_hla_per_event": 2,
    }}
    result = prefilter_splice_peptides(path, profile, tmp_path)
    selected = {row["event_id"]: row for row in read_tsv(path)}
    assert set(selected) == {"V1", "PASS", "REVIEW"}
    assert selected["PASS"]["splice_prefilter_status"] == "PASS"
    assert selected["REVIEW"]["splice_prefilter_status"] == "REVIEW"
    assert result["raw_splice_events"] == 3
    assert result["selected_events"] == 2
    decisions = {row["event_id"]: row for row in read_tsv(tmp_path / "splice_prefilter_decisions.tsv")}
    assert decisions["REJECT"]["prefilter_status"] == "REJECT"
    assert "UNIQUE_JUNCTION_READS" in decisions["REJECT"]["failed_stages"]
    funnel = {row["stage"]: row for row in read_tsv(tmp_path / "splice_prefilter_funnel.tsv")}
    assert funnel["UNIQUE_JUNCTION_READS"]["passed_events"] == "1"
    assert funnel["UNIQUE_JUNCTION_READS"]["failed_events"] == "1"
    assert funnel["UNIQUE_JUNCTION_READS"]["unassessed_events"] == "1"
    assert funnel["SELECTED_FOR_PRESENTATION"]["passed_events"] == "2"


def test_splice_prefilter_rejects_known_normal_isoform_and_unstranded_mapping(tmp_path):
    path = tmp_path / "raw_peptides.tsv"
    write_tsv(path, [
        splice_row("KNOWN", annotated_normal_isoform_status="KNOWN_NORMAL"),
        splice_row("UNSTRANDED", junction_resolution_status="RESOLVED_UNSTRANDED"),
    ], PEPTIDE_FIELDS)
    prefilter_splice_peptides(path, {"gates": {}}, tmp_path)
    assert read_tsv(path) == []
    decisions = {row["event_id"]: row for row in read_tsv(tmp_path / "splice_prefilter_decisions.tsv")}
    assert "ANNOTATED_NORMAL_ISOFORM" in decisions["KNOWN"]["failed_stages"]
    assert "ALIGNMENT_COORDINATE_QC" in decisions["UNSTRANDED"]["failed_stages"]
