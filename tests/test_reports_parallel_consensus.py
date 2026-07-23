from pathlib import Path

from neoag.reports_dual import ReportBundle, make_dual_reports


def test_parallel_consensus_is_technical_only(tmp_path: Path) -> None:
    bundle = ReportBundle(
        profile={"_profile_name": "test"},
        events=[],
        peptides=[
            {
                "peptide_id": "P_WEIGHTED",
                "event_id": "E1",
                "sample_id": "S1",
                "peptide": "AAAAAAAAA",
                "hla_allele": "HLA-A*02:01",
                "final_priority": "B",
            }
        ],
        sample_id="S1",
        provenance={
            "parallel_rankings": {
                "legacy_weighted": "scoring/ranked_peptides.tsv",
                "evidence_consensus": "scoring/ranked_peptides.evidence_consensus.tsv",
                "event_consensus": "scoring/ranked_events.evidence_consensus.tsv",
                "comparison": "scoring/ranking_compare_weighted_vs_consensus.tsv",
                "rules_name": "sarcoma_evidence_consensus_v1",
                "rules_version": "1.0",
                "rules_status": "PROVISIONAL_RESEARCH_ONLY",
            }
        },
    )

    paths = make_dual_reports(tmp_path, bundle)
    patient = Path(paths["evidence_report_patient"]).read_text(encoding="utf-8")
    technical = Path(paths["evidence_report_technical"]).read_text(encoding="utf-8")

    assert "Experimental parallel evidence-consensus ranking" not in patient
    assert "ranked_peptides.evidence_consensus.tsv" not in patient
    assert "Experimental parallel evidence-consensus ranking" in technical
    assert "does not replace the current primary weighted ranking" in technical
    assert "has not been experimentally calibrated" in technical
    assert "only for algorithm comparison and candidate review" in technical
    assert "ranked_peptides.evidence_consensus.tsv" in technical
