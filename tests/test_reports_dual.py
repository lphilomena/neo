from pathlib import Path

from neoag.reports_dual import ReportBundle, make_dual_reports, make_patient_report, make_technical_report


def _bundle():
    events = [{
        "event_id": "B2M|chr15:1A>G", "event_name": "B2M variant", "event_type": "SNV",
        "gene": "B2M", "event_score": "0.6", "ccf_status": "clonal", "safety_status": "PASS",
    }]
    peptides = [{
        "peptide_id": "P1", "event_id": "B2M|chr15:1A>G", "sample_id": "S1",
        "gene": "B2M", "peptide": "AAAAAAAAA", "wildtype_peptide": "AAABAAAAA",
        "hla_allele": "HLA-A*02:01", "event_type": "SNV", "peptide_consequence": "missense",
        "final_priority": "B", "efficacy_score": "0.72", "safety_status": "PASS",
        "presentation_evidence_grade": "A", "recommended_use": "MHC-I short peptide ELISpot/tetramer (MT + WT pair)",
    }]
    validation = [{
        "peptide_id": "P1", "validation_mode": "missense_short_pair",
        "recommended_assay": "MHC-I short peptide ELISpot/tetramer (MT + WT pair)",
        "validation_strategy": "Mutant short peptide with WT control",
    }]
    return ReportBundle(
        profile={"_profile_name": "test", "gates": {"max_el_rank": 2.0}},
        events=events,
        peptides=peptides,
        appm_summary={"mhc_i_integrity_status": "MHC_I_INTACT", "mhc_i_integrity_score": "0.95"},
        validation_rows=validation,
        sample_id="S1",
        entry_mode="snv_indel",
        provenance={"sample_id": "S1", "tools": {"netmhcpan": {"status": "real", "version": "4.2", "file": "/data/netmhcpan.xls"}}},
    )


def test_patient_report_is_plain_language(tmp_path):
    out = tmp_path / "patient.html"
    make_patient_report(out, _bundle())
    text = out.read_text(encoding="utf-8")
    assert "患者沟通版" in text or "新抗原计算分析报告" in text
    assert "不能替代临床诊断" in text
    assert "netmhcpan.xls" not in text
    assert "fastq" not in text.lower()
    assert "移码变异" not in text or "点突变" in text


def test_technical_report_has_provenance_and_thresholds(tmp_path):
    out = tmp_path / "technical.html"
    make_technical_report(out, _bundle())
    text = out.read_text(encoding="utf-8")
    assert "Tool provenance" in text
    assert "netmhcpan.xls" in text
    assert "Profile Thresholds" in text
    assert "Field Glossary" in text
    assert "Ranked Peptides (full)" in text


def test_dual_reports_writes_three_files(tmp_path):
    paths = make_dual_reports(tmp_path, _bundle())
    assert Path(paths["evidence_report_patient"]).exists()
    assert Path(paths["evidence_report_technical"]).exists()
    assert Path(paths["evidence_report"]).exists()
    assert "Technical Report" in Path(paths["evidence_report_technical"]).read_text(encoding="utf-8")
