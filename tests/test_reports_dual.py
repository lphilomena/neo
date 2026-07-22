from pathlib import Path

from neoag.reports_dual import ReportBundle, load_report_bundle, make_dual_reports, make_patient_report, make_technical_report
from neoag.utils import write_tsv


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


def test_technical_report_loads_independent_wes_qc(tmp_path):
    qc_dir = tmp_path / "qc" / "wes"
    qc_dir.mkdir(parents=True)
    write_tsv(qc_dir / "wes_qc.tsv", [{
        "sample_id": "WES1", "qc_status": "PASS_WITH_CAPTURE_RATE_UNASSESSED",
        "total_reads": "1000", "primary_mapping_rate_pct": "99.9",
        "properly_paired_rate_pct": "98.0", "duplicate_rate_pct": "20.0",
        "target_definition": "GENCODE_CDS_PROXY_NOT_ASSAY_CAPTURE_BED",
        "mean_target_coverage": "80", "pct_target_bases_20x": "95",
        "pct_target_bases_30x": "90", "on_target_rate_pct": "70",
        "capture_rate_status": "UNASSESSED_CAPTURE_BED_MISSING",
        "formal_capture_rate_pct": "",
    }])
    base = _bundle()
    bundle = load_report_bundle(
        profile=base.profile, events=base.events, peptides=base.peptides,
        appm_summary=base.appm_summary, validation_rows=base.validation_rows,
        outdir=tmp_path, sample_id="S1",
    )
    out = tmp_path / "technical_wes.html"
    make_technical_report(out, bundle)
    text = out.read_text(encoding="utf-8")
    assert "Independent WES QC" in text
    assert "PASS_WITH_CAPTURE_RATE_UNASSESSED" in text
    assert "assay-specific capture BED" in text
