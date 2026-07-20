from pathlib import Path
from neoag.pipeline import run
from neoag.cli import main

ROOT = Path(__file__).resolve().parents[1]

def test_validation_and_report_cli(tmp_path):
    out = run(
        outdir=tmp_path/"run",
        profile_name_or_path="leukemia",
        sample_id="DEMO_RC",
        pvac_paths=[ROOT/"data/fixtures/pvacseq_aggregated.tsv", ROOT/"data/fixtures/pvacfuse_aggregated.tsv"],
        netmhcpan=ROOT/"data/fixtures/netmhcpan_example.xls",
        mhcflurry=ROOT/"data/fixtures/mhcflurry_predictions.csv",
        vep_appm=ROOT/"data/fixtures/vep_appm.tsv",
        expression=ROOT/"data/fixtures/gene_expression.tsv",
        hla_loh=ROOT/"data/fixtures/hla_loh.tsv",
        purity=ROOT/"data/fixtures/purity.tsv",
        cnv=ROOT/"data/fixtures/cnv_segments.tsv",
        normal_expression=ROOT/"resources/normal_expression.example.tsv",
        normal_hla_ligands=ROOT/"resources/normal_hla_ligands.example.tsv",
    )
    val = tmp_path/"validation.tsv"
    main(["validation-plan", "--ranked-peptides", out["ranked_peptides"], "--out", str(val)])
    assert val.exists()
    rep = tmp_path/"report.html"
    main(["report", "--profile", "leukemia", "--ranked-events", out["ranked_events"], "--ranked-peptides", out["ranked_peptides"], "--appm-summary", out["appm_summary"], "--validation-plan", str(val), "--outdir", str(tmp_path/"run"), "--out", str(rep)])
    assert rep.exists()
    patient = tmp_path/"evidence_report.patient.html"
    technical = tmp_path/"evidence_report.technical.html"
    assert patient.exists()
    assert technical.exists()
    assert "新抗原" in patient.read_text(encoding="utf-8")
    assert "Technical" in technical.read_text(encoding="utf-8")
    assert out.get("evidence_report_patient") and Path(out["evidence_report_patient"]).exists()
