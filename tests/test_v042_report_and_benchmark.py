
from pathlib import Path
from neoag_v03.reports_v041 import make_report_v041
from neoag_v03.benchmark_system import run_system_benchmark
from neoag_v03.utils import read_tsv


def test_report_v041_contains_appm_cards(tmp_path):
    out = tmp_path / "report.html"
    make_report_v041(
        out, {"_profile_name":"test"},
        events=[{"event_id":"E1","event_name":"E1","event_type":"SNV","gene":"B2M","event_score":"0.5"}],
        peptides=[{"peptide_id":"P1","event_id":"E1","peptide":"AAAAAAAAA","hla_allele":"HLA-A*02:01","gene":"B2M","final_priority":"B","efficacy_score":"0.8"}],
        appm_summary={"mhc_i_integrity_score":"0.05","mhc_i_integrity_status":"MHC_I_DEFECTIVE","appm_call_confidence":"high","appm_call_confidence_score":"0.9","confidence_reason":"B2M_LOSS_HIGH_CONFIDENCE"},
        appm_submodule_scores=[{"parent_module":"MHC-I","submodule":"MHC_I_CORE","score":"0.05","status":"MHC_I_CORE_DEFECTIVE","defect_severity":"lethal","appm_call_confidence":"high","driver_defects":"B2M_biallelic_loss","action_hint":"hard_cap_mhc_i","confidence_reason":"B2M_LOSS_HIGH_CONFIDENCE"}],
        appm_gene_status=[{"gene":"B2M","pathway":"MHC-I","biallelic_status":"BIALLELIC_LOSS","functional_status":"defective","reason":"test"}],
    )
    text = out.read_text(encoding="utf-8")
    assert "APPM call confidence" in text
    assert "MHC_I_CORE" in text
    assert "Peptide Mechanism Cards" in text
    assert "v0.4.2" not in text
    assert "v0.4.3" in text
    assert "<title>NeoAg v0.4.3 Evidence Report</title>" in text


def test_benchmark_v042_external_required_outputs(tmp_path):
    out = run_system_benchmark(outdir=tmp_path, mode="ligandome-ms")
    assert Path(out["appm_ms_stratified_validation"]).exists()
    row = read_tsv(out["appm_ms_stratified_validation"])[0]
    assert row["benchmark_status"] == "external_required"
