from pathlib import Path
from neoag_v03.adapters.pvactools_parser import parse_pvactools_outputs
from neoag_v03.adapters.netmhcpan import parse_netmhcpan, write_netmhcpan_evidence
from neoag_v03.adapters.mhcflurry import parse_mhcflurry, write_mhcflurry_evidence
from neoag_v03.presentation import build_presentation_evidence
from neoag_v03.config import load_profile
from neoag_v03.appm_lite import build_appm_lite
from neoag_v03.ccf_lite import build_ccf_lite
from neoag_v03.pipeline_v03 import run_v03
from neoag_v03.utils import write_tsv

ROOT = Path(__file__).resolve().parents[1]

def test_parse_pvac():
    e,p = parse_pvactools_outputs([ROOT/"data/fixtures/pvacseq_aggregated.tsv", ROOT/"data/fixtures/pvacfuse_aggregated.tsv"], "S1", "default")
    assert any(x["event_type"] == "Fusion" for x in e)
    assert any(x["peptide"] == "SPQKQWTRV" for x in p)

def test_presentation(tmp_path):
    profile = load_profile("default")
    e,p = parse_pvactools_outputs([ROOT/"data/fixtures/pvacseq_aggregated.tsv"], "S1", "default")
    pep = tmp_path/"peptides.tsv"; write_tsv(pep, p)
    net = tmp_path/"net.tsv"; write_netmhcpan_evidence(net, parse_netmhcpan(ROOT/"data/fixtures/netmhcpan_example.xls", "S1"))
    mhc = tmp_path/"mhc.tsv"; write_mhcflurry_evidence(mhc, parse_mhcflurry(ROOT/"data/fixtures/mhcflurry_predictions.csv", "S1"))
    rows = build_presentation_evidence(pep, net, mhc, profile, tmp_path/"presentation.tsv")
    assert any(r["presentation_evidence_grade"] in {"A","B"} for r in rows)

def test_appm_and_ccf(tmp_path):
    profile = load_profile("leukemia")
    rows, summary = build_appm_lite("S1", ROOT/"data/fixtures/vep_appm.tsv", ROOT/"data/fixtures/gene_expression.tsv", ROOT/"data/fixtures/hla_loh.tsv", profile, tmp_path/"appm")
    assert float(summary["mhc_ii_integrity_score"]) < 1.0
    e,p = parse_pvactools_outputs([ROOT/"data/fixtures/pvacseq_aggregated.tsv"], "DEMO_V03", "default")
    ev = tmp_path/"events.tsv"; write_tsv(ev, e)
    ccf = build_ccf_lite(ev, ROOT/"data/fixtures/purity.tsv", ROOT/"data/fixtures/cnv_segments.tsv", profile, tmp_path/"ccf.tsv")
    assert ccf and ccf[0]["ccf_status"]

def test_run_v03(tmp_path):
    out = run_v03(
        outdir=tmp_path/"out", profile_name_or_path="leukemia", sample_id="DEMO_V03",
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
    for v in out.values():
        assert Path(v).exists()
    assert "immune_escape_summary" in out
    assert "ccf_2" in out
    assert Path(out["ccf_2"]).exists()
    assert Path(out["ccf_lite"]).exists()
    ccf_header = Path(out["ccf_2"]).read_text().splitlines()[0]
    assert "clonality_confidence" in ccf_header
    assert "ccf_resolution" in ccf_header
    txt = Path(out["ranked_peptides"]).read_text()
    assert "presentation_evidence_grade" in txt
    assert "appm_multiplier" in txt
    assert "efficacy_score" in txt
    assert "presentation_gate_status" in txt
    assert "immunogenicity_resolved" in txt
    assert "immunogenicity_composite_score" in txt
    assert "immunogenicity_source" in txt
    assert "mutation_source" in txt
    assert "peptide_consequence" in txt
    assert "immunology_composite_score" in txt
    assert "l3_hla_binding_score" in txt
