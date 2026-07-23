from pathlib import Path
from neoag.adapters.pvactools_parser import parse_pvactools_outputs
from neoag.adapters.netmhcpan import parse_netmhcpan, write_netmhcpan_evidence
from neoag.adapters.mhcflurry import parse_mhcflurry, write_mhcflurry_evidence
from neoag.presentation import build_presentation_evidence
from neoag.config import load_profile
from neoag.appm_lite import build_appm_lite
from neoag.ccf_lite import build_ccf_lite
from neoag.pipeline import run
from neoag.utils import write_tsv

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
    e,p = parse_pvactools_outputs([ROOT/"data/fixtures/pvacseq_aggregated.tsv"], "DEMO", "default")
    ev = tmp_path/"events.tsv"; write_tsv(ev, e)
    ccf = build_ccf_lite(ev, ROOT/"data/fixtures/purity.tsv", ROOT/"data/fixtures/cnv_segments.tsv", profile, tmp_path/"ccf.tsv")
    assert ccf and ccf[0]["ccf_status"]

def test_run(tmp_path):
    out = run(
        outdir=tmp_path/"out", profile_name_or_path="leukemia", sample_id="DEMO",
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
    assert "ranked_peptides_evidence_consensus" in out
    assert Path(out["ranked_peptides_evidence_consensus"]).exists()
    assert Path(out["ranked_events_evidence_consensus"]).exists()
    assert Path(out["evidence_states"]).exists()
    assert Path(out["evidence_conflicts"]).exists()
    assert Path(out["evidence_source_conflicts"]).exists()
    assert Path(out["weighted_vs_consensus_comparison"]).exists()
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
    consensus_txt = Path(out["ranked_peptides_evidence_consensus"]).read_text()
    assert "evidence_consensus_score" in consensus_txt
    assert "pareto_front" in consensus_txt
    assert "evidence_grade" in consensus_txt
    assert "evidence_completeness_score" in consensus_txt
    assert "evidence_missing_layers" in consensus_txt


def test_run_accepts_raw_inputs_already_in_outdir(tmp_path):
    outdir = tmp_path / "sample"
    parsed = outdir / "parsed"
    parsed.mkdir(parents=True)
    events, peptides = parse_pvactools_outputs(
        [ROOT / "data/fixtures/pvacseq_aggregated.tsv"],
        "SAMPLE001",
        "default",
    )
    raw_events = parsed / "raw_events.tsv"
    raw_peptides = parsed / "raw_peptides.tsv"
    write_tsv(raw_events, events)
    write_tsv(raw_peptides, peptides)

    out = run(
        outdir=outdir,
        profile_name_or_path="default",
        sample_id="SAMPLE001",
        raw_events=raw_events,
        raw_peptides=raw_peptides,
        netmhcpan=ROOT / "data/fixtures/netmhcpan_example.xls",
        mhcflurry=ROOT / "data/fixtures/mhcflurry_predictions.csv",
        expression=ROOT / "data/fixtures/gene_expression.tsv",
        hla_loh=ROOT / "data/fixtures/hla_loh.tsv",
        purity=ROOT / "data/fixtures/purity.tsv",
        cnv=ROOT / "data/fixtures/cnv_segments.tsv",
        normal_expression=ROOT / "resources/normal_expression.example.tsv",
        normal_hla_ligands=ROOT / "resources/normal_hla_ligands.example.tsv",
    )

    assert Path(out["raw_events"]) == raw_events
    assert Path(out["ranked_peptides"]).exists()
