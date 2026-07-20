import json
from pathlib import Path
import pytest

from neoag.tools import check_tool, run_tool, run_upstream, load_run_config
from neoag.tools.registry import RunContext, ROOT
from neoag.tools.runner import RUNNERS
from neoag.utils import read_tsv
from neoag.tools.prep import unique_peptide_hla_pairs, netmhcpan_allele_string
from neoag.pipeline import run
from neoag.cli import main

def test_tool_registry_covers_runners():
    assert "netmhcpan" in RUNNERS
    assert "pvacseq" in RUNNERS

def test_netmhcpan_allele_format():
    assert "HLA-A0201" in netmhcpan_allele_string(["HLA-A*02:01"])

def test_unique_peptides_from_fixture():
    pairs = unique_peptide_hla_pairs(ROOT / "data/fixtures/pvacseq_aggregated.tsv")
    assert ("VVVGADGVGK", "HLA-A*11:01") in pairs

def test_run_upstream_stub(tmp_path):
    cfg = ROOT / "conf/run.stub.toml"
    outs = run_upstream(cfg, tmp_path / "up")
    assert Path(outs["netmhcpan"]).exists()
    assert Path(outs["raw_peptides"]).exists()

def test_run_tool_stub_netmhcpan(tmp_path):
    pep = ROOT / "data/fixtures/pvacseq_aggregated.tsv"
    from neoag.adapters.pvactools_parser import parse_pvactools_outputs
    raw = tmp_path / "raw_peptides.tsv"
    parse_pvactools_outputs([pep], "S1", "default", None, raw)
    ctx = RunContext(sample_id="S1", outdir=tmp_path, stub=True, raw_peptides=raw, hla_alleles=["HLA-A*02:01"])
    out = run_tool("netmhcpan", ctx, tmp_path / "net.xls")
    assert out.exists()

def test_run_full_cli_stub(tmp_path):
    outdir = tmp_path / "full"
    main(["run-full", "--config", str(ROOT / "conf/run.stub.toml"), "--outdir", str(outdir)])
    assert (outdir / "scoring/ranked_peptides.tsv").exists()
    assert (outdir / "scoring/comprehensive_peptide_evidence.tsv").exists()


def test_multisource_peptides_all_enter_presentation_prediction(tmp_path):
    cfg = tmp_path / "multisource.toml"
    cfg.write_text(
        f'''[sample]
id = "MULTISOURCE"
profile = "default"

[tools]
stub = true
enabled = ["netmhcpan", "mhcflurry"]

[inputs]
entry_mode = "e2e"
hla_alleles = ["HLA-A*02:01"]
pvac_files = [
  "{ROOT / 'data/fixtures/pvacseq_aggregated.tsv'}",
  "{ROOT / 'data/fixtures/pvacfuse_aggregated.tsv'}",
  "{ROOT / 'data/fixtures/pvacsplice_aggregated.tsv'}",
]
expected_peptide_sources = ["pVACseq", "pVACfuse", "pVACsplice"]
required_presentation_predictors = ["netmhcpan", "mhcflurry"]
extract_appm_from_vcf = false
''',
        encoding="utf-8",
    )

    outputs = run_upstream(cfg, tmp_path / "upstream")
    peptides = read_tsv(outputs["raw_peptides"])

    assert {row["source_tool"] for row in peptides} == {"pVACseq", "pVACfuse", "pVACsplice"}
    assert outputs["peptide_sources"] == "pVACfuse,pVACseq,pVACsplice"
    coverage = read_tsv(outputs["peptide_source_coverage"])[0]
    assert coverage["status"] == "COMPLETE"
    assert Path(outputs["netmhcpan"]).is_file()
    assert Path(outputs["mhcflurry"]).is_file()


def test_multisource_missing_fusion_reports_low_confidence(tmp_path, capsys):
    cfg = tmp_path / "missing_fusion.toml"
    cfg.write_text(
        f'''[sample]
id = "INCOMPLETE"
profile = "default"

[tools]
stub = true
enabled = []

[inputs]
entry_mode = "e2e"
pvac_files = ["{ROOT / 'data/fixtures/pvacseq_aggregated.tsv'}"]
expected_peptide_sources = ["pVACseq", "pVACfuse"]
extract_appm_from_vcf = false
''',
        encoding="utf-8",
    )

    outputs = run_upstream(cfg, tmp_path / "upstream")

    assert outputs["peptide_source_completeness"] == "LOW_CONFIDENCE"
    assert outputs["missing_peptide_sources"] == "pVACfuse"
    assert "confidence is LOW" in capsys.readouterr().out
    coverage = read_tsv(outputs["peptide_source_coverage"])[0]
    assert coverage["missing_sources"] == "pVACfuse"


def test_upstream_prefers_purity_recommendation_over_facets(tmp_path):
    recommendation = tmp_path / "purity_recommendation.json"
    recommendation.write_text(json.dumps({
        "status": "CONCORDANT", "recommended_purity": 0.67, "range": "0.6500-0.6900",
        "n_tools": 2, "tool_values": {"FACETS": 0.65, "PURPLE": 0.69},
    }), encoding="utf-8")
    cfg = tmp_path / "purity.toml"
    cfg.write_text(
        f'''[sample]
id = "PURITY_CONSENSUS"
profile = "default"

[tools]
stub = true
enabled = ["facets"]

[inputs]
purity_recommendation = "{recommendation}"
extract_appm_from_vcf = false
''',
        encoding="utf-8",
    )

    outputs = run_upstream(cfg, tmp_path / "upstream")

    assert outputs["purity"] == str(recommendation)
    assert outputs["purity_recommendation"] == str(recommendation)
    assert outputs["facets_purity"].endswith("facets_purity.tsv")
    assert Path(outputs["facets_purity"]).is_file()

def test_fusion_tools_in_registry():
    from neoag.tools.registry import TOOL_REGISTRY

    for name in ("star_fusion", "arriba", "fusioncatcher", "easyfuse"):
        assert name in TOOL_REGISTRY
        assert TOOL_REGISTRY[name].category == "fusion"


def test_check_tools_runs():
    st = check_tool("netmhcpan")
    assert st.name == "netmhcpan"
