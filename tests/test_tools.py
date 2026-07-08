from pathlib import Path
import pytest

from neoag_v03.tools import check_tool, run_tool, run_upstream, load_run_config
from neoag_v03.tools.registry import RunContext, ROOT
from neoag_v03.tools.runner import RUNNERS
from neoag_v03.tools.prep import unique_peptide_hla_pairs, netmhcpan_allele_string
from neoag_v03.pipeline_v03 import run_v03
from neoag_v03.cli import main

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
    from neoag_v03.adapters.pvactools_parser import parse_pvactools_outputs
    raw = tmp_path / "raw_peptides.tsv"
    parse_pvactools_outputs([pep], "S1", "default", None, raw)
    ctx = RunContext(sample_id="S1", outdir=tmp_path, stub=True, raw_peptides=raw, hla_alleles=["HLA-A*02:01"])
    out = run_tool("netmhcpan", ctx, tmp_path / "net.xls")
    assert out.exists()

def test_run_full_cli_stub(tmp_path):
    outdir = tmp_path / "full"
    main(["run-full", "--config", str(ROOT / "conf/run.stub.toml"), "--outdir", str(outdir)])
    assert (outdir / "scoring/ranked_peptides.v03.tsv").exists()

def test_fusion_tools_in_registry():
    from neoag_v03.tools.registry import TOOL_REGISTRY

    for name in ("star_fusion", "arriba", "fusioncatcher", "easyfuse"):
        assert name in TOOL_REGISTRY
        assert TOOL_REGISTRY[name].category == "fusion"


def test_check_tools_runs():
    st = check_tool("netmhcpan")
    assert st.name == "netmhcpan"
