from __future__ import annotations

from pathlib import Path

from neoag.adapters.netmhcpan import parse_netmhcpan, write_netmhcpan_evidence
from neoag.evidence_provenance import (
    MODE_TOOL_RUN,
    SOURCE_EXTERNAL,
    SOURCE_FIXTURE,
    SOURCE_LOCAL,
    ProvenanceRegistry,
    detect_tool_version,
    provenance_for_netmhcpan_run,
    provenance_not_used,
    provenance_stub,
)
from neoag.schemas import EVIDENCE_PROVENANCE_FIELDS
from neoag.presentation import build_presentation_evidence
from neoag.config import load_profile
from neoag.utils import read_tsv, write_tsv

ROOT = Path(__file__).resolve().parents[1]


def test_netmhcpan_evidence_has_provenance_fields(tmp_path):
    rows = parse_netmhcpan(ROOT / "data/fixtures/netmhcpan_example.xls", "S1")
    out = tmp_path / "net.tsv"
    write_netmhcpan_evidence(out, rows)
    written = read_tsv(out)
    assert written
    for field in EVIDENCE_PROVENANCE_FIELDS:
        assert field in written[0]
    assert written[0]["evidence_tool"] == "netmhcpan"
    assert written[0]["evidence_status"] == "real"
    assert written[0]["evidence_source"] == SOURCE_FIXTURE
    assert written[0]["evidence_mode"] == "passthrough"


def test_netmhcpan_tool_run_provenance_columns(tmp_path, monkeypatch):
    monkeypatch.delenv("NEOAG_NETMHCPAN_BACKEND", raising=False)
    from neoag.tools.registry import RunContext
    from neoag.tools.runner import run_netmhcpan

    pep = tmp_path / "raw_peptides.tsv"
    write_tsv(
        pep,
        [{
            "peptide_id": "p1",
            "peptide": "SIINFEKL",
            "hla_allele": "HLA-A*02:01",
            "sample_id": "T1",
        }],
    )
    ctx = RunContext(sample_id="T1", outdir=tmp_path, stub=False, raw_peptides=pep)
    xls = tmp_path / "tools" / "netmhcpan.xls"
    run_netmhcpan(ctx, xls)

    prov = ctx.tool_provenance["netmhcpan"]
    assert prov.mode == MODE_TOOL_RUN
    assert prov.status == "real"
    assert prov.source in {SOURCE_LOCAL, SOURCE_EXTERNAL}
    if prov.source == SOURCE_LOCAL:
        assert "4.2" in (prov.version or detect_tool_version("netmhcpan"))
    else:
        assert prov.version == "iedb-netmhcpan"

    out = tmp_path / "net_evidence.tsv"
    write_netmhcpan_evidence(out, parse_netmhcpan(xls, "T1"), prov)
    row = read_tsv(out)[0]
    assert row["evidence_mode"] == MODE_TOOL_RUN
    assert row["evidence_source"] == prov.source
    assert row["evidence_tool"] == "netmhcpan"
    assert row["evidence_status"] == "real"
    if prov.source == SOURCE_LOCAL:
        assert "4.2" in row["evidence_tool_version"]
    else:
        assert row["evidence_tool_version"] == "iedb-netmhcpan"

    registry = ProvenanceRegistry()
    registry.set(prov)
    summary = registry.tool_summary_fields(("netmhcpan",))
    assert summary["netmhcpan_source"] == prov.source
    assert summary["netmhcpan_status"] == "real"
    if prov.source == SOURCE_LOCAL:
        assert "4.2" in summary["netmhcpan_version"]
    else:
        assert summary["netmhcpan_version"] == "iedb-netmhcpan"


def test_provenance_for_netmhcpan_run_iedb():
    rec = provenance_for_netmhcpan_run("/tmp/net.xls", "iedb")
    assert rec.source == "external"
    assert rec.mode == MODE_TOOL_RUN
    assert rec.version == "iedb-netmhcpan"


def test_presentation_tool_summary_columns(tmp_path):
    profile = load_profile("default")
    registry = ProvenanceRegistry()
    registry.register_passthrough("netmhcpan", ROOT / "data/fixtures/netmhcpan_example.xls")
    registry.register_passthrough("mhcflurry", ROOT / "data/fixtures/mhcflurry_predictions.csv")
    registry.register_not_used("bigmhc_im")
    registry.register_stub("prime")

    from neoag.adapters.pvactools_parser import parse_pvactools_outputs

    _, peptides = parse_pvactools_outputs([ROOT / "data/fixtures/pvacseq_aggregated.tsv"], "S1", "default")
    pep = tmp_path / "peptides.tsv"
    write_tsv(pep, peptides)
    net = tmp_path / "net.tsv"
    write_netmhcpan_evidence(net, parse_netmhcpan(ROOT / "data/fixtures/netmhcpan_example.xls", "S1"))
    from neoag.adapters.mhcflurry import parse_mhcflurry, write_mhcflurry_evidence

    mhc = tmp_path / "mhc.tsv"
    write_mhcflurry_evidence(mhc, parse_mhcflurry(ROOT / "data/fixtures/mhcflurry_predictions.csv", "S1"))

    rows = build_presentation_evidence(pep, net, mhc, profile, tmp_path / "pres.tsv", provenance_registry=registry)
    row = rows[0]
    assert row["netmhcpan_status"] == "real"
    assert row["bigmhc_im_status"] == "not_used"
    assert row["prime_status"] == "invalid_for_production"
    assert row["evidence_tool"] == "presentation_composite"


def test_registry_tool_summary():
    reg = ProvenanceRegistry()
    reg.set(provenance_stub("prime", production_invalid=True))
    reg.set(provenance_not_used("bigmhc_im"))
    summary = reg.tool_summary_fields(("prime", "bigmhc_im"))
    assert summary["prime_source"] == "fixture"
    assert summary["prime_status"] == "invalid_for_production"
    assert summary["bigmhc_im_status"] == "not_used"
