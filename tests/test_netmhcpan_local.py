from neoag_v03.adapters.netmhcpan import (
    parse_netmhcpan_local_stdout,
    write_netmhcpan_standard_xls,
    _hla_with_star,
)
from neoag_v03.tools.runner import resolve_netmhcpan_backend


SAMPLE_OUT = """
   1    HLA-A*02:06          SLFTLFSV SLF-TLFSV  0  0  0  3  1          SLFTLFSV         PEPLIST 0.092698    3.136 0.565176    1.113   110.47  0.00000
"""


def test_hla_with_star():
    assert _hla_with_star("HLA-A02:06") == "HLA-A*02:06"
    assert _hla_with_star("HLA-A*02:06") == "HLA-A*02:06"
    assert _hla_with_star("HLA-A*0206") == "HLA-A*02:06"


def test_parse_netmhcpan_local_stdout():
    rows = parse_netmhcpan_local_stdout(SAMPLE_OUT)
    assert len(rows) == 1
    assert rows[0]["Peptide"] == "SLFTLFSV"
    assert rows[0]["HLA"] == "HLA-A*02:06"
    assert rows[0]["Score_BA"] == "110.47"
    assert rows[0]["%Rank_BA"] == "1.113"


def test_resolve_netmhcpan_backend_default(monkeypatch):
    monkeypatch.delenv("NEOAG_NETMHCPAN_BACKEND", raising=False)
    assert resolve_netmhcpan_backend() == "local"


def test_resolve_netmhcpan_backend_iedb(monkeypatch):
    monkeypatch.setenv("NEOAG_NETMHCPAN_BACKEND", "iedb")
    assert resolve_netmhcpan_backend() == "iedb"


def test_write_netmhcpan_standard_xls(tmp_path):
    rows = parse_netmhcpan_local_stdout(SAMPLE_OUT)
    out = tmp_path / "netmhcpan.xls"
    write_netmhcpan_standard_xls(out, rows)
    text = out.read_text(encoding="utf-8")
    assert "SLFTLFSV" in text
    assert "HLA-A*02:06" in text
