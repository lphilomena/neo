from __future__ import annotations

from pathlib import Path

from neoag_v03.adapters.spechla import parse_spechla_loh_merge
from neoag_v03.tools.postprocess import spechla_to_hla_loh_tsv

MERGE_NO_LOH = """\
Sample\tHLA\tAllele1\tAllele2\tcopyratio\tKeptHLA\tLossHLA\tFreq1\tFreq2\tPurity\tHet_num\tLOH
SAMPLE_NO_LOH\tA\tA*30:01:01:01\tA*02:06:01:01\t1:1\tA*30:01:01:01\tA*02:06:01:01\t0.509\t0.491\t0.283\t95\tN
SAMPLE_NO_LOH\tB\tB*48:01:01:01\tB*13:02:01:01\t1:1\tB*48:01:01:01\tB*13:02:01:01\t0.522\t0.478\t0.283\t61\tN
"""

MERGE_WITH_LOH = """\
Sample\tHLA\tAllele1\tAllele2\tcopyratio\tKeptHLA\tLossHLA\tFreq1\tFreq2\tPurity\tHet_num\tLOH
S1\tA\tA*02:01:01:01\tA*11:01:01:01\t1:0\tA*02:01:01:01\tA*11:01:01:01\t0.8\t0.2\t0.5\t40\tY
"""


def test_parse_spechla_no_loh(tmp_path: Path):
    path = tmp_path / "merge.hla.copy.txt"
    path.write_text(MERGE_NO_LOH, encoding="utf-8")
    rows = parse_spechla_loh_merge(path)
    by_allele = {r["hla_allele"]: r["loh_status"] for r in rows}
    assert by_allele == {
        "HLA-A*30:01": "no",
        "HLA-A*02:06": "no",
        "HLA-B*48:01": "no",
        "HLA-B*13:02": "no",
    }


def test_parse_spechla_with_loh(tmp_path: Path):
    path = tmp_path / "merge.hla.copy.txt"
    path.write_text(MERGE_WITH_LOH, encoding="utf-8")
    rows = parse_spechla_loh_merge(path)
    by_allele = {r["hla_allele"]: r["loh_status"] for r in rows}
    assert by_allele["HLA-A*11:01"] == "loh"
    assert by_allele["HLA-A*02:01"] == "no"


def test_spechla_to_hla_loh_tsv(tmp_path: Path):
    path = tmp_path / "merge.hla.copy.txt"
    path.write_text(MERGE_NO_LOH, encoding="utf-8")
    out = tmp_path / "hla_loh.tsv"
    spechla_to_hla_loh_tsv(path, out)
    text = out.read_text(encoding="utf-8")
    assert "HLA-A*30:01\tno" in text
    assert "HLA-B*13:02\tno" in text
    assert "spechla" in text


def test_upstream_spechla_merge_passthrough(tmp_path: Path):
    from neoag_v03.tools.upstream import run_upstream

    merge = tmp_path / "merge.hla.copy.txt"
    merge.write_text(MERGE_NO_LOH, encoding="utf-8")
    cfg_path = tmp_path / "run.toml"
    cfg_path.write_text(
        f"""
[sample]
id = "T1"
profile = "default"

[tools]
stub = false
enabled = []

[inputs]
spechla_merge = "{merge.as_posix()}"
""",
        encoding="utf-8",
    )
    outs = run_upstream(cfg_path, tmp_path / "upstream")
    assert "hla_loh" in outs
    text = Path(outs["hla_loh"]).read_text(encoding="utf-8")
    assert "HLA-A*02:06\tno" in text
    assert "spechla" in text
