from __future__ import annotations

from pathlib import Path

from neoag.hla_loh_crosscheck import crosscheck_hla_loh, write_hla_loh_crosscheck
from neoag.utils import read_tsv


def test_hla_loh_crosscheck_consensus_discordant_and_single_tool(tmp_path: Path):
    lohhla = tmp_path / "lohhla.hla_loh.tsv"
    lohhla.write_text(
        "hla_allele\tloh_status\tconfidence\n"
        "HLA-A*02:01\tloh\t0.01\n"
        "HLA-B*07:02\tloh\t0.02\n"
        "HLA-C*07:02\tloh\t0.03\n",
        encoding="utf-8",
    )
    spechla = tmp_path / "spechla.hla_loh.tsv"
    spechla.write_text(
        "hla_allele\tloh_status\tconfidence\n"
        "HLA-A*02:01\tloh\tY\n"
        "HLA-B*07:02\tno\tN\n",
        encoding="utf-8",
    )

    rows = crosscheck_hla_loh(lohhla_hla_loh=lohhla, spechla_hla_loh=spechla)
    by = {r["hla_allele"]: r for r in rows}
    assert by["HLA-A*02:01"]["crosscheck_status"] == "CONSENSUS_LOH"
    assert by["HLA-B*07:02"]["crosscheck_status"] == "DISCORDANT"
    assert by["HLA-C*07:02"]["crosscheck_status"] == "SINGLE_TOOL_LOH"

    out = tmp_path / "hla_loh.crosscheck.tsv"
    consensus = tmp_path / "hla_loh.consensus.tsv"
    write_hla_loh_crosscheck(out, lohhla_hla_loh=lohhla, spechla_hla_loh=spechla, consensus_out=consensus)
    consensus_rows = read_tsv(consensus)
    assert {r["hla_allele"] for r in consensus_rows} == {"HLA-A*02:01", "HLA-C*07:02"}


def test_hla_loh_crosscheck_strict_consensus_only(tmp_path: Path):
    lohhla = tmp_path / "lohhla.hla_loh.tsv"
    lohhla.write_text("hla_allele\tloh_status\nHLA-A*02:01\tloh\n", encoding="utf-8")
    spechla = tmp_path / "spechla.hla_loh.tsv"
    spechla.write_text("hla_allele\tloh_status\nHLA-B*07:02\tno\n", encoding="utf-8")
    consensus = tmp_path / "hla_loh.consensus.tsv"
    write_hla_loh_crosscheck(
        tmp_path / "cross.tsv",
        lohhla_hla_loh=lohhla,
        spechla_hla_loh=spechla,
        consensus_out=consensus,
        include_single_tool=False,
    )
    assert read_tsv(consensus) == []
