from __future__ import annotations

import importlib.util
from pathlib import Path

from neoag.utils import read_tsv


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "prepare_precomputed_splicemutr.py"
SPEC = importlib.util.spec_from_file_location("prepare_precomputed_splicemutr", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_converter_does_not_mislabel_translated_protein_as_short_peptide(tmp_path: Path):
    source = tmp_path / "data_splicemutr_all_pep.txt"
    source.write_text(
        "chr\tstart\tend\tstrand\tgene\ttx_id\tjuncs\tpeptide\n"
        "chr1\t100\t200\t+\tGENE1\tTX1\tchr1:100:200:+\tABCDEFGHIJK*\n",
        encoding="utf-8",
    )
    presentation = tmp_path / "netmhcpan.xls"
    presentation.write_text("result\n", encoding="utf-8")
    output = tmp_path / "splicemutr_candidates.tsv"

    assert MODULE.convert(source, output, "S1", presentation) == 1
    row = read_tsv(output)[0]
    assert row["source_junction_id"] == "chr1:100:200:+"
    assert row["peptide"] == ""
    assert row["crosses_junction"] == "UNASSESSED"
    assert row["translated_protein_length"] == "11"
    assert row["presentation_status"] == "ASSESSED_GLOBAL_NOT_JUNCTION_PEPTIDE_MAPPED"
