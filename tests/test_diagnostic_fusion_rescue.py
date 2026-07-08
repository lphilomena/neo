from pathlib import Path

from neoag_v03.adapters.diagnostic_fusion_rescue import (
    build_diagnostic_fusion_rescue,
    diagnostic_rescue_rows_from_easyfuse,
    infer_unfiltered_easyfuse_path,
    normalize_fusion_label,
)
from neoag_v03.utils import read_tsv


def _write_easyfuse_pair(tmp_path: Path) -> tuple[Path, Path]:
    header = (
        "BPID;Fusion_Gene;Breakpoint1;Breakpoint2;FTID;prediction_class;prediction_prob;"
        "fusioncatcher_detected;star_detected;arriba_detected;tool_count;fusioncatcher_junc;"
        "fusioncatcher_span;fusioncatcher_anch;frame;type;neo_peptide_sequence;neo_peptide_sequence_bp"
    )
    pass_csv = tmp_path / "fusions.pass.csv"
    pass_csv.write_text(
        header + "\n"
        "bp_other;BCR_ABL1;22:1:+;9:2:-;tx1;positive;0.95;1;1;0;2;12;8;18;in_frame;trans;ABCDEFGHIJK;5\n",
        encoding="utf-8",
    )
    all_csv = tmp_path / "fusions.csv"
    all_csv.write_text(
        header + "\n"
        "bp_ews;EWSR1_WT1;22:29291599:+;11:32392064:-;tx_ews;;0;1;0;0;1;37;382;20;in_frame;trans;MTEYKLVVVGAG;6\n"
        "bp_other;BCR_ABL1;22:1:+;9:2:-;tx1;positive;0.95;1;1;0;2;12;8;18;in_frame;trans;ABCDEFGHIJK;5\n",
        encoding="utf-8",
    )
    return pass_csv, all_csv


def test_normalize_fusion_label_matches_common_separators():
    assert normalize_fusion_label("EWSR1::WT1") == "EWSR1_WT1"
    assert normalize_fusion_label("EWSR1_WT1") == "EWSR1_WT1"


def test_infer_unfiltered_easyfuse_path_from_pass(tmp_path):
    pass_csv, all_csv = _write_easyfuse_pair(tmp_path)
    assert infer_unfiltered_easyfuse_path(pass_csv) == all_csv


def test_diagnostic_rescue_keeps_whitelisted_fusion_outside_pass(tmp_path):
    pass_csv, all_csv = _write_easyfuse_pair(tmp_path)
    pass_keys = {("BCR_ABL1", "22:1:+", "9:2:-")}
    rows = diagnostic_rescue_rows_from_easyfuse(
        all_csv,
        sample_id="S1",
        whitelist=["EWSR1::WT1"],
        pass_keys=pass_keys,
    )
    assert len(rows) == 1
    row = rows[0]
    assert row["fusion_gene_normalized"] == "EWSR1_WT1"
    assert row["fusioncatcher_detected"] == "1"
    assert row["star_detected"] == "0"
    assert row["easyfuse_pass_status"] == "not_in_pass_table"
    assert row["diagnostic_relevance"] == "diagnostic_fusion_evidence"
    assert row["peptide_generation_status"] == "available_not_generated_by_default"


def test_build_diagnostic_rescue_writes_sidecar_without_generating_peptides(tmp_path):
    pass_csv, _all_csv = _write_easyfuse_pair(tmp_path)
    out = tmp_path / "parsed" / "diagnostic_fusion_rescue.tsv"
    result = build_diagnostic_fusion_rescue(
        {"inputs": {"easyfuse_pass_csv": str(pass_csv)}},
        sample_id="S1",
        out_path=out,
    )
    assert result["diagnostic_fusion_rescue"] == str(out)
    assert result["diagnostic_fusion_rescue_rows"] == "1"
    assert result["diagnostic_fusion_rescue_generate_peptides"] == "false"
    rows = read_tsv(out)
    assert len(rows) == 1
    assert rows[0]["fusion_gene_normalized"] == "EWSR1_WT1"
    assert rows[0]["peptide_generation_status"].endswith("not_generated_by_default")
