"""Tests for EasyFuse fusion discovery adapter."""

from pathlib import Path

import pytest

from neoag.adapters.easyfuse_adapter import (
    EasyFuseFilterConfig,
    _anchor_size,
    _tool_junction_reads,
    build_easyfuse_candidate_union,
    easyfuse_row_to_event,
    filter_easyfuse_row,
    is_easyfuse_table,
    parse_easyfuse,
    read_easyfuse_table,
)
from neoag.evidence_layer import build_standard_evidence_layer
from neoag.input_router import build_raw_intermediates

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "data/fixtures/easyfuse_fusions.pass.tsv"
FIXTURE_V2 = ROOT / "data/fixtures/easyfuse_fusions.v2.pass.csv"
REAL_PASS = ROOT / "results/easyfuse/sample/fusions.pass.csv"


def test_is_easyfuse_table(tmp_path):
    assert is_easyfuse_table(FIXTURE)
    assert not is_easyfuse_table(ROOT / "data/fixtures/pvacseq_aggregated.tsv")

    semicolon = tmp_path / "fusions.pass.csv"
    semicolon.write_text(
        "BPID;Fusion_Gene;neo_peptide_sequence\nbp1;GENE1_GENE2;PEPTIDE\n",
        encoding="utf-8",
    )
    assert is_easyfuse_table(semicolon)


def test_filter_easyfuse_row():
    rows = __import__("neoag.utils", fromlist=["read_tsv"]).read_tsv(FIXTURE)
    by_gene = {r["Fusion_Gene"]: r for r in rows}

    ok, reason = filter_easyfuse_row(by_gene["EWSR1_FLI1"])
    assert ok and reason == "pass"

    ok, reason = filter_easyfuse_row(by_gene["CBFB_MYH11"])
    assert not ok and reason.startswith("readthrough_type:")

    ok, reason = filter_easyfuse_row(by_gene["NCF1_NCF1"])
    assert not ok and reason == "prediction_class_not_positive"


def test_easyfuse_pass_plus_internal_high_confidence_union(tmp_path):
    header = (
        "BPID;context_sequence_id;FTID;Fusion_Gene;frame;type;neo_peptide_sequence;"
        "ft_junc_cnt;ft_anch_cnt;starfusion_detected;fusioncatcher_detected;"
        "arriba_detected;prediction_class;prediction_prob\n"
    )
    pass_path = tmp_path / "fusions.pass.csv"
    all_path = tmp_path / "fusions.csv"
    pass_row = "bp1;c1;t1;A_B;in_frame;trans;ABCDEFGHIJK;3;10;1;0;0;positive;0.8\n"
    rescue_row = "bp2;c2;t2;C_D;in_frame;trans;ABCDEFGHIJK;6;20;0;0;1;;0.1\n"
    weak_row = "bp3;c3;t3;E_F;in_frame;trans;ABCDEFGHIJK;3;10;0;0;1;;0.1\n"
    pass_path.write_text(header + pass_row, encoding="utf-8")
    all_path.write_text(header + pass_row + rescue_row + weak_row, encoding="utf-8")

    rows, summary = build_easyfuse_candidate_union(pass_path, all_path)

    assert len(rows) == 2
    assert summary["internal_high_confidence_added"] == 1
    by_bpid = {row["BPID"]: row for row in rows}
    assert by_bpid["bp1"]["openneo_union_source"] == "EASYFUSE_PASS"
    assert by_bpid["bp2"]["openneo_union_source"] == "INTERNAL_CALLER_HIGH_CONFIDENCE"
    ok, reason = filter_easyfuse_row(by_bpid["bp2"])
    assert ok and reason == "pass"
    event, evidence = easyfuse_row_to_event(
        by_bpid["bp2"], sample_id="S1", profile_name="default", source_path=all_path,
    )
    assert event is not None
    assert event["candidate_union_source"] == "INTERNAL_CALLER_HIGH_CONFIDENCE"
    assert evidence["candidate_union_source"] == "INTERNAL_CALLER_HIGH_CONFIDENCE"


def test_parse_easyfuse_events_and_catalog():
    parsed = parse_easyfuse(FIXTURE, "S1", "default", lengths=(8, 9))
    assert len(parsed["events"]) == 2
    genes = {e["gene"] for e in parsed["events"]}
    assert "EWSR1::FLI1" in genes
    assert "BCR::ABL1" in genes
    assert all(e["source"].startswith("easyfuse:") for e in parsed["events"])
    assert all(e["mutation_source"] == "SV" for e in parsed["events"])
    assert parsed["peptides"] == []
    assert len(parsed["catalog_rows"]) > 0
    assert len(parsed["fusion_evidence"]) == 4
    assert len(parsed["filter_qc"]) == 4
    assert len(parsed["collapse_qc"]) == 2
    passed = [r for r in parsed["fusion_evidence"] if r["filter_status"] == "pass"]
    assert len(passed) == 2


def test_v2_requant_column_mapping():
    """EasyFuse v2 exports *_cnt_best columns; adapter must not read all-zero legacy names."""
    rows = read_easyfuse_table(FIXTURE_V2)
    assert rows
    row = rows[0]
    assert _tool_junction_reads(row) >= 3
    assert _anchor_size(row) >= 10
    ok, reason = filter_easyfuse_row(row)
    assert ok and reason == "pass"


def test_parse_easyfuse_v2_fixture():
    parsed = parse_easyfuse(FIXTURE_V2, "S1", "default")
    assert len(parsed["events"]) == 3
    assert parsed["peptides"] == []
    assert len(parsed["catalog_rows"]) > 0
    passed = [r for r in parsed["fusion_evidence"] if r["filter_status"] == "pass"]
    assert len(passed) == 3
    assert all(int(r["rna_junction_reads"]) >= 3 for r in passed)
    assert all(int(r["anchor_size"]) >= 10 for r in passed)


def test_parse_real_easyfuse_pass_catalog_not_empty():
    if not REAL_PASS.is_file():
        return
    parsed = parse_easyfuse(REAL_PASS, "FP500004780_L01_203", "default")
    assert len(parsed["events"]) > 0
    assert len(parsed["fusion_evidence"]) == 2686
    passed = [r for r in parsed["fusion_evidence"] if r["filter_status"] == "pass"]
    assert len(passed) > 100


def test_build_intermediates_easyfuse_mode(tmp_path):
    cfg = {
        "sample": {"id": "EF1", "profile": "default"},
        "inputs": {
            "entry_mode": "fusion",
            "easyfuse_pass_csv": str(FIXTURE),
        },
    }
    paths = build_raw_intermediates(cfg, tmp_path / "layer", root=ROOT)
    assert Path(paths["raw_events"]).is_file()
    assert Path(paths["raw_peptides"]).is_file()
    assert Path(paths["fusion_evidence"]).is_file()

    events = __import__("neoag.utils", fromlist=["read_tsv"]).read_tsv(paths["raw_events"])
    assert len(events) == 2


def test_build_intermediates_easyfuse_expands_hla_and_requires_peptides(tmp_path):
    hla = ["HLA-A*02:01", "HLA-B*07:02"]
    cfg = {
        "sample": {"id": "EF1", "profile": "default"},
        "inputs": {
            "entry_mode": "fusion",
            "easyfuse_pass_csv": str(FIXTURE),
            "hla_alleles": hla,
            "require_nonempty_peptides": True,
        },
    }
    paths = build_raw_intermediates(cfg, tmp_path / "with_hla", root=ROOT)
    peptides = __import__("neoag.utils", fromlist=["read_tsv"]).read_tsv(paths["raw_peptides"])
    assert peptides
    assert {row["hla_allele"] for row in peptides} == set(hla)

    cfg["inputs"]["hla_alleles"] = []
    with pytest.raises(ValueError, match="candidate event.*no peptide-HLA rows"):
        build_raw_intermediates(cfg, tmp_path / "without_hla", root=ROOT)


def test_evidence_layer_uses_fusion_evidence(tmp_path):
    layer = tmp_path / "layer"
    cfg = {
        "sample": {"id": "EF1", "profile": "default"},
        "inputs": {"entry_mode": "fusion", "easyfuse_tsv": str(FIXTURE)},
    }
    build_raw_intermediates(cfg, layer, root=ROOT)
    evidence = build_standard_evidence_layer(layer, "default", sample_id="EF1")
    rna_rows = __import__("neoag.utils", fromlist=["read_tsv"]).read_tsv(evidence["rna_junction_evidence"])
    ews = [r for r in rna_rows if "EWSR1" in r.get("gene", "")]
    assert ews and int(ews[0]["junction_reads"]) >= 12
