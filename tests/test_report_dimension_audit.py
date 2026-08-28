from __future__ import annotations

from pathlib import Path

from neoag.report_dimensions import audit_report_dimensions, load_report_dimension_map
from neoag.utils import read_tsv, write_tsv


def test_report_dimension_map_loads():
    mapping = load_report_dimension_map()
    assert mapping["metadata"]["schema_version"] == "report-dimension-map-v1"
    assert "source_chain_confidence" in mapping["dimensions"]
    assert mapping["dimensions"]["functional_validation"]["layer"] == "post_ranking_validation"


def test_report_dimension_audit_tracks_present_and_missing_fields(tmp_path: Path):
    input_tsv = tmp_path / "ranked.tsv"
    output_tsv = tmp_path / "audit.tsv"
    write_tsv(input_tsv, [{
        "peptide_id": "P1",
        "source_chain_confidence_tier": "C2",
        "source_chain_confidence_label": "STRONG_COMPUTATIONAL_COMPLETE_SOURCE_CHAIN",
        "source_chain_reason_codes": "SC_ORTHOGONAL_NOT_PERFORMED",
        "source_chain_rule_version": "source-chain-v1.0",
        "safety_status": "PASS",
    }])
    result = audit_report_dimensions(input_tsv, output_tsv)
    assert result["dimensions"] >= 10
    rows = {row["dimension_key"]: row for row in read_tsv(output_tsv)}
    assert rows["source_chain_confidence"]["status"] == "PASS"
    assert rows["functional_validation"]["status"] == "MISSING"
