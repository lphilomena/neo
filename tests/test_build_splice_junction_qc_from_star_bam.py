from pathlib import Path

from scripts.build_splice_junction_qc_from_star_bam import (
    event_junction,
    load_crossvalidated_keys,
)


def test_event_junction_prefers_canonical_identifier():
    row = {
        "canonical_junction_id": "SJ|GRCh38|chr10|101035672|101036006|+",
        "event_name": "chr10:101035672",
    }
    assert event_junction(row) == ("10", 101035672, 101036006, "+")


def test_caller_consensus_requires_crossvalidated_snaf_and_splicemutr(tmp_path: Path):
    consensus = tmp_path / "splice_consensus.tsv"
    consensus.write_text(
        "event_id\tcanonical_junction_id\tsupport_tools\tstatus\n"
        "SJ|GRCh38|chr1|10|20|+\tSJ|GRCh38|chr1|10|20|+\tSNAF;SpliceMutr\tCROSS_DOMAIN_CONFIRMED_EXACT_JUNCTION\n"
        "SJ|GRCh38|chr2|30|40|-\tSJ|GRCh38|chr2|30|40|-\tSTAR;SpliceMutr\tFAIL_ARTIFACT\n",
        encoding="utf-8",
    )
    snaf, splicemutr = load_crossvalidated_keys(consensus)
    assert "SJ|GRCh38|chr1|10|20|+" in snaf
    assert "SJ|GRCh38|chr1|10|20|+" in splicemutr
    assert "SJ|GRCh38|chr2|30|40|-" not in snaf
    assert "SJ|GRCh38|chr2|30|40|-" not in splicemutr
