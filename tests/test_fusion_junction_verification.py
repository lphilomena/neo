from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "build_fusion_caller_union", ROOT / "scripts/build_fusion_caller_union.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_event_junction_matching_is_breakpoint_specific_and_deduplicated(tmp_path: Path):
    junctions = tmp_path / "Chimeric.out.junction"
    junctions.write_text(
        "chr22\t101\t+\tchr11\t201\t-\t1\t0\t0\treadA\n"
        "chr22\t101\t+\tchr11\t201\t-\t1\t0\t0\treadA\n"
        "chr22\t101\t+\tchr11\t501\t-\t1\t0\t0\treadB\n",
        encoding="utf-8",
    )
    audit = [
        {"event_id": "E1", "left_breakpoint": "chr22:100", "right_breakpoint": "chr11:200"},
        {"event_id": "E2", "left_breakpoint": "11:500", "right_breakpoint": "22:100"},
    ]

    measurements, sidecar = MODULE.verify_event_junction_reads(
        audit, [junctions], None, tolerance=2
    )

    assert measurements["E1"]["verified_count"] == 1
    assert measurements["E2"]["verified_count"] == 1
    assert measurements["E1"]["status"] == "STAR_JUNCTION_VERIFIED"
    assert {row["event_id"] for row in sidecar} == {"E1", "E2"}


def test_verified_count_replaces_caller_count_without_losing_it():
    rows = [{"event_id": "E1", "rna_junction_reads": "7"}]
    MODULE.apply_junction_measurements(
        rows,
        {
            "E1": {
                "verified_count": 1,
                "status": "BAM_VERIFIED",
                "method": "star_chimeric_breakpoint_plus_bam_qname",
                "source": "junctions;rna.bam",
            }
        },
    )

    assert rows[0]["provided_rna_junction_reads"] == "7"
    assert rows[0]["rna_junction_reads"] == "1"
    assert rows[0]["junction_match_status"] == "BAM_VERIFIED"
