from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from neoag.schemas import EVENT_FIELDS, PEPTIDE_FIELDS
from neoag.utils import read_tsv, write_tsv


ROOT = Path(__file__).resolve().parents[1]


def row(fields, **values):
    return {field: str(values.get(field, "")) for field in fields}


def test_filter_keeps_only_complete_exact_supported_snaf_candidates(tmp_path):
    events = tmp_path / "raw_events.tsv"
    peptides = tmp_path / "raw_peptides.tsv"
    consensus = tmp_path / "splice_consensus.tsv"
    write_tsv(events, [row(EVENT_FIELDS, event_id="E1"), row(EVENT_FIELDS, event_id="E2")], EVENT_FIELDS)
    write_tsv(
        peptides,
        [
            row(PEPTIDE_FIELDS, peptide_id="P1", event_id="E1", peptide="ACDEFGHIK", hla_allele="HLA-A*02:01", source_tool="SNAF", binding_rank="0.8"),
            row(PEPTIDE_FIELDS, peptide_id="P2", event_id="E2", peptide="ACDEFGHIK", hla_allele="HLA-A*02:01", source_tool="SNAF", binding_rank="0.5"),
            row(PEPTIDE_FIELDS, peptide_id="P3", event_id="E1", peptide="ACD", hla_allele="HLA-A*02:01", source_tool="SNAF", binding_rank="0.2"),
        ],
        PEPTIDE_FIELDS,
    )
    write_tsv(
        consensus,
        [
            {"event_id": "E1", "status": "CROSS_DOMAIN_CONFIRMED_EXACT_JUNCTION"},
            {"event_id": "E2", "status": "RESOLVED_SOURCE_ONLY"},
        ],
        ["event_id", "status"],
    )
    outdir = tmp_path / "selected"
    subprocess.run(
        [sys.executable, str(ROOT / "scripts/filter_splice_production_candidates.py"), "--events", str(events), "--peptides", str(peptides), "--consensus", str(consensus), "--outdir", str(outdir)],
        cwd=ROOT,
        check=True,
    )
    selected = read_tsv(outdir / "raw_peptides.tsv")
    rejected = read_tsv(outdir / "rejected_peptides.tsv")
    summary = json.loads((outdir / "production_filter_summary.json").read_text())
    assert [record["peptide_id"] for record in selected] == ["P1"]
    assert {record["peptide_id"] for record in rejected} == {"P2", "P3"}
    assert summary["selected_events"] == 1
    assert summary["selected_peptides"] == 1
