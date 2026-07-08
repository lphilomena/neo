"""Parse FACETS purity and copy-number segment outputs."""

from __future__ import annotations

from pathlib import Path

from ..utils import first, read_tsv
from ..evidence_provenance import ProvenanceRecord, provenance_from_file, without_provenance, write_evidence_tsv
from ..schemas import CNV_SEGMENT_FIELDS, PURITY_EVIDENCE_FIELDS


def parse_facets_purity(path: str | Path, sample_id: str = "") -> list[dict[str, str]]:
    rows = read_tsv(path)
    purity = "0.5"
    for row in rows:
        val = first(row, ["purity", "Purity", "purity_fraction", "TumorFraction", "cf"], "")
        if val:
            purity = str(val)
            break
    sid = sample_id or (first(rows[0], ["sample_id"], "") if rows else "") or "SAMPLE001"
    return [{"sample_id": sid, "purity": purity}]


def parse_facets_cncf(path: str | Path) -> list[dict[str, str]]:
    rows_out: list[dict[str, str]] = []
    for row in read_tsv(path):
        chrom = first(row, ["chrom", "Chromosome", "chr", "chromosome"], "")
        if chrom and not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        start = first(row, ["start", "Start", "startpos"], "")
        end = first(row, ["end", "End", "endpos"], "")
        cn = first(row, ["total_cn", "tcn.em", "tcn", "copy_number", "total_copy_number"], "")
        if chrom and start and end and cn:
            rows_out.append(
                {
                    "chrom": chrom,
                    "start": str(int(float(start))),
                    "end": str(int(float(end))),
                    "total_cn": f"{float(cn):.4f}",
                }
            )
    return rows_out


def write_purity_evidence(
    path: str | Path,
    rows: list[dict[str, str]],
    provenance: ProvenanceRecord | None = None,
) -> None:
    prov = provenance or provenance_from_file("facets", path, mode="converted")
    write_evidence_tsv(path, rows, without_provenance(PURITY_EVIDENCE_FIELDS), prov)


def write_cnv_segments(
    path: str | Path,
    rows: list[dict[str, str]],
    provenance: ProvenanceRecord | None = None,
) -> None:
    prov = provenance or provenance_from_file("facets", path, mode="converted")
    write_evidence_tsv(path, rows, without_provenance(CNV_SEGMENT_FIELDS), prov)
