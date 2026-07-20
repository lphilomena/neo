"""Parse ASCAT purity/ploidy summary and segment outputs."""

from __future__ import annotations

from pathlib import Path

from ..utils import first, read_tsv
from ..evidence_provenance import ProvenanceRecord, provenance_from_file, without_provenance, write_evidence_tsv
from ..schemas import CNV_SEGMENT_FIELDS, PURITY_EVIDENCE_FIELDS


def parse_ascat_purity(path: str | Path, sample_id: str = "") -> list[dict[str, str]]:
    rows = read_tsv(path)
    purity = "0.5"
    for row in rows:
        val = first(row, ["purity", "rho", "aberrantcellfraction"], "")
        if val and ";" not in str(val):
            purity = str(val)
            break
        if val and ";" in str(val):
            purity = str(val).split(";")[0]
            break
    sid = sample_id or (first(rows[0], ["sample_id"], "") if rows else "") or "SAMPLE001"
    return [{"sample_id": sid, "purity": purity}]


def parse_ascat_segments(path: str | Path) -> list[dict[str, str]]:
    rows_out: list[dict[str, str]] = []
    for row in read_tsv(path):
        chrom = first(row, ["chr", "chrom", "Chromosome", "chromosome"], "")
        if chrom and not str(chrom).startswith("chr"):
            chrom = f"chr{chrom}"
        start = first(row, ["startpos", "start", "Start"], "")
        end = first(row, ["endpos", "end", "End"], "")
        cn = first(row, ["nMajor", "nMinor", "copy_number", "total_cn", "tcn"], "")
        if cn and ";" in str(cn):
            parts = [p for p in str(cn).split(";") if p]
            try:
                cn = str(sum(float(p) for p in parts))
            except ValueError:
                cn = parts[0]
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
    prov = provenance or provenance_from_file("ascat", path, mode="converted")
    write_evidence_tsv(path, rows, without_provenance(PURITY_EVIDENCE_FIELDS), prov)


def write_cnv_segments(
    path: str | Path,
    rows: list[dict[str, str]],
    provenance: ProvenanceRecord | None = None,
) -> None:
    prov = provenance or provenance_from_file("ascat", path, mode="converted")
    write_evidence_tsv(path, rows, without_provenance(CNV_SEGMENT_FIELDS), prov)
