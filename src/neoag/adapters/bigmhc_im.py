"""BigMHC_IM immunogenicity adapter."""

from __future__ import annotations

import csv
from pathlib import Path

from ..utils import read_csv
from ..evidence_provenance import ProvenanceRecord, provenance_from_file, provenance_stub, without_provenance, write_evidence_tsv
from ..schemas import BIGMHC_IM_EVIDENCE_FIELDS


def _pair_key(peptide: str, hla: str) -> tuple[str, str]:
    return peptide.strip().upper(), hla.strip()


def predict_pair_stub(peptide: str, hla: str) -> str:
    h = abs(hash(f"BigMHC_IM_{peptide}_{hla}"))
    return f"{(h % 10000) / 10000.0:.6f}"


def predict_pairs(pairs: list[tuple[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for peptide, hla in pairs:
        rows.append({
            "peptide": peptide,
            "hla_allele": hla,
            "bigmhc_im_score": predict_pair_stub(peptide, hla),
        })
    return rows


def _score_column(fieldnames: list[str] | None) -> str | None:
    if not fieldnames:
        return None
    for name in fieldnames:
        low = name.lower()
        if "bigmhc" in low and "im" in low:
            return name
    for name in fieldnames:
        low = name.lower()
        if low in {"bigmhc_im", "bigmhc_im_score", "prediction", "pred", "score"}:
            return name
    return None


def parse_bigmhc_im(path: str | Path, sample_id: str = "") -> list[dict[str, str]]:
    p = Path(path)
    if p.suffix.lower() == ".csv" or p.suffix.lower() == ".prd":
        return _parse_bigmhc_csv(p, sample_id)
    return _parse_bigmhc_tsv(p, sample_id)


def _parse_bigmhc_csv(path: Path, sample_id: str) -> list[dict[str, str]]:
    raw = read_csv(path)
    if not raw:
        return []
    score_col = _score_column(list(raw[0].keys()))
    rows: list[dict[str, str]] = []
    for row in raw:
        peptide = (row.get("pep") or row.get("peptide") or row.get("Peptide") or "").strip().upper()
        hla = (row.get("mhc") or row.get("hla_allele") or row.get("HLA") or "").strip()
        if not peptide:
            continue
        score = row.get(score_col or "", "") if score_col else ""
        rows.append({
            "sample_id": sample_id,
            "peptide": peptide,
            "hla_allele": hla,
            "bigmhc_im_score": str(score),
            "source_file": str(path),
        })
    return rows


def _parse_bigmhc_tsv(path: Path, sample_id: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        score_col = _score_column(reader.fieldnames)
        for row in reader:
            peptide = (row.get("peptide") or row.get("pep") or "").strip().upper()
            hla = (row.get("hla_allele") or row.get("mhc") or "").strip()
            if not peptide:
                continue
            score = row.get(score_col or "bigmhc_im_score", "")
            rows.append({
                "sample_id": sample_id,
                "peptide": peptide,
                "hla_allele": hla,
                "bigmhc_im_score": str(score),
                "source_file": str(path),
            })
    return rows


def write_bigmhc_im_evidence(
    path: str | Path,
    rows: list[dict[str, str]],
    provenance: ProvenanceRecord | None = None,
) -> None:
    src = rows[0].get("source_file") if rows else path
    if src == "stub":
        prov = provenance or provenance_stub("bigmhc_im", production_invalid=True)
    else:
        prov = provenance or provenance_from_file("bigmhc_im", src, mode="tool_run")
    write_evidence_tsv(path, rows, without_provenance(BIGMHC_IM_EVIDENCE_FIELDS), prov)


def bigmhc_by_pair(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    out: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        key = _pair_key(row.get("peptide", ""), row.get("hla_allele", ""))
        out[key] = row
    return out
