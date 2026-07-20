"""PRIME class I immunogenicity adapter."""

from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import Mapping

from ..utils import read_tsv
from ..evidence_provenance import ProvenanceRecord, provenance_from_file, without_provenance, write_evidence_tsv
from ..schemas import PRIME_EVIDENCE_FIELDS


def _pair_key(peptide: str, hla: str) -> tuple[str, str]:
    return peptide.strip().upper(), hla.strip()


def prime_parallel_jobs(default: int = 4) -> int:
    """Worker count for peptide-level PRIME batching (NEOAG_PRIME_JOBS)."""
    raw = os.environ.get("NEOAG_PRIME_JOBS", str(default)).strip()
    try:
        n = int(raw)
    except ValueError:
        n = default
    return max(1, n)


def split_peptide_chunks(peptides: list[str], n_jobs: int) -> list[list[str]]:
    """Split sorted unique peptides into up to n_jobs contiguous batches."""
    if not peptides:
        return []
    n_jobs = max(1, min(n_jobs, len(peptides)))
    if n_jobs == 1:
        return [list(peptides)]
    size = len(peptides)
    base, extra = divmod(size, n_jobs)
    chunks: list[list[str]] = []
    start = 0
    for i in range(n_jobs):
        end = start + base + (1 if i < extra else 0)
        chunks.append(peptides[start:end])
        start = end
    return chunks


def prime_allele_tag(hla: str) -> str:
    """Convert HLA-A*02:01 to PRIME column tag A0201."""
    s = hla.upper().strip()
    if s.startswith("HLA-"):
        s = s[4:]
    return s.replace("*", "").replace(":", "")


def extract_prime_pair(row: Mapping[str, str], peptide: str, hla: str) -> tuple[str, str]:
    """Pull per-allele Score/%Rank from PRIME wide-format output."""
    tag = prime_allele_tag(hla)
    score = row.get(f"Score_{tag}") or row.get(f"score_{tag}") or ""
    rank = row.get(f"%Rank_{tag}") or row.get(f"Rank_{tag}") or ""
    if not score and tag:
        for key, val in row.items():
            low = key.lower()
            if low == f"score_{tag.lower()}" or low.endswith(f"score_{tag.lower()}"):
                score = val
            if low == f"%rank_{tag.lower()}" or low.endswith(f"rank_{tag.lower()}"):
                rank = val
    if not score:
        score = row.get("Score_bestAllele") or row.get("score_bestallele") or row.get("prime_score") or ""
    if not rank:
        rank = row.get("%Rank_bestAllele") or row.get("prime_rank") or row.get("%Rank") or ""
    return str(score), str(rank)


def read_prime_wide_rows(path: str | Path) -> list[dict[str, str]]:
    """Read PRIME wide TSV, skipping leading `#` comment lines before the header."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    header_idx = next(
        (i for i, line in enumerate(lines) if line.startswith("Peptide\t")),
        None,
    )
    if header_idx is None:
        return []
    reader = csv.DictReader(lines[header_idx:], delimiter="\t")
    rows: list[dict[str, str]] = []
    for row in reader:
        peptide = (row.get("Peptide") or row.get("peptide") or "").strip().upper()
        if peptide and not peptide.startswith("#"):
            rows.append(row)
    return rows


def predict_pair_stub(peptide: str, hla: str) -> tuple[str, str]:
    h = abs(hash(f"PRIME_{peptide}_{hla}"))
    score = f"{(h % 10000) / 10000.0:.6f}"
    rank = f"{0.1 + (h % 990) / 10.0:.2f}"
    return score, rank


def predict_pairs(pairs: list[tuple[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for peptide, hla in pairs:
        score, rank = predict_pair_stub(peptide, hla)
        rows.append({
            "peptide": peptide,
            "hla_allele": hla,
            "prime_score": score,
            "prime_rank": rank,
        })
    return rows


def parse_prime(path: str | Path, sample_id: str = "") -> list[dict[str, str]]:
    """Parse PRIME output or normalized evidence TSV."""
    p = Path(path)
    if p.suffix.lower() in {".tsv", ".txt"}:
        return _parse_prime_tsv(p, sample_id)
    return _parse_prime_tsv(p, sample_id)


def _parse_prime_tsv(path: Path, sample_id: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if any(key in path.read_text(encoding="utf-8", errors="ignore")[:500] for key in ("Score_A0201", "Score_bestAllele", "Peptide\t%Rank")):
        for row in read_prime_wide_rows(path):
            peptide = (row.get("Peptide") or row.get("peptide") or "").strip().upper()
            hla = (
                row.get("hla_allele")
                or row.get("HLA")
                or row.get("HLA Allele")
                or row.get("BestAllele")
                or row.get("Best allele")
                or ""
            ).strip()
            score, rank = extract_prime_pair(row, peptide, hla)
            rows.append({
                "sample_id": sample_id,
                "peptide": peptide,
                "hla_allele": hla,
                "prime_score": str(score),
                "prime_rank": str(rank),
                "source_file": str(path),
            })
        if rows:
            return rows

    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("Peptide", "").startswith("#") or row.get("Peptide", "").startswith("Output"):
                continue
            peptide = (
                row.get("peptide")
                or row.get("Peptide")
                or row.get("MT Epitope Seq")
                or ""
            ).strip().upper()
            if not peptide or peptide.startswith("#"):
                continue
            hla = (
                row.get("hla_allele")
                or row.get("HLA")
                or row.get("HLA Allele")
                or row.get("BestAllele")
                or row.get("Best allele")
                or ""
            ).strip()
            score, rank = extract_prime_pair(row, peptide, hla)
            if not score and not rank:
                score = (
                    row.get("prime_score")
                    or row.get("Score_bestAllele")
                    or row.get("PRIME Score")
                    or row.get("score")
                    or ""
                )
                rank = (
                    row.get("prime_rank")
                    or row.get("%Rank_bestAllele")
                    or row.get("PRIME %Rank")
                    or row.get("%Rank")
                    or row.get("Lowest %Rank")
                    or ""
                )
            rows.append({
                "sample_id": sample_id,
                "peptide": peptide,
                "hla_allele": hla,
                "prime_score": str(score),
                "prime_rank": str(rank),
                "source_file": str(path),
            })
    return rows


def write_prime_evidence(
    path: str | Path,
    rows: list[dict[str, str]],
    provenance: ProvenanceRecord | None = None,
) -> None:
    src = rows[0].get("source_file") if rows else path
    if src == "stub":
        from ..evidence_provenance import provenance_stub
        prov = provenance or provenance_stub("prime", production_invalid=True)
    else:
        prov = provenance or provenance_from_file("prime", src, mode="tool_run")
    write_evidence_tsv(path, rows, without_provenance(PRIME_EVIDENCE_FIELDS), prov)


def prime_by_pair(rows: list[dict[str, str]]) -> dict[tuple[str, str], dict[str, str]]:
    out: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        key = _pair_key(row.get("peptide", ""), row.get("hla_allele", ""))
        out[key] = row
    return out
