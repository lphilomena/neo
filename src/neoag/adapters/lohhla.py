"""Parse LOHHLA HLA loss prediction outputs into neoag hla_loh.tsv."""

from __future__ import annotations

import csv
import re
from pathlib import Path

from .peptide_input import normalize_hla_allele
from ..utils import read_tsv
from ..evidence_provenance import ProvenanceRecord, provenance_from_file, without_provenance, write_evidence_tsv
from ..schemas import HLA_LOH_EVIDENCE_FIELDS


def normalize_lohhla_allele(raw: str) -> str:
    """Convert LOHHLA token hla_a_24_02_01_01 → HLA-A*24:02:01."""
    token = raw.strip().lower()
    if not token or token in {"na", "nan", "none", "-"}:
        return ""
    if token.startswith("hla-"):
        return normalize_hla_allele(token.upper())
    m = re.match(r"^hla_([a-z])_(\d{2})_(\d{2})(?:_\d+)*$", token)
    if not m:
        return normalize_hla_allele(token.upper().replace("_", "*"))
    gene, d1, d2 = m.group(1).upper(), m.group(2), m.group(3)
    return normalize_hla_allele(f"HLA-{gene}*{d1}:{d2}")


def _read_rows(path: Path) -> list[dict[str, str]]:
    suffix = path.suffix.lower()
    if suffix in {".xls", ".xlsx", ".txt", ".tsv"} or "HLAlossPrediction" in path.name:
        with path.open(encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            return [dict(row) for row in reader]
    return read_tsv(path)


def parse_lohhla_prediction(path: str | Path) -> list[dict[str, str]]:
    """Parse LOHHLA *HLAlossPrediction_CI* output to per-allele LOH status."""
    path = Path(path)
    rows_in = _read_rows(path)
    loh_alleles: set[str] = set()
    all_alleles: set[str] = set()

    for row in rows_in:
        type_cols = [
            row.get("HLA_A_type1", ""),
            row.get("HLA_A_type2", ""),
            row.get("HLA_B_type1", ""),
            row.get("HLA_B_type2", ""),
            row.get("HLA_C_type1", ""),
            row.get("HLA_C_type2", ""),
        ]
        for col in ("HLA_type1", "HLA_type2", "type1", "type2", "allele1", "allele2"):
            type_cols.extend([row.get(col, "")])
        for raw in type_cols:
            allele = normalize_lohhla_allele(raw)
            if allele:
                all_alleles.add(allele)

        for col in ("LossAllele", "loss_allele", "HLA_loss", "lost_allele"):
            lost = normalize_lohhla_allele(row.get(col, ""))
            if lost:
                loh_alleles.add(lost)

    if not all_alleles and loh_alleles:
        all_alleles = set(loh_alleles)

    if not all_alleles:
        return []

    out: list[dict[str, str]] = []
    for allele in sorted(all_alleles):
        status = "loh" if allele in loh_alleles else "no"
        out.append({"hla_allele": allele, "loh_status": status})
    return out


def write_hla_loh_evidence(
    path: str | Path,
    rows: list[dict[str, str]],
    provenance: ProvenanceRecord | None = None,
) -> None:
    prov = provenance or provenance_from_file("lohhla", path, mode="converted")
    write_evidence_tsv(path, rows, without_provenance(HLA_LOH_EVIDENCE_FIELDS), prov)
