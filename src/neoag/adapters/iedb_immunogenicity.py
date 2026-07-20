"""IEDB Class I pMHC immunogenicity predictor (v3.0 algorithm, local)."""

from __future__ import annotations

from pathlib import Path

from ..evidence_provenance import ProvenanceRecord, provenance_derived, without_provenance, write_evidence_tsv
from ..schemas import IEDB_IMMUNOGENICITY_FIELDS

# Allele-specific anchor masks from IEDB Immunogenicity tool 3.0.
ALLELE_MASK: dict[str, str] = {
    "H-2-Db": "2,5,9",
    "H-2-Dd": "2,3,5",
    "H-2-Kb": "2,3,9",
    "H-2-Kd": "2,5,9",
    "H-2-Kk": "2,8,9",
    "H-2-Ld": "2,5,9",
    "HLA-A0101": "2,3,9",
    "HLA-A0201": "1,2,9",
    "HLA-A0202": "1,2,9",
    "HLA-A0203": "1,2,9",
    "HLA-A0206": "1,2,9",
    "HLA-A0211": "1,2,9",
    "HLA-A0301": "1,2,9",
    "HLA-A1101": "1,2,9",
    "HLA-A2301": "2,7,9",
    "HLA-A2402": "2,7,9",
    "HLA-A2601": "1,2,9",
    "HLA-A2902": "2,7,9",
    "HLA-A3001": "1,3,9",
    "HLA-A3002": "2,7,9",
    "HLA-A3101": "1,2,9",
    "HLA-A3201": "1,2,9",
    "HLA-A3301": "1,2,9",
    "HLA-A6801": "1,2,9",
    "HLA-A6802": "1,2,9",
    "HLA-A6901": "1,2,9",
    "HLA-B0702": "1,2,9",
    "HLA-B0801": "2,5,9",
    "HLA-B1501": "1,2,9",
    "HLA-B1502": "1,2,9",
    "HLA-B1801": "1,2,9",
    "HLA-B2705": "2,3,9",
    "HLA-B3501": "1,2,9",
    "HLA-B3901": "1,2,9",
    "HLA-B4001": "1,2,9",
    "HLA-B4002": "1,2,9",
    "HLA-B4402": "2,3,9",
    "HLA-B4403": "2,3,9",
    "HLA-B4501": "1,2,9",
    "HLA-B4601": "1,2,9",
    "HLA-B5101": "1,2,9",
    "HLA-B5301": "1,2,9",
    "HLA-B5401": "1,2,9",
    "HLA-B5701": "1,2,9",
    "HLA-B5801": "1,2,9",
}

IMMUNOSCALE: dict[str, float] = {
    "A": 0.127,
    "C": -0.175,
    "D": 0.072,
    "E": 0.325,
    "F": 0.380,
    "G": 0.110,
    "H": 0.105,
    "I": 0.432,
    "K": -0.700,
    "L": -0.036,
    "M": -0.570,
    "N": -0.021,
    "P": -0.036,
    "Q": -0.376,
    "R": 0.168,
    "S": -0.537,
    "T": 0.126,
    "V": 0.134,
    "W": 0.719,
    "Y": -0.012,
}
IMMUNOWEIGHT_BASE = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]


def iedb_allele_key(hla: str) -> str:
    s = hla.strip().upper()
    if s.startswith("HLA-"):
        s = s[4:]
    elif s.startswith("HLA"):
        s = s[3:]
    return "HLA-" + s.replace("*", "").replace(":", "")


def _mask_positions(hla: str) -> list[int]:
    key = iedb_allele_key(hla)
    if key in ALLELE_MASK:
        return [int(x) - 1 for x in ALLELE_MASK[key].split(",")]
    return [0, 1, -1]


def predict_immunogenicity(peptide: str, hla: str) -> float:
    """Return IEDB immunogenicity score (>0 favors response)."""
    peptide = peptide.strip().upper()
    mask_num = _mask_positions(hla)
    if -1 in mask_num:
        mask_num = [m if m >= 0 else len(peptide) - 1 for m in mask_num]
    peplen = len(peptide)
    if peplen > 9:
        pepweight = IMMUNOWEIGHT_BASE[:5] + ((peplen - 9) * [0.30]) + IMMUNOWEIGHT_BASE[5:]
    else:
        pepweight = IMMUNOWEIGHT_BASE
    score = 0.0
    count = 0
    for aa in peptide:
        if count in mask_num:
            count += 1
            continue
        score += pepweight[count] * IMMUNOSCALE[aa]
        count += 1
    return round(score, 5)


def predict_pairs(pairs: list[tuple[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for peptide, hla in pairs:
        rows.append({
            "peptide": peptide,
            "hla_allele": hla,
            "iedb_immunogenicity_score": f"{predict_immunogenicity(peptide, hla):.5f}",
            "iedb_allele_key": iedb_allele_key(hla),
            "iedb_mask_known": "yes" if iedb_allele_key(hla) in ALLELE_MASK else "no",
        })
    return rows


def write_immunogenicity_evidence(
    path: str | Path,
    rows: list[dict[str, str]],
    provenance: ProvenanceRecord | None = None,
) -> None:
    prov = provenance or provenance_derived("iedb", path, upstream="iedb_immunogenicity_model")
    write_evidence_tsv(path, rows, without_provenance(IEDB_IMMUNOGENICITY_FIELDS), prov)
