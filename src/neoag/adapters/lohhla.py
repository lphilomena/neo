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


def _float(row: dict[str, str], *keys: str) -> float | None:
    for key in keys:
        value = str(row.get(key, "") or "").strip()
        if not value:
            continue
        try:
            return float(value)
        except ValueError:
            continue
    return None


def _raw(row: dict[str, str], *keys: str) -> str:
    return next((str(row.get(key, "") or "") for key in keys if str(row.get(key, "") or "")), "")


def parse_lohhla_prediction(
    path: str | Path,
    *,
    max_copy_number: float = 0.5,
    max_pval: float = 0.01,
) -> list[dict[str, str]]:
    """Parse LOHHLA output using its paired-P and allele-copy-number rule.

    ``LossAllele`` identifies the lower-copy candidate; it is not by itself a
    positive LOH call. A candidate is LOST only when its BAF-informed copy
    number is below ``max_copy_number`` and the paired P value is below
    ``max_pval``. Raw statistics are retained for audit and cross-tool review.
    """
    path = Path(path)
    rows_in = _read_rows(path)
    evidence: dict[str, dict[str, str]] = {}

    for row in rows_in:
        pairs: list[tuple[str, str]] = []
        for locus in ("A", "B", "C"):
            first = normalize_lohhla_allele(_raw(row, f"HLA_{locus}_type1"))
            second = normalize_lohhla_allele(_raw(row, f"HLA_{locus}_type2"))
            if first or second:
                pairs.append((first, second))
        if not pairs:
            pairs.append((
                normalize_lohhla_allele(_raw(row, "HLA_type1", "type1", "allele1")),
                normalize_lohhla_allele(_raw(row, "HLA_type2", "type2", "allele2")),
            ))

        lost_raw = _raw(row, "LossAllele", "loss_allele", "HLA_loss", "lost_allele")
        kept_raw = _raw(row, "KeptAllele", "kept_allele", "HLA_kept")
        lost = normalize_lohhla_allele(lost_raw)
        pval = _float(row, "PVal", "Pval", "pval")
        common = {
            "call_rule": f"candidate_loss_copy_number_with_baf<{max_copy_number:g} AND paired_pval<{max_pval:g}",
            "lohhla_pval": _raw(row, "PVal", "Pval", "pval"),
            "lohhla_unpaired_pval": _raw(row, "UnPairedPval", "UnpairedPval", "unpaired_pval"),
            "lohhla_pval_unique": _raw(row, "PVal_unique", "Pval_unique", "pval_unique"),
            "lohhla_unpaired_pval_unique": _raw(row, "UnPairedPval_unique", "UnpairedPval_unique"),
            "lohhla_mismatch_sites": _raw(row, "numMisMatchSitesCov", "mismatch_sites"),
            "lohhla_prop_supportive_sites": _raw(row, "propSupportiveSites", "supportive_sites"),
            "lohhla_loss_allele_raw": lost_raw,
            "lohhla_kept_allele_raw": kept_raw,
        }
        for allele1, allele2 in pairs:
            for index, allele in enumerate((allele1, allele2), 1):
                if not allele:
                    continue
                cn_raw = _raw(row, f"HLA_type{index}copyNum_withBAF")
                cn = _float(row, f"HLA_type{index}copyNum_withBAF")
                qc: list[str] = []
                if cn is None or pval is None:
                    status = "unassessed"
                    qc.append("MISSING_PVAL_OR_COPY_NUMBER")
                elif allele == lost and cn < max_copy_number and pval < max_pval:
                    status = "loh"
                else:
                    status = "no"
                if cn is not None and cn < 0:
                    qc.append("NEGATIVE_COPY_NUMBER_LOW_PURITY_FIT")
                evidence[allele] = {
                    "hla_allele": allele,
                    "loh_status": status,
                    **common,
                    "call_qc": ";".join(qc) or "PASS",
                    "lohhla_copy_number_with_baf": cn_raw,
                    "lohhla_cn_lower": _raw(row, f"HLA_type{index}copyNum_withBAF_lower"),
                    "lohhla_cn_upper": _raw(row, f"HLA_type{index}copyNum_withBAF_upper"),
                }
    return [evidence[allele] for allele in sorted(evidence)]


def write_hla_loh_evidence(
    path: str | Path,
    rows: list[dict[str, str]],
    provenance: ProvenanceRecord | None = None,
) -> None:
    prov = provenance or provenance_from_file("lohhla", path, mode="converted")
    write_evidence_tsv(path, rows, without_provenance(HLA_LOH_EVIDENCE_FIELDS), prov)
