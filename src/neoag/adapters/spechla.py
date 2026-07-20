"""Parse SpecHLA merge.hla.copy.txt into neoag hla_loh.tsv."""

from __future__ import annotations

from pathlib import Path

from .lohhla import write_hla_loh_evidence
from .peptide_input import normalize_hla_allele
from ..evidence_provenance import provenance_from_file
from ..utils import first, read_tsv


def normalize_spechla_allele(locus: str, raw: str) -> str:
    token = raw.strip()
    if not token or token in {".", "-", "NA", "NAN", "NONE"}:
        return ""
    if "*" in token:
        gene_part, fields = token.split("*", 1)
    else:
        gene_part, fields = locus.strip(), token
    gene_part = gene_part.upper().removeprefix("HLA-")
    field_parts = [p for p in fields.split(":") if p]
    if not field_parts:
        return ""
    allele = f"HLA-{gene_part}*{field_parts[0]}"
    if len(field_parts) >= 2:
        allele = f"{allele}:{field_parts[1]}"
    if gene_part in {"A", "B", "C"}:
        return normalize_hla_allele(allele)
    return allele


def _loh_positive(raw: str) -> bool:
    return raw.strip().upper() in {"Y", "YES", "LOH", "LOSS", "LOST", "1", "TRUE"}


def parse_spechla_loh_merge(path: str | Path) -> list[dict[str, str]]:
    """Parse SpecHLA merge.hla.copy.txt to per-allele LOH status."""
    allele_status: dict[str, str] = {}
    for row in read_tsv(path):
        locus = first(row, ["HLA", "hla", "locus"], "")
        allele1 = normalize_spechla_allele(locus, first(row, ["Allele1", "allele1"], ""))
        allele2 = normalize_spechla_allele(locus, first(row, ["Allele2", "allele2"], ""))
        lost = normalize_spechla_allele(locus, first(row, ["LossHLA", "loss_hla", "LossAllele"], ""))
        kept = normalize_spechla_allele(locus, first(row, ["KeptHLA", "kept_hla"], ""))
        alleles = [a for a in (allele1, allele2) if a]
        if _loh_positive(first(row, ["LOH", "loh", "loh_status"], "")):
            if lost:
                allele_status[lost] = "loh"
            for allele in alleles:
                if allele != lost:
                    allele_status.setdefault(allele, "no")
            if kept and kept != lost:
                allele_status[kept] = "no"
        else:
            for allele in alleles:
                allele_status.setdefault(allele, "no")
    return [{"hla_allele": allele, "loh_status": status} for allele, status in sorted(allele_status.items())]


def write_spechla_hla_loh_evidence(
    out_path: str | Path,
    rows: list[dict[str, str]],
    *,
    source_path: str | Path,
) -> None:
    write_hla_loh_evidence(
        out_path,
        rows,
        provenance=provenance_from_file("spechla", source_path, mode="converted"),
    )
