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


def parse_spechla_loh_merge(path: str | Path, *, min_het: float = 5) -> list[dict[str, str]]:
    """Parse SpecHLA merge.hla.copy.txt to per-allele LOH status."""
    evidence: dict[str, dict[str, str]] = {}
    for row in read_tsv(path):
        locus = first(row, ["HLA", "hla", "locus"], "")
        allele1 = normalize_spechla_allele(locus, first(row, ["Allele1", "allele1"], ""))
        allele2 = normalize_spechla_allele(locus, first(row, ["Allele2", "allele2"], ""))
        lost = normalize_spechla_allele(locus, first(row, ["LossHLA", "loss_hla", "LossAllele"], ""))
        kept = normalize_spechla_allele(locus, first(row, ["KeptHLA", "kept_hla"], ""))
        alleles = [a for a in (allele1, allele2) if a]
        het_raw = first(row, ["Het_num", "het_num", "heterozygous_snp_count"], "")
        try:
            informative = float(het_raw) >= min_het
        except (TypeError, ValueError):
            informative = True
        homogeneous = first(row, ["LossHLA", "loss_hla", "LossAllele"], "").strip().lower() == "homogeneous"
        loh_raw = first(row, ["LOH", "loh", "loh_status"], "")
        positive = _loh_positive(loh_raw)
        common = {
            "call_rule": f"SpecHLA LOH flag with informative Het_num>={min_het:g}",
            "call_qc": "HOMOZYGOUS_OR_UNINFORMATIVE" if homogeneous or not informative else "PASS",
            "spechla_loh_raw": loh_raw,
            "spechla_copyratio": first(row, ["copyratio", "CopyRatio", "copy_ratio"], ""),
            "spechla_purity": first(row, ["purity", "Purity"], ""),
            "spechla_ploidy": first(row, ["ploidy", "Ploidy"], ""),
            "spechla_het_num": str(het_raw),
            "spechla_loss_hla_raw": first(row, ["LossHLA", "loss_hla", "LossAllele"], ""),
            "spechla_kept_hla_raw": first(row, ["KeptHLA", "kept_hla"], ""),
        }
        if homogeneous or not informative:
            statuses = {allele: "unassessed" for allele in alleles}
        elif positive:
            statuses = {allele: "loh" if allele == lost else "no" for allele in alleles}
            if lost:
                statuses[lost] = "loh"
            if kept and kept != lost:
                statuses[kept] = "no"
        else:
            statuses = {allele: "no" for allele in alleles}
        fallback_frequency = first(row, ["AlleleFreq", "allele_frequency", "frequency"], "")
        frequencies = {
            allele1: first(row, ["Freq1", "freq1", "Allele1Freq"], fallback_frequency),
            allele2: first(row, ["Freq2", "freq2", "Allele2Freq"], fallback_frequency),
        }
        for allele in alleles:
            evidence[allele] = {
                "hla_allele": allele,
                "loh_status": statuses.get(allele, "unassessed"),
                **common,
                "spechla_allele_frequency": frequencies.get(allele, fallback_frequency),
            }
    return [evidence[allele] for allele in sorted(evidence)]


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
