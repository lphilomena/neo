from __future__ import annotations
from pathlib import Path
from ..utils import read_tsv, first, write_tsv
from ..adapters.peptide_input import (
    pair_key,
    unique_peptide_hla_pairs_from_table,
)


def unique_peptide_hla_pairs(raw_peptides_tsv: str | Path) -> list[tuple[str, str]]:
    """Return unique (peptide, HLA) pairs from raw_peptides.tsv or flexible input tables."""
    path = Path(raw_peptides_tsv)
    if not path.is_file():
        return []

    # Standard pipeline output already uses peptide/hla_allele columns.
    rows = read_tsv(path)
    if rows and first(rows[0], ["peptide"], "") and first(rows[0], ["hla_allele", "HLA Allele"], ""):
        seen: set[tuple[str, str]] = set()
        pairs: list[tuple[str, str]] = []
        for row in rows:
            hla = first(row, ["HLA Allele", "hla_allele", "Allele", "allele"], "")
            if not hla:
                continue
            for pep in (
                first(row, ["MT Epitope Seq", "MT Epitope", "peptide", "Peptide"], ""),
                first(row, ["WT Epitope Seq", "WT Epitope", "wildtype_peptide"], ""),
            ):
                if not pep or str(pep).strip().upper() == "NA":
                    continue
                key = pair_key(pep, hla)
                if key in seen:
                    continue
                seen.add(key)
                pairs.append(key)
        return pairs

    return unique_peptide_hla_pairs_from_table(path)


def write_peptide_fasta(pairs: list[tuple[str, str]], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for i, (pep, hla) in enumerate(pairs):
        header = f">{i}_{pep}_{hla.replace('*', '').replace(':', '')}"
        lines.append(header)
        lines.append(pep)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_mhcflurry_peptides_csv(pairs: list[tuple[str, str]], path: str | Path) -> None:
    import csv
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["peptide", "allele"])
        for pep, hla in pairs:
            writer.writerow([pep, hla])


def netmhcpan_pmhc_allele(hla: str) -> str:
    """Format allele for NetMHCpan PEPTIDEMHC input (HLA-A02:06)."""
    s = str(hla or "").strip().upper()
    if s.startswith("HLA-"):
        s = s[4:]
    s = s.replace("*", "")
    return f"HLA-{s}"


def write_netmhcpan_pmhc_input(pairs: list[tuple[str, str]], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{pep} {netmhcpan_pmhc_allele(hla)} 0" for pep, hla in pairs]
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def netmhcpan_allele_string(alleles: list[str]) -> str:
    """Convert HLA-A*02:01 style alleles to NetMHCpan -a format (HLA-A0201)."""
    out: list[str] = []
    for allele in alleles:
        s = str(allele or "").strip().upper()
        if s.startswith("HLA-"):
            s = s[4:]
        s = s.replace("*", "").replace(":", "")
        out.append(f"HLA-{s}")
    return ",".join(out)


def mhcflurry_allele_list(alleles: list[str]) -> str:
    return ",".join(alleles)
