"""Parse OptiType HLA typing output CSV and normalise allele names."""

from __future__ import annotations
from pathlib import Path
from .peptide_input import normalize_hla_allele


def parse_optitype_result(csv_path: str | Path) -> list[str]:
    """Parse OptiType _result.tsv and return normalised HLA alleles.

    OptiType output columns: A1, A2, B1, B2, C1, C2, nof_reads, obj
    We extract the best (first) row and collect unique non-empty alleles.
    """
    import csv

    p = Path(csv_path)
    if not p.is_file():
        raise FileNotFoundError(f"OptiType result not found: {p}")

    with open(p, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            alleles = []
            for col in ("A1", "A2", "B1", "B2", "C1", "C2"):
                val = (row.get(col) or "").strip()
                if val and val not in alleles:
                    norm = normalize_hla_allele(val)
                    if norm and norm not in alleles:
                        alleles.append(norm)
            return sorted(alleles)

    raise ValueError(f"No HLA typing results in {p}")


def write_optitype_hla_alleles(
    csv_path: str | Path,
    out_path: str | Path,
    sample_id: str = "",
) -> None:
    """Parse OptiType result CSV and write a plain-text HLA alleles file.

    One allele per line, compatible with tools that accept --hla-alleles lists.
    """
    alleles = parse_optitype_result(csv_path)
    Path(out_path).write_text("\n".join(alleles) + "\n")
