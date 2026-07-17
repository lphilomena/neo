"""Tumor-specificity scoring via GTEx median-TPM percentile rank (scoring
audit fix #3).

Previously ``tumor_specificity`` was a hardcoded ``"0.7"`` default in every
adapter with zero computation behind it anywhere in the codebase, despite
carrying weight in both the L3 composite (``l3_weights.tumor_specificity``)
and the event_score formula (``event_weights.tumor_specificity``).

This module implements the approach requested for the first pass: for a
given gene, compare the tumor's own expression (TPM) for that gene against
the distribution of GTEx median TPM values for that same gene across normal
tissues, and take the percentile rank of the tumor value within that
distribution.

    tumor_specificity = (# GTEx tissues with median TPM <= tumor TPM)
                         / (# GTEx tissues with data for this gene)

Interpretation: if the tumor's expression of a gene is higher than (nearly)
every normal tissue's baseline expression of that same gene, the tumor
looks like an outlier relative to normal biology for this gene -> high
percentile -> treated as more "tumor-specific" in the sense that few/no
normal tissues already present this gene's products at a comparable level.
If several normal tissues already express the gene at or above the tumor's
level, the percentile -- and therefore the specificity score -- is low,
flagging a higher on-target/off-tumor risk.

The bundled reference matrix (``resources/gtex_median_tpm.tsv``) is a small,
ILLUSTRATIVE, hand-assembled approximation covering a representative subset
of GTEx tissues and a limited gene set (see the file header for the
production-swap instructions: replace it with a pivoted export of GTEx's
official "gene median TPM" bulk-tissue-expression download, keeping the same
``gene`` + arbitrary tissue-columns schema). Genes not present in the
bundled/production matrix fall back to a neutral, profile-configurable
default rather than a value that implies the gene is confirmed
tissue-restricted or confirmed ubiquitous.
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Mapping

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_GTEX_PATH = ROOT / "resources" / "gtex_median_tpm.tsv"

# Returned when the gene isn't present in the GTEx reference matrix at all --
# i.e. genuinely unknown, not "confirmed broadly expressed" (which would
# argue for a low score) or "confirmed tissue-restricted" (which would argue
# for a high score). 0.5 is the neutral midpoint.
DEFAULT_UNKNOWN_GENE_SPECIFICITY = 0.5

_CACHE: dict[str, dict[str, list[float]]] = {}


def load_gtex_median_tpm(path: str | Path | None = None) -> dict[str, list[float]]:
    """Load {gene_symbol: [median_tpm_per_tissue, ...]} from a GTEx-style TSV.

    Lines starting with ``#`` are treated as documentation/comments and
    skipped, so the bundled resource file can carry a human-readable header
    explaining its provenance without breaking the parser.
    """
    p = Path(path) if path else DEFAULT_GTEX_PATH
    key = str(p)
    if key in _CACHE:
        return _CACHE[key]
    table: dict[str, list[float]] = {}
    if p.exists():
        with p.open("r", encoding="utf-8", newline="") as fh:
            lines = [ln for ln in fh if not ln.lstrip().startswith("#") and ln.strip()]
        reader = csv.DictReader(lines, delimiter="\t")
        tissue_cols = [c for c in (reader.fieldnames or []) if c and c != "gene"]
        for row in reader:
            gene = (row.get("gene") or "").strip().upper()
            if not gene:
                continue
            values: list[float] = []
            for col in tissue_cols:
                raw = row.get(col)
                if raw in (None, ""):
                    continue
                try:
                    values.append(float(raw))
                except ValueError:
                    continue
            if values:
                table[gene] = values
    _CACHE[key] = table
    return table


def _percentile_rank(value: float, distribution: list[float]) -> float:
    """Fraction of ``distribution`` that is <= value, in [0, 1]."""
    if not distribution:
        return DEFAULT_UNKNOWN_GENE_SPECIFICITY
    n_le = sum(1 for x in distribution if x <= value)
    return n_le / len(distribution)


def compute_tumor_specificity(
    gene: str,
    tumor_tpm: float,
    profile: Mapping | None = None,
    table: Mapping[str, list[float]] | None = None,
) -> float:
    """Return a tumor_specificity score in [0, 1] for ``gene``.

    ``tumor_tpm`` should be the tumor/event-level expression already
    computed for this gene (``event_expression``), so this needs no separate
    RNA-seq quantification step -- it only needs the GTEx comparison table.
    """
    cfg = dict((profile or {}).get("tumor_specificity_gtex", {}))
    if not cfg.get("enabled", True):
        return float(cfg.get("unknown_gene_specificity", DEFAULT_UNKNOWN_GENE_SPECIFICITY))
    if table is None:
        table = load_gtex_median_tpm(cfg.get("reference_path"))
    unknown_default = float(cfg.get("unknown_gene_specificity", DEFAULT_UNKNOWN_GENE_SPECIFICITY))
    distribution = table.get((gene or "").strip().upper())
    if not distribution:
        return unknown_default
    return _percentile_rank(float(tumor_tpm), distribution)
