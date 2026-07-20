"""Driver-gene relevance lookup (scoring audit fix #2).

Previously ``driver_relevance`` was hardcoded to ``"0.0"`` in every adapter
unless the upstream pVACseq/VCF row happened to carry a non-standard
``driver_relevance`` / ``Driver Relevance`` / ``driver`` column (which almost
never exists in practice). Because ``driver_relevance`` carries the
second-largest weight in the event_score formula (0.18, see
``profiles/default.toml`` -> ``[event_weights]``), this made ~1/5 of the
event-level ranking score a constant that never actually discriminated
between candidates.

This module ships a small, self-contained reference table
(``resources/driver_genes.tsv``) distilled from widely-cited, publicly
published pan-cancer driver gene compendia (Vogelstein et al. 2013 Science
"Cancer Genome Landscapes"; Bailey et al. 2018 Cell "Comprehensive
Characterization of Cancer Driver Genes and Mutations"; COSMIC Cancer Gene
Census Tier 1 genes that are consistently reported across both). It is
intentionally a *simplified* illustrative list (~100 genes), not a
substitute for a maintained, licensed source.

For production use, swap ``resources/driver_genes.tsv`` for an export from
an actively curated database such as OncoKB's cancerGeneList, COSMIC Cancer
Gene Census, or CIViC -- the lookup function and TSV schema
(``gene\tcategory\tdriver_relevance\tsource_note``) are deliberately kept
generic so any of those exports can be dropped in without code changes, as
long as they are reduced to the same four columns.
"""
from __future__ import annotations

from pathlib import Path
from typing import Mapping

from .utils import read_tsv, to_float

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DRIVER_GENE_PATH = ROOT / "resources" / "driver_genes.tsv"

# Default score assigned to a gene that is *not* found in the reference
# table. This is deliberately a neutral mid-range value rather than 0.0:
# "not in our simplified ~100-gene list" is evidence of absence, not
# evidence of a confirmed non-driver -- collapsing "unknown" and "confirmed
# passenger" into the same 0.0 score is exactly the bug this fix removes.
DEFAULT_UNKNOWN_GENE_RELEVANCE = 0.3

# OncoKB cancerGeneList.tsv boolean/curation columns used to derive a
# consensus-strength score. "OncoKB Annotated" is included alongside the six
# underlying source databases because OncoKB's own curation counts as a
# resource in its own right in the exported "# of occurrence" column.
_ONCOKB_SOURCE_COLUMNS = [
    "OncoKB Annotated", "MSK-IMPACT", "MSK-HEME",
    "FOUNDATION ONE", "FOUNDATION ONE HEME", "Vogelstein",
]
_ONCOKB_COSMIC_PREFIX = "COSMIC CGC"  # column name carries a version suffix, e.g. "COSMIC CGC (v99)"

_CACHE: dict[str, dict[str, float]] = {}


def _parse_simple_driver_tsv(rows: list[dict[str, str]]) -> dict[str, float]:
    """Bundled-resource format: gene\tcategory\tdriver_relevance\tsource_note."""
    table: dict[str, float] = {}
    for row in rows:
        gene = (row.get("gene") or "").strip().upper()
        if not gene:
            continue
        table[gene] = to_float(row.get("driver_relevance"), DEFAULT_UNKNOWN_GENE_RELEVANCE)
    return table


def _parse_oncokb_gene_list(rows: list[dict[str, str]]) -> dict[str, float]:
    """Official OncoKB ``cancerGeneList.tsv`` export.

    Columns of interest: ``Hugo Symbol``, ``Gene Type`` (ONCOGENE/TSG),
    ``OncoKB Annotated`` plus the five other curated-resource Yes/No columns
    (MSK-IMPACT, MSK-HEME, FOUNDATION ONE, FOUNDATION ONE HEME, Vogelstein,
    COSMIC CGC <version>). ``driver_relevance`` is derived from how many of
    those independent, curated resources list the gene as cancer-relevant
    (consensus strength), rather than trusted as a single binary flag:

        driver_relevance = 0.3 + 0.65 * (n_sources_agreeing / n_source_columns)

    A gene absent from every source column but still present in the file
    (edge case) lands at 0.3 -- the same neutral value used for genes not in
    the table at all -- while a gene confirmed across every source (e.g.
    TP53, KRAS) lands at 0.95. The occurrence count is recomputed directly
    from the Yes/No columns rather than trusted from the file's own
    "# of occurrence" column, so this stays correct even if that column's
    exact definition changes between OncoKB export versions.
    """
    table: dict[str, float] = {}
    for row in rows:
        gene = (row.get("Hugo Symbol") or "").strip().upper()
        if not gene:
            continue
        cosmic_col = next((k for k in row if k.startswith(_ONCOKB_COSMIC_PREFIX)), None)
        source_cols = list(_ONCOKB_SOURCE_COLUMNS) + ([cosmic_col] if cosmic_col else [])
        if not source_cols:
            continue
        hits = sum(1 for c in source_cols if str(row.get(c, "")).strip().lower() == "yes")
        score = 0.3 + 0.65 * (hits / len(source_cols))
        table[gene] = max(0.0, min(1.0, score))
    return table


def load_driver_gene_table(path: str | Path | None = None) -> dict[str, float]:
    """Load {gene_symbol: driver_relevance} from a driver-gene reference file.

    Auto-detects two formats so an official OncoKB ``cancerGeneList.tsv``
    export can be dropped in without any manual conversion step:
      - the bundled simplified resource format (header starts with "gene");
      - the official OncoKB gene list (header includes "Hugo Symbol" and
        "Gene Type" -- see https://www.oncokb.org/cancer-genes).

    Results are cached per resolved path so repeated per-row lookups during
    adapter parsing don't re-read the file from disk.
    """
    p = Path(path) if path else DEFAULT_DRIVER_GENE_PATH
    key = str(p)
    if key in _CACHE:
        return _CACHE[key]
    table: dict[str, float] = {}
    if p.exists():
        rows = read_tsv(p)
        if rows and "Hugo Symbol" in rows[0] and "Gene Type" in rows[0]:
            table = _parse_oncokb_gene_list(rows)
        else:
            table = _parse_simple_driver_tsv(rows)
    _CACHE[key] = table
    return table


def lookup_driver_relevance(
    gene: str,
    profile: Mapping | None = None,
    table: Mapping[str, float] | None = None,
) -> float:
    """Return a driver_relevance score in [0, 1] for ``gene``.

    Resolution order:
      1. If the profile disables driver-gene lookup
         (``[driver_genes] enabled = false``), always return the neutral
         unknown-gene default without consulting the table.
      2. Look the gene up (case-insensitively) in the reference table
         (bundled ``resources/driver_genes.tsv`` by default, or a path
         supplied via ``profile["driver_genes"]["reference_path"]``).
      3. If not found, return the profile-configurable
         ``unknown_gene_relevance`` (default 0.3).
    """
    cfg = dict((profile or {}).get("driver_genes", {}))
    if not cfg.get("enabled", True):
        return float(cfg.get("unknown_gene_relevance", DEFAULT_UNKNOWN_GENE_RELEVANCE))
    if table is None:
        table = load_driver_gene_table(cfg.get("reference_path"))
    unknown_default = float(cfg.get("unknown_gene_relevance", DEFAULT_UNKNOWN_GENE_RELEVANCE))

    # Fusion gene names use the HGVS "::" convention (e.g. "EWSR1::WT1", see
    # resources/normal_expression.example.tsv) -> look up each partner
    # independently and take the higher relevance; a fusion driven by a known
    # driver on either side is itself a driver event.
    raw = (gene or "").strip().upper()
    if "::" in raw:
        parts = [p for p in raw.split("::") if p]
        if parts:
            return max(table.get(p, unknown_default) for p in parts)
    return float(table.get(raw, unknown_default))
