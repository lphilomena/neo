"""Cancer-gene context annotation without changing neoantigen scores."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping

from .utils import read_tsv


CANCER_GENE_FIELDS = [
    "cancer_gene_list_status",
    "cancer_gene_symbols",
    "cancer_gene_types",
    "cancer_driver_context",
    "oncokb_annotated",
    "cosmic_cgc_flag",
    "cancer_gene_source_count",
    "cancer_gene_sources",
    "cancer_gene_match_basis",
    "cancer_gene_context",
]

DRIVER_TYPES = {"ONCOGENE", "TSG", "ONCOGENE_AND_TSG"}
SOURCE_COLUMNS = (
    "OncoKB Annotated",
    "MSK-IMPACT",
    "MSK-HEME",
    "FOUNDATION ONE",
    "FOUNDATION ONE HEME",
    "Vogelstein",
    "COSMIC CGC (v99)",
)
GENE_SPLIT = re.compile(r"::|--|[;,/|]")


def _yes(value: Any) -> bool:
    return str(value or "").strip().lower() in {"yes", "true", "1", "y"}


def _source_count(value: Any) -> int:
    try:
        return int(float(str(value or "0").strip()))
    except ValueError:
        return 0


def _tokens(value: Any) -> list[str]:
    return [token.strip().upper() for token in GENE_SPLIT.split(str(value or "")) if token.strip()]


class CancerGeneIndex:
    def __init__(self, rows: Iterable[Mapping[str, Any]]):
        self.by_symbol: dict[str, dict[str, str]] = {}
        alias_candidates: dict[str, set[str]] = defaultdict(set)
        for raw in rows:
            row = {str(key): str(value or "").strip() for key, value in raw.items()}
            symbol = row.get("Hugo Symbol", "").upper()
            if not symbol:
                continue
            self.by_symbol[symbol] = row
            for alias in re.split(r"\s*,\s*", row.get("Gene Aliases", "")):
                if alias.strip():
                    alias_candidates[alias.strip().upper()].add(symbol)
        self.by_alias = {
            alias: next(iter(symbols))
            for alias, symbols in alias_candidates.items()
            if len(symbols) == 1 and alias not in self.by_symbol
        }

    @classmethod
    def from_tsv(cls, path: str | Path) -> "CancerGeneIndex":
        rows = read_tsv(path)
        if rows and "Hugo Symbol" not in rows[0]:
            raise ValueError("Cancer gene list must contain a 'Hugo Symbol' column")
        return cls(rows)

    def resolve(self, gene: str) -> tuple[list[dict[str, str]], str]:
        matched: list[dict[str, str]] = []
        bases: set[str] = set()
        seen: set[str] = set()
        for token in _tokens(gene):
            symbol = token if token in self.by_symbol else self.by_alias.get(token, "")
            if not symbol or symbol in seen:
                continue
            seen.add(symbol)
            matched.append(self.by_symbol[symbol])
            bases.add("symbol" if token == symbol else "alias")
        return matched, "+".join(sorted(bases))


def load_cancer_gene_index(path: str | Path | None) -> CancerGeneIndex | None:
    if not path:
        return None
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    return CancerGeneIndex.from_tsv(source)


def annotate_cancer_gene_context(row: Mapping[str, Any], index: CancerGeneIndex | None) -> dict[str, Any]:
    annotated = dict(row)
    if index is None:
        for field in CANCER_GENE_FIELDS:
            annotated.setdefault(field, "")
        annotated["cancer_gene_list_status"] = "UNASSESSED"
        annotated["cancer_driver_context"] = "UNASSESSED"
        return annotated

    matches, basis = index.resolve(str(row.get("gene") or row.get("event_name") or ""))
    if not matches:
        annotated.update({field: "" for field in CANCER_GENE_FIELDS})
        annotated["cancer_gene_list_status"] = "NOT_LISTED"
        annotated["cancer_driver_context"] = "NOT_LISTED"
        annotated["cancer_gene_context"] = "not present in supplied cancer-gene list"
        return annotated

    symbols = [record.get("Hugo Symbol", "") for record in matches]
    typed = [f"{record.get('Hugo Symbol', '')}:{record.get('Gene Type', '') or 'UNCLASSIFIED'}" for record in matches]
    driver_records = [record for record in matches if record.get("Gene Type", "").upper() in DRIVER_TYPES]
    sources = sorted({column for record in matches for column in SOURCE_COLUMNS if _yes(record.get(column))})
    source_counts = [_source_count(record.get("# of occurrence within resources (Column K-P)")) for record in matches]
    driver_status = "DRIVER_CONTEXT" if driver_records else "LISTED_NO_DRIVER_CLASS"
    annotated.update({
        "cancer_gene_list_status": "ANNOTATED",
        "cancer_gene_symbols": ";".join(symbols),
        "cancer_gene_types": ";".join(typed),
        "cancer_driver_context": driver_status,
        "oncokb_annotated": "yes" if any(_yes(record.get("OncoKB Annotated")) for record in matches) else "no",
        "cosmic_cgc_flag": "yes" if any(_yes(record.get("COSMIC CGC (v99)")) for record in matches) else "no",
        "cancer_gene_source_count": str(max(source_counts, default=0)),
        "cancer_gene_sources": ";".join(sources),
        "cancer_gene_match_basis": basis,
        "cancer_gene_context": (
            f"driver-context cancer gene: {';'.join(typed)}"
            if driver_records else
            f"listed without oncogene/TSG classification: {';'.join(typed)}"
        ),
    })
    return annotated


def annotate_events(rows: Iterable[Mapping[str, Any]], index: CancerGeneIndex | None) -> list[dict[str, Any]]:
    return [annotate_cancer_gene_context(row, index) for row in rows]


def propagate_to_peptide(peptide: Mapping[str, Any], event: Mapping[str, Any]) -> dict[str, Any]:
    annotated = dict(peptide)
    for field in CANCER_GENE_FIELDS:
        annotated[field] = event.get(field, "")
    return annotated
