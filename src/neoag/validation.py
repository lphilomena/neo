from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .validation_design import (
    VALIDATION_PLAN_FIELDS,
    design_validation_row,
    load_peptide_catalog_index,
    lookup_peptide_catalog,
    resolve_peptide_catalog,
)


def make_validation_plan(
    peptides: list[Mapping[str, Any]],
    *,
    peptide_catalog_tsv: str | Path | None = None,
    outdir: str | Path | None = None,
) -> list[dict[str, str]]:
    catalog_path = resolve_peptide_catalog(peptide_catalog_tsv, outdir=outdir)
    catalog_index = load_peptide_catalog_index(catalog_path)
    rows = []
    for p in peptides:
        catalog_row = lookup_peptide_catalog(p, catalog_index)
        rows.append(design_validation_row(p, catalog_row=catalog_row))
    return rows


def validation_plan_fieldnames() -> list[str]:
    return list(VALIDATION_PLAN_FIELDS)
