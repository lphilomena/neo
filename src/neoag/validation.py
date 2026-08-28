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
    group_counts: dict[str, int] = {}
    for p in peptides:
        catalog_row = lookup_peptide_catalog(p, catalog_index)
        row = design_validation_row(p, catalog_row=catalog_row)
        group = str(p.get("redundancy_group") or p.get("phase_group_id") or "").strip()
        if row["validation_mode"] in {"do_not_advance", "phasing_required"}:
            row["shortlist_status"] = "NOT_ELIGIBLE"
        elif group:
            rank = group_counts.get(group, 0) + 1
            group_counts[group] = rank
            row["shortlist_rank"] = str(rank)
            if rank <= 2:
                row["shortlist_status"] = "SHORTLISTED"
            else:
                row["shortlist_status"] = "REDUNDANT_NOT_SHORTLISTED"
                row["validation_notes"] = "; ".join(x for x in [row["validation_notes"], "overlapping phased-haplotype candidate; retain only top 2 peptide-HLA combinations"] if x)
        else:
            row["shortlist_status"] = "SHORTLISTED"
        rows.append(row)
    return rows


def validation_plan_fieldnames() -> list[str]:
    return list(VALIDATION_PLAN_FIELDS)
