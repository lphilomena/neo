from __future__ import annotations

from pathlib import Path
from typing import Mapping, Any

from ..utils import read_tsv, first, to_float


def load_expression(path: str | Path | None) -> dict[str, float]:
    if not path or not Path(path).exists():
        return {}
    out: dict[str, float] = {}
    for r in read_tsv(path):
        gene = first(r, ["gene", "Gene", "gene_name", "symbol", "SYMBOL"], "")
        if gene:
            out[gene] = to_float(first(r, ["TPM", "tpm", "expression", "expr", "gene_tpm"], "0"), 0.0)
    return out


def load_normal_expression(path: str | Path | None) -> dict[str, dict[str, Any]]:
    if not path or not Path(path).exists():
        return {}
    out: dict[str, dict[str, Any]] = {}
    for r in read_tsv(path):
        gene = first(r, ["gene", "Gene", "gene_name", "symbol", "SYMBOL"], "")
        if not gene:
            continue
        out[gene] = {
            "normal_tissue_max_tpm": to_float(first(r, ["normal_tissue_max_tpm", "max_tpm", "TPM"], "0"), 0.0),
            "normal_hspc_tpm": to_float(first(r, ["normal_hspc_tpm", "hspc_tpm"], "0"), 0.0),
            "critical_tissue_hit": first(r, ["critical_tissue_hit", "critical"], "no"),
        }
    return out


def normal_expr_for_genes(genes: list[str], normal_expr: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    nt = 0.0
    nh = 0.0
    crit = False
    for g in genes:
        e = normal_expr.get(g, {})
        nt = max(nt, to_float(e.get("normal_tissue_max_tpm"), 0.0))
        nh = max(nh, to_float(e.get("normal_hspc_tpm"), 0.0))
        crit = crit or str(e.get("critical_tissue_hit", "")).strip().lower() in {"1", "true", "yes", "y"}
    return {"normal_tissue_max_tpm": nt, "normal_hspc_tpm": nh, "critical_tissue_hit": "yes" if crit else "no"}


def load_normal_ligands(path: str | Path | None) -> set[str]:
    if not path or not Path(path).exists():
        return set()
    return {
        first(r, ["peptide", "Peptide", "sequence"], "").strip().upper()
        for r in read_tsv(path)
        if first(r, ["peptide", "Peptide", "sequence"], "").strip()
    }


def load_junction_reads(path: str | Path | None) -> dict[str, int]:
    """Load optional RNA junction support.

    Accepted keys: event_id, sv_event_id, gene_pair (GENE1::GENE2), or gene1/gene2.
    Value columns: junction_reads, split_reads, spanning_reads, reads.
    """
    if not path or not Path(path).exists():
        return {}
    out: dict[str, int] = {}
    for r in read_tsv(path):
        val = int(to_float(first(r, ["junction_reads", "split_reads", "spanning_reads", "reads"], "0"), 0.0))
        keys = [first(r, ["event_id"], ""), first(r, ["sv_event_id"], ""), first(r, ["gene_pair", "fusion"], "")]
        g1 = first(r, ["gene1", "left_gene"], "")
        g2 = first(r, ["gene2", "right_gene"], "")
        if g1 and g2:
            keys.extend([f"{g1}::{g2}", f"{g2}::{g1}"])
        for k in keys:
            if k:
                out[k] = max(out.get(k, 0), val)
    return out
