#!/usr/bin/env python3
"""Build a full-gene GTEx normal tissue + HPA HSPC safety reference."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import zipfile
from collections import defaultdict
from pathlib import Path


HSPC_CELL_TYPES = {
    "hematopoietic stem cells",
    "erythrocyte progenitors",
    "megakaryocyte progenitors",
    "megakaryocyte-erythroid progenitors",
    "monocyte progenitors",
    "neutrophil progenitors",
}

CRITICAL_TISSUE_PREFIXES = (
    "Brain_", "Heart_", "Kidney_", "Liver", "Lung", "Muscle_Skeletal",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtex-median-gct", type=Path, required=True)
    parser.add_argument("--hpa-cell-type-zip", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--qc-out", type=Path, required=True)
    return parser.parse_args()


def update_max(store: dict[str, tuple[float, str]], gene: str, value: float, label: str) -> None:
    if gene not in store or value > store[gene][0]:
        store[gene] = (value, label)


def main() -> int:
    args = parse_args()
    tissue_max: dict[str, tuple[float, str]] = {}
    critical_max: dict[str, tuple[float, str]] = {}
    ensembl: dict[str, set[str]] = defaultdict(set)
    with gzip.open(args.gtex_median_gct, "rt", encoding="utf-8") as handle:
        handle.readline()
        handle.readline()
        reader = csv.DictReader(handle, delimiter="\t")
        tissue_columns = [name for name in reader.fieldnames or [] if name not in {"Name", "Description"}]
        critical_columns = [name for name in tissue_columns if name.startswith(CRITICAL_TISSUE_PREFIXES)]
        for row in reader:
            gene = (row.get("Description") or "").strip()
            if not gene:
                continue
            ensembl[gene].add((row.get("Name") or "").split(".")[0])
            values = [(float(row[name] or 0), name) for name in tissue_columns]
            critical_values = [(float(row[name] or 0), name) for name in critical_columns]
            if values:
                value, tissue = max(values)
                update_max(tissue_max, gene, value, tissue)
            if critical_values:
                value, tissue = max(critical_values)
                update_max(critical_max, gene, value, tissue)

    hspc_max: dict[str, tuple[float, str]] = {}
    with zipfile.ZipFile(args.hpa_cell_type_zip) as archive:
        member = archive.namelist()[0]
        with archive.open(member) as raw:
            reader = csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8"), delimiter="\t")
            for row in reader:
                cell_type = (row.get("Cell type") or "").strip().lower()
                if cell_type not in HSPC_CELL_TYPES:
                    continue
                gene = (row.get("Gene name") or "").strip()
                if not gene:
                    continue
                ensembl[gene].add((row.get("Gene") or "").split(".")[0])
                update_max(hspc_max, gene, float(row.get("nCPM") or 0), cell_type)

    genes = sorted(set(tissue_max) | set(hspc_max))
    fields = [
        "gene", "ensembl_gene_id", "normal_tissue_max_tpm", "normal_tissue_max_tissue",
        "critical_tissue_max_tpm", "critical_tissue_name", "critical_tissue_hit",
        "normal_hspc_tpm", "normal_hspc_cell_type", "normal_hspc_unit",
        "normal_expression_status", "normal_hspc_status", "gtex_version", "hpa_version",
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gene in genes:
            normal_value, normal_tissue = tissue_max.get(gene, (0.0, ""))
            critical_value, critical_tissue = critical_max.get(gene, (0.0, ""))
            hspc_value, hspc_type = hspc_max.get(gene, (0.0, ""))
            writer.writerow({
                "gene": gene,
                "ensembl_gene_id": ";".join(sorted(x for x in ensembl.get(gene, set()) if x)),
                "normal_tissue_max_tpm": f"{normal_value:.6f}",
                "normal_tissue_max_tissue": normal_tissue,
                "critical_tissue_max_tpm": f"{critical_value:.6f}",
                "critical_tissue_name": critical_tissue,
                "critical_tissue_hit": "yes" if critical_value > 0 else "no",
                "normal_hspc_tpm": f"{hspc_value:.6f}",
                "normal_hspc_cell_type": hspc_type,
                "normal_hspc_unit": "HPA_nCPM",
                "normal_expression_status": "ASSESSED" if gene in tissue_max else "UNASSESSED",
                "normal_hspc_status": "ASSESSED" if gene in hspc_max else "UNASSESSED",
                "gtex_version": "GTEx_v11_2025-08-22_GENCODE47",
                "hpa_version": "HPA_25.1_2025-12-12_Ensembl109",
            })

    qc = {
        "genes_total": len(genes),
        "genes_with_gtex": len(tissue_max),
        "genes_with_hspc": len(hspc_max),
        "gtex_tissues": len(tissue_columns),
        "critical_tissues": len(critical_columns),
        "hspc_cell_types": len(HSPC_CELL_TYPES),
    }
    with args.qc_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(qc.items())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
