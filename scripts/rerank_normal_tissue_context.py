#!/usr/bin/env python3
"""Annotate and rerank peptide candidates with GTEx tissue context.

The original model score is preserved. This creates a review-oriented score
and a diversity-limited experimental shortlist; it does not claim that normal
gene expression proves presentation of the candidate peptide in normal tissue.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


CRITICAL_ORGANS = {"Brain", "Heart", "Kidney", "Liver", "Lung", "Muscle"}
ORGAN_PREFIXES = (
    ("Small_Intestine", "Small intestine"),
    ("Minor_Salivary_Gland", "Salivary gland"),
    ("Esophagus", "Esophagus"),
    ("Adipose", "Adipose"),
    ("Adrenal_Gland", "Adrenal gland"),
    ("Blood", "Blood"),
    ("Brain", "Brain"),
    ("Breast", "Breast"),
    ("Colon", "Colon"),
    ("Heart", "Heart"),
    ("Kidney", "Kidney"),
    ("Liver", "Liver"),
    ("Lung", "Lung"),
    ("Muscle", "Muscle"),
    ("Nerve", "Nerve"),
    ("Pancreas", "Pancreas"),
    ("Pituitary", "Pituitary"),
    ("Skin", "Skin"),
    ("Stomach", "Stomach"),
    ("Testis", "Testis"),
    ("Thyroid", "Thyroid"),
    ("Uterus", "Uterus"),
    ("Vagina", "Vagina"),
)


def organ_for_tissue(tissue: str) -> str:
    for prefix, organ in ORGAN_PREFIXES:
        if tissue == prefix or tissue.startswith(prefix + "_"):
            return organ
    return tissue.split("_", 1)[0].replace("-", " ")


def as_float(value: object, default: float = 0.0) -> float:
    try:
        return float(str(value or "").strip())
    except ValueError:
        return default


def read_rows(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def load_gtex(path: Path, genes: set[str]) -> dict[str, dict[str, object]]:
    result: dict[str, dict[str, object]] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        handle.readline()
        handle.readline()
        reader = csv.DictReader(handle, delimiter="\t")
        tissues = [x for x in reader.fieldnames or [] if x not in {"Name", "Description"}]
        for row in reader:
            gene = (row.get("Description") or "").strip()
            if gene not in genes:
                continue
            tissue_values = sorted(
                ((as_float(row.get(tissue)), tissue) for tissue in tissues), reverse=True
            )
            organ_values: dict[str, float] = defaultdict(float)
            for value, tissue in tissue_values:
                organ = organ_for_tissue(tissue)
                organ_values[organ] = max(organ_values[organ], value)
            organs = sorted(((value, organ) for organ, value in organ_values.items()), reverse=True)
            top_tpm, top_tissue = tissue_values[0]
            top_organ_tpm, top_organ = organs[0]
            second_tpm, second_organ = organs[1] if len(organs) > 1 else (0.0, "")
            result[gene] = {
                "max_tpm": top_tpm,
                "max_tissue": top_tissue,
                "max_organ_tpm": top_organ_tpm,
                "max_organ": top_organ,
                "second_organ_tpm": second_tpm,
                "second_organ": second_organ,
                "specificity_ratio": top_organ_tpm / max(second_tpm, 0.01),
            }
    return result


def classify(
    gtex: dict[str, object], biopsy_tissue: str, min_tpm: float, min_ratio: float
) -> tuple[str, float, str, str]:
    if not gtex:
        return "UNASSESSED", 0.75, "GTEx gene not found", "REVIEW_GTEX_MISSING"
    max_tpm = as_float(gtex.get("max_organ_tpm"))
    ratio = as_float(gtex.get("specificity_ratio"))
    organ = str(gtex.get("max_organ") or "")
    enriched = max_tpm >= min_tpm and ratio >= min_ratio
    biopsy_match = organ.casefold() == biopsy_tissue.casefold() if biopsy_tissue else False
    if enriched and biopsy_match:
        return (
            "HIGH_CONTAMINATION_OR_OFF_TUMOR_RISK",
            0.25,
            f"GTEx tissue-enriched in biopsy organ ({organ}); tumor-cell origin unresolved",
            "HOLD_NORMAL_TISSUE_CONTEXT",
        )
    if max_tpm >= min_tpm and organ in CRITICAL_ORGANS:
        return (
            "CRITICAL_TISSUE_REVIEW",
            0.50,
            f"high GTEx expression in critical organ ({organ})",
            "HOLD_NORMAL_TISSUE_CONTEXT",
        )
    if enriched:
        return (
            "TISSUE_ENRICHED_REVIEW",
            0.75,
            f"GTEx tissue-enriched in {organ}",
            "REVIEW_NORMAL_TISSUE_CONTEXT",
        )
    return "LOW_CONTEXT_RISK", 1.0, "no strong GTEx tissue-enrichment signal", "ELIGIBLE"


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranked-peptides", type=Path, required=True)
    parser.add_argument("--gtex-median-gct", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--shortlist-out", type=Path, required=True)
    parser.add_argument("--biopsy-tissue", default="")
    parser.add_argument("--min-tissue-tpm", type=float, default=10.0)
    parser.add_argument("--min-organ-specificity-ratio", type=float, default=4.0)
    parser.add_argument("--max-per-gene", type=int, default=2)
    parser.add_argument("--max-per-event", type=int, default=1)
    parser.add_argument("--shortlist-size", type=int, default=50)
    args = parser.parse_args()

    rows, fields = read_rows(args.ranked_peptides)
    genes = {(row.get("gene") or "").split("::", 1)[0].strip() for row in rows}
    gtex_by_gene = load_gtex(args.gtex_median_gct, genes)
    extra = [
        "gtex_max_tpm", "gtex_max_tissue", "gtex_max_organ", "gtex_second_organ_tpm",
        "gtex_second_organ", "gtex_organ_specificity_ratio", "gtex_tissue_enriched",
        "biopsy_tissue_match", "normal_tissue_context_risk", "normal_tissue_context_multiplier",
        "normal_tissue_context_reason", "context_adjusted_score", "context_rank",
        "experimental_shortlist_status",
    ]
    for row in rows:
        gene = (row.get("gene") or "").split("::", 1)[0].strip()
        gtex = gtex_by_gene.get(gene, {})
        risk, multiplier, reason, shortlist_status = classify(
            gtex, args.biopsy_tissue, args.min_tissue_tpm, args.min_organ_specificity_ratio
        )
        ratio = as_float(gtex.get("specificity_ratio"))
        enriched = as_float(gtex.get("max_organ_tpm")) >= args.min_tissue_tpm and ratio >= args.min_organ_specificity_ratio
        row.update({
            "gtex_max_tpm": f"{as_float(gtex.get('max_tpm')):.4f}" if gtex else "",
            "gtex_max_tissue": gtex.get("max_tissue", ""),
            "gtex_max_organ": gtex.get("max_organ", ""),
            "gtex_second_organ_tpm": f"{as_float(gtex.get('second_organ_tpm')):.4f}" if gtex else "",
            "gtex_second_organ": gtex.get("second_organ", ""),
            "gtex_organ_specificity_ratio": f"{ratio:.4f}" if gtex else "",
            "gtex_tissue_enriched": "yes" if enriched else "no" if gtex else "not_assessed",
            "biopsy_tissue_match": "yes" if gtex and str(gtex.get("max_organ", "")).casefold() == args.biopsy_tissue.casefold() else "no",
            "normal_tissue_context_risk": risk,
            "normal_tissue_context_multiplier": f"{multiplier:.4f}",
            "normal_tissue_context_reason": reason,
            "context_adjusted_score": f"{as_float(row.get('efficacy_score')) * multiplier:.6f}",
            "experimental_shortlist_status": shortlist_status,
        })
    rows.sort(key=lambda row: as_float(row.get("context_adjusted_score")), reverse=True)
    for rank, row in enumerate(rows, 1):
        row["context_rank"] = str(rank)

    gene_counts: dict[str, int] = defaultdict(int)
    event_counts: dict[str, int] = defaultdict(int)
    shortlist: list[dict[str, object]] = []
    for row in rows:
        if row.get("experimental_shortlist_status") != "ELIGIBLE":
            continue
        gene = str(row.get("gene") or "")
        event = str(row.get("event_id") or "")
        if gene_counts[gene] >= args.max_per_gene or event_counts[event] >= args.max_per_event:
            continue
        row["experimental_shortlist_status"] = "SELECTED"
        shortlist.append(row)
        gene_counts[gene] += 1
        event_counts[event] += 1
        if len(shortlist) >= args.shortlist_size:
            break

    output_fields = fields + [field for field in extra if field not in fields]
    write_tsv(args.out, rows, output_fields)
    write_tsv(args.shortlist_out, shortlist, output_fields)
    print(f"rows={len(rows)} shortlist={len(shortlist)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
