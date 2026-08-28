#!/usr/bin/env python3
"""Integrate batch SNAF splice WT reconstruction and predictions into a branch."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--branch-run", required=True, type=Path)
    parser.add_argument("--pairs", required=True, type=Path)
    parser.add_argument("--prediction", required=True, type=Path)
    args = parser.parse_args()

    _, pairs = read_tsv(args.pairs)
    _, predictions = read_tsv(args.prediction)
    prediction_index = {
        (row.get("peptide", ""), row.get("hla_allele", "")): row for row in predictions
    }
    pair_index = {
        (row["event_id"], row["mutant_peptide"], row["hla_allele"]): row for row in pairs
    }

    raw_paths = [
        args.branch_run / "inputs" / "combined_raw_peptides.tsv",
        args.branch_run / "splice_snaf" / "snaf_raw_peptides.tsv",
    ]
    raw_updates: dict[str, int] = {}
    peptide_id_to_pair: dict[str, dict[str, str]] = {}
    for path in raw_paths:
        fields, rows = read_tsv(path)
        count = 0
        for row in rows:
            key = (row.get("event_id", ""), row.get("peptide", ""), row.get("hla_allele", ""))
            pair = pair_index.get(key)
            if not pair or pair.get("status") not in {"PASS", "NON_MUTANT_SEQUENCE"}:
                continue
            row["wildtype_peptide"] = pair["wildtype_peptide"]
            peptide_id_to_pair[row["peptide_id"]] = pair
            count += 1
        write_tsv(path, fields, rows)
        raw_updates[str(path)] = count

    presentation_path = args.branch_run / "presentation" / "presentation_evidence.tsv"
    presentation_fields, presentation_rows = read_tsv(presentation_path)
    wt_fields = [
        "netmhcpan_wt_rank_ba",
        "netmhcpan_wt_rank_el",
        "mhcflurry_wt_affinity_percentile",
        "mhcflurry_wt_processing_score",
        "mhcflurry_wt_presentation_score",
        "prime_wt_score",
        "prime_wt_rank",
        "bigmhc_im_wt_score",
    ]
    for field in wt_fields:
        if field not in presentation_fields:
            presentation_fields.append(field)
    prediction_updates = identical_updates = missing_predictions = 0
    missing_examples: list[str] = []
    for row in presentation_rows:
        pair = peptide_id_to_pair.get(row.get("peptide_id", ""))
        if not pair:
            continue
        if pair["status"] == "NON_MUTANT_SEQUENCE":
            identical_updates += 1
            continue
        wt = prediction_index.get((pair["wildtype_peptide"], pair["hla_allele"]))
        if not wt:
            missing_predictions += 1
            if len(missing_examples) < 10:
                missing_examples.append(
                    f"{pair['event_id']}|{pair['wildtype_peptide']}|{pair['hla_allele']}"
                )
            continue
        row.update(
            {
                "netmhcpan_wt_rank_ba": wt.get("netmhcpan_ba_rank", ""),
                "netmhcpan_wt_rank_el": wt.get("netmhcpan_el_rank", ""),
                "mhcflurry_wt_affinity_percentile": wt.get("mhcflurry_affinity_percentile", ""),
                "mhcflurry_wt_processing_score": wt.get("mhcflurry_processing_score", ""),
                "mhcflurry_wt_presentation_score": wt.get("mhcflurry_presentation_score", ""),
                "prime_wt_score": wt.get("prime_score", ""),
                "prime_wt_rank": wt.get("prime_rank", ""),
                "bigmhc_im_wt_score": wt.get("bigmhc_im_score", ""),
            }
        )
        prediction_updates += 1
    if missing_predictions:
        raise ValueError(
            f"Missing WT predictions for {missing_predictions} rows: {missing_examples}"
        )
    write_tsv(presentation_path, presentation_fields, presentation_rows)

    status_counts = Counter(pair["status"] for pair in pairs)
    summary = {
        "status": "PASS",
        "pair_status_counts": dict(sorted(status_counts.items())),
        "raw_updates": raw_updates,
        "presentation_prediction_updates": prediction_updates,
        "identical_sequence_updates": identical_updates,
        "unresolved_rows": status_counts.get("WT_NOT_RECONSTRUCTABLE", 0),
        "prediction_unique_pairs": len(prediction_index),
        "presentation_path": str(presentation_path),
    }
    out = args.pairs.parent / "splice_mt_wt_integration_summary.json"
    out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
