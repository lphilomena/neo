#!/usr/bin/env python3
"""Attach an independently predicted WT control to one splice peptide-HLA row."""

from __future__ import annotations

import argparse
import csv
import json
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


def update_rows(
    path: Path,
    event_id: str,
    mutant_peptide: str,
    hla: str,
    values: dict[str, str],
) -> int:
    fields, rows = read_tsv(path)
    missing_fields = [field for field in values if field not in fields]
    if missing_fields:
        fields.extend(missing_fields)
    count = 0
    for row in rows:
        if (
            row.get("event_id") == event_id
            and row.get("peptide") == mutant_peptide
            and row.get("hla_allele") == hla
        ):
            row.update(values)
            count += 1
    if count != 1:
        raise ValueError(f"Expected one row in {path}, found {count}")
    write_tsv(path, fields, rows)
    return count


def sync_safety(path: Path, safety_path: Path) -> int:
    fields, rows = read_tsv(path)
    _, safety_rows = read_tsv(safety_path)
    safety_index = {row.get("peptide_id", ""): row for row in safety_rows}
    safety_fields = [
        "safety_tier",
        "safety_status",
        "safety_reason",
        "safety_multiplier",
        "review_required",
        "reference_proteome_exact_match",
        "normal_ligand_tissue",
        "mutation_anchor_only",
    ]
    for field in safety_fields:
        if field not in fields:
            fields.append(field)
    count = 0
    for row in rows:
        safety = safety_index.get(row.get("peptide_id", ""))
        if not safety:
            continue
        for field in safety_fields:
            value = str(safety.get(field, "")).strip()
            if value:
                row[field] = value
        count += 1
    write_tsv(path, fields, rows)
    return count


def sync_event_safety(path: Path, safety_path: Path) -> int:
    fields, rows = read_tsv(path)
    _, safety_rows = read_tsv(safety_path)
    safety_index = {row.get("event_id", ""): row for row in safety_rows}
    mappings = {
        "event_safety_status": "safety_status",
        "event_safety_reason": "safety_reason",
        "safety_evidence_completeness": "safety_evidence_completeness",
        "safety_missing_layers": "safety_missing_layers",
    }
    for destination in mappings.values():
        if destination not in fields:
            fields.append(destination)
    count = 0
    for row in rows:
        safety = safety_index.get(row.get("event_id", ""))
        if not safety:
            continue
        for source, destination in mappings.items():
            value = str(safety.get(source, "")).strip()
            if value:
                row[destination] = value
        count += 1
    write_tsv(path, fields, rows)
    return count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--branch-run", required=True, type=Path)
    parser.add_argument("--prediction", required=True, type=Path)
    parser.add_argument("--event-id", required=True)
    parser.add_argument("--mutant-peptide", required=True)
    parser.add_argument("--wildtype-peptide", required=True)
    parser.add_argument("--hla", required=True)
    args = parser.parse_args()

    _, predictions = read_tsv(args.prediction)
    index = {(row.get("peptide"), row.get("hla_allele")): row for row in predictions}
    mt = index.get((args.mutant_peptide, args.hla))
    wt = index.get((args.wildtype_peptide, args.hla))
    if mt is None or wt is None:
        raise ValueError("The prediction table must contain both MT and WT peptide-HLA rows")

    raw_values = {"wildtype_peptide": args.wildtype_peptide}
    presentation_values = {
        "netmhcpan_mt_rank_ba": mt.get("netmhcpan_ba_rank", ""),
        "netmhcpan_mt_rank_el": mt.get("netmhcpan_el_rank", ""),
        "netmhcpan_wt_rank_ba": wt.get("netmhcpan_ba_rank", ""),
        "netmhcpan_wt_rank_el": wt.get("netmhcpan_el_rank", ""),
        "mhcflurry_wt_affinity_percentile": wt.get("mhcflurry_affinity_percentile", ""),
        "mhcflurry_wt_processing_score": wt.get("mhcflurry_processing_score", ""),
        "mhcflurry_wt_presentation_score": wt.get("mhcflurry_presentation_score", ""),
        "prime_wt_score": wt.get("prime_score", ""),
        "prime_wt_rank": wt.get("prime_rank", ""),
        "bigmhc_im_wt_score": wt.get("bigmhc_im_score", ""),
    }

    raw_paths = [
        args.branch_run / "inputs" / "combined_raw_peptides.tsv",
        args.branch_run / "splice_snaf" / "snaf_raw_peptides.tsv",
    ]
    updated = {
        str(path): update_rows(
            path,
            args.event_id,
            args.mutant_peptide,
            args.hla,
            raw_values,
        )
        for path in raw_paths
    }
    presentation_path = args.branch_run / "presentation" / "presentation_evidence.tsv"
    updated[str(presentation_path)] = update_rows(
        presentation_path,
        args.event_id,
        args.mutant_peptide,
        args.hla,
        presentation_values,
    )
    safety_path = args.branch_run / "safety" / "peptide_safety.tsv"
    updated[f"{raw_paths[0]}#safety_sync"] = sync_safety(raw_paths[0], safety_path)
    updated[f"{raw_paths[1]}#safety_sync"] = sync_safety(raw_paths[1], safety_path)
    event_safety_path = args.branch_run / "safety" / "event_safety.tsv"
    event_paths = [
        args.branch_run / "inputs" / "combined_raw_events.tsv",
        args.branch_run / "splice_snaf" / "snaf_raw_events.tsv",
    ]
    updated[f"{event_paths[0]}#safety_sync"] = sync_event_safety(
        event_paths[0], event_safety_path
    )
    updated[f"{event_paths[1]}#safety_sync"] = sync_event_safety(
        event_paths[1], event_safety_path
    )

    mt_el = float(mt["netmhcpan_el_rank"])
    wt_el = float(wt["netmhcpan_el_rank"])
    summary = {
        "status": "PASS",
        "event_id": args.event_id,
        "hla_allele": args.hla,
        "mutant_peptide": args.mutant_peptide,
        "wildtype_peptide": args.wildtype_peptide,
        "netmhcpan_mt_el_rank": mt_el,
        "netmhcpan_wt_el_rank": wt_el,
        "agretopicity_wt_over_mt": wt_el / mt_el if mt_el else None,
        "mhcflurry_mt_presentation_score": mt.get("mhcflurry_presentation_score", ""),
        "mhcflurry_wt_presentation_score": wt.get("mhcflurry_presentation_score", ""),
        "prime_mt_score": mt.get("prime_score", ""),
        "prime_wt_score": wt.get("prime_score", ""),
        "bigmhc_mt_score": mt.get("bigmhc_im_score", ""),
        "bigmhc_wt_score": wt.get("bigmhc_im_score", ""),
        "updated_files": updated,
    }
    out = args.prediction.parent.parent / "mt_wt_comparison.json"
    out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
