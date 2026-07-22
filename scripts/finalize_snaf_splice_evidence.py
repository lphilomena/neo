#!/usr/bin/env python3
"""Merge SNAF splice predictions with an existing NeoAg evidence branch."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def file_record(path: Path) -> dict[str, object]:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(path), "size": path.stat().st_size, "sha256": digest.hexdigest()}


def clone_by_pair(
    source_rows: list[dict[str, str]],
    target_rows: list[dict[str, str]],
    overrides: dict[str, str] | None = None,
) -> list[dict[str, str]]:
    pair_index: dict[tuple[str, str], dict[str, str]] = {}
    allele_index: dict[str, dict[str, str]] = {}
    for row in source_rows:
        pair_index[(row.get("peptide", ""), row.get("hla_allele", ""))] = row
        allele_index.setdefault(row.get("hla_allele", ""), row)

    output: list[dict[str, str]] = []
    for target in target_rows:
        source = pair_index.get((target.get("peptide", ""), target.get("hla_allele", "")))
        if source is None:
            source = allele_index.get(target.get("hla_allele", ""))
        if source is None:
            raise ValueError(
                f"No evidence template for {target.get('peptide')} / {target.get('hla_allele')}"
            )
        row = dict(source)
        for key in ("peptide_id", "event_id", "sample_id", "peptide", "hla_allele", "mhc_class"):
            if key in target:
                row[key] = target[key]
        if overrides:
            row.update(overrides)
        output.append(row)
    return output


def make_ccf_rows(
    fields: list[str], events: list[dict[str, str]]
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for event in events:
        row = {field: "" for field in fields}
        for key in ("event_id", "sample_id", "gene", "event_type", "mutation_source", "chrom", "pos"):
            row[key] = event.get(key, "")
        row.update(
            {
                "ccf_status": "RNA_ONLY_UNRESOLVED",
                "clonality_status": "unresolved",
                "clonality_confidence": "unresolved",
                "clonality_multiplier": "0.8500",
                "ccf_method": "RNA_ONLY_UNRESOLVED",
                "ccf_confidence": "unresolved",
                "ccf_warning": "RNA-only splice event; DNA clonality is not estimated",
                "ccf_resolution": "unresolved",
                "ccf_resolution_reason": "RNA_ONLY_UNRESOLVED",
            }
        )
        rows.append(row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-run", required=True, type=Path)
    parser.add_argument("--branch-run", required=True, type=Path)
    args = parser.parse_args()

    base = args.base_run.resolve()
    branch = args.branch_run.resolve()
    snaf_dir = branch / "splice_snaf"

    raw_peptide_path = snaf_dir / "snaf_raw_peptides.tsv"
    raw_event_path = snaf_dir / "snaf_raw_events.tsv"
    predicted_path = snaf_dir / "predictors" / "presentation" / "presentation_evidence.tsv"
    base_presentation_path = base / "presentation" / "presentation_evidence.tsv"
    base_appm_path = base / "appm" / "appm_peptide_modifiers.tsv"
    base_escape_path = base / "immune_escape" / "peptide_escape_flags.tsv"
    base_ccf_path = base / "clonality" / "ccf_2.tsv"

    raw_peptide_fields, raw_peptides = read_tsv(raw_peptide_path)
    _, raw_events = read_tsv(raw_event_path)
    prediction_fields, predictions = read_tsv(predicted_path)
    base_prediction_fields, base_predictions = read_tsv(base_presentation_path)
    appm_fields, base_appm = read_tsv(base_appm_path)
    escape_fields, base_escape = read_tsv(base_escape_path)
    ccf_fields, base_ccf = read_tsv(base_ccf_path)

    prediction_index: dict[tuple[str, str], dict[str, str]] = {}
    for row in predictions:
        prediction_index[(row.get("peptide", ""), row.get("hla_allele", ""))] = row

    mapped_predictions: list[dict[str, str]] = []
    missing_pairs: list[tuple[str, str]] = []
    for peptide in raw_peptides:
        key = (peptide.get("peptide", ""), peptide.get("hla_allele", ""))
        source = prediction_index.get(key)
        if source is None:
            missing_pairs.append(key)
            continue
        row = dict(source)
        for field in ("peptide_id", "event_id", "sample_id", "peptide", "hla_allele", "mhc_class"):
            row[field] = peptide.get(field, row.get(field, ""))
        mapped_predictions.append(row)

    if missing_pairs:
        examples = ", ".join(f"{p}/{h}" for p, h in missing_pairs[:5])
        raise ValueError(f"Missing predictions for {len(missing_pairs)} SNAF pairs: {examples}")

    if prediction_fields != base_prediction_fields:
        missing = sorted(set(base_prediction_fields) - set(prediction_fields))
        if missing:
            raise ValueError(f"Prediction schema is missing base fields: {missing}")

    appm_rows = clone_by_pair(
        base_appm,
        raw_peptides,
        overrides={
            "priority_cap": "C_CAUTION",
            "appm_review_required": "yes",
        },
    )
    escape_rows = clone_by_pair(base_escape, raw_peptides)
    ccf_rows = make_ccf_rows(ccf_fields, raw_events)

    outputs = {
        "presentation": branch / "presentation" / "presentation_evidence.tsv",
        "appm_modifiers": branch / "appm" / "appm_peptide_modifiers.tsv",
        "appm_summary": branch / "appm" / "appm_summary.tsv",
        "escape": branch / "immune_escape" / "peptide_escape_flags.tsv",
        "ccf": branch / "clonality" / "ccf_2.tsv",
    }
    write_tsv(outputs["presentation"], base_prediction_fields, base_predictions + mapped_predictions)
    write_tsv(outputs["appm_modifiers"], appm_fields, base_appm + appm_rows)
    write_tsv(outputs["escape"], escape_fields, base_escape + escape_rows)
    write_tsv(outputs["ccf"], ccf_fields, base_ccf + ccf_rows)

    summary_fields, summary_rows = read_tsv(base / "appm" / "appm_summary.tsv")
    write_tsv(outputs["appm_summary"], summary_fields, summary_rows)

    duplicate_counts: dict[tuple[str, str], int] = defaultdict(int)
    for row in raw_peptides:
        duplicate_counts[(row.get("peptide", ""), row.get("hla_allele", ""))] += 1

    manifest = {
        "status": "PASS",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "base_run": str(base),
        "branch_run": str(branch),
        "snaf_events": len(raw_events),
        "snaf_peptide_hla_rows": len(raw_peptides),
        "unique_prediction_pairs": len(prediction_index),
        "mapped_prediction_rows": len(mapped_predictions),
        "pairs_reused_across_events": sum(1 for value in duplicate_counts.values() if value > 1),
        "ccf_policy": "RNA_ONLY_UNRESOLVED; ccf_estimate remains blank; multiplier=0.85",
        "priority_policy": "SNAF splice candidates remain capped at C_CAUTION pending orthogonal validation",
        "outputs": {name: file_record(path) for name, path in outputs.items()},
    }
    manifest_path = branch / "qc" / "snaf_evidence_merge_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
