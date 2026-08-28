#!/usr/bin/env python3
"""Join long-read splice provenance to the project's presentation evidence."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--provenance", required=True, type=Path)
    parser.add_argument("--presentation", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    provenance_fields, provenance = rows(args.provenance)
    prediction_fields, predictions = rows(args.presentation)
    prediction_fields = [f for f in prediction_fields if f not in {"peptide_id", "event_id", "sample_id", "peptide", "hla_allele"}]
    index = {(row.get("peptide", ""), row.get("hla_allele", "")): row for row in predictions}
    output_fields = provenance_fields + [f for f in prediction_fields if f not in provenance_fields]
    output = []
    missing = []
    for row in provenance:
        key = (row.get("peptide", ""), row.get("hla", ""))
        pred = index.get(key)
        merged = dict(row)
        if pred is None:
            merged["prediction_status"] = "missing"
            missing.append(key)
        else:
            merged.update({field: pred.get(field, "") for field in prediction_fields})
            merged["prediction_status"] = "complete"
        output.append(merged)
    if "prediction_status" not in output_fields:
        output_fields.append("prediction_status")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output)
    summary = {
        "status": "PASS" if not missing else "WARN",
        "provenance_rows": len(provenance),
        "presentation_rows": len(predictions),
        "merged_rows": len(output),
        "complete_rows": sum(row.get("prediction_status") == "complete" for row in output),
        "missing_prediction_rows": len(missing),
        "longrna_supported_rows": sum(row.get("longrna_support") == "true" for row in output),
        "presentation_grade_counts": dict(Counter(row.get("presentation_evidence_grade", "") for row in output)),
        "netmhcstabpan_policy": "required_for_production",
        "output": str(args.output),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())
