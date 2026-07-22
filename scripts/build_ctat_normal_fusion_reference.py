#!/usr/bin/env python3
"""Convert the CTAT GTEx recurrent fusion matrix into NeoAg junction evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path


VALUE_RE = re.compile(r"^(?P<pct>[0-9.]+)(?:;n=(?P<count>[0-9]+))?$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--qc-out", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    args = parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    rows_written = 0
    tissues_seen: set[str] = set()
    fields = [
        "gene_pair", "gene1", "gene2", "normal_sample_count",
        "normal_reads", "max_tissue_prevalence_pct", "tissue",
        "tissue_count", "source", "dataset", "junction_class",
    ]
    with args.input.open(encoding="utf-8") as source, args.out.open(
        "w", encoding="utf-8", newline=""
    ) as target:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        tissue_columns = [column for column in reader.fieldnames or [] if column != "Fusion"]
        for row in reader:
            fusion = (row.get("Fusion") or "").strip()
            if "--" not in fusion:
                continue
            gene1, gene2 = fusion.split("--", 1)
            total_samples = 0
            max_pct = 0.0
            max_tissue = ""
            max_tissue_count = 0
            for tissue in tissue_columns:
                raw = (row.get(tissue) or "0").strip()
                match = VALUE_RE.match(raw)
                if not match:
                    continue
                pct = float(match.group("pct") or 0)
                count = int(match.group("count") or 0)
                total_samples += count
                if pct > max_pct or (pct == max_pct and count > max_tissue_count):
                    max_pct = pct
                    max_tissue = tissue.removeprefix("GTEx-")
                    max_tissue_count = count
                if count:
                    tissues_seen.add(tissue)
            writer.writerow({
                "gene_pair": f"{gene1}::{gene2}",
                "gene1": gene1,
                "gene2": gene2,
                "normal_sample_count": total_samples,
                # Compatibility with the existing safety output field. The unit
                # is normal samples, not junction-supporting reads.
                "normal_reads": total_samples,
                "max_tissue_prevalence_pct": f"{max_pct:.3f}",
                "tissue": max_tissue,
                "tissue_count": max_tissue_count,
                "source": "CTAT_HumanFusionLib",
                "dataset": "GTEx_recurrent_StarF2019",
                "junction_class": "normal_recurrent_fusion",
            })
            rows_written += 1

    qc = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "input": str(args.input),
        "input_sha256": sha256(args.input),
        "output": str(args.out),
        "output_sha256": sha256(args.out),
        "fusion_pairs": rows_written,
        "tissues_with_observations": len(tissues_seen),
        "source": "FusionAnnotator/CTAT_HumanFusionLib",
        "dataset": "GTEx_recurrent_StarF2019",
        "scope": "gene-pair recurrent normal fusion background; not exact breakpoint evidence",
    }
    args.qc_out.parent.mkdir(parents=True, exist_ok=True)
    args.qc_out.write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
