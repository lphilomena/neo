#!/usr/bin/env python3
"""Build a compact normal-junction catalog from the GTEx v11 junction GCT.

The source GCT is sample-level and very wide. This builder streams it and
retains the complete junction coordinate set without loading the matrix into
memory. A row in the GCT establishes presence in at least one GTEx sample; the
catalog deliberately records that as a lower bound instead of inventing a
sample frequency.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path


FIELDS = [
    "junction_id",
    "gene_id",
    "gene_pair",
    "normal_sample_count",
    "count_metric",
    "source",
    "tissue",
    "junction_class",
]


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def build_catalog(gct: Path, output: Path) -> int:
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    with _open_text(gct) as source, output.open("w", encoding="utf-8", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for line_number, line in enumerate(source, start=1):
            if line_number <= 2:
                continue
            if line_number == 3:
                if not line.startswith("Name\tDescription\t"):
                    raise ValueError("Unexpected GTEx GCT header")
                continue
            first_tab = line.find("\t")
            second_tab = line.find("\t", first_tab + 1)
            if first_tab < 1 or second_tab < 0:
                continue
            junction_id = line[:first_tab]
            gene_id = line[first_tab + 1:second_tab]
            writer.writerow({
                "junction_id": junction_id,
                "gene_id": gene_id,
                "gene_pair": "",
                "normal_sample_count": "1",
                "count_metric": "presence_lower_bound",
                "source": "GTEx_v11_exon_junctions",
                "tissue": "GTEx_all_tissues",
                "junction_class": "normal_splice_junction",
            })
            rows += 1
    return rows


def combine_catalogs(splice_catalog: Path, fusion_catalog: Path | None, output: Path) -> int:
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    with output.open("w", encoding="utf-8", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for path in [splice_catalog, fusion_catalog]:
            if not path:
                continue
            with path.open("r", encoding="utf-8", newline="") as source:
                for row in csv.DictReader(source, delimiter="\t"):
                    writer.writerow({field: row.get(field, "") for field in FIELDS})
                    rows += 1
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gct", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--fusion-catalog", type=Path)
    parser.add_argument("--combined-out", type=Path)
    parser.add_argument("--qc-json", type=Path)
    args = parser.parse_args()

    splice_rows = build_catalog(args.gct, args.out)
    combined_rows = None
    if args.combined_out:
        combined_rows = combine_catalogs(args.out, args.fusion_catalog, args.combined_out)
    qc = {
        "source": str(args.gct),
        "junction_class": "normal_splice_junction",
        "splice_junction_rows": splice_rows,
        "count_metric": "presence_lower_bound",
        "fusion_catalog": str(args.fusion_catalog) if args.fusion_catalog else None,
        "combined_rows": combined_rows,
    }
    if args.qc_json:
        args.qc_json.parent.mkdir(parents=True, exist_ok=True)
        args.qc_json.write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
