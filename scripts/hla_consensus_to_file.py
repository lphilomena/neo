#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert HLA typing consensus TSV to one allele per line")
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    rows = list(csv.DictReader(Path(args.consensus).open(encoding="utf-8"), delimiter="\t"))
    alleles: list[str] = []
    for row in rows:
        for value in str(row.get("consensus_lowres") or "").split(" / "):
            value = value.strip()
            if value and "*" in value:
                allele = value if value.upper().startswith("HLA-") else "HLA-" + value
                if allele not in alleles:
                    alleles.append(allele)
    if not any(value.startswith("HLA-A*") for value in alleles) or not any(value.startswith("HLA-B*") for value in alleles):
        raise SystemExit("HLA consensus lacks sufficient class-I calls")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(alleles) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
