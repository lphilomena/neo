#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_gtf_attributes(attr: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for raw in attr.strip().rstrip(";").split(";"):
        raw = raw.strip()
        if not raw:
            continue
        if " " in raw:
            k, v = raw.split(" ", 1)
            out[k] = v.strip().strip('"')
        elif "=" in raw:
            k, v = raw.split("=", 1)
            out[k] = v.strip().strip('"')
    return out


def merge_len(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    intervals.sort()
    total = 0
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
    total += cur_e - cur_s + 1
    return total


def load_gene_lengths(gtf_path: Path) -> dict[str, int]:
    by_gene: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with gtf_path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            start, end, attrs = int(parts[3]), int(parts[4]), parts[8]
            ad = parse_gtf_attributes(attrs)
            gene = (ad.get("gene_name") or ad.get("gene_id") or "").upper()
            if gene:
                by_gene[gene].append((start, end))
    return {gene: merge_len(iv) for gene, iv in by_gene.items()}


def sniff_delimiter(path: Path) -> str:
    if path.suffix.lower() in {".tsv", ".txt"}:
        sample = path.read_text(errors="ignore")[:4096]
        return "\t" if "\t" in sample else ","
    return ","


def load_counts(path: Path) -> list[tuple[str, float]]:
    delim = sniff_delimiter(path)
    rows: list[tuple[str, float]] = []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        for row in reader:
            gene = (
                row.get("gene")
                or row.get("Gene")
                or row.get("gene_name")
                or row.get("symbol")
                or row.get("SYMBOL")
                or ""
            ).strip().upper()
            if not gene:
                continue
            raw = (
                row.get("count")
                or row.get("counts")
                or row.get("raw_count")
                or row.get("raw_counts")
                or row.get("expression")
                or row.get("expr")
                or row.get("TPM")
                or "0"
            )
            try:
                count = float(raw or 0)
            except Exception:
                count = 0.0
            rows.append((gene, max(count, 0.0)))
    return rows


def counts_to_tpm(rows: list[tuple[str, float]], gene_lengths: dict[str, int]) -> tuple[list[tuple[str, float]], int]:
    rpk: list[tuple[str, float]] = []
    matched = 0
    for gene, count in rows:
        glen = gene_lengths.get(gene, 0)
        if glen > 0:
            matched += 1
            rpk.append((gene, count / (glen / 1000.0)))
        else:
            rpk.append((gene, 0.0))
    scale = sum(v for _, v in rpk) / 1_000_000.0
    if scale <= 0:
        return [(gene, 0.0) for gene, _ in rows], matched
    return [(gene, v / scale) for gene, v in rpk], matched


def main() -> None:
    ap = argparse.ArgumentParser(description="Convert gene raw counts to TPM using exon-union gene lengths from GTF.")
    ap.add_argument("--counts", required=True, type=Path)
    ap.add_argument("--gtf", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    gene_lengths = load_gene_lengths(args.gtf)
    rows = load_counts(args.counts)
    tpm_rows, matched = counts_to_tpm(rows, gene_lengths)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["gene", "TPM"])
        for gene, tpm in tpm_rows:
            writer.writerow([gene, f"{tpm:.6f}"])

    print(f"genes_in_counts={len(rows)}")
    print(f"genes_with_gtf_length={matched}")
    print(f"out={args.out}")


if __name__ == "__main__":
    main()
