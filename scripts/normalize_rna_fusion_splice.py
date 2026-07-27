#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


EVENT_FIELDS = ["sample_id", "event_id", "event_type", "gene", "chrom", "start", "end", "rna_junction_reads", "source_tool", "normal_junction_status", "confidence"]
PEPTIDE_FIELDS = ["sample_id", "event_id", "peptide_id", "peptide", "hla_allele", "source_type", "crosses_junction", "generation_status"]


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        first = handle.readline()
        handle.seek(0)
        delimiter = "\t" if "\t" in first else ","
        return [{str(k): str(v or "") for k, v in row.items()} for row in csv.DictReader(handle, delimiter=delimiter)]


def first(row: dict[str, str], names: list[str]) -> str:
    lower = {str(k).lower(): str(v or "") for k, v in row.items()}
    for name in names:
        value = lower.get(name.lower(), "").strip()
        if value:
            return value
    return ""


def write_rows(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def regtools_rows(path: Path) -> list[dict[str, str]]:
    rows = read_rows(path)
    if rows and any(key.lower() in {"chrom", "chromosome", "junction"} for key in rows[0]):
        return rows
    parsed: list[dict[str, str]] = []
    if not path.is_file():
        return parsed
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) >= 5:
            parsed.append({"chrom": parts[0], "start": parts[1], "end": parts[2], "junction_id": parts[3], "read_count": parts[4]})
    return parsed


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--junctions", required=True, type=Path)
    ap.add_argument("--snaf", type=Path)
    ap.add_argument("--splicemutr", type=Path)
    ap.add_argument("--normal-junctions", type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()

    normal_keys: set[str] = set()
    if args.normal_junctions and args.normal_junctions.is_file():
        for row in read_rows(args.normal_junctions):
            chrom = first(row, ["chrom", "chromosome", "chr"])
            start = first(row, ["start", "junction_start", "donor"])
            end = first(row, ["end", "junction_end", "acceptor"])
            if chrom and start and end:
                normal_keys.add(f"{chrom}:{start}-{end}")

    sources = [("RegTools", regtools_rows(args.junctions))]
    if args.snaf:
        sources.append(("SNAF", read_rows(args.snaf)))
    if args.splicemutr:
        sources.append(("SpliceMutr", read_rows(args.splicemutr)))
    events: list[dict[str, str]] = []
    peptides: list[dict[str, str]] = []
    by_key: dict[str, set[str]] = defaultdict(set)
    event_seen: set[tuple[str, str]] = set()
    peptide_seen: set[tuple[str, str, str]] = set()
    for tool, rows in sources:
        for index, row in enumerate(rows, 1):
            chrom = first(row, ["chrom", "chromosome", "chr", "seqnames"])
            start = first(row, ["start", "junction_start", "intron_start", "donor"])
            end = first(row, ["end", "junction_end", "intron_end", "acceptor"])
            gene = first(row, ["gene", "gene_name", "symbol"])
            key = first(row, ["event_id", "junction_id", "uid"]) or (f"{chrom}:{start}-{end}" if chrom and start and end else f"{tool}:{index}")
            event_id = f"SPLICE|{key}"
            by_key[key].add(tool)
            if (event_id, tool) not in event_seen:
                event_seen.add((event_id, tool))
                normal_status = "UNASSESSED" if not normal_keys else ("DETECTED" if f"{chrom}:{start}-{end}" in normal_keys else "NOT_DETECTED")
                events.append({
                    "sample_id": args.sample_id, "event_id": event_id, "event_type": "Splice", "gene": gene,
                    "chrom": chrom, "start": start, "end": end,
                    "rna_junction_reads": first(row, ["junction_reads", "read_count", "reads", "score"]),
                    "source_tool": tool, "normal_junction_status": normal_status,
                    "confidence": "normal_junction_review" if normal_status != "NOT_DETECTED" else "rna_junction_candidate",
                })
            peptide = first(row, ["peptide", "junction_peptide", "mutant_peptide", "neoepitope"])
            hla = first(row, ["hla_allele", "hla", "allele"])
            if peptide and (event_id, peptide, hla) not in peptide_seen:
                peptide_seen.add((event_id, peptide, hla))
                peptides.append({
                    "sample_id": args.sample_id, "event_id": event_id,
                    "peptide_id": f"{event_id}|{peptide}|{hla or 'HLA_UNASSESSED'}",
                    "peptide": peptide, "hla_allele": hla, "source_type": tool.lower(),
                    "crosses_junction": "true", "generation_status": "provided_by_splice_caller",
                })

    consensus = [{"junction_key": key, "support_tools": ",".join(sorted(tools)), "n_tools": str(len(tools)), "status": "CROSS_VALIDATED" if len(tools) >= 2 else "SINGLE_TOOL"} for key, tools in sorted(by_key.items())]
    evidence = [{"event_id": row["event_id"], "gene": row["gene"], "rna_junction_reads": row["rna_junction_reads"], "source_tool": row["source_tool"], "normal_junction_status": row["normal_junction_status"]} for row in events]
    write_rows(args.outdir / "raw_events.tsv", EVENT_FIELDS, events)
    write_rows(args.outdir / "raw_peptides.tsv", PEPTIDE_FIELDS, peptides)
    write_rows(args.outdir / "rna_junction_evidence.tsv", ["event_id", "gene", "rna_junction_reads", "source_tool", "normal_junction_status"], evidence)
    write_rows(args.outdir / "splice_consensus.tsv", ["junction_key", "support_tools", "n_tools", "status"], consensus)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
