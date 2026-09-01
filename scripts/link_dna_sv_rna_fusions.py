#!/usr/bin/env python3
"""Attach DNA-SV evidence to RNA fusion entities without gene-pair evidence leakage."""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path


def read_table(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        header = handle.readline()
        handle.seek(0)
        delimiter = "\t" if "\t" in header else ","
        return [{str(k): str(v or "") for k, v in row.items()} for row in csv.DictReader(handle, delimiter=delimiter)]


def write_table(path: Path, rows: list[dict[str, str]], preferred: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(preferred or [])
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean_gene(value: str) -> str:
    return re.split(r"[\^|]", str(value or ""), maxsplit=1)[0].strip().upper()


def gene_pair(value: str) -> tuple[str, str] | None:
    parts = re.split(r"::|--|/", str(value or ""))
    genes = tuple(clean_gene(part) for part in parts if clean_gene(part))
    return genes[:2] if len(genes) >= 2 else None


def breakpoint(value: str) -> tuple[str, int | None, str]:
    parts = str(value or "").replace(";", ":").split(":")
    if len(parts) < 2:
        return "", None, ""
    chrom = parts[0] if parts[0].startswith("chr") else f"chr{parts[0]}"
    try:
        pos = int(float(parts[1]))
    except ValueError:
        pos = None
    return chrom, pos, parts[2] if len(parts) > 2 else ""


def compatible_strand(left: str, right: str) -> bool:
    return not left or left == "." or not right or right == "." or left == right


@dataclass(frozen=True)
class DnaEvent:
    row: dict[str, str]
    genes: tuple[str, str]
    chroms: tuple[str, str]
    positions: tuple[int, int]
    strands: tuple[str, str]


def dna_events(rows: list[dict[str, str]]) -> list[DnaEvent]:
    events: list[DnaEvent] = []
    for row in rows:
        genes = (clean_gene(row.get("gene1", "")), clean_gene(row.get("gene2", "")))
        try:
            positions = (int(float(row.get("pos1", ""))), int(float(row.get("pos2", ""))))
        except ValueError:
            continue
        chroms = tuple(c if c.startswith("chr") else f"chr{c}" for c in (row.get("chrom1", ""), row.get("chrom2", "")))
        if not all(genes) or not all(chroms):
            continue
        events.append(DnaEvent(row, genes, chroms, positions, (row.get("strand1", ""), row.get("strand2", ""))))
    return events


def ordered_candidate(dna: DnaEvent, genes: tuple[str, str], chroms: tuple[str, str]) -> tuple[DnaEvent, bool] | None:
    if dna.genes == genes and dna.chroms == chroms:
        return dna, False
    if dna.genes == genes[::-1] and dna.chroms == chroms[::-1]:
        return dna, True
    return None


def match_event(
    fusion: dict[str, str],
    candidates: list[DnaEvent],
    *,
    exact_tolerance: int,
) -> tuple[DnaEvent | None, str, str]:
    genes = gene_pair(fusion.get("gene_pair") or fusion.get("gene") or fusion.get("event_name") or "")
    left = breakpoint(fusion.get("left_breakpoint") or fusion.get("breakpoint1") or "")
    right = breakpoint(fusion.get("right_breakpoint") or fusion.get("breakpoint2") or "")
    if not genes or not left[0] or not right[0]:
        return None, "UNASSESSED", "RNA fusion lacks a two-sided gene pair or genomic breakpoint"
    ordered: list[tuple[DnaEvent, bool]] = []
    for item in candidates:
        matched = ordered_candidate(item, genes, (left[0], right[0]))
        if not matched:
            continue
        _, reversed_order = matched
        strands = item.strands[::-1] if reversed_order else item.strands
        if compatible_strand(strands[0], left[2]) and compatible_strand(strands[1], right[2]):
            ordered.append(matched)
    exact: list[DnaEvent] = []
    for item, reversed_order in ordered:
        pos = item.positions[::-1] if reversed_order else item.positions
        strand = item.strands[::-1] if reversed_order else item.strands
        if left[1] is None or right[1] is None:
            continue
        if (
            abs(pos[0] - left[1]) <= exact_tolerance
            and abs(pos[1] - right[1]) <= exact_tolerance
            and compatible_strand(strand[0], left[2])
            and compatible_strand(strand[1], right[2])
        ):
            exact.append(item)
    if len(exact) == 1:
        return exact[0], "EXACT_ADJACENCY", f"both breakends matched within {exact_tolerance} bp"
    if len(exact) > 1:
        return None, "AMBIGUOUS", "multiple DNA events match the exact RNA adjacency"
    if len(ordered) == 1:
        return ordered[0][0], "TRANSCRIPT_PROJECTION_UNIQUE", "unique PASS DNA adjacency for the ordered gene/chromosome pair; RNA exon and DNA intron coordinates are not equated"
    if len(ordered) > 1:
        return None, "AMBIGUOUS", "multiple DNA adjacencies share the RNA gene/chromosome pair; no evidence was borrowed"
    return None, "NOT_DETECTED", "no compatible PASS DNA adjacency"


DNA_FIELDS = [
    "dna_sv_confirmation_status", "dna_sv_event_id", "dna_sv_adjacency_key",
    "dna_sv_callers", "dna_sv_caller_count", "dna_sv_support_reads",
    "dna_sv_split_reads", "dna_sv_discordant_pairs", "dna_sv_match_method",
    "dna_sv_match_reason",
]


def annotation(dna: DnaEvent | None, method: str, reason: str) -> dict[str, str]:
    if dna is None:
        status = "AMBIGUOUS" if method == "AMBIGUOUS" else "UNASSESSED" if method == "UNASSESSED" else "NOT_DETECTED"
        return {**{field: "" for field in DNA_FIELDS}, "dna_sv_confirmation_status": status, "dna_sv_match_method": method, "dna_sv_match_reason": reason}
    row = dna.row
    return {
        "dna_sv_confirmation_status": "SUPPORTED",
        "dna_sv_event_id": row.get("sv_event_id") or row.get("event_id", ""),
        "dna_sv_adjacency_key": row.get("adjacency_key", ""),
        "dna_sv_callers": row.get("callers", ""),
        "dna_sv_caller_count": row.get("caller_count", ""),
        "dna_sv_support_reads": row.get("tumor_alt_support", ""),
        "dna_sv_split_reads": row.get("tumor_sr", ""),
        "dna_sv_discordant_pairs": row.get("tumor_pe", ""),
        "dna_sv_match_method": method,
        "dna_sv_match_reason": reason,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fusion-events", required=True, type=Path)
    parser.add_argument("--fusion-peptides", required=True, type=Path)
    parser.add_argument("--fusion-union", required=True, type=Path)
    parser.add_argument("--sv-events", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--exact-tolerance-bp", type=int, default=3)
    args = parser.parse_args()

    events = read_table(args.fusion_events)
    peptides = read_table(args.fusion_peptides)
    union = read_table(args.fusion_union)
    sv_rows = read_table(args.sv_events) if args.sv_events else []
    sv_candidates = dna_events(sv_rows)
    union_by_event = {row.get("event_id", ""): row for row in union if row.get("event_id")}
    annotations: dict[str, dict[str, str]] = {}
    links: list[dict[str, str]] = []
    for event in events:
        event_id = event.get("event_id", "")
        source = dict(event)
        source.update(union_by_event.get(event_id, {}))
        dna, method, reason = match_event(source, sv_candidates, exact_tolerance=args.exact_tolerance_bp)
        values = annotation(dna, method, reason)
        annotations[event_id] = values
        event.update(values)
        links.append({"fusion_event_id": event_id, **values})
    for peptide in peptides:
        peptide.update(annotations.get(peptide.get("event_id", ""), annotation(None, "UNASSESSED", "fusion event not present in normalized event table")))

    write_table(args.outdir / "raw_events.tsv", events)
    write_table(args.outdir / "raw_peptides.tsv", peptides)
    write_table(args.outdir / "dna_sv_rna_fusion_links.tsv", links, ["fusion_event_id", *DNA_FIELDS])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
