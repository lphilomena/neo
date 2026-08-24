"""Resolve unstranded canonical junctions from an exactly matching GTF.

Only transcript-specific pairs of exon boundaries are accepted.  Gene-only,
nearest-coordinate, and single-boundary matches never assign a strand.
"""
from __future__ import annotations

import gzip
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from .coordinates import CanonicalJunction, JunctionSourceRecord, normalize_chromosome

_ATTR = re.compile(r'(\S+)\s+"([^"]*)"')


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.suffix == ".gz" else path.open(
        "r", encoding="utf-8", errors="replace"
    )


def _attrs(value: str) -> dict[str, str]:
    return {key: item for key, item in _ATTR.findall(value)}


def resolve_gtf_junction_strands(
    records: Iterable[JunctionSourceRecord], annotation_gtf: str | Path
) -> list[JunctionSourceRecord]:
    """Assign strand/transcript metadata to exact unstranded introns.

    A junction is resolved only when the GTF contains an exon ending at
    ``intron_start - 1`` and another exon from the same transcript starting at
    ``intron_end + 1``.  All matching transcripts must agree on strand.
    """
    rows = list(records)
    targets = {
        record.junction.unstranded_key
        for record in rows
        if record.junction is not None and record.junction.strand == "."
    }
    if not targets:
        return rows
    gtf = Path(annotation_gtf)
    if not gtf.is_file():
        raise FileNotFoundError(f"Missing matched annotation GTF: {gtf}")

    left: dict[tuple[str, str, int, int], dict[str, tuple[str, str, str]]] = defaultdict(dict)
    right: dict[tuple[str, str, int, int], dict[str, tuple[str, str, str]]] = defaultdict(dict)
    by_left: dict[tuple[str, int], list[tuple[str, str, int, int]]] = defaultdict(list)
    by_right: dict[tuple[str, int], list[tuple[str, str, int, int]]] = defaultdict(list)
    for key in targets:
        _, chrom, start, end = key
        by_left[(chrom, start - 1)].append(key)
        by_right[(chrom, end + 1)].append(key)

    with _open_text(gtf) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon" or fields[6] not in {"+", "-"}:
                continue
            chrom = normalize_chromosome(fields[0])
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                continue
            left_keys = by_left.get((chrom, end), ())
            right_keys = by_right.get((chrom, start), ())
            if not left_keys and not right_keys:
                continue
            attrs = _attrs(fields[8])
            transcript = attrs.get("transcript_id", "")
            if not transcript:
                continue
            meta = (fields[6], attrs.get("gene_id", ""), attrs.get("gene_name", ""))
            for key in left_keys:
                left[key][transcript] = meta
            for key in right_keys:
                right[key][transcript] = meta

    resolved: dict[tuple[str, str, int, int], tuple[str, str, str, str]] = {}
    for key in targets:
        transcripts = sorted(set(left.get(key, {})) & set(right.get(key, {})))
        strands = {left[key][tx][0] for tx in transcripts if left[key][tx] == right[key][tx]}
        transcripts = [tx for tx in transcripts if left[key][tx] == right[key][tx]]
        if transcripts and len(strands) == 1:
            genes = sorted({left[key][tx][1] for tx in transcripts if left[key][tx][1]})
            names = sorted({left[key][tx][2] for tx in transcripts if left[key][tx][2]})
            resolved[key] = (next(iter(strands)), ";".join(transcripts), ";".join(genes), ";".join(names))

    for record in rows:
        junction = record.junction
        if junction is None or junction.strand != ".":
            continue
        match = resolved.get(junction.unstranded_key)
        if not match:
            if junction.unstranded_key in targets:
                record.coordinate_warning = "Matched GTF did not yield one unambiguous transcript-specific strand."
            continue
        strand, transcripts, gene_ids, gene_names = match
        record.junction = CanonicalJunction(
            junction.genome_build, junction.chrom, junction.intron_start_1based,
            junction.intron_end_1based, strand,
        )
        record.transcript_ids = transcripts
        if not record.gene_id:
            record.gene_id = gene_ids
        if not record.gene:
            record.gene = gene_names
        record.resolution_status = "RESOLVED"
        record.resolution_method = "GTF_EXACT_TRANSCRIPT_EXON_BOUNDARIES"
        record.coordinate_warning = ""
    return rows
