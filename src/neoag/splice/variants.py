"""Exact variant and junction identifiers used by the DNA-causal branch."""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Mapping

from .coordinates import CanonicalJunction, normalize_chromosome, normalize_genome_build, normalize_strand

_DNA = re.compile(r"^[ACGTN]+$", re.I)
_JUNCTION = re.compile(
    r"^(?P<chrom>(?:chr)?(?:[0-9]+|X|Y|M|MT))[:_](?P<start>\d+)[-:](?P<end>\d+)[:_](?P<strand>[+\-.])$",
    re.I,
)
_VARIANT = re.compile(
    r"^(?P<chrom>(?:chr)?(?:[0-9]+|X|Y|M|MT))[:_](?P<pos>\d+)[:_](?P<ref>[ACGTN.-]+)[:_](?P<alt>[ACGTN.-]+)$",
    re.I,
)


def clean_allele(value: Any) -> str:
    allele = str(value or "").strip().upper()
    return "-" if allele in {"", ".", "-"} else allele


def canonical_variant_id(genome_build: str, chrom: str, pos_1based: int, ref: str, alt: str) -> str:
    build = normalize_genome_build(genome_build)
    c = normalize_chromosome(chrom)
    pos = int(pos_1based)
    r, a = clean_allele(ref), clean_allele(alt)
    if not build or not c or pos < 1 or r == a:
        raise ValueError(f"Invalid variant identity: {genome_build!r} {chrom!r} {pos!r} {ref!r}>{alt!r}")
    if r != "-" and not _DNA.fullmatch(r):
        raise ValueError(f"Invalid REF allele: {r}")
    if a != "-" and not _DNA.fullmatch(a):
        raise ValueError(f"Invalid ALT allele: {a}")
    return f"VAR|{build}|{c}|{pos}|{r}|{a}"


def variant_type(ref: str, alt: str) -> str:
    r, a = clean_allele(ref), clean_allele(alt)
    if r == "-" or len(a) > len(r):
        return "INSERTION" if r == "-" else "INDEL"
    if a == "-" or len(r) > len(a):
        return "DELETION" if a == "-" else "INDEL"
    return "SNV" if len(r) == len(a) == 1 else "MNV"


def parse_variant_token(value: str, *, genome_build: str = "GRCh38") -> tuple[str, str, int, str, str] | None:
    text = str(value or "").strip()
    if text.startswith("VAR|"):
        parts = text.split("|")
        if len(parts) == 6:
            return parts[1], parts[2], int(parts[3]), parts[4], parts[5]
    match = _VARIANT.match(text)
    if not match:
        return None
    return (
        normalize_genome_build(genome_build),
        normalize_chromosome(match.group("chrom")),
        int(match.group("pos")),
        clean_allele(match.group("ref")),
        clean_allele(match.group("alt")),
    )


def parse_junction_token(value: str, *, genome_build: str = "GRCh38") -> CanonicalJunction | None:
    text = str(value or "").strip()
    if text.startswith("SJ|"):
        parts = text.split("|")
        if len(parts) == 6:
            return CanonicalJunction(parts[1], parts[2], int(parts[3]), int(parts[4]), parts[5])
    match = _JUNCTION.match(text)
    if not match:
        return None
    return CanonicalJunction(
        normalize_genome_build(genome_build),
        normalize_chromosome(match.group("chrom")),
        min(int(match.group("start")), int(match.group("end"))),
        max(int(match.group("start")), int(match.group("end"))),
        normalize_strand(match.group("strand")),
    )


def variant_from_row(
    row: Mapping[str, Any],
    *,
    genome_build: str = "GRCh38",
    zero_based_half_open: bool = False,
) -> tuple[str, str, int, str, str] | None:
    def first(*names: str) -> str:
        for name in names:
            if name in row and str(row.get(name, "")).strip():
                return str(row.get(name, "")).strip()
        return ""

    token = first("variant_id", "canonical_variant_id", "variant_key")
    parsed = parse_variant_token(token, genome_build=genome_build) if token else None
    if parsed:
        return parsed
    chrom = first("Chromosome", "chromosome", "chrom", "chr")
    pos_text = first("Start", "start", "pos", "position", "POS")
    ref = first("Reference", "reference", "ref", "REF")
    alt = first("Variant", "variant", "alt", "ALT")
    if not all((chrom, pos_text, ref, alt)):
        return None
    try:
        pos = int(float(pos_text)) + (1 if zero_based_half_open else 0)
    except Exception:
        return None
    return normalize_genome_build(genome_build), normalize_chromosome(chrom), pos, clean_allele(ref), clean_allele(alt)


@dataclass(frozen=True)
class ExactEventIndex:
    """Resolve event identity only from exact junction sets.

    Gene names are used only as a tie-breaker among already identical junction
    sets. They are never used to transfer reads or link unrelated events.
    """

    by_junction_set: Mapping[frozenset[str], tuple[str, ...]]
    event_gene: Mapping[str, str]

    @classmethod
    def from_tables(cls, tables: Mapping[str, list[Mapping[str, str]]]) -> "ExactEventIndex":
        grouped: dict[frozenset[str], list[str]] = {}
        genes: dict[str, str] = {}
        for event in tables.get("events", []):
            event_id = str(event.get("splice_event_id", ""))
            jids = frozenset(x for x in str(event.get("junction_ids", "")).split(";") if x)
            if event_id and jids:
                grouped.setdefault(jids, []).append(event_id)
                genes[event_id] = str(event.get("gene", ""))
        return cls({key: tuple(sorted(set(value))) for key, value in grouped.items()}, genes)

    def resolve(self, junction_ids: list[str], *, gene: str = "") -> tuple[str, str]:
        key = frozenset(x for x in junction_ids if x)
        candidates = list(self.by_junction_set.get(key, ()))
        if len(candidates) == 1:
            return candidates[0], "RESOLVED_EXACT_JUNCTION_SET"
        if len(candidates) > 1 and gene:
            exact_gene = [event_id for event_id in candidates if self.event_gene.get(event_id) == gene]
            if len(exact_gene) == 1:
                return exact_gene[0], "RESOLVED_EXACT_JUNCTION_SET_AND_GENE"
        return "", "AMBIGUOUS_EXACT_JUNCTION_SET" if candidates else "NO_EXACT_EVENT"
