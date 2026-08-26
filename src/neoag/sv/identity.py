"""Canonical identity primitives for structural-variant adjacencies.

The key deliberately excludes gene names.  Gene annotations can change between
releases and two distinct breakpoints may involve the same gene pair.  Evidence
may only be joined through this genomic identity (or an event_id derived from it).
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass


def normalize_contig(chrom: str) -> str:
    value = str(chrom or "").strip()
    return value[3:] if value.lower().startswith("chr") else value


@dataclass(frozen=True, order=True)
class Breakend:
    chrom: str
    pos: int
    strand: str = "."

    @classmethod
    def create(cls, chrom: str, pos: int, strand: str = ".") -> "Breakend":
        strand = strand if strand in {"+", "-"} else "."
        return cls(normalize_contig(chrom), int(pos), strand)

    def token(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.strand}"


def canonical_breakpoint_key(
    genome_build: str,
    chrom1: str,
    pos1: int,
    strand1: str,
    chrom2: str,
    pos2: int,
    strand2: str,
) -> str:
    """Return an order-independent key for one physical adjacency."""

    build = str(genome_build or "").strip().lower()
    if not build:
        raise ValueError("genome_build is required for exact SV evidence identity")
    left = Breakend.create(chrom1, pos1, strand1)
    right = Breakend.create(chrom2, pos2, strand2)
    first, second = sorted((left, right))
    return f"{build}|{first.token()}|{second.token()}"


def stable_identifier(prefix: str, *parts: object, length: int = 16) -> str:
    payload = "\x1f".join(str(part) for part in parts)
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:length]
    return f"{prefix}_{digest}"

