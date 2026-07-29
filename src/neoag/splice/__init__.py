"""Splice-junction normalization and provenance primitives (v0.4.4)."""

from .coordinates import (
    CANONICAL_JUNCTION_FIELDS,
    CANONICAL_JUNCTION_SCHEMA_VERSION,
    CanonicalJunction,
    JunctionNormalizationError,
    JunctionSourceRecord,
    canonical_junction_id,
    iter_junction_records,
    junction_record_from_row,
    normalize_chromosome,
    normalize_genome_build,
    parse_canonical_junction_id,
    read_source_rows,
)
from .normalization import normalize_splice_sources
from .registry import (
    JunctionRegistry,
    JunctionResolution,
    JunctionSupportIndex,
    JunctionSupportMatch,
    build_support_index,
    resolve_junction_support,
)

__all__ = [
    "CANONICAL_JUNCTION_FIELDS",
    "CANONICAL_JUNCTION_SCHEMA_VERSION",
    "CanonicalJunction",
    "JunctionNormalizationError",
    "JunctionSourceRecord",
    "JunctionRegistry",
    "JunctionResolution",
    "JunctionSupportIndex",
    "JunctionSupportMatch",
    "build_support_index",
    "canonical_junction_id",
    "iter_junction_records",
    "junction_record_from_row",
    "normalize_chromosome",
    "normalize_genome_build",
    "normalize_splice_sources",
    "parse_canonical_junction_id",
    "read_source_rows",
    "resolve_junction_support",
]
