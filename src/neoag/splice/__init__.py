"""Splice-junction normalization and formal provenance primitives (v0.5.0)."""

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
from .identifiers import (
    orf_id,
    peptide_id,
    peptide_origin_id,
    splice_event_id,
    transcript_hypothesis_id,
)
from .normalization import normalize_splice_sources
from .pipeline import SpliceLayer, build_splice_provenance_layer
from .registry import (
    JunctionRegistry,
    JunctionResolution,
    JunctionSupportIndex,
    JunctionSupportMatch,
    build_support_index,
    resolve_junction_support,
)
from .schemas import SPLICE_PROVENANCE_SCHEMA_VERSION

__all__ = [
    "CANONICAL_JUNCTION_FIELDS", "CANONICAL_JUNCTION_SCHEMA_VERSION",
    "SPLICE_PROVENANCE_SCHEMA_VERSION", "CanonicalJunction", "JunctionNormalizationError",
    "JunctionSourceRecord", "JunctionRegistry", "JunctionResolution", "JunctionSupportIndex",
    "JunctionSupportMatch", "SpliceLayer", "build_support_index",
    "build_splice_provenance_layer", "canonical_junction_id", "iter_junction_records",
    "junction_record_from_row", "normalize_chromosome", "normalize_genome_build",
    "normalize_splice_sources", "orf_id", "parse_canonical_junction_id", "peptide_id",
    "peptide_origin_id", "read_source_rows", "resolve_junction_support", "splice_event_id",
    "transcript_hypothesis_id",
]
