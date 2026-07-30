"""Stable identifiers for the v0.5.0 Splice Provenance Layer.

Identifiers are content-derived and independent of input row order.  The readable
prefix exposes the entity class while the SHA-256 digest protects against unsafe
string concatenation and keeps identifiers compact enough for TSV workflows.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any, Iterable, Mapping


def _normalise(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, Mapping):
        return {str(k): _normalise(value[k]) for k in sorted(value, key=str)}
    if isinstance(value, (set, frozenset)):
        return sorted((_normalise(v) for v in value), key=lambda x: json.dumps(x, sort_keys=True))
    if isinstance(value, (list, tuple)):
        return [_normalise(v) for v in value]
    if isinstance(value, bool):
        return bool(value)
    return str(value).strip()


def stable_digest(*parts: Any, length: int = 24) -> str:
    """Return a deterministic SHA-256 prefix for structured content."""
    payload = json.dumps(_normalise(parts), sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:length]


def stable_id(prefix: str, *parts: Any, length: int = 24) -> str:
    return f"{prefix}|{stable_digest(*parts, length=length)}"


def splice_event_id(
    *,
    genome_build: str,
    event_type: str,
    strand: str,
    junction_ids: Iterable[str],
    gene: str = "",
    affected_exons: Iterable[str] = (),
) -> str:
    return stable_id(
        "SEV",
        genome_build,
        event_type.upper(),
        strand,
        sorted({x for x in junction_ids if x}),
        gene,
        sorted({x for x in affected_exons if x}),
    )


def unresolved_splice_event_id(*, source_tool: str, source_record_id: str, event_type: str = "UNRESOLVED") -> str:
    return stable_id("USEV", source_tool, source_record_id, event_type)


def transcript_hypothesis_id(
    *,
    splice_event_id_value: str,
    reference_transcript_id: str = "",
    exon_chain: Iterable[str] = (),
    junction_chain: Iterable[str] = (),
    path_role: str = "",
    sequence_sha256: str = "",
    source_generator: str = "",
) -> str:
    return stable_id(
        "STH",
        splice_event_id_value,
        reference_transcript_id,
        list(exon_chain),
        list(junction_chain),
        path_role,
        sequence_sha256,
        source_generator,
    )


def orf_id(
    *,
    transcript_hypothesis_id_value: str,
    protein_sequence_sha256: str,
    orf_start: Any = "",
    orf_stop: Any = "",
    frame_status: str = "",
) -> str:
    return stable_id(
        "ORF",
        transcript_hypothesis_id_value,
        protein_sequence_sha256,
        orf_start,
        orf_stop,
        frame_status,
    )


def peptide_id(sequence: str) -> str:
    return stable_id("PEP", sequence.upper())


def peptide_origin_id(
    *,
    orf_id_value: str,
    splice_event_id_value: str,
    peptide_sequence: str,
    protein_start: Any = "",
    protein_end: Any = "",
    junction_offset: Any = "",
) -> str:
    return stable_id(
        "POR",
        orf_id_value,
        splice_event_id_value,
        peptide_sequence.upper(),
        protein_start,
        protein_end,
        junction_offset,
    )


def link_id(prefix: str, *parts: Any) -> str:
    return stable_id(prefix, *parts)


def sequence_sha256(sequence: str) -> str:
    cleaned = "".join(str(sequence or "").split()).upper()
    return hashlib.sha256(cleaned.encode("ascii", errors="ignore")).hexdigest() if cleaned else ""
