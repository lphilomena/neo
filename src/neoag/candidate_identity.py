"""Stable candidate identities used by ranking and report de-duplication."""

from __future__ import annotations

import hashlib
import re
from typing import Any, Mapping


IDENTITY_FIELDS = (
    "event_identity_id",
    "protein_change_identity_id",
    "peptide_sequence_id",
    "peptide_hla_id",
)


def _text(row: Mapping[str, Any], *fields: str) -> str:
    for field in fields:
        value = str(row.get(field) or "").strip()
        if value and value.upper() not in {"NA", "N/A", "NONE", ".", "UNASSESSED"}:
            return value
    return ""


def _stable_id(prefix: str, *parts: str) -> str:
    payload = "|".join(part.strip() for part in parts)
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]
    return f"{prefix}_{digest}"


def normalize_peptide(value: Any) -> str:
    return re.sub(r"\s+", "", str(value or "")).upper()


def normalize_hla(value: Any) -> str:
    text = str(value or "").strip().upper().replace("_", "-")
    text = re.sub(r"^HLA-?", "", text)
    match = re.match(r"^([A-Z0-9]+)\*?([0-9]{2,3}):?([0-9]{2,3})(?::?([0-9]{2,3}))?", text)
    if not match:
        return text
    locus, first, second, third = match.groups()
    suffix = f":{third}" if third else ""
    return f"HLA-{locus}*{first}:{second}{suffix}"


def candidate_identity(row: Mapping[str, Any]) -> dict[str, str]:
    """Return event, protein-change, peptide and peptide-HLA identities.

    The event identity remains event-specific. A shared peptide-HLA identity can
    therefore collapse duplicate ranking positions while preserving links to
    distinct breakpoints, transcript hypotheses or ORFs.
    """
    event_source = _text(
        row,
        "canonical_event_id",
        "event_group_id",
        "event_id",
        "source_event_id",
        "event_name",
    )
    if not event_source:
        event_source = "|".join(
            _text(row, field)
            for field in ("event_type", "gene", "chrom", "pos", "ref", "alt")
        )
    event_id = _stable_id("EVT", event_source or "UNRESOLVED_EVENT")

    protein_change = _text(
        row,
        "protein_change_id",
        "orf_id",
        "transcript_hypothesis_id",
        "combined_protein_change",
        "protein_change",
        "hgvsp",
        "amino_acid_change",
    )
    protein_id = _stable_id("PROT", event_id, protein_change or "UNRESOLVED_PROTEIN_CHANGE")

    peptide = normalize_peptide(_text(row, "peptide", "mutant_peptide", "mt_peptide"))
    peptide_id = _stable_id("PEP", peptide or "UNRESOLVED_PEPTIDE")
    hla = normalize_hla(_text(row, "hla_allele", "hla", "allele", "restricting_hla"))
    peptide_hla_id = _stable_id("PHLA", peptide or "UNRESOLVED_PEPTIDE", hla or "UNRESOLVED_HLA")
    return {
        "event_identity_id": event_id,
        "protein_change_identity_id": protein_id,
        "peptide_sequence_id": peptide_id,
        "peptide_hla_id": peptide_hla_id,
    }


def identity_value(row: Mapping[str, Any], field: str) -> str:
    value = str(row.get(field) or "").strip()
    return value or candidate_identity(row)[field]
