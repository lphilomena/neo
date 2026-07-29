"""Canonical splice-junction registry and non-leaking evidence resolution.

v0.4.4 removes all gene-level and nearest-locus read transfer. Junction support
may be transferred only through an exact canonical junction, an exact unique
source junction identifier, or an explicit unique variant-to-junction relation
recorded by the source caller (for example RegTools ``variant_info``).
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping

from neoag.utils import MISSING, first, to_float

from .coordinates import (
    CanonicalJunction,
    JunctionSourceRecord,
    iter_junction_records,
    normalize_chromosome,
    normalize_genome_build,
)

JUNCTION_ENTITY_FIELDS = [
    "junction_id",
    "sample_id",
    "genome_build",
    "chrom",
    "intron_start_1based",
    "intron_end_1based",
    "strand",
    "donor_1based",
    "acceptor_1based",
    "gene",
    "gene_ids",
    "transcript_ids",
    "rna_junction_reads",
    "read_count_values",
    "read_count_selection_rule",
    "source_junction_ids",
    "source_tools",
    "source_files",
    "source_record_ids",
    "provenance_record_count",
    "junction_resolution_status",
    "evidence_conflict_status",
]

SPLICE_TOOL_EVIDENCE_FIELDS = [
    "evidence_id",
    "entity_type",
    "evidence_role",
    "sample_id",
    "junction_id",
    "source_tool",
    "source_tool_version",
    "source_file",
    "source_row_number",
    "source_record_id",
    "source_row_sha256",
    "source_junction_id",
    "source_coordinate_system",
    "source_chrom",
    "source_start",
    "source_end",
    "genome_build",
    "chrom",
    "intron_start_1based",
    "intron_end_1based",
    "strand",
    "donor_1based",
    "acceptor_1based",
    "resolution_status",
    "resolution_method",
    "coordinate_warning",
    "gene",
    "gene_id",
    "transcript_ids",
    "variant_info",
    "variant_key",
    "rna_junction_reads",
    "peptide",
    "hla_allele",
    "raw_event_id",
]

SPLICE_PEPTIDE_PROVENANCE_FIELDS = [
    "evidence_id",
    "entity_type",
    "sample_id",
    "event_id",
    "peptide_id",
    "peptide",
    "hla_allele",
    "source_tool",
    "source_file",
    "source_record_id",
    "junction_id",
    "source_junction_id",
    "junction_match_status",
    "junction_match_method",
    "junction_support_status",
    "provided_rna_junction_reads",
    "resolved_rna_junction_reads",
    "junction_support_conflict",
    "source_row_sha256",
]

SPLICE_CONFLICT_FIELDS = [
    "evidence_domain",
    "record_id",
    "conflict_type",
    "details",
    "source_tool",
    "source_file",
    "source_row_number",
]


def _split_values(value: Any) -> set[str]:
    result: set[str] = set()
    raw = str(value or "").replace(",", ";")
    for part in raw.split(";"):
        token = part.strip()
        if token and token not in MISSING:
            result.add(token)
    return result


def _variant_key(value: Any, *, genome_build: str) -> str:
    raw = str(value or "").strip()
    if not raw:
        return ""
    # RegTools commonly emits chr:start-end. Keep the exact interval; do not
    # turn it into a nearby-locus search.
    if ":" not in raw:
        return ""
    chrom, interval = raw.split(":", 1)
    if "-" in interval:
        start, end = interval.split("-", 1)
    else:
        start = end = interval
    try:
        start_i, end_i = int(start), int(end)
    except ValueError:
        return ""
    if end_i < start_i:
        start_i, end_i = end_i, start_i
    return (
        f"VAR|{normalize_genome_build(genome_build)}|{normalize_chromosome(chrom)}|"
        f"{start_i}|{end_i}"
    )


def _row_variant_candidates(row: Mapping[str, Any], *, genome_build: str) -> list[str]:
    values: list[str] = []
    explicit = str(row.get("variant_key") or "").strip()
    if explicit:
        values.append(explicit)
    info = first(row, ["variant_info", "Variant Info", "variant_locus"], "")
    key = _variant_key(info, genome_build=genome_build)
    if key and key not in values:
        values.append(key)
    chrom = normalize_chromosome(row.get("chrom"))
    pos = int(to_float(row.get("pos"), 0.0))
    if chrom and pos > 0:
        values.append(f"VARPOS|{normalize_genome_build(genome_build)}|{chrom}|{pos}")
    return values


def unresolved_event_id(record: JunctionSourceRecord) -> str:
    digest = hashlib.sha256(
        f"{record.source_tool}|{record.source_file}|{record.source_row_number}|{record.record_sha256}".encode("utf-8")
    ).hexdigest()[:20]
    return f"UNRESOLVED_SJ|{record.source_tool}|{digest}"


@dataclass(frozen=True)
class JunctionResolution:
    junction: CanonicalJunction | None
    status: str
    method: str
    warning: str = ""

    @property
    def junction_id(self) -> str:
        return self.junction.junction_id if self.junction else ""


@dataclass
class JunctionEntity:
    junction: CanonicalJunction
    records: list[JunctionSourceRecord] = field(default_factory=list)
    source_tools: set[str] = field(default_factory=set)
    source_record_ids: set[str] = field(default_factory=set)
    source_junction_ids: set[str] = field(default_factory=set)
    genes: set[str] = field(default_factory=set)
    gene_ids: set[str] = field(default_factory=set)
    transcript_ids: set[str] = field(default_factory=set)

    def add(self, record: JunctionSourceRecord) -> None:
        self.records.append(record)
        if record.source_tool:
            self.source_tools.add(record.source_tool)
        if record.source_record_id:
            self.source_record_ids.add(record.source_record_id)
        if record.source_junction_id:
            self.source_junction_ids.add(record.source_junction_id)
        if record.gene and record.gene.upper() not in {"NA", "N/A", "."}:
            self.genes.add(record.gene)
        self.gene_ids.update(_split_values(record.gene_id))
        self.transcript_ids.update(_split_values(record.transcript_ids))

    def maximum_reads(self, source_tools: set[str] | None = None) -> int:
        allowed = {tool.casefold() for tool in source_tools} if source_tools else None
        return max(
            (
                record.total_split_reads
                for record in self.records
                if allowed is None or record.source_tool.casefold() in allowed
            ),
            default=0,
        )

    def as_dict(self, *, sample_id: str, primary_tools: set[str] | None = None) -> dict[str, str]:
        reads = self.maximum_reads(primary_tools)
        read_values = sorted(
            {
                record.total_split_reads
                for record in self.records
                if primary_tools is None or record.source_tool.casefold() in {tool.casefold() for tool in primary_tools}
            }
        )
        return {
            "junction_id": self.junction.junction_id,
            "sample_id": sample_id,
            "genome_build": self.junction.genome_build,
            "chrom": self.junction.chrom,
            "intron_start_1based": str(self.junction.intron_start_1based),
            "intron_end_1based": str(self.junction.intron_end_1based),
            "strand": self.junction.strand,
            "donor_1based": str(self.junction.donor_1based),
            "acceptor_1based": str(self.junction.acceptor_1based),
            "gene": ";".join(sorted(self.genes)),
            "gene_ids": ";".join(sorted(self.gene_ids)),
            "transcript_ids": ";".join(sorted(self.transcript_ids)),
            "rna_junction_reads": str(reads),
            "read_count_values": ";".join(str(value) for value in read_values),
            "read_count_selection_rule": "max_within_exact_canonical_junction",
            "source_junction_ids": ";".join(sorted(self.source_junction_ids)),
            "source_tools": ";".join(sorted(self.source_tools)),
            "source_files": ";".join(sorted({record.source_file for record in self.records})),
            "source_record_ids": ";".join(sorted(self.source_record_ids)),
            "provenance_record_count": str(len(self.records)),
            "junction_resolution_status": "RESOLVED",
            "evidence_conflict_status": "READ_COUNT_VARIATION" if len(read_values) > 1 else "NONE",
        }


class JunctionRegistry:
    """Registry of resolved junctions and exact aliases."""

    def __init__(self) -> None:
        self.entities: dict[str, JunctionEntity] = {}
        self.alias_to_ids: dict[str, set[str]] = defaultdict(set)
        self.unstranded_to_ids: dict[tuple[str, str, int, int], set[str]] = defaultdict(set)
        self.variant_to_ids: dict[str, set[str]] = defaultdict(set)
        self.conflicts: list[dict[str, str]] = []

    def add(self, record: JunctionSourceRecord, *, junction: CanonicalJunction | None = None) -> JunctionResolution:
        target = junction or record.junction
        if target is None:
            return JunctionResolution(None, "UNRESOLVED", "coordinates_missing", record.coordinate_warning)
        junction_id = target.junction_id
        entity = self.entities.setdefault(junction_id, JunctionEntity(target))
        entity.add(record)
        self.unstranded_to_ids[target.unstranded_key].add(junction_id)
        if record.source_junction_id:
            previous = set(self.alias_to_ids[record.source_junction_id])
            self.alias_to_ids[record.source_junction_id].add(junction_id)
            if previous and junction_id not in previous:
                self._record_conflict(
                    record,
                    "AMBIGUOUS_SOURCE_JUNCTION_ALIAS",
                    f"alias maps to {sorted(previous | {junction_id})}",
                )
        variant_key = _variant_key(
            first(record.row, ["variant_info", "Variant Info", "variant_locus"], ""),
            genome_build=target.genome_build,
        )
        if variant_key:
            self.variant_to_ids[variant_key].add(junction_id)
            _, build, chrom, start, end = variant_key.split("|", 4)
            self.variant_to_ids[f"VARPOS|{build}|{chrom}|{start}"].add(junction_id)
            self.variant_to_ids[f"VARPOS|{build}|{chrom}|{end}"].add(junction_id)
        return JunctionResolution(target, record.resolution_status, record.resolution_method, record.coordinate_warning)

    def add_many(self, records: Iterable[JunctionSourceRecord]) -> None:
        for record in records:
            self.add(record)

    def resolve(self, record: JunctionSourceRecord) -> JunctionResolution:
        direct = record.junction
        alias = record.source_junction_id
        alias_ids = set(self.alias_to_ids.get(alias, set())) if alias else set()

        if direct is not None:
            direct_id = direct.junction_id
            if direct_id in self.entities:
                if alias_ids and direct_id not in alias_ids:
                    warning = f"alias {alias!r} maps to {sorted(alias_ids)} but coordinates map to {direct_id}"
                    self._record_conflict(record, "ALIAS_COORDINATE_CONFLICT", warning)
                    return JunctionResolution(None, "CONFLICT", "alias_coordinate_conflict", warning)
                return JunctionResolution(direct, "RESOLVED", "EXACT_CANONICAL_COORDINATES")

            candidates = set(self.unstranded_to_ids.get(direct.unstranded_key, set()))
            if direct.strand == "." and len(candidates) == 1:
                target = self.entities[next(iter(candidates))].junction
                if alias_ids and target.junction_id not in alias_ids:
                    warning = f"alias {alias!r} conflicts with unique unstranded coordinate {target.junction_id}"
                    self._record_conflict(record, "ALIAS_COORDINATE_CONFLICT", warning)
                    return JunctionResolution(None, "CONFLICT", "alias_coordinate_conflict", warning)
                return JunctionResolution(
                    target,
                    "RESOLVED_WITH_CAUTION",
                    "UNIQUE_COORDINATE_WITHOUT_STRAND",
                    "source strand unavailable; unique stranded junction adopted",
                )
            if direct.strand == "." and len(candidates) > 1:
                warning = f"unstranded interval maps to multiple junctions: {sorted(candidates)}"
                self._record_conflict(record, "AMBIGUOUS_UNSTRANDED_COORDINATE", warning)
                return JunctionResolution(None, "AMBIGUOUS", "unstranded_coordinate_ambiguous", warning)

            if alias_ids and direct_id not in alias_ids:
                warning = f"alias {alias!r} maps to {sorted(alias_ids)} but row coordinates map to {direct_id}"
                self._record_conflict(record, "ALIAS_COORDINATE_CONFLICT", warning)
                return JunctionResolution(None, "CONFLICT", "alias_coordinate_conflict", warning)
            return JunctionResolution(direct, "RESOLVED_SOURCE_ONLY", "SOURCE_CANONICAL_COORDINATES")

        if alias:
            if len(alias_ids) == 1:
                target_id = next(iter(alias_ids))
                return JunctionResolution(self.entities[target_id].junction, "RESOLVED", "EXACT_SOURCE_JUNCTION_ALIAS")
            if len(alias_ids) > 1:
                warning = f"source alias {alias!r} maps to multiple canonical junctions: {sorted(alias_ids)}"
                self._record_conflict(record, "AMBIGUOUS_SOURCE_JUNCTION_ALIAS", warning)
                return JunctionResolution(None, "AMBIGUOUS", "source_alias_ambiguous", warning)

        return JunctionResolution(None, "UNRESOLVED", "NO_EXACT_JUNCTION_KEY", record.coordinate_warning)

    def entity_rows(self, *, sample_id: str, primary_tools: set[str] | None = None) -> list[dict[str, str]]:
        return [
            self.entities[junction_id].as_dict(sample_id=sample_id, primary_tools=primary_tools)
            for junction_id in sorted(self.entities)
        ]

    def _record_conflict(self, record: JunctionSourceRecord, conflict_type: str, details: str) -> None:
        self.conflicts.append(
            {
                "evidence_domain": "splice_junction",
                "record_id": record.source_record_id,
                "conflict_type": conflict_type,
                "details": details,
                "source_tool": record.source_tool,
                "source_file": record.source_file,
                "source_row_number": str(record.source_row_number),
            }
        )


@dataclass
class JunctionSupportIndex:
    registry: JunctionRegistry
    primary_tools: set[str] = field(default_factory=lambda: {"RegTools"})

    @property
    def by_canonical_id(self) -> dict[str, JunctionEntity]:
        return self.registry.entities

    @property
    def by_source_junction_id(self) -> dict[str, set[str]]:
        return self.registry.alias_to_ids

    @property
    def by_variant_key(self) -> dict[str, set[str]]:
        return self.registry.variant_to_ids


@dataclass(frozen=True)
class JunctionSupportMatch:
    junction_id: str = ""
    source_junction_id: str = ""
    selected_reads: int = 0
    provided_reads: int = 0
    match_status: str = "UNRESOLVED"
    match_method: str = "none"
    support_status: str = "UNRESOLVED"
    conflict: str = "NONE"
    reason: str = "No exact junction support link was found."


def build_support_index(
    splice_path: str | Path,
    *,
    sample_id: str = "",
    genome_build: str = "GRCh38",
    coordinate_system: str = "auto",
    source_tool: str = "RegTools",
) -> JunctionSupportIndex:
    registry = JunctionRegistry()
    for record in iter_junction_records(
        splice_path,
        sample_id=sample_id,
        source_tool=source_tool,
        genome_build=genome_build,
        coordinate_system=coordinate_system,
    ):
        registry.add(record)
    return JunctionSupportIndex(registry=registry, primary_tools={source_tool})


def _provided_reads(peptide: Mapping[str, Any], event: Mapping[str, Any] | None) -> int:
    for source in (peptide, event or {}):
        raw = first(source, ["provided_rna_junction_reads", "rna_junction_reads", "RNA Junction Reads"], "")
        if raw and raw not in MISSING:
            return max(0, int(to_float(raw, 0.0)))
    return 0


def _candidate_canonical_ids(peptide: Mapping[str, Any], event: Mapping[str, Any] | None) -> list[str]:
    values: list[str] = []
    for source in (peptide, event or {}):
        for field in ("canonical_junction_id", "junction_id", "event_id"):
            value = str(source.get(field) or "").strip()
            if value.startswith("SJ|") and value not in values:
                values.append(value)
    return values


def _candidate_source_ids(peptide: Mapping[str, Any], event: Mapping[str, Any] | None) -> list[str]:
    values: list[str] = []
    for source in (peptide, event or {}):
        for field in ("source_junction_id", "Splice Junction", "splice_junction", "event_name", "junction_id"):
            value = str(source.get(field) or "").strip()
            if value and not value.startswith("SJ|") and value not in values:
                values.append(value)
    return values


def _coordinate_id(peptide: Mapping[str, Any], event: Mapping[str, Any] | None) -> str:
    merged = {**(event or {}), **peptide}
    chrom = normalize_chromosome(
        merged.get("junction_chrom") or merged.get("chrom")
    )
    start = int(
        to_float(
            merged.get("junction_start") or merged.get("intron_start_1based"),
            0.0,
        )
    )
    end = int(
        to_float(
            merged.get("junction_end") or merged.get("intron_end_1based"),
            0.0,
        )
    )
    strand = str(
        merged.get("junction_strand") or merged.get("strand") or "."
    ).strip()
    if not chrom or start < 1 or end < start or strand not in {"+", "-"}:
        return ""
    return CanonicalJunction(
        normalize_genome_build(merged.get("genome_build")), chrom, start, end, strand
    ).junction_id


def _match_entity(
    index: JunctionSupportIndex,
    entity_ids: set[str],
    *,
    provided: int,
    source_id: str,
    method: str,
    status: str = "EXACT",
) -> JunctionSupportMatch:
    if len(entity_ids) != 1:
        ambiguous = bool(entity_ids)
        return JunctionSupportMatch(
            source_junction_id=source_id,
            selected_reads=0,
            provided_reads=provided,
            match_status="AMBIGUOUS" if ambiguous else "UNRESOLVED",
            match_method=method,
            support_status="AMBIGUOUS" if ambiguous else ("PROVIDED_UNVERIFIED" if provided else "UNRESOLVED"),
            conflict="AMBIGUOUS_LINK" if ambiguous else "NONE",
            reason=(
                f"{len(entity_ids)} exact candidates matched; no read count was transferred."
                if ambiguous
                else "No exact candidate matched; any upstream count remains provenance-only."
            ),
        )
    junction_id = next(iter(entity_ids))
    entity = index.registry.entities[junction_id]
    stranded = entity.junction.strand in {"+", "-"}
    resolved = entity.maximum_reads(index.primary_tools) if stranded else 0
    conflict = "NONE" if not provided or provided == resolved else f"PROVIDED_{provided}_NE_RESOLVED_{resolved}"
    return JunctionSupportMatch(
        junction_id=junction_id,
        source_junction_id=source_id,
        selected_reads=resolved,
        provided_reads=provided,
        match_status=status,
        match_method=method,
        support_status=(
            "SUPPORTED_EXACT_JUNCTION"
            if resolved > 0
            else "MATCHED_ZERO_READS"
            if stranded
            else "UNSTRANDED_PRIMARY_UNVERIFIED"
        ),
        conflict=conflict,
        reason=(
            "Reads were transferred only from the uniquely matched canonical junction entity."
            if stranded
            else "The exact interval is unstranded; its count is retained as provenance and not transferred as verified support."
        ),
    )


def resolve_junction_support(
    peptide: Mapping[str, Any],
    index: JunctionSupportIndex,
    *,
    event: Mapping[str, Any] | None = None,
) -> JunctionSupportMatch:
    provided = _provided_reads(peptide, event)

    for junction_id in _candidate_canonical_ids(peptide, event):
        if junction_id in index.registry.entities:
            return _match_entity(
                index,
                {junction_id},
                provided=provided,
                source_id="",
                method="exact_canonical_junction_id",
            )

    coordinate = _coordinate_id(peptide, event)
    if coordinate and coordinate in index.registry.entities:
        return _match_entity(
            index,
            {coordinate},
            provided=provided,
            source_id="",
            method="exact_build_chrom_intron_strand",
        )

    for source_id in _candidate_source_ids(peptide, event):
        if source_id in index.registry.alias_to_ids:
            return _match_entity(
                index,
                set(index.registry.alias_to_ids[source_id]),
                provided=provided,
                source_id=source_id,
                method="exact_source_junction_id",
            )

    merged = {**(event or {}), **peptide}
    for variant_key in _row_variant_candidates(
        merged,
        genome_build=normalize_genome_build(merged.get("genome_build")),
    ):
        if variant_key in index.registry.variant_to_ids:
            return _match_entity(
                index,
                set(index.registry.variant_to_ids[variant_key]),
                provided=provided,
                source_id="",
                method="unique_explicit_variant_to_junction_link",
                status="UNIQUE_VARIANT_LINK",
            )

    if provided:
        return JunctionSupportMatch(
            selected_reads=0,
            provided_reads=provided,
            match_status="PROVIDED_UNVERIFIED",
            match_method="upstream_provided_only",
            support_status="PROVIDED_UNVERIFIED",
            reason="The upstream value is retained only in provided_rna_junction_reads; it is not used as verified RNA support.",
        )
    return JunctionSupportMatch()


def peptide_provenance_row(
    peptide: Mapping[str, Any],
    match: JunctionSupportMatch,
) -> dict[str, str]:
    payload = json.dumps(
        {str(key): "" if value is None else str(value) for key, value in sorted(peptide.items())},
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    material = f"{peptide.get('peptide_id')}|{peptide.get('source_tool')}|{digest}"
    return {
        "evidence_id": f"EVID|PEP|{hashlib.sha256(material.encode('utf-8')).hexdigest()[:20]}",
        "entity_type": "splice_peptide",
        "sample_id": str(peptide.get("sample_id") or ""),
        "event_id": str(peptide.get("event_id") or ""),
        "peptide_id": str(peptide.get("peptide_id") or ""),
        "peptide": str(peptide.get("peptide") or ""),
        "hla_allele": str(peptide.get("hla_allele") or ""),
        "source_tool": str(peptide.get("source_tool") or ""),
        "source_file": str(peptide.get("source_file") or ""),
        "source_record_id": str(peptide.get("source_record_id") or ""),
        "junction_id": match.junction_id,
        "source_junction_id": match.source_junction_id,
        "junction_match_status": match.match_status,
        "junction_match_method": match.match_method,
        "junction_support_status": match.support_status,
        "provided_rna_junction_reads": str(match.provided_reads),
        "resolved_rna_junction_reads": str(match.selected_reads),
        "junction_support_conflict": match.conflict,
        "source_row_sha256": digest,
    }
