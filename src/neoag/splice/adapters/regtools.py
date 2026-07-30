"""Adapter for exact junction sources (RegTools, STAR SJ, and strict generic TSV)."""
from __future__ import annotations

from pathlib import Path
from typing import Any

from neoag.splice.coordinates import iter_junction_records

from .base import row_hash


def parse_junction_source(
    path: str | Path,
    *,
    sample_id: str,
    source_tool: str = "RegTools",
    source_tool_version: str = "UNASSESSED",
    genome_build: str = "GRCh38",
    coordinate_system: str = "auto",
    strict: bool = False,
) -> dict[str, list[dict[str, str]]]:
    junctions: list[dict[str, str]] = []
    evidence: list[dict[str, str]] = []
    conflicts: list[dict[str, str]] = []
    records = list(
        iter_junction_records(
            path,
            sample_id=sample_id,
            source_tool=source_tool,
            genome_build=genome_build,
            coordinate_system=coordinate_system,
            strict=strict,
        )
    )
    for record in records:
        raw_hash = row_hash(record.row)
        if record.junction is None:
            entity_id = record.source_record_id
            conflicts.append({
                "entity_type": "JUNCTION",
                "entity_id": entity_id,
                "sample_id": sample_id,
                "conflict_type": "JUNCTION_COORDINATE_UNRESOLVED",
                "field_name": "coordinates",
                "observed_values": str(record.row),
                "source_tools": source_tool,
                "source_record_ids": record.source_record_id,
                "severity": "ERROR" if strict else "WARNING",
                "resolution_status": "UNRESOLVED",
                "resolution_reason": record.coordinate_warning or "No exact canonical junction could be derived.",
            })
            evidence.append({
                "entity_type": "JUNCTION", "entity_id": entity_id, "sample_id": sample_id,
                "evidence_group": "RNA_JUNCTION", "evidence_type": "UNRESOLVED_SOURCE_RECORD",
                "source_tool": source_tool, "source_tool_version": source_tool_version,
                "source_file": str(Path(path)), "source_row_number": str(record.source_row_number),
                "source_record_id": record.source_record_id, "provided_value": str(record.total_split_reads),
                "verified_value": "", "resolution_status": "UNRESOLVED",
                "resolution_reason": record.coordinate_warning or "coordinate resolution failed",
                "raw_payload_sha256": raw_hash,
            })
            continue
        j = record.junction
        exact_resolution = record.resolution_status == "RESOLVED"
        resolution_status = "RESOLVED_EXACT" if exact_resolution else record.resolution_status
        if not exact_resolution:
            conflicts.append({
                "entity_type": "JUNCTION", "entity_id": j.junction_id, "sample_id": sample_id,
                "conflict_type": "JUNCTION_STRAND_UNRESOLVED", "field_name": "strand",
                "observed_values": j.strand, "source_tools": source_tool,
                "source_record_ids": record.source_record_id, "severity": "WARNING",
                "resolution_status": "UNRESOLVED",
                "resolution_reason": "Unstranded coordinates are retained for review but cannot contribute exact junction support.",
            })
        junctions.append({
            "junction_id": j.junction_id,
            "sample_id": sample_id,
            "genome_build": j.genome_build,
            "chrom": j.chrom,
            "intron_start_1based": str(j.intron_start_1based),
            "intron_end_1based": str(j.intron_end_1based),
            "strand": j.strand,
            "donor_1based": str(j.donor_1based),
            "acceptor_1based": str(j.acceptor_1based),
            "splice_motif": record.splice_motif,
            "annotation_status": record.known_junction or "UNASSESSED",
            "unique_split_reads": str(record.unique_split_reads),
            "multi_split_reads": str(record.multi_split_reads),
            "total_split_reads": str(record.total_split_reads),
            "max_overhang": "",
            "source_coordinate_systems": record.source_coordinate_system,
            "source_tools": source_tool,
            "source_tool_versions": source_tool_version,
            "source_files": str(Path(path)),
            "source_record_ids": record.source_record_id,
            "provenance_record_count": "1",
            "junction_resolution_status": resolution_status,
            "evidence_conflict_status": "NONE" if exact_resolution else "JUNCTION_STRAND_UNRESOLVED",
        })
        evidence.append({
            "entity_type": "JUNCTION", "entity_id": j.junction_id, "sample_id": sample_id,
            "evidence_group": "RNA_JUNCTION",
            "evidence_type": "EXACT_SPLIT_READ_SUPPORT" if exact_resolution else "UNSTRANDED_SPLIT_READ_SUPPORT_UNVERIFIED",
            "source_tool": source_tool, "source_tool_version": source_tool_version,
            "source_file": str(Path(path)), "source_row_number": str(record.source_row_number),
            "source_record_id": record.source_record_id, "provided_value": str(record.total_split_reads),
            "verified_value": str(record.total_split_reads) if exact_resolution else "",
            "resolution_status": resolution_status,
            "resolution_reason": record.resolution_method if exact_resolution else "strand unavailable; exact support withheld",
            "raw_payload_sha256": raw_hash,
        })
    return {"junctions": junctions, "tool_evidence": evidence, "conflicts": conflicts}
