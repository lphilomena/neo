"""v0.4.4 exact splice-junction normalization and provenance pipeline.

The production invariant is deliberately narrow: verified RNA junction support
may cross a source boundary only after an exact canonical junction has been
resolved.  A canonical junction is identified by genome build, chromosome,
1-based closed intron start/end, and strand.  Gene-level, nearest-locus, and
"largest junction in the gene" fallbacks are forbidden.

The module emits backward-compatible ``raw_events.tsv``/``raw_peptides.tsv``
plus auditable entity, provenance, conflict, consensus, and QC tables.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from neoag.model_layers import enrich_event_layers, enrich_peptide_layers
from neoag.provenance import (
    CONFLICT_FIELDS,
    PROVENANCE_FIELDS,
    merge_rows_preserving_provenance,
)
from neoag.schemas import EVENT_FIELDS, PEPTIDE_FIELDS, RNA_JUNCTION_EVIDENCE_FIELDS
from neoag.utils import first, safe_id, to_float, write_json, write_tsv

from .coordinates import (
    CANONICAL_JUNCTION_SCHEMA_VERSION,
    JunctionSourceRecord,
    file_sha256,
    iter_junction_records,
    peptide_metadata,
)
from .registry import (
    JUNCTION_ENTITY_FIELDS,
    SPLICE_CONFLICT_FIELDS,
    SPLICE_PEPTIDE_PROVENANCE_FIELDS,
    SPLICE_TOOL_EVIDENCE_FIELDS,
    JunctionRegistry,
    JunctionResolution,
    peptide_provenance_row,
    unresolved_event_id,
)


@dataclass(frozen=True)
class SpliceSource:
    tool: str
    path: Path
    role: str
    coordinate_system: str = "auto"
    version: str = "UNASSESSED"


@dataclass(frozen=True)
class NormalizedRecord:
    record: JunctionSourceRecord
    resolution: JunctionResolution
    role: str
    event_id: str


PRIMARY_JUNCTION_TOOLS = {"RegTools"}

CONSENSUS_PROVENANCE_FIELDS = [
    "event_id",
    "canonical_junction_id",
    "evidence_domain",
    "tool",
    "source_file",
    "source_row_number",
    "source_record_id",
    "source_junction_id",
    "resolution_status",
    "resolution_method",
    "coordinate_warning",
    "peptide_present",
]


def _tokens(value: Any) -> list[str]:
    result: list[str] = []
    for token in str(value or "").replace(",", ";").split(";"):
        item = token.strip()
        if item and item not in result:
            result.append(item)
    return result


def _join(values: Iterable[str]) -> str:
    return ";".join(sorted({value for value in values if value}))


def _source_record_fields(record: JunctionSourceRecord) -> dict[str, str]:
    return {
        "source_file": record.source_file,
        "source_row_number": str(record.source_row_number),
        "source_record_id": record.source_record_id,
        "source_tools": record.source_tool,
        "source_records": record.source_record_id,
        "provenance_record_count": "1",
        "evidence_conflict_status": "NONE",
    }


def _junction_fields(item: NormalizedRecord) -> dict[str, str]:
    junction = item.resolution.junction
    if junction is None:
        return {
            "genome_build": "",
            "canonical_junction_id": "",
            "source_junction_id": item.record.source_junction_id,
            "junction_chrom": "",
            "junction_start": "",
            "junction_end": "",
            "junction_strand": "",
            "junction_donor": "",
            "junction_acceptor": "",
            "junction_coordinate_system": item.record.source_coordinate_system,
            "junction_resolution_status": item.resolution.status,
            "junction_resolution_reason": item.resolution.method
            + (f": {item.resolution.warning}" if item.resolution.warning else ""),
        }
    return {
        "genome_build": junction.genome_build,
        "canonical_junction_id": junction.junction_id,
        "source_junction_id": item.record.source_junction_id,
        "junction_chrom": junction.chrom,
        "junction_start": str(junction.intron_start_1based),
        "junction_end": str(junction.intron_end_1based),
        "junction_strand": junction.strand,
        "junction_donor": str(junction.donor_1based),
        "junction_acceptor": str(junction.acceptor_1based),
        "junction_coordinate_system": "intron_1based_closed",
        "junction_resolution_status": item.resolution.status,
        "junction_resolution_reason": item.resolution.method
        + (f": {item.resolution.warning}" if item.resolution.warning else ""),
    }


def _entity_primary_support(
    item: NormalizedRecord,
    registry: JunctionRegistry,
    primary_tools: set[str],
) -> tuple[int, set[str]]:
    junction = item.resolution.junction
    if junction is None or junction.junction_id not in registry.entities:
        return 0, set()
    folded = {tool.casefold() for tool in primary_tools}
    entity = registry.entities[junction.junction_id]
    # Strand is part of the canonical identity. An unstranded primary row is
    # retained as provenance but cannot donate verified support across tools.
    if entity.junction.strand not in {"+", "-"}:
        return 0, set()
    tools = {
        record.source_tool
        for record in entity.records
        if record.source_tool.casefold() in folded
    }
    if not tools:
        return 0, set()
    return entity.maximum_reads(primary_tools), tools


def _support_fields(
    item: NormalizedRecord,
    registry: JunctionRegistry,
    primary_tools: set[str],
) -> dict[str, str]:
    verified, primary_sources = _entity_primary_support(item, registry, primary_tools)
    provided = item.record.total_split_reads
    resolved = item.resolution.junction is not None
    if primary_sources:
        status = "SUPPORTED_EXACT_JUNCTION" if verified > 0 else "MATCHED_ZERO_READS"
        match_status = "EXACT"
        reason = "Verified reads were transferred only from the exact canonical primary junction entity."
    elif resolved:
        status = "RESOLVED_WITHOUT_PRIMARY_SUPPORT"
        match_status = "RESOLVED_SOURCE_ONLY"
        reason = "Coordinates were normalized, but no exact primary RNA-junction record supports this entity."
    elif provided > 0:
        status = "PROVIDED_UNVERIFIED"
        match_status = item.resolution.status or "UNRESOLVED"
        reason = "The caller-provided count is retained only as provenance and is excluded from verified RNA support."
    else:
        status = "UNRESOLVED"
        match_status = item.resolution.status or "UNRESOLVED"
        reason = "No exact canonical junction or verified primary RNA support was found."

    conflict = "NONE"
    if primary_sources and provided > 0 and provided != verified:
        conflict = f"PROVIDED_{provided}_NE_RESOLVED_{verified}"
    return {
        "junction_match_status": match_status,
        "junction_match_method": item.resolution.method,
        "junction_support_status": status,
        "junction_support_conflict": conflict,
        "junction_support_reason": reason,
        "provided_rna_junction_reads": str(provided),
        "rna_junction_reads": str(verified),
        "rna_junction_source": _join(primary_sources),
    }


def _normal_status(
    event_id: str,
    *,
    normal_file_declared: bool,
    normal_resolvable_rows: int,
    normal_detected: set[str],
) -> str:
    if not normal_file_declared:
        return "UNASSESSED"
    if event_id in normal_detected:
        return "DETECTED"
    if normal_resolvable_rows > 0:
        # A panel-level absence is not proof of adequate locus coverage.
        return "NOT_DETECTED_COVERAGE_UNASSESSED"
    return "UNASSESSED"


def _normal_specificity(status: str) -> str:
    if status == "DETECTED":
        return "0.1"
    if status == "NOT_DETECTED_COVERAGE_UNASSESSED":
        return "0.6"
    return "0.5"


def _event_source_row(
    item: NormalizedRecord,
    *,
    sample_id: str,
    profile_name: str,
    registry: JunctionRegistry,
    primary_tools: set[str],
    normal_status: str,
) -> dict[str, str]:
    record = item.record
    junction = item.resolution.junction
    gene = record.gene.strip() or (
        f"{junction.chrom}:{junction.intron_start_1based}" if junction else f"{record.source_tool}:row{record.source_row_number}"
    )
    verified, primary_sources = _entity_primary_support(item, registry, primary_tools)
    if junction is not None and primary_sources:
        confidence = "0.85" if len(registry.entities[junction.junction_id].source_tools) >= 2 else "0.70"
    elif junction is not None:
        confidence = "0.45"
    else:
        confidence = "0.20"

    row: dict[str, str] = {
        "event_id": item.event_id,
        "sample_id": sample_id,
        "disease_profile": profile_name,
        "event_type": "Splice",
        "mutation_source": "Other",
        "peptide_consequence": "splice_junction",
        "gene": gene,
        "event_name": record.source_junction_id or item.event_id,
        # Caller coordinate is retained for compatibility. Canonical coordinate
        # fields below are authoritative for any cross-source join.
        "chrom": record.source_chrom or (junction.chrom if junction else ""),
        "pos": record.source_start or (str(junction.intron_start_1based) if junction else ""),
        "transcript_id": _tokens(record.transcript_ids)[0] if _tokens(record.transcript_ids) else "",
        "consequence": "splice_junction",
        "normal_junction_status": normal_status,
        "event_confidence": confidence,
        "event_expression": "0.0",
        "driver_relevance": "0.0",
        "tumor_vaf": "0.0",
        "clonality": "0.5",
        "persistence": "0.5",
        "tumor_specificity": _normal_specificity(normal_status),
        "source": f"splice_source:{record.source_tool}",
        **_junction_fields(item),
        **_support_fields(item, registry, primary_tools),
        **_source_record_fields(record),
    }
    # Do not accidentally signal an assessed negative in the evidence layer.
    if not primary_sources:
        row["rna_junction_source"] = ""
        row["rna_junction_reads"] = "0"
    return enrich_event_layers(row)


def _mhc_class(hla: str) -> str:
    if not hla:
        return ""
    return "II" if any(token in hla.upper() for token in ("DR", "DQ", "DP")) else "I"


def _peptide_source_row(
    item: NormalizedRecord,
    *,
    sample_id: str,
    registry: JunctionRegistry,
    primary_tools: set[str],
    event: dict[str, str],
) -> tuple[dict[str, str] | None, dict[str, str] | None, dict[str, str] | None]:
    metadata = peptide_metadata(item.record)
    peptide = metadata["peptide"].strip().upper()
    if not peptide:
        return None, None, None
    hla = metadata["hla_allele"].strip()
    support = _support_fields(item, registry, primary_tools)
    peptide_id = safe_id(f"{item.event_id}_{hla or 'HLA_UNASSESSED'}_{peptide}")
    row: dict[str, str] = {
        "peptide_id": peptide_id,
        "event_id": item.event_id,
        "sample_id": sample_id,
        "event_type": "Splice",
        "mutation_source": "Other",
        "peptide_consequence": "splice_junction",
        "gene": event.get("gene", item.record.gene),
        "peptide": peptide,
        "wildtype_peptide": first(
            item.record.row,
            ["wildtype_peptide", "WT Epitope Seq", "WT Epitope", "wt_peptide"],
            "",
        ),
        "crosses_junction": (
            "yes"
            if str(metadata["crosses_junction"]).strip().lower() in {"1", "true", "yes", "y", "pass"}
            else str(metadata["crosses_junction"] or "")
        ),
        "contains_novel_aa": first(
            item.record.row,
            ["contains_novel_aa", "Contains Novel AA"],
            "",
        ),
        "hla_allele": hla,
        "mhc_class": _mhc_class(hla),
        "source_tool": item.record.source_tool,
        "generation_status": metadata["generation_status"],
        "binding_rank": str(
            to_float(
                first(
                    item.record.row,
                    ["binding_rank", "Best MT Score", "Median MT Score", "MT %Rank", "%ile MT"],
                    "99",
                ),
                99.0,
            )
        ),
        "el_rank": str(
            to_float(
                first(item.record.row, ["el_rank", "EL Rank", "Best MT EL Score"], "99"),
                99.0,
            )
        ),
        "presentation_score": first(item.record.row, ["presentation_score", "Presentation Score"], "0.0"),
        "immunogenicity_score": first(
            item.record.row,
            ["immunogenicity_score", "Immunogenicity Score"],
            "0.5",
        ),
        **_junction_fields(item),
        **support,
        **_source_record_fields(item.record),
    }
    # Verified evidence is zero unless an exact primary entity was resolved.
    if support["junction_support_status"] not in {"SUPPORTED_EXACT_JUNCTION", "MATCHED_ZERO_READS"}:
        row["rna_junction_reads"] = "0"
        row["rna_junction_source"] = ""

    enriched = enrich_peptide_layers(row, event)
    match_like = type(
        "JunctionMatch",
        (),
        {
            "junction_id": enriched.get("canonical_junction_id", ""),
            "source_junction_id": enriched.get("source_junction_id", ""),
            "match_status": enriched.get("junction_match_status", ""),
            "match_method": enriched.get("junction_match_method", ""),
            "support_status": enriched.get("junction_support_status", ""),
            "provided_reads": int(to_float(enriched.get("provided_rna_junction_reads"), 0.0)),
            "selected_reads": int(to_float(enriched.get("rna_junction_reads"), 0.0)),
            "conflict": enriched.get("junction_support_conflict", "NONE"),
        },
    )()
    provenance = peptide_provenance_row(enriched, match_like)
    conflict: dict[str, str] | None = None
    if match_like.conflict and match_like.conflict != "NONE":
        conflict = {
            "evidence_domain": "splice_peptide",
            "record_id": peptide_id,
            "conflict_type": "JUNCTION_READ_COUNT_CONFLICT",
            "details": str(match_like.conflict),
            "source_tool": item.record.source_tool,
            "source_file": item.record.source_file,
            "source_row_number": str(item.record.source_row_number),
        }
    return enriched, provenance, conflict


def _tool_evidence_row(item: NormalizedRecord) -> dict[str, str]:
    base = item.record.as_dict()
    junction = item.resolution.junction
    if junction is not None:
        base.update(
            {
                "junction_id": junction.junction_id,
                "genome_build": junction.genome_build,
                "chrom": junction.chrom,
                "intron_start_1based": str(junction.intron_start_1based),
                "intron_end_1based": str(junction.intron_end_1based),
                "strand": junction.strand,
                "donor_1based": str(junction.donor_1based),
                "acceptor_1based": str(junction.acceptor_1based),
            }
        )
    metadata = peptide_metadata(item.record)
    material = f"{item.record.record_sha256}|{item.event_id}|{item.role}"
    base.update(
        {
            "evidence_id": f"EVID|SJ|{hashlib.sha256(material.encode('utf-8')).hexdigest()[:20]}",
            "entity_type": "splice_junction",
            "evidence_role": item.role,
            "resolution_status": item.resolution.status,
            "resolution_method": item.resolution.method,
            "coordinate_warning": item.resolution.warning or item.record.coordinate_warning,
            "variant_info": first(item.record.row, ["variant_info", "Variant Info", "variant_locus"], ""),
            "variant_key": "",
            # This is the source-row count, not necessarily verified support.
            "rna_junction_reads": str(item.record.total_split_reads),
            "peptide": metadata["peptide"],
            "hla_allele": metadata["hla_allele"],
            "raw_event_id": metadata["raw_event_id"],
        }
    )
    return {field: str(base.get(field) or "") for field in SPLICE_TOOL_EVIDENCE_FIELDS}


def _merge_conflict_to_public(row: dict[str, str], domain: str) -> dict[str, str]:
    return {
        "evidence_domain": domain,
        "record_id": row.get("merge_key", ""),
        "conflict_type": row.get("conflict_type", "FIELD_CONFLICT"),
        "details": f"{row.get('field', '')}: selected={row.get('selected_value', '')}; observed={row.get('observed_values', '')}",
        "source_tool": row.get("source_tools", ""),
        "source_file": "",
        "source_row_number": "",
    }


def normalize_splice_sources(
    *,
    sample_id: str,
    junctions: str | Path,
    outdir: str | Path,
    profile_name: str = "default",
    snaf: str | Path | None = None,
    splicemutr: str | Path | None = None,
    normal_junctions: str | Path | None = None,
    genome_build: str = "GRCh38",
    junction_coordinate_system: str = "auto",
    snaf_coordinate_system: str = "auto",
    splicemutr_coordinate_system: str = "auto",
    normal_coordinate_system: str = "auto",
    strict: bool = False,
) -> dict[str, str]:
    """Normalize splice sources and emit canonical entities plus full provenance."""

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    primary_path = Path(junctions)
    if not primary_path.is_file():
        raise FileNotFoundError(f"Missing primary splice junction table: {primary_path}")

    sources = [
        SpliceSource("RegTools", primary_path, "rna_junction", junction_coordinate_system),
    ]
    if snaf and Path(snaf).is_file():
        sources.append(SpliceSource("SNAF", Path(snaf), "neoantigen", snaf_coordinate_system))
    if splicemutr and Path(splicemutr).is_file():
        sources.append(
            SpliceSource("SpliceMutr", Path(splicemutr), "neoantigen", splicemutr_coordinate_system)
        )

    primary_tools = set(PRIMARY_JUNCTION_TOOLS)
    registry = JunctionRegistry()
    normalized: list[NormalizedRecord] = []

    # Primary RNA evidence is registered first.  Secondary callers can then be
    # linked only by exact coordinate/strand, exact unique source alias, or an
    # explicit unique variant-to-junction relation.
    primary = sources[0]
    for record in iter_junction_records(
        primary.path,
        sample_id=sample_id,
        source_tool=primary.tool,
        genome_build=genome_build,
        coordinate_system=primary.coordinate_system,
        source_tool_version=primary.version,
        strict=strict,
    ):
        resolution = registry.add(record)
        normalized.append(
            NormalizedRecord(
                record,
                resolution,
                primary.role,
                resolution.junction_id or unresolved_event_id(record),
            )
        )

    for source in sources[1:]:
        for record in iter_junction_records(
            source.path,
            sample_id=sample_id,
            source_tool=source.tool,
            genome_build=genome_build,
            coordinate_system=source.coordinate_system,
            source_tool_version=source.version,
            strict=False,
        ):
            resolution = registry.resolve(record)
            if resolution.junction is not None:
                registry.add(record, junction=resolution.junction)
            normalized.append(
                NormalizedRecord(
                    record,
                    resolution,
                    source.role,
                    resolution.junction_id or unresolved_event_id(record),
                )
            )

    normal_declared = bool(normal_junctions and Path(normal_junctions).is_file())
    normal_detected: set[str] = set()
    normal_resolvable_rows = 0
    if normal_declared:
        for record in iter_junction_records(
            Path(normal_junctions),
            sample_id=sample_id,
            source_tool="NormalJunctionPanel",
            genome_build=genome_build,
            coordinate_system=normal_coordinate_system,
            strict=False,
        ):
            if record.junction is not None:
                normal_resolvable_rows += 1
            resolution = registry.resolve(record)
            if resolution.junction is not None and resolution.junction.junction_id in registry.entities:
                normal_detected.add(resolution.junction.junction_id)
            normalized.append(
                NormalizedRecord(
                    record,
                    resolution,
                    "normal_background",
                    resolution.junction_id or unresolved_event_id(record),
                )
            )

    tumor_items = [item for item in normalized if item.role != "normal_background"]
    event_source_rows: list[dict[str, str]] = []
    for item in tumor_items:
        status = _normal_status(
            item.event_id,
            normal_file_declared=normal_declared,
            normal_resolvable_rows=normal_resolvable_rows,
            normal_detected=normal_detected,
        )
        event_source_rows.append(
            _event_source_row(
                item,
                sample_id=sample_id,
                profile_name=profile_name,
                registry=registry,
                primary_tools=primary_tools,
                normal_status=status,
            )
        )

    events, event_merge_provenance, event_merge_conflicts = merge_rows_preserving_provenance(
        event_source_rows,
        EVENT_FIELDS,
        ("event_id",),
        entity_type="splice_event",
    )
    event_by_id = {row["event_id"]: row for row in events}

    peptide_source_rows: list[dict[str, str]] = []
    peptide_support_provenance: list[dict[str, str]] = []
    public_conflicts: list[dict[str, str]] = list(registry.conflicts)
    for item in tumor_items:
        peptide, provenance, conflict = _peptide_source_row(
            item,
            sample_id=sample_id,
            registry=registry,
            primary_tools=primary_tools,
            event=event_by_id[item.event_id],
        )
        if peptide is not None:
            peptide_source_rows.append(peptide)
        if provenance is not None:
            peptide_support_provenance.append(provenance)
        if conflict is not None:
            public_conflicts.append(conflict)

    peptides, peptide_merge_provenance, peptide_merge_conflicts = merge_rows_preserving_provenance(
        peptide_source_rows,
        PEPTIDE_FIELDS,
        ("event_id", "peptide", "hla_allele"),
        entity_type="splice_peptide",
    )

    public_conflicts.extend(
        _merge_conflict_to_public(row, "splice_event_merge") for row in event_merge_conflicts
    )
    public_conflicts.extend(
        _merge_conflict_to_public(row, "splice_peptide_merge") for row in peptide_merge_conflicts
    )

    tool_evidence = [_tool_evidence_row(item) for item in normalized]
    entity_rows = registry.entity_rows(sample_id=sample_id, primary_tools=primary_tools)

    consensus_rows: list[dict[str, str]] = []
    for event_id in sorted(registry.entities):
        entity = registry.entities[event_id]
        tools = sorted(entity.source_tools)
        primary_present = any(tool.casefold() in {x.casefold() for x in primary_tools} for tool in tools)
        peptide_tools = sorted(
            {
                item.record.source_tool
                for item in tumor_items
                if item.event_id == event_id and peptide_metadata(item.record)["peptide"]
            }
        )
        if primary_present and peptide_tools:
            status = "CROSS_DOMAIN_CONFIRMED_EXACT_JUNCTION"
        elif primary_present and len(tools) >= 2:
            status = "MULTI_TOOL_EXACT_JUNCTION"
        elif primary_present:
            status = "SINGLE_TOOL_EXACT_JUNCTION"
        else:
            status = "RESOLVED_SOURCE_ONLY"
        consensus_rows.append(
            {
                "event_id": event_id,
                "canonical_junction_id": event_id,
                "support_tools": ";".join(tools),
                "peptide_tools": ";".join(peptide_tools),
                "n_tools": str(len(tools)),
                "rna_junction_reads": str(entity.maximum_reads(primary_tools)),
                "normal_junction_status": _normal_status(
                    event_id,
                    normal_file_declared=normal_declared,
                    normal_resolvable_rows=normal_resolvable_rows,
                    normal_detected=normal_detected,
                ),
                "coordinate_resolution": "RESOLVED",
                "status": status,
            }
        )
    for item in tumor_items:
        if item.resolution.junction is not None:
            continue
        consensus_rows.append(
            {
                "event_id": item.event_id,
                "canonical_junction_id": "",
                "support_tools": item.record.source_tool,
                "peptide_tools": item.record.source_tool if peptide_metadata(item.record)["peptide"] else "",
                "n_tools": "1",
                "rna_junction_reads": "0",
                "normal_junction_status": "UNASSESSED",
                "coordinate_resolution": item.resolution.status,
                "status": "UNRESOLVED_SOURCE_RECORD",
            }
        )

    alias_rows = [
        {
            "source_junction_id": alias,
            "canonical_junction_ids": ";".join(sorted(event_ids)),
            "n_junctions": str(len(event_ids)),
            "alias_status": "UNIQUE" if len(event_ids) == 1 else "AMBIGUOUS",
        }
        for alias, event_ids in sorted(registry.alias_to_ids.items())
    ]

    domain_by_role = {
        "rna_junction": "splice_rna",
        "neoantigen": "splice_neoantigen",
        "normal_background": "normal_background",
    }
    consensus_provenance_rows = [
        {
            "event_id": item.event_id,
            "canonical_junction_id": item.resolution.junction_id,
            "evidence_domain": domain_by_role.get(item.role, item.role),
            "tool": item.record.source_tool,
            "source_file": item.record.source_file,
            "source_row_number": str(item.record.source_row_number),
            "source_record_id": item.record.source_record_id,
            "source_junction_id": item.record.source_junction_id,
            "resolution_status": item.resolution.status,
            "resolution_method": item.resolution.method,
            "coordinate_warning": item.resolution.warning or item.record.coordinate_warning,
            "peptide_present": "yes" if peptide_metadata(item.record)["peptide"] else "no",
        }
        for item in normalized
    ]

    rna_evidence_rows: list[dict[str, str]] = []
    for event in events:
        event_id = event.get("event_id", "")
        rna_evidence_rows.append(
            {
                "evidence_id": f"RNAJ|{hashlib.sha256(event_id.encode('utf-8')).hexdigest()[:20]}",
                "event_id": event_id,
                "peptide_id": "",
                "sample_id": sample_id,
                "gene": event.get("gene", ""),
                "gene_pair": "",
                "junction_reads": event.get("rna_junction_reads", "0"),
                "junction_source": event.get("rna_junction_source", ""),
                "mutation_source": event.get("mutation_source", "Other"),
                "peptide_consequence": "splice_junction",
                "rna_frame_status": event.get("rna_frame_status", ""),
                "rna_support_status": event.get("rna_support_status", ""),
                "rna_evidence_completeness": event.get("rna_evidence_completeness", ""),
                "rna_evidence_score": event.get("rna_evidence_score", ""),
                "targeted_validation_status": "UNASSESSED",
                "targeted_validation_source": "",
                "targeted_validation_method": "",
            }
        )

    unresolved_tumor = [item for item in tumor_items if item.resolution.junction is None]
    exact_cross_domain = sum(
        row.get("status") == "CROSS_DOMAIN_CONFIRMED_EXACT_JUNCTION" for row in consensus_rows
    )
    qc_rows = [
        {"metric": "primary_junction_records", "value": str(sum(item.role == "rna_junction" for item in normalized))},
        {"metric": "resolved_junction_entities", "value": str(len(registry.entities))},
        {"metric": "unresolved_tumor_records", "value": str(len(unresolved_tumor))},
        {"metric": "peptide_rows_before_merge", "value": str(len(peptide_source_rows))},
        {"metric": "peptide_rows_after_merge", "value": str(len(peptides))},
        {"metric": "event_provenance_records", "value": str(len(event_merge_provenance))},
        {"metric": "peptide_provenance_records", "value": str(len(peptide_merge_provenance))},
        {"metric": "junction_support_provenance_records", "value": str(len(peptide_support_provenance))},
        {"metric": "conflicts", "value": str(len(public_conflicts))},
        {"metric": "exact_cross_domain_confirmed", "value": str(exact_cross_domain)},
        {"metric": "normal_junction_file_declared", "value": str(normal_declared).lower()},
        {"metric": "normal_resolvable_rows", "value": str(normal_resolvable_rows)},
        {"metric": "normal_exact_tumor_junction_hits", "value": str(len(normal_detected))},
        {"metric": "gene_or_nearest_locus_fallbacks_used", "value": "0"},
    ]

    paths = {
        "raw_events": out / "raw_events.tsv",
        "raw_peptides": out / "raw_peptides.tsv",
        "splice_junctions": out / "splice_junctions.tsv",
        "splice_tool_evidence": out / "splice_tool_evidence.long.tsv",
        "splice_peptide_provenance": out / "splice_peptide_provenance.tsv",
        "event_merge_provenance": out / "splice_event_merge_provenance.tsv",
        "peptide_merge_provenance": out / "splice_peptide_merge_provenance.tsv",
        "splice_merge_conflicts": out / "splice_merge_conflicts.tsv",
        "rna_junction_evidence": out / "rna_junction_evidence.tsv",
        "splice_consensus": out / "splice_consensus.tsv",
        "splice_consensus_provenance": out / "splice_consensus_provenance.tsv",
        "splice_consensus_conflicts": out / "splice_consensus_conflicts.tsv",
        "junction_aliases": out / "junction_aliases.tsv",
        "evidence_conflicts": out / "evidence_conflicts.tsv",
        "splice_qc": out / "splice_qc.tsv",
        "provenance_manifest": out / "provenance_manifest.json",
    }

    write_tsv(paths["raw_events"], events, EVENT_FIELDS)
    write_tsv(paths["raw_peptides"], peptides, PEPTIDE_FIELDS)
    write_tsv(paths["splice_junctions"], entity_rows, JUNCTION_ENTITY_FIELDS)
    write_tsv(paths["splice_tool_evidence"], tool_evidence, SPLICE_TOOL_EVIDENCE_FIELDS)
    write_tsv(
        paths["splice_peptide_provenance"],
        peptide_support_provenance,
        SPLICE_PEPTIDE_PROVENANCE_FIELDS,
    )
    write_tsv(paths["event_merge_provenance"], event_merge_provenance, PROVENANCE_FIELDS)
    write_tsv(paths["peptide_merge_provenance"], peptide_merge_provenance, PROVENANCE_FIELDS)
    write_tsv(
        paths["splice_merge_conflicts"],
        event_merge_conflicts + peptide_merge_conflicts,
        CONFLICT_FIELDS,
    )
    write_tsv(paths["rna_junction_evidence"], rna_evidence_rows, RNA_JUNCTION_EVIDENCE_FIELDS)
    write_tsv(
        paths["splice_consensus"],
        consensus_rows,
        [
            "event_id",
            "canonical_junction_id",
            "support_tools",
            "peptide_tools",
            "n_tools",
            "rna_junction_reads",
            "normal_junction_status",
            "coordinate_resolution",
            "status",
        ],
    )
    write_tsv(
        paths["splice_consensus_provenance"],
        consensus_provenance_rows,
        CONSENSUS_PROVENANCE_FIELDS,
    )
    write_tsv(
        paths["splice_consensus_conflicts"],
        public_conflicts,
        SPLICE_CONFLICT_FIELDS,
    )
    write_tsv(
        paths["junction_aliases"],
        alias_rows,
        ["source_junction_id", "canonical_junction_ids", "n_junctions", "alias_status"],
    )
    write_tsv(paths["evidence_conflicts"], public_conflicts, SPLICE_CONFLICT_FIELDS)
    write_tsv(paths["splice_qc"], qc_rows, ["metric", "value"])

    input_paths = [source.path for source in sources]
    if normal_declared:
        input_paths.append(Path(normal_junctions))
    write_json(
        paths["provenance_manifest"],
        {
            "schema_version": CANONICAL_JUNCTION_SCHEMA_VERSION,
            "software_version": "0.4.4",
            "sample_id": sample_id,
            "profile_name": profile_name,
            "genome_build": genome_build,
            "canonical_coordinate_system": "1-based closed intron interval",
            "matching_policy": [
                "exact_canonical_junction_id",
                "exact_build_chrom_intron_start_intron_end_strand",
                "exact_unique_source_junction_id",
                "unique_unstranded_coordinate_with_explicit_caution",
                "unique_explicit_variant_to_junction_link",
                "no_gene_level_fallback",
                "no_nearest_locus_fallback",
            ],
            "evidence_policy": {
                "verified_rna_junction_reads": "primary exact canonical junction only",
                "caller_provided_unverified_reads": "provided_rna_junction_reads only",
                "normal_panel_absence": "NOT_DETECTED_COVERAGE_UNASSESSED unless per-locus coverage is supplied",
            },
            "inputs": [
                {"path": str(path), "sha256": file_sha256(path)}
                for path in input_paths
                if path.is_file()
            ],
            "outputs": {key: str(path) for key, path in paths.items()},
        },
    )
    return {key: str(path) for key, path in paths.items()}
