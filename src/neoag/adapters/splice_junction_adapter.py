"""Mode C splice-junction adapter for NeoAg v0.4.4.

The adapter transfers RNA junction reads only through an exact canonical
junction, an exact unique caller junction ID, or an explicit unique
variant-to-junction relation. Gene-level and nearest-locus fallbacks are
forbidden because they can assign a highly expressed canonical junction to an
unrelated weak/novel junction in the same gene.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..model_layers import enrich_peptide_layers
from ..provenance import (
    CONFLICT_FIELDS,
    PROVENANCE_FIELDS,
    merge_rows_preserving_provenance,
)
from ..schemas import EVENT_FIELDS, PEPTIDE_FIELDS
from ..splice.registry import (
    JUNCTION_ENTITY_FIELDS,
    SPLICE_PEPTIDE_PROVENANCE_FIELDS,
    JunctionSupportIndex,
    build_support_index,
    peptide_provenance_row,
    resolve_junction_support,
)
from ..utils import read_tsv, to_float, write_tsv
from .event_catalog import parse_splice_catalog

SPLICE_CONSEQUENCES = frozenset(
    {
        "splice_donor_variant",
        "splice_acceptor_variant",
        "splice_donor_region_variant",
        "splice_acceptor_region_variant",
        "splice_region_variant",
        "splice_polypyrimidine_tract_variant",
    }
)

_VARIANT_INFO_RE = re.compile(r"^(?P<chrom>[^:]+):(?P<start>\d+)(?:-(?P<end>\d+))?$")


@dataclass
class SpliceJunctionFilterConfig:
    min_junction_reads: int = 3
    require_verified_junction_support: bool = True
    genome_build: str = "GRCh38"
    coordinate_system: str = "auto"


def _consequences(cons_field: str) -> set[str]:
    return {value.strip() for value in str(cons_field or "").split("&") if value.strip()}


def is_splice_consequence(cons_field: str) -> bool:
    return bool(_consequences(cons_field) & SPLICE_CONSEQUENCES)


def parse_regtools_variant_info(value: str) -> tuple[str, str, str] | None:
    match = _VARIANT_INFO_RE.match(str(value or "").strip())
    if not match:
        return None
    chrom = match.group("chrom")
    start = match.group("start")
    end = match.group("end") or start
    return chrom, start, end


def build_junction_support_index(
    splice_path: str | Path,
    *,
    sample_id: str = "",
    genome_build: str = "GRCh38",
    coordinate_system: str = "auto",
) -> JunctionSupportIndex:
    """Build an exact RegTools junction index without gene/locus fallbacks."""

    return build_support_index(
        splice_path,
        sample_id=sample_id,
        genome_build=genome_build,
        coordinate_system=coordinate_system,
        source_tool="RegTools",
    )


def junction_reads_for_peptide(
    peptide: dict[str, str],
    index: JunctionSupportIndex,
    *,
    event: dict[str, str] | None = None,
) -> int:
    """Return only verified reads; unresolved upstream counts return zero."""

    return resolve_junction_support(peptide, index, event=event).selected_reads


def _apply_match_fields(
    row: dict[str, str],
    index: JunctionSupportIndex,
    *,
    event: dict[str, str] | None,
) -> tuple[dict[str, str], Any]:
    match = resolve_junction_support(row, index, event=event)
    out = dict(row)
    out.update(
        {
            "provided_rna_junction_reads": str(match.provided_reads),
            "rna_junction_reads": str(match.selected_reads),
            "junction_match_status": match.match_status,
            "junction_match_method": match.match_method,
            "junction_support_status": match.support_status,
            "junction_support_conflict": match.conflict,
            "junction_support_reason": match.reason,
        }
    )
    if match.source_junction_id and not out.get("source_junction_id"):
        out["source_junction_id"] = match.source_junction_id

    if match.junction_id:
        entity = index.registry.entities[match.junction_id]
        junction = entity.junction
        out.update(
            {
                "genome_build": junction.genome_build,
                "canonical_junction_id": junction.junction_id,
                "junction_chrom": junction.chrom,
                "junction_start": str(junction.intron_start_1based),
                "junction_end": str(junction.intron_end_1based),
                "junction_strand": junction.strand,
                "junction_donor": str(junction.donor_1based),
                "junction_acceptor": str(junction.acceptor_1based),
                "junction_coordinate_system": "intron_1based_closed",
                "junction_resolution_status": "RESOLVED",
                "junction_resolution_reason": match.match_method,
                "rna_junction_source": ";".join(
                    sorted(
                        {
                            record.source_tool
                            for record in entity.records
                            if record.source_tool.casefold()
                            in {tool.casefold() for tool in index.primary_tools}
                        }
                    )
                ),
            }
        )
    else:
        out.setdefault("junction_resolution_status", match.match_status)
        out.setdefault("junction_resolution_reason", match.reason)
    return out, match


def enrich_splice_peptide_layers(
    peptide: dict[str, str],
    *,
    index: JunctionSupportIndex | None = None,
    event: dict[str, str] | None = None,
) -> dict[str, str]:
    out = dict(peptide)
    out["peptide_consequence"] = "splice_junction"
    if not out.get("crosses_junction"):
        out["crosses_junction"] = "yes"
    if not out.get("source_tool"):
        out["source_tool"] = "splice-junction-adapter"
    out.setdefault("source_tools", out["source_tool"])
    out.setdefault("source_record_id", out.get("peptide_id", ""))
    out.setdefault("source_records", out.get("source_record_id", ""))
    out.setdefault("provenance_record_count", "1")
    out.setdefault("evidence_conflict_status", "NONE")

    if index is not None:
        out, _ = _apply_match_fields(out, index, event=event)
    return enrich_peptide_layers(out, event)


def filter_splice_variant_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in rows if is_splice_consequence(row.get("consequence", ""))]


def _canonicalize_splice_events(
    events: list[dict[str, str]],
    index: JunctionSupportIndex,
) -> tuple[list[dict[str, str]], dict[str, str]]:
    canonicalized: list[dict[str, str]] = []
    id_map: dict[str, str] = {}

    for raw in events:
        event = dict(raw)
        old_id = str(event.get("event_id") or "")
        is_splice = (
            str(event.get("event_type") or "").lower() == "splice"
            or str(event.get("peptide_consequence") or "").lower() == "splice_junction"
            or str(event.get("source") or "").lower() == "pvacsplice"
        )
        if not is_splice:
            canonicalized.append(event)
            continue

        event.setdefault("source_tools", event.get("source") or "pVACsplice")
        event.setdefault("source_record_id", old_id)
        event.setdefault("source_records", event.get("source_record_id") or old_id)
        event.setdefault("provenance_record_count", "1")
        event.setdefault("evidence_conflict_status", "NONE")
        event, match = _apply_match_fields(event, index, event=None)

        if match.junction_id:
            event["event_id"] = match.junction_id
            # Keep caller-level chrom/pos as variant provenance. Canonical
            # junction coordinates are in junction_* fields.
            if old_id:
                id_map[old_id] = match.junction_id
        canonicalized.append(event)

    return canonicalized, id_map


def _verified_for_filter(peptide: dict[str, str]) -> bool:
    status = str(peptide.get("junction_support_status") or "").upper()
    return status in {"SUPPORTED_EXACT_JUNCTION", "MATCHED_ZERO_READS"}


def build_splice_peptides_from_vcf(
    variants_vcf: str | Path,
    splice_path: str | Path,
    *,
    sample_id: str,
    profile_name: str,
    hla_alleles: list[str],
    cfg: dict[str, Any] | None = None,
    tools_dir: Path | None = None,
    min_junction_reads: int = 0,
) -> dict[str, Any]:
    """Extract splice-affecting variant peptides and attach exact RNA support."""

    from ..vep.extract_peptides import extract_variant_peptides_from_vcf
    from .variant_peptide_adapter import (
        resolve_variant_peptide_options,
        variant_peptide_rows_to_raw_tables,
    )

    if not hla_alleles:
        raise ValueError("build_splice_peptides_from_vcf requires inputs.hla_alleles")

    cfg = cfg or {}
    inputs = cfg.get("inputs") or {}
    opts = resolve_variant_peptide_options(cfg)
    out_tsv = (tools_dir or Path(".")) / "splice_variant_peptides.tsv"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    summary = extract_variant_peptides_from_vcf(
        variants_vcf,
        out_tsv,
        lengths=opts["lengths"],
        sample_id=sample_id,
        exclude_multi_aa=opts["exclude_multi_aa"],
        single_aa_only=opts["single_aa_only"],
        mini_len=opts["mini_len"],
        normal_proteome_fasta=opts["normal_proteome_fasta"],
        filter_normal_proteome=opts["filter_normal_proteome"],
        hla_alleles=hla_alleles,
        tumor_sample_name=str(inputs.get("tumor_sample_name") or sample_id) or None,
        rna_sample_name=str(inputs.get("rna_sample_name") or "") or None,
        consequence_filter="splice",
    )
    splice_rows = filter_splice_variant_rows(read_tsv(out_tsv))
    events, peptides = variant_peptide_rows_to_raw_tables(
        splice_rows,
        sample_id=sample_id,
        profile_name=profile_name,
        hla_alleles=hla_alleles,
    )
    index = build_junction_support_index(
        splice_path,
        sample_id=sample_id,
        genome_build=str(inputs.get("genome_build") or "GRCh38"),
        coordinate_system=str(inputs.get("splice_coordinate_system") or "auto"),
    )
    events, id_map = _canonicalize_splice_events(events, index)
    events_by_id = {event["event_id"]: event for event in events}
    enriched: list[dict[str, str]] = []
    provenance: list[dict[str, str]] = []

    for raw in peptides:
        peptide = dict(raw)
        peptide["event_id"] = id_map.get(
            peptide.get("event_id", ""), peptide.get("event_id", "")
        )
        event = events_by_id.get(peptide["event_id"])
        peptide = enrich_splice_peptide_layers(peptide, index=index, event=event)
        match = resolve_junction_support(peptide, index, event=event)
        provenance.append(peptide_provenance_row(peptide, match))
        if min_junction_reads:
            if not _verified_for_filter(peptide):
                continue
            if int(to_float(peptide.get("rna_junction_reads"), 0.0)) < min_junction_reads:
                continue
        enriched.append(peptide)

    return {
        "events": events,
        "peptides": enriched,
        "peptide_provenance": provenance,
        "junction_index": index,
        "splice_variant_peptides": str(out_tsv),
        "summary": summary,
        "splice_variant_rows": len(splice_rows),
    }


def merge_splice_into_catalog(
    splice_path: str | Path,
    sample_id: str,
    profile_name: str,
    events: list[dict[str, str]],
    peptides: list[dict[str, str]],
    *,
    variants_vcf: str | Path | None = None,
    hla_alleles: list[str] | None = None,
    cfg: dict[str, Any] | None = None,
    tools_dir: Path | None = None,
    filter_cfg: SpliceJunctionFilterConfig | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Merge RegTools events and peptides while retaining all provenance."""

    filter_cfg = filter_cfg or SpliceJunctionFilterConfig()
    cfg = cfg or {}
    inputs = cfg.get("inputs") or {}
    build = str(inputs.get("genome_build") or filter_cfg.genome_build)
    coordinate_system = str(
        inputs.get("splice_coordinate_system") or filter_cfg.coordinate_system
    )
    index = build_junction_support_index(
        splice_path,
        sample_id=sample_id,
        genome_build=build,
        coordinate_system=coordinate_system,
    )

    canonical_existing, id_map = _canonicalize_splice_events(events, index)
    splice_events = parse_splice_catalog(
        splice_path,
        sample_id,
        profile_name,
        genome_build=build,
        coordinate_system=coordinate_system,
    )
    merged_events, event_provenance, event_conflicts = merge_rows_preserving_provenance(
        canonical_existing + splice_events,
        EVENT_FIELDS,
        ("event_id",),
        entity_type="event",
    )
    event_by_id = {event["event_id"]: event for event in merged_events}

    peptide_rows: list[dict[str, str]] = []
    peptide_support_provenance: list[dict[str, str]] = []
    for raw in peptides:
        peptide = dict(raw)
        peptide["event_id"] = id_map.get(
            peptide.get("event_id", ""), peptide.get("event_id", "")
        )
        event = event_by_id.get(peptide.get("event_id", ""))
        peptide = enrich_splice_peptide_layers(peptide, index=index, event=event)
        match = resolve_junction_support(peptide, index, event=event)
        peptide_support_provenance.append(peptide_provenance_row(peptide, match))
        peptide_rows.append(peptide)

    if not peptide_rows and variants_vcf and hla_alleles:
        built = build_splice_peptides_from_vcf(
            variants_vcf,
            splice_path,
            sample_id=sample_id,
            profile_name=profile_name,
            hla_alleles=list(hla_alleles),
            cfg=cfg,
            tools_dir=tools_dir,
            min_junction_reads=filter_cfg.min_junction_reads,
        )
        merged_events, extra_event_provenance, extra_event_conflicts = (
            merge_rows_preserving_provenance(
                merged_events + built["events"],
                EVENT_FIELDS,
                ("event_id",),
                entity_type="event",
            )
        )
        event_provenance.extend(extra_event_provenance)
        event_conflicts.extend(extra_event_conflicts)
        peptide_rows = built["peptides"]
        peptide_support_provenance.extend(built["peptide_provenance"])

    if filter_cfg.min_junction_reads:
        filtered: list[dict[str, str]] = []
        for peptide in peptide_rows:
            if filter_cfg.require_verified_junction_support and not _verified_for_filter(peptide):
                continue
            if int(to_float(peptide.get("rna_junction_reads"), 0.0)) >= filter_cfg.min_junction_reads:
                filtered.append(peptide)
        peptide_rows = filtered

    merged_peptides, peptide_merge_provenance, peptide_conflicts = (
        merge_rows_preserving_provenance(
            peptide_rows,
            PEPTIDE_FIELDS,
            ("event_id", "peptide", "hla_allele"),
            entity_type="peptide",
        )
    )

    if tools_dir:
        tools_dir.mkdir(parents=True, exist_ok=True)
        write_tsv(
            tools_dir / "splice_junctions.tsv",
            index.registry.entity_rows(sample_id=sample_id, primary_tools={"RegTools"}),
            JUNCTION_ENTITY_FIELDS,
        )
        write_tsv(
            tools_dir / "splice_peptide_provenance.tsv",
            peptide_support_provenance,
            SPLICE_PEPTIDE_PROVENANCE_FIELDS,
        )
        write_tsv(
            tools_dir / "splice_event_merge_provenance.tsv",
            event_provenance,
            PROVENANCE_FIELDS,
        )
        write_tsv(
            tools_dir / "splice_peptide_merge_provenance.tsv",
            peptide_merge_provenance,
            PROVENANCE_FIELDS,
        )
        write_tsv(
            tools_dir / "splice_merge_conflicts.tsv",
            event_conflicts + peptide_conflicts,
            CONFLICT_FIELDS,
        )

    return merged_events, merged_peptides


def run_splice_junction_upstream(
    cfg: dict[str, Any],
    *,
    splice_path: Path,
    variants_vcf: Path | None,
    parsed_dir: Path,
    tools_dir: Path,
    sample_id: str,
    profile_name: str,
    hla_alleles: list[str],
    pvacsplice_tsv: Path | None = None,
) -> dict[str, str]:
    """Build parsed/raw tables for splice_junction entry mode."""

    from .pvactools_parser import parse_pvactools_outputs

    events: list[dict[str, str]] = []
    peptides: list[dict[str, str]] = []
    outputs: dict[str, str] = {"splice_junction_tsv": str(splice_path)}

    if pvacsplice_tsv and pvacsplice_tsv.is_file():
        events, peptides = parse_pvactools_outputs(
            [pvacsplice_tsv], sample_id, profile_name
        )
        outputs["pvacsplice"] = str(pvacsplice_tsv)
        outputs["peptide_source"] = "pvacsplice"

    events, peptides = merge_splice_into_catalog(
        splice_path,
        sample_id,
        profile_name,
        events,
        peptides,
        variants_vcf=variants_vcf,
        hla_alleles=hla_alleles,
        cfg=cfg,
        tools_dir=tools_dir,
    )

    if not peptides:
        raise ValueError(
            "splice_junction mode produced no peptides. Provide pVACsplice output or "
            "variants_vcf + hla_alleles with splice-affecting PASS variants."
        )

    parsed_dir.mkdir(parents=True, exist_ok=True)
    raw_events = parsed_dir / "raw_events.tsv"
    raw_peptides = parsed_dir / "raw_peptides.tsv"
    write_tsv(raw_events, events, EVENT_FIELDS)
    write_tsv(raw_peptides, peptides, PEPTIDE_FIELDS)
    outputs.update(
        {
            "raw_events": str(raw_events),
            "raw_peptides": str(raw_peptides),
            "splice_junctions": str(tools_dir / "splice_junctions.tsv"),
            "splice_peptide_provenance": str(
                tools_dir / "splice_peptide_provenance.tsv"
            ),
            "splice_event_merge_provenance": str(
                tools_dir / "splice_event_merge_provenance.tsv"
            ),
            "splice_peptide_merge_provenance": str(
                tools_dir / "splice_peptide_merge_provenance.tsv"
            ),
            "splice_merge_conflicts": str(tools_dir / "splice_merge_conflicts.tsv"),
        }
    )
    if outputs.get("peptide_source") is None:
        outputs["peptide_source"] = "splice-variant-peptides"
    return outputs
