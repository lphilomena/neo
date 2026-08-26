"""Catalog adapters: fusion / splice event stubs for multi-entry ingestion (Layer 1 events)."""

from __future__ import annotations

from pathlib import Path
import re

from ..model_layers import enrich_event_layers, infer_mutation_source, infer_peptide_consequence
from ..schemas import EVENT_FIELDS
from ..utils import first, read_tsv, safe_id, to_float, write_tsv


def _base_event(
    *,
    sample_id: str,
    profile_name: str,
    event_type: str,
    gene: str,
    event_name: str,
    source: str,
    mutation_source: str,
    peptide_consequence: str,
    rna_reads: int = 0,
    expression: float = 0.0,
) -> dict[str, str]:
    eid = safe_id(f"{sample_id}_{event_type}_{gene}_{event_name}")
    base = {
        "event_id": eid,
        "sample_id": sample_id,
        "disease_profile": profile_name,
        "event_type": event_type,
        "mutation_source": mutation_source,
        "peptide_consequence": peptide_consequence,
        "gene": gene,
        "event_name": event_name,
        "chrom": "",
        "pos": "",
        "ref": "",
        "alt": "",
        "transcript_id": "",
        "consequence": peptide_consequence.replace("_", " "),
        "rna_junction_reads": str(rna_reads),
        "event_confidence": "0.7",
        "event_expression": f"{expression:.4f}",
        "driver_relevance": "0.0",
        "tumor_vaf": "0.0",
        "tumor_depth": "",
        "tumor_alt_count": "",
        "rna_vaf": "",
        "rna_alt_reads": "",
        "rna_depth": "",
        "clonality": "0.5",
        "persistence": "0.5",
        "tumor_specificity": "0.7",
        "source": source,
    }
    return enrich_event_layers(base)


def parse_fusion_catalog(
    fusion_path: str | Path,
    sample_id: str,
    profile_name: str,
) -> list[dict[str, str]]:
    """Ingest STAR-Fusion / Arriba / AGFusion-style fusion tables as Layer-1 events."""
    events: dict[str, dict[str, str]] = {}
    for row in read_tsv(fusion_path):
        g1 = first(row, ["LeftGene", "left_gene", "gene1", "Gene1", "#FusionName"], "")
        g2 = first(row, ["RightGene", "right_gene", "gene2", "Gene2"], "")
        if "#FusionName" in row and not g2:
            parts = str(row["#FusionName"]).split("--")
            if len(parts) >= 2:
                g1, g2 = parts[0], parts[1]
        if not g1:
            continue
        gene = f"{g1}::{g2}" if g2 and g2 != g1 else g1
        reads = int(to_float(first(row, ["JunctionReadCount", "junction_reads", "split_reads", "reads"], "0"), 0.0))
        ev = _base_event(
            sample_id=sample_id,
            profile_name=profile_name,
            event_type="Fusion",
            gene=gene,
            event_name=gene,
            source=f"fusion_catalog:{Path(fusion_path).name}",
            mutation_source=infer_mutation_source(event_type="Fusion", tool="pVACfuse"),
            peptide_consequence="fusion",
            rna_reads=reads,
        )
        events[ev["event_id"]] = ev
    return list(events.values())


def _splice_gene_name_from_record(record) -> str:
    gene = str(record.gene or "").strip()
    if gene.upper() in {"", "NA", "N/A", "."}:
        gene = ""
    source_id = str(record.source_junction_id or "").strip()
    ensg_match = re.search(r"ENSG\d+", source_id)
    if (not gene or re.fullmatch(r"(?:chr)?[0-9XYM]+:\d+(?:[-:]\d+)?", gene, re.IGNORECASE)) and ensg_match:
        gene = ensg_match.group(0)
    if not gene and record.source_chrom:
        gene = f"{record.source_chrom}:{record.source_start}"
    return gene


def parse_splice_catalog(
    splice_path: str | Path,
    sample_id: str,
    profile_name: str,
    *,
    genome_build: str = "GRCh38",
    coordinate_system: str = "auto",
    source_tool: str = "RegTools",
) -> list[dict[str, str]]:
    """Ingest a splice-junction catalog with canonical, auditable coordinates.

    Resolved rows use ``SJ|build|chrom|intron_start|intron_end|strand`` as the
    event ID. Unresolved rows receive a source-scoped deterministic ID and are
    never merged by gene or nearby genomic position.
    """

    from ..provenance import merge_rows_preserving_provenance
    from ..splice.coordinates import iter_junction_records
    from ..splice.registry import unresolved_event_id

    rows: list[dict[str, str]] = []
    for record in iter_junction_records(
        splice_path,
        sample_id=sample_id,
        source_tool=source_tool,
        genome_build=genome_build,
        coordinate_system=coordinate_system,
        strict=False,
    ):
        junction = record.junction
        event_id = junction.junction_id if junction else unresolved_event_id(record)
        gene = _splice_gene_name_from_record(record)
        if not gene:
            gene = f"{source_tool}:row{record.source_row_number}"

        ev = _base_event(
            sample_id=sample_id,
            profile_name=profile_name,
            event_type="Splice",
            gene=gene,
            event_name=record.source_junction_id or event_id,
            source=f"splice_catalog:{Path(splice_path).name}",
            mutation_source=infer_mutation_source(event_type="Splice", tool="pVACsplice"),
            peptide_consequence="splice_junction",
            rna_reads=record.total_split_reads,
        )
        ev["event_id"] = event_id
        ev["chrom"] = record.source_chrom or (junction.chrom if junction else "")
        # Preserve the caller's coordinate for backward compatibility; the
        # exact normalized intron coordinate is stored separately below.
        ev["pos"] = record.source_start
        if record.transcript_ids:
            ev["transcript_id"] = record.transcript_ids.replace(",", ";").split(";")[0]

        ev.update(
            {
                "source_junction_id": record.source_junction_id,
                "junction_coordinate_system": record.source_coordinate_system,
                "junction_resolution_status": record.resolution_status,
                "junction_resolution_reason": record.resolution_method
                + (f": {record.coordinate_warning}" if record.coordinate_warning else ""),
                "junction_match_status": "PRIMARY_SOURCE_RECORD",
                "junction_match_method": "source_row_coordinates",
                "junction_support_status": (
                    "SUPPORTED_EXACT_JUNCTION"
                    if junction is not None and record.total_split_reads > 0
                    else "MATCHED_ZERO_READS"
                    if junction is not None
                    else "PROVIDED_SOURCE_RECORD_UNRESOLVED_COORDINATES"
                ),
                "junction_support_conflict": "NONE",
                "junction_support_reason": (
                    "Read count belongs to this exact normalized source junction."
                    if junction is not None
                    else "Read count belongs to this source row but cannot be transferred to another event until coordinates are resolved."
                ),
                "provided_rna_junction_reads": str(record.total_split_reads),
                "unique_junction_reads": str(record.unique_split_reads),
                "junction_total_coverage": str(record.total_split_reads),
                "splice_psi": first(record.row, ["splice_psi", "junction_psi", "tumor_psi", "psi", "PSI"], ""),
                "splice_alignment_qc_status": "PASS" if junction is not None else "FAIL_UNRESOLVED_COORDINATES",
                "matched_normal_junction_status": first(record.row, ["matched_normal_junction_status", "patient_normal_junction_status"], ""),
                "matched_normal_junction_reads": first(record.row, ["matched_normal_junction_reads", "normal_junction_reads"], ""),
                "normal_cohort_junction_status": first(record.row, ["normal_cohort_junction_status", "gtex_junction_status"], ""),
                "annotated_normal_isoform_status": first(record.row, ["annotated_normal_isoform_status", "normal_isoform_status"], ""),
                "splice_annotation_status": first(record.row, ["splice_annotation_status", "junction_annotation_status", "known_junction"], ""),
                "splice_orf_status": first(record.row, ["splice_orf_status", "orf_status", "frame_status"], ""),
                "splice_nmd_status": first(record.row, ["splice_nmd_status", "nmd_status", "predicted_NMD", "predicted_nmd"], ""),
                "rna_junction_source": source_tool,
                "source_file": record.source_file,
                "source_row_number": str(record.source_row_number),
                "source_record_id": record.source_record_id,
                "source_tools": source_tool,
                "source_records": record.source_record_id,
                "provenance_record_count": "1",
                "evidence_conflict_status": "NONE",
            }
        )
        if junction is not None:
            ev.update(
                {
                    "genome_build": junction.genome_build,
                    "canonical_junction_id": junction.junction_id,
                    "junction_chrom": junction.chrom,
                    "junction_start": str(junction.intron_start_1based),
                    "junction_end": str(junction.intron_end_1based),
                    "junction_strand": junction.strand,
                    "junction_donor": str(junction.donor_1based),
                    "junction_acceptor": str(junction.acceptor_1based),
                }
            )
        rows.append(enrich_event_layers(ev))

    merged, _, _ = merge_rows_preserving_provenance(
        rows,
        EVENT_FIELDS,
        ("event_id",),
        entity_type="splice_event",
    )
    return merged


def write_event_catalog(events: list[dict[str, str]], out_path: str | Path) -> None:
    write_tsv(out_path, events, EVENT_FIELDS)
