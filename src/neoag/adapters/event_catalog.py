"""Catalog adapters: fusion / splice event stubs for multi-entry ingestion (Layer 1 events)."""

from __future__ import annotations

from pathlib import Path

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


def _splice_gene_name(row: dict[str, str]) -> str:
    gene = first(row, ["gene", "Gene", "gene_name", "gene_names", "symbol"], "")
    if gene.upper() in {"", "NA", "N/A", "."}:
        gene = ""
    chrom = first(row, ["chrom1", "chrom", "Chromosome"], "")
    pos = first(row, ["start1", "start", "Start"], "")
    if not gene and chrom:
        gene = f"{chrom}:{pos}"
    return gene


def parse_splice_catalog(
    splice_path: str | Path,
    sample_id: str,
    profile_name: str,
) -> list[dict[str, str]]:
    """Ingest RegTools / junction annotation tables as splice/junction Layer-1 events."""
    events: dict[str, dict[str, str]] = {}
    for row in read_tsv(splice_path):
        gene = _splice_gene_name(row)
        if not gene:
            continue
        chrom = first(row, ["chrom1", "chrom", "Chromosome"], "")
        pos = first(row, ["start1", "start", "Start"], "")
        reads = int(
            to_float(
                first(
                    row,
                    [
                        "counts",
                        "junction_reads",
                        "reads",
                        "split_reads",
                        "score",  # RegTools junction read support
                        "JunctionReadCount",
                    ],
                    "0",
                ),
                0.0,
            )
        )
        ev = _base_event(
            sample_id=sample_id,
            profile_name=profile_name,
            event_type="Splice",
            gene=gene,
            event_name=first(row, ["event_id", "name", "junction_id"], gene),
            source=f"splice_catalog:{Path(splice_path).name}",
            mutation_source=infer_mutation_source(event_type="Splice", tool="pVACsplice"),
            peptide_consequence="splice_junction",
            rna_reads=reads,
        )
        if chrom:
            ev["chrom"] = chrom
        if pos:
            ev["pos"] = pos
        transcript = first(row, ["transcript_id", "transcripts", "Transcript"], "")
        if transcript and transcript.upper() not in {"NA", "N/A", "."}:
            ev["transcript_id"] = transcript.split(",")[0]
        events[ev["event_id"]] = ev
    return list(events.values())


def write_event_catalog(events: list[dict[str, str]], out_path: str | Path) -> None:
    write_tsv(out_path, events, EVENT_FIELDS)
