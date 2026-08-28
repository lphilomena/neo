"""Strict RNA-junction and expressed-product inputs for the SV pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..utils import first, read_tsv, to_float
from .identity import canonical_breakpoint_key
from .protein_reconstruct import ProteinReconstruction


def _required(row: dict[str, str], names: list[str], label: str) -> str:
    value = first(row, names, "").strip()
    if not value:
        raise ValueError(f"Exact SV evidence row is missing required column/value: {label}")
    return value


def _key(row: dict[str, str], default_build: str) -> str:
    build = first(row, ["genome_build", "build", "assembly"], default_build).strip()
    return canonical_breakpoint_key(
        build,
        _required(row, ["chrom1", "chr1"], "chrom1"),
        int(float(_required(row, ["pos1", "breakpoint1"], "pos1"))),
        _required(row, ["strand1", "orientation1"], "strand1"),
        _required(row, ["chrom2", "chr2"], "chrom2"),
        int(float(_required(row, ["pos2", "breakpoint2"], "pos2"))),
        _required(row, ["strand2", "orientation2"], "strand2"),
    )


@dataclass(frozen=True)
class ExactJunctionEvidence:
    adjacency_key: str
    split_reads: int
    spanning_reads: int
    unique_start_count: int
    min_anchor_bp: int
    min_mapq: int
    source_tool: str
    source_record_id: str

    @property
    def support_reads(self) -> int:
        # Split reads establish the exact junction; spanning reads are auxiliary.
        return self.split_reads

    @property
    def qc_pass(self) -> bool:
        return (
            self.split_reads >= 3
            and self.unique_start_count >= 2
            and self.min_anchor_bp >= 10
            and self.min_mapq >= 20
        )


def load_exact_junction_evidence(
    path: str | Path | None,
    *,
    default_build: str,
) -> dict[str, ExactJunctionEvidence]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Missing RNA junction evidence: {p}")
    out: dict[str, ExactJunctionEvidence] = {}
    for row in read_tsv(p):
        key = _key(row, default_build)
        evidence = ExactJunctionEvidence(
            adjacency_key=key,
            split_reads=int(to_float(first(row, ["split_reads", "junction_reads"], "0"), 0)),
            spanning_reads=int(to_float(first(row, ["spanning_reads"], "0"), 0)),
            unique_start_count=int(to_float(first(row, ["unique_start_count", "unique_starts"], "0"), 0)),
            min_anchor_bp=int(to_float(first(row, ["min_anchor_bp", "anchor_bp"], "0"), 0)),
            min_mapq=int(to_float(first(row, ["min_mapq", "mapq"], "0"), 0)),
            source_tool=first(row, ["source_tool", "caller"], ""),
            source_record_id=first(row, ["source_record_id", "record_id"], ""),
        )
        previous = out.get(key)
        if previous is None or (evidence.qc_pass, evidence.support_reads) > (previous.qc_pass, previous.support_reads):
            out[key] = evidence
    return out


@dataclass(frozen=True)
class ExpressedProduct:
    adjacency_key: str
    gene1: str
    gene2: str
    transcript1: str
    transcript2: str
    protein_sequence: str
    wildtype_protein_sequence: str
    junction_aa_position: int
    in_frame: str
    orf_status: str
    source_tool: str
    source_record_id: str

    def to_reconstruction(self, event_id: str, sample_id: str) -> ProteinReconstruction:
        sequence = self.protein_sequence.upper().replace("*", "")
        junction = self.junction_aa_position
        return ProteinReconstruction(
            protein_sequence_id=f"PROT_{event_id}",
            event_id=event_id,
            sample_id=sample_id,
            gene=f"{self.gene1}::{self.gene2}",
            transcript_id=f"{self.transcript1}::{self.transcript2}",
            protein_type="SV_Fusion",
            protein_sequence=sequence,
            wt_protein_sequence=self.wildtype_protein_sequence.upper().replace("*", ""),
            wt_prefix_aa=self.wildtype_protein_sequence.upper().replace("*", "")[:junction],
            novel_aa=sequence[junction : junction + 80],
            junction_aa_position=str(junction),
            novel_start_aa=str(junction),
            frameshift_start_aa="" if self.in_frame == "yes" else str(junction),
            in_frame=self.in_frame,
            reconstruction_method=f"external_expressed_transcript:{self.source_tool or 'unspecified'}",
            reconstruction_confidence="high",
            reconstruction_reason=(
                f"exact_adjacency={self.adjacency_key};orf_status={self.orf_status};"
                f"source_record_id={self.source_record_id}"
            ),
        )


def load_expressed_products(
    path: str | Path | None,
    *,
    default_build: str,
) -> dict[str, ExpressedProduct]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Missing expressed-product table: {p}")
    out: dict[str, ExpressedProduct] = {}
    for row in read_tsv(p):
        key = _key(row, default_build)
        sequence = _required(row, ["protein_sequence", "fusion_protein"], "protein_sequence").upper().replace("*", "")
        junction = int(float(_required(row, ["junction_aa_position", "junction_aa"], "junction_aa_position")))
        if not (0 < junction < len(sequence)):
            raise ValueError(f"Invalid junction_aa_position={junction} for product {key}")
        orf_status = _required(row, ["orf_status", "translation_status"], "orf_status").upper()
        if orf_status not in {"CONFIRMED", "COMPLETE", "TRANSLATABLE"}:
            raise ValueError(f"Expressed product {key} is not a confirmed translatable ORF: {orf_status}")
        out[key] = ExpressedProduct(
            adjacency_key=key,
            gene1=_required(row, ["gene1", "left_gene"], "gene1"),
            gene2=_required(row, ["gene2", "right_gene"], "gene2"),
            transcript1=_required(row, ["transcript1", "left_transcript"], "transcript1"),
            transcript2=_required(row, ["transcript2", "right_transcript"], "transcript2"),
            protein_sequence=sequence,
            wildtype_protein_sequence=first(row, ["wildtype_protein_sequence", "wt_protein_sequence"], ""),
            junction_aa_position=junction,
            in_frame="yes" if first(row, ["in_frame", "frame_status"], "").lower() in {"yes", "true", "1", "in_frame", "inframe"} else "no",
            orf_status=orf_status,
            source_tool=first(row, ["source_tool", "caller"], ""),
            source_record_id=first(row, ["source_record_id", "record_id"], ""),
        )
    return out


def product_meta(product: ExpressedProduct) -> dict[str, Any]:
    return {
        "gene1": product.gene1,
        "gene2": product.gene2,
        "transcript1": product.transcript1,
        "transcript2": product.transcript2,
        "effect_class": "SV_Fusion",
        "fusion_in_frame": product.in_frame,
        "frameshift": "no" if product.in_frame == "yes" else "yes",
        "junction_aa_position": str(product.junction_aa_position),
        "reconstruction_status": "confirmed_expressed_product",
        "reconstruction_reason": (
            f"exact RNA product from {product.source_tool or 'external tool'}; "
            f"record={product.source_record_id or 'not supplied'}"
        ),
    }
