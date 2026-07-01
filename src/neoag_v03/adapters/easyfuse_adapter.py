"""EasyFuse fusion discovery adapter — RNA fusion evidence, not full neoantigen scoring.

EasyFuse is a fusion metacaller (STAR-Fusion + FusionCatcher + Arriba + ML filter).
Project B uses it as Mode B fusion/RNA evidence input only; candidates still pass
HLA binding, safety, and score_v03 like any other source.
"""

from __future__ import annotations

import csv
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..model_layers import enrich_event_layers, enrich_peptide_layers, infer_mutation_source
from ..schemas import EVENT_FIELDS, FUSION_EVIDENCE_FIELDS, PEPTIDE_FIELDS
from ..utils import first, safe_id, to_float
from ..evidence_provenance import ProvenanceRecord, provenance_from_file, without_provenance, write_evidence_tsv

# Common read-through / low-confidence fusion types in normal tissue (EasyFuse taxonomy).
READTHROUGH_TYPES = frozenset({"cis_near"})


@dataclass
class EasyFuseFilterConfig:
    """Project B biological filters on top of EasyFuse prediction_class."""

    require_positive_class: bool = True
    min_prediction_prob: float = 0.5
    min_junction_reads: int = 3
    min_anchor_size: int = 10
    exclude_no_frame: bool = True
    require_neo_peptide: bool = True
    exclude_cis_near: bool = True
    exclude_readthrough_types: frozenset[str] = field(default_factory=lambda: READTHROUGH_TYPES)


def _easyfuse_delimiter(header: str) -> str | None:
    """Return delimiter when header looks like EasyFuse (BPID + Fusion_Gene)."""
    for delim in (";", "\t", ","):
        keys = {c.strip() for c in header.strip().split(delim) if c and c.strip()}
        if "BPID" in keys and "Fusion_Gene" in keys:
            return delim
    return None


def is_easyfuse_table(path: str | Path) -> bool:
    """Detect EasyFuse fusions.csv / fusions.pass.csv by header."""
    p = Path(path)
    with p.open("r", encoding="utf-8", errors="ignore") as fh:
        header = fh.readline()
    return _easyfuse_delimiter(header) is not None


def read_easyfuse_table(path: str | Path) -> list[dict[str, str]]:
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"Missing EasyFuse table: {p}")
    with p.open("r", encoding="utf-8", newline="") as fh:
        header = fh.readline()
        delim = _easyfuse_delimiter(header)
        if delim is None:
            delim = "\t"
        fh.seek(0)
        return [dict(row) for row in csv.DictReader(fh, delimiter=delim)]


def _gene_pair(row: dict[str, str]) -> tuple[str, str]:
    fusion_gene = first(row, ["Fusion_Gene", "fusion_gene"], "")
    if fusion_gene:
        if "::" in fusion_gene:
            parts = fusion_gene.split("::", 1)
            return parts[0], parts[1]
        if "_" in fusion_gene:
            parts = fusion_gene.split("_", 1)
            return parts[0], parts[1] if len(parts) > 1 else ""
        return fusion_gene, ""
    bpid = first(row, ["BPID", "bpid"], "")
    if bpid and "_" in bpid:
        return bpid.split("_", 1)[0], ""
    return "", ""


def _gene_label(g1: str, g2: str) -> str:
    if g1 and g2:
        return f"{g1}::{g2}"
    return g1 or g2 or "UNKNOWN"


def easyfuse_event_id(row: dict[str, str]) -> str:
    bpid = first(row, ["BPID", "bpid"], "")
    if bpid:
        return safe_id(f"EF_{bpid}")
    ftid = first(row, ["FTID", "ftid"], "")
    if ftid:
        return safe_id(f"EF_{ftid}")
    g1, g2 = _gene_pair(row)
    return safe_id(f"EF_{g1}_{g2}_{first(row, ['Breakpoint1', 'breakpoint1'], '')}")


def _neo_peptide_sequence(row: dict[str, str]) -> str:
    neo = str(first(row, ["neo_peptide_sequence", "Neo_peptide_sequence"], "")).strip()
    if neo in {"", "0", "NA", "na"}:
        return ""
    return neo


def _requant_count_max(row: dict[str, str], metric: str) -> int:
    """Max requant read count across ft/wt1/wt2 (EasyFuse v2 *_cnt_best) with v1 fallback."""
    vals: list[int] = []
    for key in (
        f"ft_{metric}_cnt",
        f"FT_{metric}_cnt",
        f"ft_{metric}_cnt_best",
        f"wt1_{metric}_cnt_best",
        f"wt2_{metric}_cnt_best",
    ):
        val = int(to_float(first(row, [key], "0"), 0.0))
        if val > 0:
            vals.append(val)
    if vals:
        return max(vals)
    # v1 fixtures: integer junction/span counts reported per caller tool.
    tool_keys = [f"{tool}_{metric}" for tool in (
        "star-fusion", "STAR-Fusion", "star", "starfusion",
        "arriba", "fusioncatcher", "FusionCatcher",
    )]
    total = 0
    for key in tool_keys:
        val = to_float(first(row, [key], "0"), 0.0)
        if val >= 1.0 and abs(val - round(val)) < 1e-6:
            total += int(round(val))
    return total


def _tool_junction_reads(row: dict[str, str]) -> int:
    return _requant_count_max(row, "junc")


def _tool_spanning_reads(row: dict[str, str]) -> int:
    return _requant_count_max(row, "span")


def _anchor_size(row: dict[str, str]) -> int:
    return _requant_count_max(row, "anch")


_TOOL_DETECTED_ALIASES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("star-fusion", (
        "star-fusion_detected", "STAR-Fusion_detected",
        "star_detected", "starfusion_detected",
    )),
    ("fusioncatcher", ("fusioncatcher_detected", "FusionCatcher_detected")),
    ("arriba", ("arriba_detected",)),
)


def _tools_detected(row: dict[str, str]) -> str:
    hits: list[str] = []
    for label, keys in _TOOL_DETECTED_ALIASES:
        if any(str(first(row, [k], "0")).strip() in {"1", "True", "true"} for k in keys):
            hits.append(label)
    return ",".join(sorted(set(hits)))


def filter_easyfuse_row(
    row: dict[str, str],
    cfg: EasyFuseFilterConfig | None = None,
) -> tuple[bool, str]:
    """Apply Project B filters; do not trust prediction_class alone."""
    cfg = cfg or EasyFuseFilterConfig()
    pred_class = str(first(row, ["prediction_class", "Prediction_class"], "")).strip().lower()
    pred_prob = to_float(first(row, ["prediction_prob", "Prediction_prob"], "0"), 0.0)
    frame = str(first(row, ["frame", "Frame"], "")).strip().lower()
    fusion_type = str(first(row, ["type", "Type"], "")).strip().lower()
    neo_pep = _neo_peptide_sequence(row)

    if cfg.require_positive_class and pred_class and pred_class != "positive":
        return False, "prediction_class_not_positive"
    if pred_prob < cfg.min_prediction_prob:
        return False, f"prediction_prob<{cfg.min_prediction_prob}"
    if cfg.exclude_no_frame and frame == "no_frame":
        return False, "frame_no_frame"
    if cfg.require_neo_peptide and not neo_pep:
        return False, "missing_neo_peptide"
    if _tool_junction_reads(row) < cfg.min_junction_reads:
        return False, f"junction_reads<{cfg.min_junction_reads}"
    if _anchor_size(row) < cfg.min_anchor_size:
        return False, f"anchor_size<{cfg.min_anchor_size}"
    if cfg.exclude_cis_near and fusion_type in cfg.exclude_readthrough_types:
        return False, f"readthrough_type:{fusion_type}"
    return True, "pass"


def easyfuse_row_to_event(
    row: dict[str, str],
    *,
    sample_id: str,
    profile_name: str,
    source_path: Path,
    filter_cfg: EasyFuseFilterConfig | None = None,
) -> tuple[dict[str, str] | None, dict[str, str]]:
    """Map one EasyFuse row to Layer-1 event + fusion evidence sidecar."""
    passed, reason = filter_easyfuse_row(row, filter_cfg)
    g1, g2 = _gene_pair(row)
    gene = _gene_label(g1, g2)
    eid = easyfuse_event_id(row)
    junction = _tool_junction_reads(row)
    pred_prob = to_float(first(row, ["prediction_prob"], "0"), 0.0)
    pred_class = first(row, ["prediction_class"], "")

    evidence = {
        "evidence_id": safe_id(f"FUSE_{eid}"),
        "event_id": eid,
        "sample_id": sample_id,
        "bpid": first(row, ["BPID", "bpid"], ""),
        "ftid": first(row, ["FTID", "ftid"], ""),
        "fusion_gene": first(row, ["Fusion_Gene", "fusion_gene"], gene.replace("::", "_")),
        "breakpoint1": first(row, ["Breakpoint1", "breakpoint1"], ""),
        "breakpoint2": first(row, ["Breakpoint2", "breakpoint2"], ""),
        "fusion_type": first(row, ["type", "Type"], ""),
        "frame_status": first(row, ["frame", "Frame"], ""),
        "bp1_frame": first(row, ["bp1_frame"], ""),
        "bp2_frame": first(row, ["bp2_frame"], ""),
        "exon_boundary": first(row, ["exon_boundary"], ""),
        "neo_peptide_sequence": first(row, ["neo_peptide_sequence"], ""),
        "fusion_protein_sequence": first(row, ["fusion_protein_sequence"], ""),
        "rna_junction_reads": str(junction),
        "rna_spanning_reads": str(_tool_spanning_reads(row)),
        "anchor_size": str(_anchor_size(row)),
        "caller_support_frac": str(to_float(first(row, ["tool_frac"], "0"), 0.0)),
        "caller_prob": f"{pred_prob:.4f}",
        "caller_pass": pred_class,
        "tools_detected": _tools_detected(row),
        "filter_status": "pass" if passed else "fail",
        "filter_reason": reason,
        "source_file": str(source_path),
    }

    if not passed:
        return None, evidence

    bp1 = first(row, ["Breakpoint1", "breakpoint1"], "")
    chrom = bp1.split(":", 1)[0] if bp1 else ""
    pos = bp1.split(":", 2)[1] if bp1.count(":") >= 1 else ""
    frame = str(first(row, ["frame", "Frame"], "")).lower()
    confidence = min(0.95, max(0.5, pred_prob))

    base = {
        "event_id": eid,
        "sample_id": sample_id,
        "disease_profile": profile_name,
        "event_type": "Fusion",
        "mutation_source": infer_mutation_source(event_type="Fusion", tool="EasyFuse"),
        "peptide_consequence": "fusion",
        "gene": gene,
        "event_name": first(row, ["Fusion_Gene", "fusion_gene"], gene.replace("::", "_")),
        "chrom": chrom,
        "pos": pos,
        "ref": "",
        "alt": "",
        "transcript_id": first(row, ["FTID", "ftid"], ""),
        "consequence": frame or "fusion",
        "rna_junction_reads": str(junction),
        "event_confidence": f"{confidence:.3f}",
        "event_expression": "0.0",
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
        "source": f"easyfuse:{source_path.name}",
    }
    return enrich_event_layers(base), evidence


def easyfuse_row_to_peptide(
    row: dict[str, str],
    event: dict[str, str],
    *,
    sample_id: str,
) -> dict[str, str] | None:
    """Deprecated: use build_easyfuse_catalog sliding-window catalog rows instead."""
    warnings.warn(
        "easyfuse_row_to_peptide is deprecated; use build_easyfuse_catalog",
        DeprecationWarning,
        stacklevel=2,
    )
    neo = _neo_peptide_sequence(row)
    if not neo or len(neo) < 8:
        return None
    # Prefer an 8–11 mer around the breakpoint when sequence is long.
    bp_pos = int(to_float(first(row, ["neo_peptide_sequence_bp"], "0"), 0.0))
    if len(neo) > 11 and bp_pos > 0:
        start = max(0, min(bp_pos - 6, len(neo) - 9))
        peptide = neo[start : start + 9]
    else:
        peptide = neo[:11] if len(neo) > 11 else neo

    frame = str(first(row, ["frame", "Frame"], "")).lower()
    base = {
        "peptide_id": safe_id(f"{event['event_id']}_{peptide}"),
        "event_id": event["event_id"],
        "sample_id": sample_id,
        "event_type": "Fusion",
        "mutation_source": event.get("mutation_source", ""),
        "peptide_consequence": "fusion",
        "gene": event.get("gene", ""),
        "peptide": peptide,
        "wildtype_peptide": "",
        "crosses_junction": "yes",
        "contains_novel_aa": "yes" if frame in {"neo_frame", "out_frame"} else "no",
        "rna_junction_reads": event.get("rna_junction_reads", "0"),
        "hla_allele": "",
        "mhc_class": "I",
        "source_tool": "EasyFuse",
        "binding_rank": "99",
        "el_rank": "99",
        "presentation_score": "0.0",
        "immunogenicity_score": "0.5",
        "wildtype_binding_rank": "99",
        "self_similarity_score": "0.0",
        "normal_hla_ligand_overlap": "no",
    }
    return enrich_peptide_layers(base, event)


def parse_easyfuse(
    path: str | Path,
    sample_id: str,
    profile_name: str,
    *,
    filter_cfg: EasyFuseFilterConfig | None = None,
    include_failed_evidence: bool = True,
    lengths: tuple[int, ...] = (8, 9, 10, 11),
    peptide_cfg: Any | None = None,
) -> dict[str, list[dict[str, str]]]:
    """Parse EasyFuse via unified catalog builder (legacy one-peptide path removed)."""
    from .easyfuse_variant_peptide import (
        EasyFusePeptideConfig,
        build_easyfuse_catalog,
    )

    result = build_easyfuse_catalog(
        path,
        sample_id,
        profile_name,
        lengths=lengths,
        filter_cfg=filter_cfg,
        peptide_cfg=peptide_cfg or EasyFusePeptideConfig(),
    )
    fusion_evidence = result.fusion_evidence
    if not include_failed_evidence:
        fusion_evidence = [r for r in fusion_evidence if r.get("filter_status") == "pass"]

    return {
        "events": result.events,
        "peptides": [],
        "catalog_rows": result.catalog_rows,
        "fusion_evidence": fusion_evidence,
        "filter_qc": result.filter_qc,
        "collapse_qc": result.collapse_qc,
        "summary_qc": [result.summary_qc],
    }


def write_fusion_evidence(
    rows: list[dict[str, str]],
    out_path: str | Path,
    provenance: ProvenanceRecord | None = None,
) -> None:
    src = rows[0].get("source_file") if rows else out_path
    prov = provenance or provenance_from_file("easyfuse", src, mode="passthrough")
    write_evidence_tsv(out_path, rows, without_provenance(FUSION_EVIDENCE_FIELDS), prov)


def merge_easyfuse_into_catalog(
    path: str | Path,
    sample_id: str,
    profile_name: str,
    *,
    filter_cfg: EasyFuseFilterConfig | None = None,
) -> dict[str, list[dict[str, str]]]:
    """Entry point for input_router: parse EasyFuse table with Project B filters."""
    return parse_easyfuse(path, sample_id, profile_name, filter_cfg=filter_cfg)
