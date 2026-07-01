"""Unified EasyFuse catalog builder: filter → isoform collapse → sliding-window peptides."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

from ..utils import first, to_float, write_tsv
from .easyfuse_adapter import (
    EasyFuseFilterConfig,
    easyfuse_event_id,
    easyfuse_row_to_event,
    filter_easyfuse_row,
    read_easyfuse_table,
    _anchor_size,
    _gene_label,
    _gene_pair,
    _neo_peptide_sequence,
    _tool_junction_reads,
)

EasyFuseIsoformStrategy = Literal[
    "shortest_neo",
    "longest_neo",
    "max_junction_reads",
    "max_prediction_prob",
    "first_in_file",
]

ISOFORM_STRATEGIES: tuple[str, ...] = (
    "shortest_neo",
    "longest_neo",
    "max_junction_reads",
    "max_prediction_prob",
    "first_in_file",
)

FILTER_QC_FIELDS = [
    "sample_id",
    "row_index",
    "bpid",
    "ftid",
    "fusion_gene",
    "event_id",
    "filter_status",
    "filter_reason",
    "neo_peptide_length",
    "junction_reads",
    "anchor_size",
    "prediction_prob",
    "prediction_class",
    "frame",
]

COLLAPSE_QC_FIELDS = [
    "sample_id",
    "event_id",
    "fusion_gene",
    "isoform_strategy",
    "dedup_per_event",
    "candidate_count",
    "selected_row_index",
    "selected_bpid",
    "selected_ftid",
    "selected_neo_length",
    "selected_junction_reads",
    "selected_prediction_prob",
    "rejected_row_indices",
    "rejected_ftids",
]

SUMMARY_QC_FIELDS = [
    "sample_id",
    "source_file",
    "isoform_strategy",
    "junction_only",
    "dedup_per_event",
    "peptide_lengths",
    "input_rows",
    "pipeline_pass_rows",
    "unique_events",
    "selected_isoforms",
    "catalog_peptide_rows",
    "generation_method",
]


@dataclass
class EasyFusePeptideConfig:
    """Production defaults: strict sliding (junction-only + per-event collapse)."""

    junction_only: bool = True
    dedup_per_event: bool = True
    isoform_strategy: EasyFuseIsoformStrategy = "max_junction_reads"
    mini_len: int = 10
    minigene_total_len: int | None = None


@dataclass
class EasyFuseCatalogResult:
    catalog_rows: list[dict[str, str]] = field(default_factory=list)
    events: list[dict[str, str]] = field(default_factory=list)
    fusion_evidence: list[dict[str, str]] = field(default_factory=list)
    filter_qc: list[dict[str, str]] = field(default_factory=list)
    collapse_qc: list[dict[str, str]] = field(default_factory=list)
    summary_qc: dict[str, str] = field(default_factory=dict)


def fusion_generation_method(
    lengths: tuple[int, ...],
    *,
    junction_only: bool,
    dedup_per_event: bool,
    isoform_strategy: str = "shortest_neo",
) -> str:
    k = ",".join(str(x) for x in lengths)
    parts = [f"easyfuse_neo_peptide_sliding_window;k={k}"]
    parts.append("fusion_minigene=peptide_centered")
    if junction_only:
        parts.append("junction_covering_only")
    if dedup_per_event:
        parts.append("dedup_per_event")
    parts.append(f"isoform={isoform_strategy}")
    return ";".join(parts)


def build_fusion_centered_minigene(
    neo_sequence: str,
    *,
    peptide_start_aa: int | str,
    peptide_end_aa: int | str,
    mini_len: int = 10,
    minigene_total_len: int | None = None,
) -> str:
    neo = str(neo_sequence or "").strip().upper()
    if not neo or mini_len < 0:
        return ""
    try:
        start0 = max(0, int(peptide_start_aa) - 1)
        end0 = min(len(neo), int(peptide_end_aa))
    except (TypeError, ValueError):
        return neo
    if start0 >= end0:
        return neo
    if minigene_total_len is not None:
        from ..vep.extract_peptides import _centered_minigene_segments

        seq_f, center, seq_b = _centered_minigene_segments(
            neo, start0, end0, int(minigene_total_len)
        )
    else:
        seq_f = neo[max(0, start0 - mini_len):start0]
        center = neo[start0:end0]
        seq_b = neo[end0:end0 + mini_len]
    return f"{seq_f}|{center}|{seq_b}"


def sliding_fusion_neo_peptides(
    neo_sequence: str,
    lengths: tuple[int, ...],
    *,
    bp_pos: int = 0,
) -> list[dict[str, str]]:
    """8–11 aa windows on EasyFuse neo_peptide_sequence."""
    neo = str(neo_sequence or "").strip().upper()
    if not neo:
        return []
    min_len = min(lengths) if lengths else 8
    if len(neo) < min_len:
        return []

    windows: list[dict[str, str]] = []
    seen: set[str] = set()
    for length in lengths:
        if length > len(neo):
            continue
        for start in range(0, len(neo) - length + 1):
            peptide = neo[start : start + length]
            if peptide in seen:
                continue
            seen.add(peptide)
            mut_pos = ""
            if bp_pos > 0:
                rel = bp_pos - start
                if 1 <= rel <= length:
                    mut_pos = str(rel)
            windows.append({
                "mutant_peptide": peptide,
                "wildtype_peptide": "",
                "peptide_length": str(length),
                "peptide_start_aa": str(start + 1),
                "peptide_end_aa": str(start + length),
                "mutation_position_in_peptide": mut_pos,
            })
    return windows


def _breakpoint_chrom_pos(bp: str) -> tuple[str, str]:
    bp = str(bp or "").strip()
    if not bp:
        return "", ""
    chrom = bp.split(":", 1)[0]
    pos = bp.split(":", 2)[1] if bp.count(":") >= 1 else ""
    return chrom, pos


def _isoform_sort_key(
    strategy: str,
    row_index: int,
    row: dict[str, str],
) -> tuple:
    neo_len = len(_neo_peptide_sequence(row))
    junc = _tool_junction_reads(row)
    prob = to_float(first(row, ["prediction_prob", "Prediction_prob"], "0"), 0.0)
    ftid = first(row, ["FTID", "ftid"], "")
    if strategy == "longest_neo":
        return (-neo_len, -junc, -prob, row_index, ftid)
    if strategy == "max_junction_reads":
        return (-junc, -prob, -neo_len, row_index, ftid)
    if strategy == "max_prediction_prob":
        return (-prob, -junc, -neo_len, row_index, ftid)
    if strategy == "first_in_file":
        return (row_index, ftid)
    return (neo_len, row_index, ftid)


def _filter_qc_row(
    row: dict[str, str],
    *,
    sample_id: str,
    row_index: int,
    filter_status: str,
    filter_reason: str,
) -> dict[str, str]:
    fusion_gene = first(row, ["Fusion_Gene", "fusion_gene"], "")
    return {
        "sample_id": sample_id,
        "row_index": str(row_index),
        "bpid": first(row, ["BPID", "bpid"], ""),
        "ftid": first(row, ["FTID", "ftid"], ""),
        "fusion_gene": fusion_gene,
        "event_id": easyfuse_event_id(row),
        "filter_status": filter_status,
        "filter_reason": filter_reason,
        "neo_peptide_length": str(len(_neo_peptide_sequence(row))),
        "junction_reads": str(_tool_junction_reads(row)),
        "anchor_size": str(_anchor_size(row)),
        "prediction_prob": f"{to_float(first(row, ['prediction_prob', 'Prediction_prob'], '0'), 0.0):.4f}",
        "prediction_class": first(row, ["prediction_class", "Prediction_class"], ""),
        "frame": first(row, ["frame", "Frame"], ""),
    }


def _collapse_qc_row(
    *,
    sample_id: str,
    event_id: str,
    fusion_gene: str,
    peptide_cfg: EasyFusePeptideConfig,
    candidates: list[tuple[int, dict[str, str], dict[str, str]]],
    selected: tuple[int, dict[str, str], dict[str, str]],
) -> dict[str, str]:
    sel_idx, sel_row, _sel_event = selected
    rejected = [item for item in candidates if item[0] != sel_idx]
    return {
        "sample_id": sample_id,
        "event_id": event_id,
        "fusion_gene": fusion_gene,
        "isoform_strategy": peptide_cfg.isoform_strategy,
        "dedup_per_event": "yes" if peptide_cfg.dedup_per_event else "no",
        "candidate_count": str(len(candidates)),
        "selected_row_index": str(sel_idx),
        "selected_bpid": first(sel_row, ["BPID", "bpid"], ""),
        "selected_ftid": first(sel_row, ["FTID", "ftid"], ""),
        "selected_neo_length": str(len(_neo_peptide_sequence(sel_row))),
        "selected_junction_reads": str(_tool_junction_reads(sel_row)),
        "selected_prediction_prob": (
            f"{to_float(first(sel_row, ['prediction_prob', 'Prediction_prob'], '0'), 0.0):.4f}"
        ),
        "rejected_row_indices": ",".join(str(item[0]) for item in rejected),
        "rejected_ftids": ",".join(first(item[1], ["FTID", "ftid"], "") for item in rejected),
    }


def easyfuse_row_to_variant_peptide_rows(
    row: dict[str, str],
    event: dict[str, str],
    *,
    sample_id: str,
    lengths: tuple[int, ...],
    row_index: int,
    peptide_cfg: EasyFusePeptideConfig | None = None,
) -> list[dict[str, str]]:
    peptide_cfg = peptide_cfg or EasyFusePeptideConfig()
    neo = _neo_peptide_sequence(row)
    if not neo:
        return []
    bp_pos = int(to_float(first(row, ["neo_peptide_sequence_bp"], "0"), 0.0))
    g1, g2 = _gene_pair(row)
    gene = _gene_label(g1, g2)
    chrom, pos = _breakpoint_chrom_pos(first(row, ["Breakpoint1", "breakpoint1"], ""))
    frame = str(first(row, ["frame", "Frame"], "")).lower()
    fusion_gene = first(row, ["Fusion_Gene", "fusion_gene"], gene.replace("::", "_"))
    eid = event["event_id"]
    variant_key = eid

    generation_method = fusion_generation_method(
        lengths,
        junction_only=peptide_cfg.junction_only,
        dedup_per_event=peptide_cfg.dedup_per_event,
        isoform_strategy=peptide_cfg.isoform_strategy,
    )
    windows = sliding_fusion_neo_peptides(neo, lengths, bp_pos=bp_pos)
    if peptide_cfg.junction_only:
        windows = [w for w in windows if w.get("mutation_position_in_peptide")]

    out: list[dict[str, str]] = []
    for win_idx, win in enumerate(windows, start=1):
        pid = f"{sample_id}.EF{row_index}.{win_idx}"
        out.append({
            "peptide_id": pid,
            "gene": gene,
            "ensembl_gene_id": "",
            "transcript_id": first(row, ["FTID", "ftid"], ""),
            "hgvsc": "",
            "hgvsp": fusion_gene,
            "chrom": chrom,
            "pos": pos,
            "ref": "",
            "alt": "",
            "vaf": "",
            "consequence": frame or "fusion",
            "protein_position": "",
            "amino_acids": "",
            "multi_aa_flag": "fusion_neo",
            "minigene": build_fusion_centered_minigene(
                neo,
                peptide_start_aa=win["peptide_start_aa"],
                peptide_end_aa=win["peptide_end_aa"],
                mini_len=peptide_cfg.mini_len,
                minigene_total_len=peptide_cfg.minigene_total_len,
            ),
            "minigene_nt": "",
            "in_normal_proteome": "no",
            "variant_key": variant_key,
            "peptide_source": "easyfuse",
            "peptide_label": (
                f"GENE={gene}|FUSION={fusion_gene}|EF={eid}"
                f"|NEO={win['mutant_peptide']}|BP={bp_pos}"
            ),
            "generation_method": generation_method,
            "fusion_generation_method": generation_method,
            **win,
        })
    return out


def build_easyfuse_catalog(
    path: str | Path,
    sample_id: str,
    profile_name: str,
    *,
    lengths: tuple[int, ...] = (8, 9, 10, 11),
    filter_cfg: EasyFuseFilterConfig | None = None,
    peptide_cfg: EasyFusePeptideConfig | None = None,
) -> EasyFuseCatalogResult:
    """Single production path: filter QC → isoform collapse QC → sliding-window catalog."""
    source_path = Path(path)
    filter_cfg = filter_cfg or EasyFuseFilterConfig()
    peptide_cfg = peptide_cfg or EasyFusePeptideConfig()
    if peptide_cfg.isoform_strategy not in ISOFORM_STRATEGIES:
        raise ValueError(
            f"Unknown easyfuse_isoform_strategy: {peptide_cfg.isoform_strategy!r}. "
            f"Choose from: {', '.join(ISOFORM_STRATEGIES)}"
        )

    generation_method = fusion_generation_method(
        lengths,
        junction_only=peptide_cfg.junction_only,
        dedup_per_event=peptide_cfg.dedup_per_event,
        isoform_strategy=peptide_cfg.isoform_strategy,
    )

    catalog_rows: list[dict[str, str]] = []
    events: dict[str, dict[str, str]] = {}
    fusion_evidence: list[dict[str, str]] = []
    filter_qc: list[dict[str, str]] = []
    collapse_qc: list[dict[str, str]] = []
    by_event: dict[str, list[tuple[int, dict[str, str], dict[str, str]]]] = {}

    table_rows = read_easyfuse_table(source_path)
    pipeline_pass_rows = 0

    for row_index, row in enumerate(table_rows, start=1):
        passed, reason = filter_easyfuse_row(row, filter_cfg)
        event, evidence = easyfuse_row_to_event(
            row,
            sample_id=sample_id,
            profile_name=profile_name,
            source_path=source_path,
            filter_cfg=filter_cfg,
        )
        fusion_evidence.append(evidence)
        filter_qc.append(
            _filter_qc_row(
                row,
                sample_id=sample_id,
                row_index=row_index,
                filter_status=evidence["filter_status"],
                filter_reason=reason,
            )
        )
        if event is None:
            continue
        pipeline_pass_rows += 1
        by_event.setdefault(event["event_id"], []).append((row_index, row, event))

    selected_rows: list[tuple[int, dict[str, str], dict[str, str]]] = []
    for event_id, candidates in by_event.items():
        fusion_gene = first(candidates[0][1], ["Fusion_Gene", "fusion_gene"], "")
        if peptide_cfg.dedup_per_event:
            best = (
                min(
                    candidates,
                    key=lambda item: _isoform_sort_key(peptide_cfg.isoform_strategy, item[0], item[1]),
                )
                if len(candidates) > 1
                else candidates[0]
            )
            selected_rows.append(best)
            collapse_qc.append(
                _collapse_qc_row(
                    sample_id=sample_id,
                    event_id=event_id,
                    fusion_gene=fusion_gene,
                    peptide_cfg=peptide_cfg,
                    candidates=candidates,
                    selected=best,
                )
            )
        else:
            selected_rows.extend(candidates)
            for candidate in candidates:
                collapse_qc.append(
                    _collapse_qc_row(
                        sample_id=sample_id,
                        event_id=event_id,
                        fusion_gene=fusion_gene,
                        peptide_cfg=peptide_cfg,
                        candidates=[candidate],
                        selected=candidate,
                    )
                )

    for row_index, row, event in selected_rows:
        events[event["event_id"]] = event
        catalog_rows.extend(
            easyfuse_row_to_variant_peptide_rows(
                row,
                event,
                sample_id=sample_id,
                lengths=lengths,
                row_index=row_index,
                peptide_cfg=peptide_cfg,
            )
        )

    summary_qc = {
        "sample_id": sample_id,
        "source_file": str(source_path),
        "isoform_strategy": peptide_cfg.isoform_strategy,
        "junction_only": "yes" if peptide_cfg.junction_only else "no",
        "dedup_per_event": "yes" if peptide_cfg.dedup_per_event else "no",
        "peptide_lengths": ",".join(str(x) for x in lengths),
        "input_rows": str(len(table_rows)),
        "pipeline_pass_rows": str(pipeline_pass_rows),
        "unique_events": str(len(by_event)),
        "selected_isoforms": str(len(selected_rows)),
        "catalog_peptide_rows": str(len(catalog_rows)),
        "generation_method": generation_method,
    }

    return EasyFuseCatalogResult(
        catalog_rows=catalog_rows,
        events=list(events.values()),
        fusion_evidence=fusion_evidence,
        filter_qc=filter_qc,
        collapse_qc=collapse_qc,
        summary_qc=summary_qc,
    )


def write_easyfuse_qc_tables(
    result: EasyFuseCatalogResult,
    out_dir: str | Path,
    *,
    prefix: str = "easyfuse",
) -> dict[str, str]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        "filter_qc": str(out / f"{prefix}_filter_qc.tsv"),
        "collapse_qc": str(out / f"{prefix}_collapse_qc.tsv"),
        "summary_qc": str(out / f"{prefix}_summary_qc.tsv"),
    }
    write_tsv(paths["filter_qc"], result.filter_qc, FILTER_QC_FIELDS)
    write_tsv(paths["collapse_qc"], result.collapse_qc, COLLAPSE_QC_FIELDS)
    write_tsv(paths["summary_qc"], [result.summary_qc], SUMMARY_QC_FIELDS)
    return paths


def easyfuse_to_variant_peptide_rows(
    path: str | Path,
    sample_id: str,
    profile_name: str,
    *,
    lengths: tuple[int, ...] = (8, 9, 10, 11),
    filter_cfg: EasyFuseFilterConfig | None = None,
    peptide_cfg: EasyFusePeptideConfig | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    """Backward-compatible wrapper; legacy peptide list is always empty."""
    result = build_easyfuse_catalog(
        path,
        sample_id,
        profile_name,
        lengths=lengths,
        filter_cfg=filter_cfg,
        peptide_cfg=peptide_cfg,
    )
    return result.catalog_rows, result.events, []


def resolve_easyfuse_filter_config(cfg: dict[str, Any]) -> EasyFuseFilterConfig:
    inputs = cfg.get("inputs") or {}
    return EasyFuseFilterConfig(
        require_positive_class=bool(inputs.get("easyfuse_require_positive_class", True)),
        min_prediction_prob=float(inputs.get("easyfuse_min_prediction_prob", 0.5)),
        min_junction_reads=int(inputs.get("easyfuse_min_junction_reads", 3)),
        min_anchor_size=int(inputs.get("easyfuse_min_anchor_size", 10)),
        exclude_no_frame=bool(inputs.get("easyfuse_exclude_no_frame", True)),
        require_neo_peptide=bool(inputs.get("easyfuse_require_neo_peptide", True)),
        exclude_cis_near=bool(inputs.get("easyfuse_exclude_cis_near", True)),
    )


def resolve_easyfuse_peptide_config(cfg: dict[str, Any]) -> EasyFusePeptideConfig:
    inputs = cfg.get("inputs") or {}
    strategy = str(inputs.get("easyfuse_isoform_strategy", "max_junction_reads")).strip().lower()
    return EasyFusePeptideConfig(
        junction_only=bool(inputs.get("easyfuse_junction_only", True)),
        dedup_per_event=bool(inputs.get("easyfuse_dedup_per_event", True)),
        isoform_strategy=strategy,  # type: ignore[arg-type]
        mini_len=int(
            inputs.get("easyfuse_mini_len")
            or inputs.get("variant_peptide_mini_len")
            or inputs.get("mini_len")
            or 10
        ),
        minigene_total_len=(
            int(inputs.get("easyfuse_minigene_total_len") or inputs.get("variant_peptide_minigene_total_len"))
            if (inputs.get("easyfuse_minigene_total_len") or inputs.get("variant_peptide_minigene_total_len"))
            else None
        ),
    )
