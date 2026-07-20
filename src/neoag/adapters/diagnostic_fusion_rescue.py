"""Diagnostic fusion whitelist rescue layer.

This layer preserves diagnostically important fusions observed in raw/unfiltered
fusion caller outputs even when they are absent from EasyFuse pass tables. It is
explicitly evidence-only by default: rescued rows document diagnostic support and
validation needs, but they are not automatically converted into ranked peptides.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..evidence_provenance import provenance_from_file
from ..schemas import DIAGNOSTIC_FUSION_RESCUE_FIELDS, EVIDENCE_PROVENANCE_FIELDS
from ..utils import first, safe_id, to_float, write_tsv
from .easyfuse_adapter import (
    _anchor_size,
    _gene_label,
    _gene_pair,
    _neo_peptide_sequence,
    _tool_junction_reads,
    _tool_spanning_reads,
    _tools_detected,
    is_easyfuse_table,
    read_easyfuse_table,
)

DEFAULT_DIAGNOSTIC_FUSION_WHITELIST = ("EWSR1_WT1", "EWSR1::WT1")


def normalize_fusion_label(value: Any) -> str:
    """Canonicalize fusion labels for whitelist matching."""
    label = str(value or "").strip().upper().replace("::", "_")
    return "_".join(x for x in label.split("_") if x)


def _input_list(value: Any, default: tuple[str, ...] = ()) -> list[str]:
    if value is None:
        return list(default)
    if isinstance(value, str):
        return [x.strip() for x in value.replace(";", ",").split(",") if x.strip()]
    return [str(x).strip() for x in value if str(x).strip()]


def _path_or_none(value: Any, root: Path | None = None) -> Path | None:
    if not value:
        return None
    p = Path(str(value))
    if not p.is_absolute() and root is not None:
        p = (root / p).resolve()
    return p if p.is_file() else None


def infer_unfiltered_easyfuse_path(pass_path: str | Path | None) -> Path | None:
    """Infer sibling fusions.csv from a fusions.pass.csv path when available."""
    if not pass_path:
        return None
    p = Path(pass_path)
    candidates = []
    if p.name == "fusions.pass.csv":
        candidates.append(p.with_name("fusions.csv"))
    if p.name == "fusions.pass.tsv":
        candidates.append(p.with_name("fusions.tsv"))
    candidates.append(p.with_name("fusions.csv"))
    try:
        resolved = p.resolve(strict=True)
    except OSError:
        resolved = p
    if resolved != p:
        if resolved.name == "fusions.pass.csv":
            candidates.append(resolved.with_name("fusions.csv"))
        candidates.append(resolved.with_name("fusions.csv"))
    for candidate in candidates:
        if candidate.is_file() and candidate != p:
            return candidate
    return None


def resolve_diagnostic_fusion_rescue_inputs(cfg: dict[str, Any], *, root: Path | None = None) -> dict[str, Any]:
    """Resolve rescue settings from run TOML inputs."""
    inputs = cfg.get("inputs") or {}
    enabled = bool(inputs.get("diagnostic_fusion_rescue_enabled", True))
    whitelist = _input_list(
        inputs.get("diagnostic_fusion_whitelist"),
        DEFAULT_DIAGNOSTIC_FUSION_WHITELIST,
    )
    whitelist_norm = sorted({normalize_fusion_label(x) for x in whitelist if normalize_fusion_label(x)})

    pass_path = _path_or_none(
        inputs.get("easyfuse_pass_csv") or inputs.get("easyfuse_tsv") or inputs.get("fusion_tsv"),
        root,
    )
    explicit_sources = _input_list(
        inputs.get("diagnostic_fusion_rescue_sources")
        or inputs.get("easyfuse_unfiltered_csv")
        or inputs.get("easyfuse_all_csv")
        or inputs.get("diagnostic_fusion_rescue_tsv")
    )
    source_paths = []
    for src in explicit_sources:
        p = _path_or_none(src, root)
        if p:
            source_paths.append(p)
    inferred = infer_unfiltered_easyfuse_path(pass_path)
    if inferred:
        source_paths.append(inferred)

    # Stable de-duplication while preserving user-provided order.
    seen: set[str] = set()
    unique_sources: list[Path] = []
    for p in source_paths:
        key = str(p.resolve())
        if key not in seen:
            seen.add(key)
            unique_sources.append(p)

    return {
        "enabled": enabled,
        "whitelist": whitelist_norm,
        "source_paths": unique_sources,
        "pass_path": pass_path,
        "generate_peptides": bool(inputs.get("diagnostic_fusion_rescue_generate_peptides", False)),
    }


def _pass_fusion_keys(pass_path: Path | None) -> set[tuple[str, str, str]]:
    if not pass_path or not pass_path.is_file() or not is_easyfuse_table(pass_path):
        return set()
    keys: set[tuple[str, str, str]] = set()
    for row in read_easyfuse_table(pass_path):
        fusion = normalize_fusion_label(first(row, ["Fusion_Gene", "fusion_gene"], ""))
        bp1 = first(row, ["Breakpoint1", "breakpoint1"], "")
        bp2 = first(row, ["Breakpoint2", "breakpoint2"], "")
        if fusion:
            keys.add((fusion, bp1, bp2))
    return keys


def _row_fusion_label(row: dict[str, str]) -> tuple[str, str, str]:
    g1, g2 = _gene_pair(row)
    label = _gene_label(g1, g2)
    raw = first(row, ["Fusion_Gene", "fusion_gene"], label.replace("::", "_"))
    return raw, label, normalize_fusion_label(raw or label)


def _tool_flag(row: dict[str, str], keys: list[str]) -> str:
    return "1" if any(str(first(row, [k], "0")).strip().lower() in {"1", "true", "yes"} for k in keys) else "0"


def diagnostic_rescue_rows_from_easyfuse(
    source_path: str | Path,
    *,
    sample_id: str,
    whitelist: list[str],
    pass_keys: set[tuple[str, str, str]] | None = None,
) -> list[dict[str, str]]:
    """Extract whitelisted diagnostic fusion evidence from EasyFuse-like tables."""
    path = Path(source_path)
    if not path.is_file() or not is_easyfuse_table(path):
        return []
    pass_keys = pass_keys or set()
    whitelist_set = {normalize_fusion_label(x) for x in whitelist}
    out: list[dict[str, str]] = []
    seen: set[tuple[str, str, str, str]] = set()

    for row in read_easyfuse_table(path):
        raw_label, display_label, norm_label = _row_fusion_label(row)
        if norm_label not in whitelist_set:
            continue
        bp1 = first(row, ["Breakpoint1", "breakpoint1"], "")
        bp2 = first(row, ["Breakpoint2", "breakpoint2"], "")
        key = (norm_label, bp1, bp2, first(row, ["FTID", "ftid"], ""))
        if key in seen:
            continue
        seen.add(key)
        g1, g2 = _gene_pair(row)
        in_pass = (norm_label, bp1, bp2) in pass_keys
        pred_class = first(row, ["prediction_class", "Prediction_class"], "")
        pred_prob = first(row, ["prediction_prob", "Prediction_prob"], "")
        rescue_reason = "diagnostic_whitelist_not_in_easyfuse_pass" if not in_pass else "diagnostic_whitelist_in_easyfuse_pass"
        peptide_status = "not_generated_by_default"
        if _neo_peptide_sequence(row):
            peptide_status = "available_not_generated_by_default"
        out.append({
            "rescue_id": safe_id(f"DFR_{sample_id}_{norm_label}_{bp1}_{bp2}_{len(out) + 1}"),
            "sample_id": sample_id,
            "fusion_gene": display_label,
            "fusion_gene_raw": raw_label,
            "fusion_gene_normalized": norm_label,
            "gene5": g1,
            "gene3": g2,
            "breakpoint1": bp1,
            "breakpoint2": bp2,
            "ftid": first(row, ["FTID", "ftid"], ""),
            "fusion_type": first(row, ["type", "Type"], ""),
            "frame_status": first(row, ["frame", "Frame"], ""),
            "neo_peptide_sequence": _neo_peptide_sequence(row),
            "neo_peptide_sequence_bp": first(row, ["neo_peptide_sequence_bp"], ""),
            "fusion_protein_sequence": first(row, ["fusion_protein_sequence", "fusion_peptide"], ""),
            "rna_junction_reads": str(_tool_junction_reads(row)),
            "rna_spanning_reads": str(_tool_spanning_reads(row)),
            "anchor_size": str(_anchor_size(row)),
            "star_detected": _tool_flag(row, ["star-fusion_detected", "STAR-Fusion_detected", "star_detected", "starfusion_detected"]),
            "fusioncatcher_detected": _tool_flag(row, ["fusioncatcher_detected", "FusionCatcher_detected"]),
            "arriba_detected": _tool_flag(row, ["arriba_detected"]),
            "tools_detected": _tools_detected(row),
            "tool_count": first(row, ["tool_count", "tool_frac"], ""),
            "prediction_class": pred_class,
            "prediction_prob": f"{to_float(pred_prob, 0.0):.4f}" if pred_prob else "",
            "easyfuse_pass_status": "in_pass_table" if in_pass else "not_in_pass_table",
            "diagnostic_whitelist_status": "whitelisted",
            "diagnostic_relevance": "diagnostic_fusion_evidence",
            "rescue_reason": rescue_reason,
            "peptide_generation_status": peptide_status,
            "source_file": str(path),
            "notes": "Evidence-only rescue; not automatically included in ranked neoantigen peptides.",
        })
    collapsed: dict[tuple[str, str, str], dict[str, str]] = {}
    collapsed_counts: dict[tuple[str, str, str], int] = {}
    for row in out:
        key = (row.get("fusion_gene_normalized", ""), row.get("breakpoint1", ""), row.get("breakpoint2", ""))
        collapsed_counts[key] = collapsed_counts.get(key, 0) + 1
        score = (
            int(to_float(row.get("rna_junction_reads"), 0.0)),
            int(to_float(row.get("rna_spanning_reads"), 0.0)),
            int(to_float(row.get("anchor_size"), 0.0)),
            1 if row.get("neo_peptide_sequence") else 0,
        )
        prev = collapsed.get(key)
        prev_score = (
            int(to_float(prev.get("rna_junction_reads"), 0.0)),
            int(to_float(prev.get("rna_spanning_reads"), 0.0)),
            int(to_float(prev.get("anchor_size"), 0.0)),
            1 if prev and prev.get("neo_peptide_sequence") else 0,
        ) if prev else (-1, -1, -1, -1)
        if prev is None or score > prev_score:
            collapsed[key] = row
    rows = list(collapsed.values())
    for row in rows:
        key = (row.get("fusion_gene_normalized", ""), row.get("breakpoint1", ""), row.get("breakpoint2", ""))
        n = collapsed_counts.get(key, 1)
        if n > 1:
            row["notes"] = row.get("notes", "") + f" Collapsed from {n} unfiltered EasyFuse rows for this breakpoint."
    return rows


def write_diagnostic_fusion_rescue(
    rows: list[dict[str, str]],
    out_path: str | Path,
    *,
    source_path: str | Path | None = None,
) -> None:
    p = Path(out_path)
    fields = DIAGNOSTIC_FUSION_RESCUE_FIELDS + EVIDENCE_PROVENANCE_FIELDS
    if not rows:
        write_tsv(p, [], fields)
        return
    prov = provenance_from_file("diagnostic_fusion_rescue", source_path or rows[0].get("source_file"), mode="converted")
    prov_fields = prov.as_fields()
    write_tsv(p, [{**row, **prov_fields} for row in rows], fields)


def build_diagnostic_fusion_rescue(
    cfg: dict[str, Any],
    *,
    sample_id: str,
    out_path: str | Path,
    root: Path | None = None,
) -> dict[str, str]:
    """Build a diagnostic rescue sidecar if configured and source tables exist."""
    resolved = resolve_diagnostic_fusion_rescue_inputs(cfg, root=root)
    if not resolved["enabled"] or not resolved["whitelist"] or not resolved["source_paths"]:
        return {}
    pass_keys = _pass_fusion_keys(resolved["pass_path"])
    rows: list[dict[str, str]] = []
    for source_path in resolved["source_paths"]:
        rows.extend(
            diagnostic_rescue_rows_from_easyfuse(
                source_path,
                sample_id=sample_id,
                whitelist=resolved["whitelist"],
                pass_keys=pass_keys,
            )
        )
    write_diagnostic_fusion_rescue(rows, out_path, source_path=resolved["source_paths"][0])
    return {
        "diagnostic_fusion_rescue": str(out_path),
        "diagnostic_fusion_rescue_rows": str(len(rows)),
        "diagnostic_fusion_rescue_generate_peptides": str(resolved["generate_peptides"]).lower(),
    }
