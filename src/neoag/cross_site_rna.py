from __future__ import annotations

"""Exact, provenance-preserving comparison of RNA evidence across tumor sites."""

from collections import defaultdict
from pathlib import Path
import re
from typing import Any, Mapping

from .utils import read_tsv, write_tsv


CROSS_SITE_FIELDS = [
    "event_id", "secondary_event_id", "secondary_sample_id", "sample_identity_status",
    "cross_site_status", "cross_site_exact_support", "cross_site_match_method",
    "cross_site_review_status", "cross_site_review_reason",
    "cross_site_event_key", "secondary_gene_expression_tpm",
    "secondary_transcript_expression_tpm", "secondary_rna_ref_reads",
    "secondary_rna_alt_reads", "secondary_rna_depth", "secondary_rna_vaf",
    "secondary_rna_junction_reads", "secondary_rna_support_status",
    "secondary_source_tools", "cross_site_reason",
]


def _text(row: Mapping[str, Any], *fields: str) -> str:
    for field in fields:
        value = str(row.get(field) or "").strip()
        if value:
            return value
    return ""


def _chrom(value: str) -> str:
    value = str(value or "").strip()
    if not value:
        return ""
    return value if value.lower().startswith("chr") else f"chr{value}"


def _event_type(row: Mapping[str, Any]) -> str:
    value = _text(row, "event_type", "mutation_source", "source").upper()
    if "FUSION" in value:
        return "FUSION"
    if "SPLICE" in value or "JUNCTION" in value:
        return "SPLICE"
    if "INDEL" in value or "FRAME" in value or "INSERT" in value or "DELET" in value:
        return "INDEL"
    return "SNV" if value else "OTHER"


def _gene_pair(row: Mapping[str, Any]) -> tuple[str, str] | None:
    left = _text(row, "gene5", "gene_5prime", "left_gene")
    right = _text(row, "gene3", "gene_3prime", "right_gene")
    if not left or not right:
        candidate = _text(row, "fusion_name", "event_name", "gene", "event_id")
        match = re.search(r"([A-Za-z0-9_.-]+)\s*(?:::|--|/|_)\s*([A-Za-z0-9_.-]+)", candidate)
        if match:
            left, right = match.group(1), match.group(2)
    if not left or not right:
        return None
    return left.upper(), right.upper()


def _breakpoints(row: Mapping[str, Any]) -> tuple[str, str] | None:
    left = _text(row, "breakpoint1", "breakpoint_5prime", "left_breakpoint", "breakpoint5")
    right = _text(row, "breakpoint2", "breakpoint_3prime", "right_breakpoint", "breakpoint3")
    if left and right:
        return left, right
    return None


def event_keys(row: Mapping[str, Any]) -> list[tuple[str, str, str]]:
    """Return (specificity, method, key); exact keys always precede broad keys."""
    track = _event_type(row)
    if track in {"SNV", "INDEL"}:
        chrom = _chrom(_text(row, "chrom", "chromosome"))
        pos = _text(row, "pos", "position", "start")
        ref = _text(row, "ref", "reference")
        alt = _text(row, "alt", "alternate")
        if chrom and pos and ref and alt:
            return [("exact", "BUILD_CHROM_POS_REF_ALT", f"GRCh38|{chrom}|{pos}|{ref}|{alt}")]
    if track == "SPLICE":
        canonical = _text(row, "canonical_junction_id")
        if canonical.startswith("SJ|"):
            return [("exact", "CANONICAL_JUNCTION_ID", canonical)]
        chrom = _chrom(_text(row, "junction_chrom", "chrom"))
        start = _text(row, "junction_start", "intron_start")
        end = _text(row, "junction_end", "intron_end")
        strand = _text(row, "junction_strand", "strand")
        build = _text(row, "genome_build") or "GRCh38"
        if chrom and start and end and strand in {"+", "-"}:
            return [("exact", "BUILD_CHROM_INTRON_STRAND", f"SJ|{build}|{chrom}|{start}|{end}|{strand}")]
    if track == "FUSION":
        pair = _gene_pair(row)
        breakpoints = _breakpoints(row)
        keys: list[tuple[str, str, str]] = []
        if pair and breakpoints:
            keys.append(("exact", "GENE_PAIR_EXACT_BREAKPOINTS", f"FUSION|{pair[0]}|{pair[1]}|{breakpoints[0]}|{breakpoints[1]}"))
        if pair:
            keys.append(("broad", "GENE_PAIR_ONLY", f"FUSION|{pair[0]}|{pair[1]}"))
        return keys
    event_id = _text(row, "event_id")
    return [("broad", "EVENT_ID_ONLY", event_id)] if event_id else []


def _number(row: Mapping[str, Any], *fields: str) -> float | None:
    value = _text(row, *fields)
    try:
        return float(value) if value else None
    except ValueError:
        return None


def _support(row: Mapping[str, Any], track: str) -> tuple[str, str]:
    if track in {"SNV", "INDEL"}:
        alt = _number(row, "rna_alt_reads")
        depth = _number(row, "rna_depth")
        if alt is not None and alt > 0:
            return "SUPPORTED", f"secondary RNA ALT reads={alt:g}; depth={depth if depth is not None else 'NA'}"
        if depth is not None and depth >= 10:
            return "NEGATIVE_ADEQUATE_COVERAGE", f"secondary RNA depth={depth:g}; ALT reads={alt or 0:g}"
        return "INDETERMINATE_LOW_POWER", f"secondary RNA depth={depth if depth is not None else 'NA'}"
    reads = _number(row, "rna_junction_reads", "junction_reads", "provided_rna_junction_reads")
    if reads is not None and reads > 0:
        return "SUPPORTED", f"secondary exact junction reads={reads:g}"
    return "INDETERMINATE_LOW_POWER", f"secondary exact junction reads={reads if reads is not None else 'NA'}"


def build_cross_site_rna_evidence(
    primary_events: str | Path,
    secondary_events: str | Path,
    output_tsv: str | Path,
    *,
    secondary_sample_id: str,
    identity_status: str = "UNASSESSED",
) -> dict[str, Any]:
    primary = read_tsv(primary_events)
    secondary = read_tsv(secondary_events)
    index: dict[str, list[tuple[str, str, dict[str, str]]]] = defaultdict(list)
    for row in secondary:
        for specificity, method, key in event_keys(row):
            index[key].append((specificity, method, row))

    identity = str(identity_status or "UNASSESSED").strip().upper()
    identity_confirmed = identity in {"CONFIRMED", "MATCH", "MATCHED", "PASS", "SAME_PATIENT"}
    output: list[dict[str, str]] = []
    counts: dict[str, int] = defaultdict(int)
    for row in primary:
        event_id = _text(row, "event_id")
        track = _event_type(row)
        match = None
        used_key = ""
        for specificity, method, key in event_keys(row):
            candidates = index.get(key, [])
            if len(candidates) == 1:
                match = (specificity, method, candidates[0][2])
                used_key = key
                break
            if len(candidates) > 1:
                match = ("ambiguous", method, candidates[0][2])
                used_key = key
                break
        if not match:
            status, exact, method, reason = "PRIMARY_ONLY", "no", "NO_EXACT_MATCH", "No matching secondary-site event"
            secondary_row: dict[str, str] = {}
        else:
            specificity, method, secondary_row = match
            support, support_reason = _support(secondary_row, track)
            exact_match = specificity == "exact" and support == "SUPPORTED"
            if specificity == "ambiguous":
                status, exact, reason = "AMBIGUOUS_SECONDARY_MATCH", "no", "Multiple secondary records share the same key"
            elif specificity != "exact":
                status, exact, reason = "GENE_PAIR_SHARED_BREAKPOINT_UNASSESSED", "no", "Gene pair matches but exact breakpoints are unavailable"
            elif support == "NEGATIVE_ADEQUATE_COVERAGE":
                status, exact, reason = "SECONDARY_NEGATIVE_ADEQUATE_COVERAGE", "no", support_reason
            elif support != "SUPPORTED":
                status, exact, reason = "SECONDARY_MATCH_LOW_POWER", "no", support_reason
            elif not identity_confirmed:
                status, exact, reason = "EXACT_MATCH_IDENTITY_UNASSESSED", "no", f"{support_reason}; sample identity={identity}"
            else:
                status, exact, reason = "EXACT_SHARED", "yes", support_reason
        counts[status] += 1
        review_status = "PASS" if status == "EXACT_SHARED" else (
            "UNASSESSED" if status == "PRIMARY_ONLY" else "REVIEW_REQUIRED"
        )
        review_reason = {
            "EXACT_MATCH_IDENTITY_UNASSESSED": "待核对：事件和直接RNA支持可精确匹配，但尚未完成两份样本的RNA-DNA指纹确认",
            "GENE_PAIR_SHARED_BREAKPOINT_UNASSESSED": "待核对：两部位检出相同融合基因对，但至少一侧缺少可比较的精确断点，不能判定为同一融合事件",
            "SECONDARY_NEGATIVE_ADEQUATE_COVERAGE": "待核对：另一部位覆盖充分但未检出直接RNA支持，应结合部位差异和样本成分解释",
            "SECONDARY_MATCH_LOW_POWER": "待核对：事件键可匹配，但另一部位的位点或junction覆盖不足，不能解释为阴性",
            "AMBIGUOUS_SECONDARY_MATCH": "待核对：同一匹配键对应多个记录，需人工确认坐标、方向和来源记录",
            "PRIMARY_ONLY": "当前仅主肿瘤部位获得证据；另一部位未匹配不代表阴性",
            "EXACT_SHARED": "两部位样本身份已确认，且事件与直接RNA支持均可精确匹配",
        }.get(status, reason)
        output.append({
            "event_id": event_id,
            "secondary_event_id": _text(secondary_row, "event_id"),
            "secondary_sample_id": secondary_sample_id,
            "sample_identity_status": identity,
            "cross_site_status": status,
            "cross_site_exact_support": exact,
            "cross_site_match_method": method,
            "cross_site_review_status": review_status,
            "cross_site_review_reason": review_reason,
            "cross_site_event_key": used_key,
            "secondary_gene_expression_tpm": _text(secondary_row, "gene_expression_tpm", "expression_tpm"),
            "secondary_transcript_expression_tpm": _text(secondary_row, "transcript_expression_tpm"),
            "secondary_rna_ref_reads": _text(secondary_row, "rna_ref_reads"),
            "secondary_rna_alt_reads": _text(secondary_row, "rna_alt_reads"),
            "secondary_rna_depth": _text(secondary_row, "rna_depth"),
            "secondary_rna_vaf": _text(secondary_row, "rna_vaf"),
            "secondary_rna_junction_reads": _text(secondary_row, "rna_junction_reads", "junction_reads", "provided_rna_junction_reads"),
            "secondary_rna_support_status": _text(secondary_row, "rna_support_status"),
            "secondary_source_tools": _text(secondary_row, "source_tools", "source_tool"),
            "cross_site_reason": reason,
        })
    write_tsv(output_tsv, output, CROSS_SITE_FIELDS)
    return {"rows": len(output), "status_counts": dict(sorted(counts.items())), "output": str(output_tsv)}


__all__ = ["CROSS_SITE_FIELDS", "build_cross_site_rna_evidence", "event_keys"]
