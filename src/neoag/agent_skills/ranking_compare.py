from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
from typing import Any, Iterable

from .common import ensure_dir, markdown_table, read_tsv, write_tsv

TOP_NS = [10, 20, 50, 100]


def candidate_id(row: dict[str, str], idx: int) -> str:
    explicit = (row.get("peptide_id") or "").strip()
    if explicit:
        return explicit
    parts = [row.get("gene", ""), row.get("peptide", ""), row.get("hla_allele", ""), row.get("event_type", "")]
    key = "|".join(part.strip() for part in parts if part and part.strip())
    return key or f"row_{idx + 1}"


def row_ids(rows: list[dict[str, str]]) -> list[str]:
    seen: dict[str, int] = {}
    ids: list[str] = []
    for index, row in enumerate(rows):
        base = candidate_id(row, index)
        occurrence = seen.get(base, 0)
        seen[base] = occurrence + 1
        ids.append(base if occurrence == 0 else f"{base}#{occurrence + 1}")
    return ids


def rank_from_ids(ids: list[str]) -> dict[str, int]:
    return {candidate: index + 1 for index, candidate in enumerate(ids)}


def spearman_from_ranks(left: dict[str, int], right: dict[str, int]) -> float | None:
    keys = sorted(set(left) & set(right))
    if len(keys) < 2:
        return None
    xs = [left[key] for key in keys]
    ys = [right[key] for key in keys]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    vx = sum((value - mx) ** 2 for value in xs)
    vy = sum((value - my) ** 2 for value in ys)
    if vx <= 0 or vy <= 0:
        return None
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(vx * vy)


def _truthy(value: Any) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y"}


def _first(row: dict[str, str], fields: Iterable[str]) -> str:
    for field in fields:
        value = str(row.get(field, "")).strip()
        if value:
            return value
    return ""


def _has_conflict(row: dict[str, str]) -> bool:
    text = " ".join(_first(row, [field]) for field in (
        "evidence_conflict_fields", "evidence_conflict_layers",
        "presentation_consensus_state", "comprehensive_evidence_status",
    )).upper()
    return bool(str(row.get("evidence_conflict_fields", "")).strip()) or any(token in text for token in ("CONFLICT", "DISCORDANT"))


def _has_missing_evidence(row: dict[str, str]) -> bool:
    explicit = " ".join(str(row.get(field, "")).strip() for field in (
        "evidence_missing_layers", "safety_missing_layers", "event_safety_missing_layers",
    ))
    if explicit.strip():
        return True
    states = " ".join(str(value).upper() for key, value in row.items() if key.endswith("_state") or key.endswith("_status"))
    return any(token in states for token in ("UNASSESSED", "MISSING", "UNRESOLVED"))


def _manual_review(row: dict[str, str]) -> bool:
    return _truthy(row.get("manual_review_required")) or "MANUAL_REVIEW" in str(row.get("evidence_track", "")).upper()


def _hard_fail(row: dict[str, str]) -> bool:
    return _truthy(row.get("hard_failure")) or bool(str(row.get("hard_failure_codes", "")).strip())


def _quality_rows(name: str, rows: list[dict[str, str]]) -> list[dict[str, str]]:
    metrics = {
        "tool_or_evidence_conflict": sum(_has_conflict(row) for row in rows),
        "missing_evidence": sum(_has_missing_evidence(row) for row in rows),
        "manual_review": sum(_manual_review(row) for row in rows),
        "hard_failure": sum(_hard_fail(row) for row in rows),
    }
    denominator = max(len(rows), 1)
    return [
        {"ranking": name, "metric": metric, "count": str(count), "fraction": f"{count / denominator:.6f}"}
        for metric, count in metrics.items()
    ]


def _composition_rows(name: str, rows: list[dict[str, str]]) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    for top_n in TOP_NS:
        subset = rows[:top_n]
        for dimension, field in (("event_type", "event_type"), ("hla_allele", "hla_allele")):
            counts = Counter(str(row.get(field, "") or "UNSPECIFIED") for row in subset)
            for value, count in sorted(counts.items()):
                output.append({
                    "ranking": name, "top_n": str(top_n), "dimension": dimension,
                    "value": value, "count": str(count),
                    "fraction": f"{count / max(len(subset), 1):.6f}",
                })
    return output


def compare_rankings(
    left_path: str | Path,
    right_path: str | Path,
    outdir: str | Path,
    *,
    left_name: str = "left",
    right_name: str = "right",
) -> dict[str, Any]:
    output_dir = ensure_dir(outdir)
    _, left_rows = read_tsv(left_path)
    _, right_rows = read_tsv(right_path)
    left_ids, right_ids = row_ids(left_rows), row_ids(right_rows)
    left_rank, right_rank = rank_from_ids(left_ids), rank_from_ids(right_ids)
    left_by_id, right_by_id = dict(zip(left_ids, left_rows)), dict(zip(right_ids, right_rows))
    common = set(left_rank) & set(right_rank)
    correlation = spearman_from_ranks(left_rank, right_rank)

    overlap_rows = []
    for top_n in TOP_NS:
        left_top, right_top = set(left_ids[:top_n]), set(right_ids[:top_n])
        overlap = len(left_top & right_top)
        denominator = max(min(top_n, len(left_ids)), min(top_n, len(right_ids)), 1)
        union = len(left_top | right_top)
        overlap_rows.append({
            "top_n": str(top_n), "overlap": str(overlap),
            "overlap_fraction": f"{overlap / denominator:.6f}",
            "jaccard": f"{overlap / max(union, 1):.6f}",
        })
    write_tsv(output_dir / "topn_overlap.tsv", overlap_rows)

    change_rows: list[dict[str, str]] = []
    for candidate in sorted(set(left_ids) | set(right_ids)):
        left_position, right_position = left_rank.get(candidate), right_rank.get(candidate)
        row = right_by_id.get(candidate) or left_by_id.get(candidate) or {}
        if left_position is None:
            direction, delta = "RIGHT_ONLY", ""
        elif right_position is None:
            direction, delta = "LEFT_ONLY", ""
        else:
            signed = left_position - right_position
            direction = "PROMOTED_IN_RIGHT" if signed > 0 else "DEMOTED_IN_RIGHT" if signed < 0 else "UNCHANGED"
            delta = str(signed)
        change_rows.append({
            "candidate_id": candidate,
            "left_rank": str(left_position or ""), "right_rank": str(right_position or ""),
            "rank_delta_left_minus_right": delta, "change": direction,
            "left_grade": _first(left_by_id.get(candidate, {}), ["evidence_grade", "final_priority"]),
            "right_grade": _first(right_by_id.get(candidate, {}), ["evidence_grade", "final_priority"]),
            "gene": str(row.get("gene", "")), "peptide": str(row.get("peptide", "")),
            "hla_allele": str(row.get("hla_allele", "")), "event_type": str(row.get("event_type", "")),
            "hard_failure": "yes" if _hard_fail(row) else "no",
            "hard_failure_codes": str(row.get("hard_failure_codes", "")),
            "manual_review_required": "yes" if _manual_review(row) else "no",
            "evidence_missing_layers": str(row.get("evidence_missing_layers", "")),
            "evidence_conflict_layers": str(row.get("evidence_conflict_layers", "")),
        })
    change_rows.sort(key=lambda row: (
        0 if row["right_rank"] else 1,
        int(row["right_rank"] or 10**12),
        int(row["left_rank"] or 10**12),
        row["candidate_id"],
    ))
    write_tsv(output_dir / "candidate_rank_changes.tsv", change_rows)
    legacy_shift_rows = [
        {"peptide_id": row["candidate_id"], **{key: value for key, value in row.items() if key != "candidate_id"}}
        for row in change_rows
    ]
    write_tsv(output_dir / "rank_shift.tsv", legacy_shift_rows)

    high_hard_fail = [
        row for row in change_rows
        if row["hard_failure"] == "yes" and (
            0 < int(row["left_rank"] or 10**12) <= 100
            or 0 < int(row["right_rank"] or 10**12) <= 100
        )
    ]
    write_tsv(output_dir / "high_rank_hard_fail.tsv", high_hard_fail)

    composition_rows = _composition_rows(left_name, left_rows) + _composition_rows(right_name, right_rows)
    write_tsv(output_dir / "top_composition.tsv", composition_rows)
    quality_rows = _quality_rows(left_name, left_rows) + _quality_rows(right_name, right_rows)
    write_tsv(output_dir / "evidence_qc_summary.tsv", quality_rows)

    manual_rows = [row for row in change_rows if row["manual_review_required"] == "yes"]
    write_tsv(output_dir / "manual_review_candidates.tsv", manual_rows)
    shared_top20 = [row for row in change_rows if row["left_rank"] and row["right_rank"] and int(row["left_rank"]) <= 20 and int(row["right_rank"]) <= 20]
    write_tsv(output_dir / "shared_top20_candidates.tsv", shared_top20)

    summary = {
        "schema_version": "2.0",
        "left": {"name": left_name, "path": str(left_path), "rows": len(left_rows)},
        "right": {"name": right_name, "path": str(right_path), "rows": len(right_rows)},
        "common_candidates": len(common),
        "spearman_correlation": correlation,
        "high_rank_hard_fail_candidates": len(high_hard_fail),
        "manual_review_candidates": len(manual_rows),
        "outputs": {
            "topn_overlap": str(output_dir / "topn_overlap.tsv"),
            "candidate_rank_changes": str(output_dir / "candidate_rank_changes.tsv"),
            "high_rank_hard_fail": str(output_dir / "high_rank_hard_fail.tsv"),
            "top_composition": str(output_dir / "top_composition.tsv"),
            "evidence_qc_summary": str(output_dir / "evidence_qc_summary.tsv"),
            "manual_review_candidates": str(output_dir / "manual_review_candidates.tsv"),
        },
    }
    (output_dir / "ranking_comparison_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")

    event20 = [row for row in composition_rows if row["top_n"] == "20" and row["dimension"] == "event_type"]
    hla20 = [row for row in composition_rows if row["top_n"] == "20" and row["dimension"] == "hla_allele"]
    md = [
        f"# Ranking comparison: {left_name} vs {right_name}", "",
        f"Rows: {left_name}={len(left_rows)}, {right_name}={len(right_rows)}",
        f"Common candidate IDs: {len(common)}",
        f"Spearman rank correlation: {correlation:.6f}" if correlation is not None else "Spearman rank correlation: NA",
        "", "## Top-N overlap", markdown_table(overlap_rows),
        "", "## Candidate promotions and demotions", markdown_table(change_rows, max_rows=50),
        "", "## High-ranking hard-fail audit", markdown_table(high_hard_fail, max_rows=100),
        "", "## Top20 event-type composition", markdown_table(event20),
        "", "## Top20 HLA coverage", markdown_table(hla20),
        "", "## Conflict and missing-evidence rates", markdown_table(quality_rows),
        "", "## Manual-review candidates", markdown_table(manual_rows, max_rows=100),
        "", "## Interpretation boundary",
        "This comparison describes ranking behavior and evidence availability. It does not replace experimental validation or establish a treatment decision.",
    ]
    report_path = output_dir / "ranking_compare_report.md"
    report_path.write_text("\n".join(md) + "\n", encoding="utf-8")
    summary["outputs"]["report"] = str(report_path)
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Compare any two peptide rankings")
    parser.add_argument("--left")
    parser.add_argument("--left-name", default="left")
    parser.add_argument("--right")
    parser.add_argument("--right-name", default="right")
    parser.add_argument("--netmhcpan42", help="Compatibility alias for --left")
    parser.add_argument("--recommendation", help="Compatibility alias for --right")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args(argv)
    left = args.left or args.netmhcpan42
    right = args.right or args.recommendation
    if not left or not right:
        parser.error("--left and --right are required (legacy --netmhcpan42/--recommendation are accepted)")
    left_name = args.left_name if args.left else "netmhcpan42"
    right_name = args.right_name if args.right else "recommendation"
    compare_rankings(left, right, args.outdir, left_name=left_name, right_name=right_name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
