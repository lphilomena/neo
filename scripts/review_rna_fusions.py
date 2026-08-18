#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_rows(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        first_line = handle.readline()
        handle.seek(0)
        delimiter = "\t" if "\t" in first_line else ","
        return [{str(key): str(value or "") for key, value in row.items()} for row in csv.DictReader(handle, delimiter=delimiter)]


def first(row: dict[str, str], names: list[str]) -> str:
    normalized = {key.lower().replace("#", ""): value for key, value in row.items()}
    for name in names:
        value = normalized.get(name.lower().replace("#", ""), "").strip()
        if value:
            return value
    return ""


def pair(row: dict[str, str]) -> str:
    combined = first(row, ["fusion", "fusion_name", "fusionname", "fusion_gene"])
    if combined:
        return combined.replace("--", "::")
    left = first(row, ["gene1", "gene5", "left_gene", "gene_1_symbol(5end_fusion_partner)", "gene_1_symbol", "5end_fusion_partner"])
    right = first(row, ["gene2", "gene3", "right_gene", "gene_2_symbol(3end_fusion_partner)", "gene_2_symbol", "3end_fusion_partner"])
    return f"{left}::{right}" if left and right else ""


def write_rows(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--easyfuse", type=Path)
    ap.add_argument("--star-fusion", type=Path)
    ap.add_argument("--arriba", type=Path)
    ap.add_argument("--fusioncatcher", type=Path)
    ap.add_argument("--normal-readthrough", type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()

    normal_pairs = {pair(row).upper() for row in read_rows(args.normal_readthrough) if pair(row)}
    evidence: dict[str, dict[str, object]] = defaultdict(lambda: {"tools": set(), "frames": set(), "junction_reads": []})
    for tool, path in (("EasyFuse", args.easyfuse), ("STAR-Fusion", args.star_fusion), ("Arriba", args.arriba), ("FusionCatcher", args.fusioncatcher)):
        for row in read_rows(path):
            fusion = pair(row)
            if not fusion:
                continue
            evidence[fusion]["tools"].add(tool)  # type: ignore[union-attr]
            frame = first(row, ["frame", "reading_frame", "in_frame"])
            if frame:
                evidence[fusion]["frames"].add(frame)  # type: ignore[union-attr]
            reads = first(row, ["junction_reads", "junctionreadcount", "split_reads", "supporting_reads", "junction_reads1", "spanning_pairs", "spanning_unique_reads"])
            if reads:
                evidence[fusion]["junction_reads"].append(reads)  # type: ignore[union-attr]

    rows: list[dict[str, str]] = []
    for fusion, item in evidence.items():
        tools = sorted(item["tools"])  # type: ignore[arg-type]
        in_normal = fusion.upper() in normal_pairs if normal_pairs else None
        rows.append({
            "fusion": fusion,
            "support_tools": ",".join(tools),
            "n_tools": str(len(tools)),
            "frame_status": ";".join(sorted(item["frames"])),  # type: ignore[arg-type]
            "junction_reads": ";".join(item["junction_reads"]),  # type: ignore[arg-type]
            "normal_readthrough_status": "UNASSESSED" if in_normal is None else ("DETECTED" if in_normal else "NOT_DETECTED"),
            "status": "NORMAL_BACKGROUND_REVIEW" if in_normal else ("CROSS_VALIDATED" if len(tools) >= 2 else "SINGLE_TOOL"),
        })
    rows.sort(key=lambda row: (-int(row["n_tools"]), row["fusion"]))
    fields = ["fusion", "support_tools", "n_tools", "frame_status", "junction_reads", "normal_readthrough_status", "status"]
    write_rows(args.outdir / "fusion_consensus.tsv", fields, rows)
    write_rows(args.outdir / "fusion_background_review.tsv", fields, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
