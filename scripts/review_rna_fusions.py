#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

STAR_FUSION_PATTERNS = (
    "**/star-fusion.fusion_predictions.abridged.tsv",
    "**/star-fusion.fusion_predictions.tsv",
    "**/*fusion_predictions.abridged.tsv",
    "**/*fusion_predictions.tsv",
)
ARRIBA_PATTERNS = ("**/*.fusions.tsv", "**/fusions.tsv")
FUSIONCATCHER_PATTERNS = (
    "**/fusioncatcher.final-list.txt",
    "**/final-list_candidate-fusion-genes*.txt",
    "**/final-list_candidate-fusion-genes*",
)
JAFFAL_PATTERNS = ("**/jaffa_results.csv", "**/jaffal_results.csv")
EASYFUSE_PATTERNS = ("**/fusions.pass.csv",)


def read_rows(path: Path | None) -> list[dict[str, str]]:
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        first_line = handle.readline()
        handle.seek(0)
        delimiter = "\t" if "\t" in first_line else (";" if ";" in first_line else ",")
        return [{str(key): str(value or "") for key, value in row.items()} for row in csv.DictReader(handle, delimiter=delimiter)]


def existing(paths: list[Path | None]) -> list[Path]:
    seen: set[Path] = set()
    result: list[Path] = []
    for path in paths:
        if path and path.is_file() and path.stat().st_size > 0 and path not in seen:
            seen.add(path)
            result.append(path)
    return result


def discover_files(roots: list[Path], patterns: tuple[str, ...]) -> list[Path]:
    found: list[Path] = []
    seen: set[Path] = set()
    for root in roots:
        if not root or not root.exists():
            continue
        search_root = root if root.is_dir() else root.parent
        for pattern in patterns:
            for path in sorted(search_root.glob(pattern)):
                if path.is_file() and path.stat().st_size > 0 and path not in seen:
                    seen.add(path)
                    found.append(path)
    return found


def first(row: dict[str, str], names: list[str]) -> str:
    normalized = {key.lower().replace("#", ""): value for key, value in row.items()}
    for name in names:
        value = normalized.get(name.lower().replace("#", ""), "").strip()
        if value:
            return value
    return ""


def pair(row: dict[str, str]) -> str:
    combined = first(row, ["fusion", "fusion_name", "fusionname", "fusion_gene", "fusion genes", "fusion_genes"])
    if combined:
        normalized = combined.replace("--", "::")
        if "::" not in normalized and "_" in normalized:
            left, right = normalized.split("_", 1)
            return f"{left}::{right}"
        return normalized
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
    ap.add_argument("--jaffal", type=Path)
    ap.add_argument("--caller-root", action="append", type=Path, default=[])
    ap.add_argument("--normal-readthrough", type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()

    normal_pairs = {pair(row).upper() for row in read_rows(args.normal_readthrough) if pair(row)}
    evidence: dict[str, dict[str, object]] = defaultdict(lambda: {"tools": set(), "frames": set(), "junction_reads": []})
    roots = [root for root in args.caller_root if root]
    caller_paths = [
        ("EasyFuse", existing([args.easyfuse]) or discover_files(roots, EASYFUSE_PATTERNS)),
        ("STAR-Fusion", existing([args.star_fusion]) + discover_files(roots, STAR_FUSION_PATTERNS)),
        ("Arriba", existing([args.arriba]) + [path for path in discover_files(roots, ARRIBA_PATTERNS) if path.name != "fusions.pass.csv"]),
        ("FusionCatcher", existing([args.fusioncatcher]) + discover_files(roots, FUSIONCATCHER_PATTERNS)),
        ("JAFFAL", existing([args.jaffal]) + discover_files(roots, JAFFAL_PATTERNS)),
    ]
    for tool, paths in caller_paths:
        for path in existing(paths):
            for row in read_rows(path):
                if tool == "JAFFAL":
                    classification = first(row, ["classification", "confidence"])
                    if classification and classification.lower() != "highconfidence":
                        continue
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
