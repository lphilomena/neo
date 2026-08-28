from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

from .common import ensure_dir, write_json, write_tsv

TOOL_NAMES = ["EasyFuse", "STAR-Fusion", "Arriba", "FusionCatcher"]


def infer_tool(path: Path) -> str:
    s = str(path).lower()
    if "easyfuse" in s:
        return "EasyFuse"
    if "star-fusion" in s or "star_fusion" in s or "starfusion" in s:
        return "STAR-Fusion"
    if "arriba" in s:
        return "Arriba"
    if "fusioncatcher" in s or "fusion-catcher" in s:
        return "FusionCatcher"
    return "Unknown"


def looks_like_fusion_file(path: Path) -> bool:
    s = str(path).lower()
    name = path.name.lower()
    if path.stat().st_size > 100_000_000:
        return False
    if any(x in name for x in ["agent.log", "ranked_peptides", "ranked_events", "raw_peptides", "task_spec", "context", "case_state", "coordinator_plan", "final_response", "recommendation"]):
        return False
    if "agent_fusion" in s or "neoag-fusion-rna-run" in s:
        return False
    known_result_name = any(x in name for x in ["fusions.tsv", "fusion_predictions.tsv", "arriba", "star-fusion", "star_fusion", "fusioncatcher", "final-list", "easyfuse"] )
    known_result_path = any(x in s for x in ["/easyfuse/", "/arriba/", "/star-fusion/", "/star_fusion/", "/fusioncatcher/"])
    return known_result_name or known_result_path


def candidate_files(paths: list[Path], sample_id: str | None = None, max_files: int = 600) -> list[Path]:
    pats = ["*.tsv", "*.csv", "*.txt", "*.json", "*.out", "*.html"]
    files: list[Path] = []
    for base in paths:
        if not base or not base.exists():
            continue
        if base.is_file():
            if looks_like_fusion_file(base):
                files.append(base)
        else:
            for pat in pats:
                files.extend(p for p in base.rglob(pat) if p.is_file() and looks_like_fusion_file(p))
    sid = sample_id.lower() if sample_id else ""
    out: list[Path] = []
    seen: set[str] = set()
    for p in sorted(files, key=lambda x: x.stat().st_mtime, reverse=True):
        key = str(p.resolve())
        if key in seen:
            continue
        if sid and sid not in str(p).lower():
            try:
                if sid not in p.read_text(encoding="utf-8", errors="replace")[:200_000].lower():
                    continue
            except Exception:
                continue
        seen.add(key)
        out.append(p)
        if len(out) >= max_files:
            break
    return out


def parse_rows(path: Path, limit: int = 5000) -> list[dict[str, str]]:
    text = path.read_text(encoding="utf-8", errors="replace")[:2_000_000]
    if not text.strip():
        return []
    lines = [x for x in text.splitlines() if x.strip() and not x.startswith("##")]
    if not lines:
        return []
    delim = "\t" if "\t" in lines[0] else ","
    try:
        reader = csv.DictReader(lines, delimiter=delim)
        rows = []
        for i, row in enumerate(reader):
            if i >= limit:
                break
            rows.append({str(k): (v or "") for k, v in row.items()})
        return rows
    except Exception:
        return []


def fusion_id(row: dict[str, str]) -> str:
    lower = {str(k).lower(): v for k, v in row.items()}
    candidates = [
        lower.get("fusion"), lower.get("fusion_name"), lower.get("#fusionname"), lower.get("fusionname"),
        lower.get("gene1") and lower.get("gene2") and f"{lower.get('gene1')}--{lower.get('gene2')}",
        lower.get("left_gene") and lower.get("right_gene") and f"{lower.get('left_gene')}--{lower.get('right_gene')}",
    ]
    for c in candidates:
        if c:
            return str(c).replace("::", "--")
    genes = []
    for key, val in lower.items():
        if "gene" in key and val and len(val) < 80:
            genes.append(val)
    if len(genes) >= 2:
        return f"{genes[0]}--{genes[1]}"
    return ""


def collect_results(paths: list[Path], sample_id: str | None = None) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    summaries: list[dict[str, Any]] = []
    calls: list[dict[str, Any]] = []
    for path in candidate_files(paths, sample_id=sample_id):
        tool = infer_tool(path)
        rows = parse_rows(path)
        fusions = []
        for row in rows:
            fid = fusion_id(row)
            if fid and fid not in fusions:
                fusions.append(fid)
        summaries.append({
            "tool": tool,
            "source_file": str(path),
            "rows": len(rows),
            "n_unique_fusions": len(fusions),
            "top_fusions": ";".join(fusions[:20]),
        })
        for fid in fusions[:200]:
            calls.append({"tool": tool, "fusion": fid, "source_file": str(path)})
    return summaries, calls


def consensus(calls: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_fusion: dict[str, set[str]] = {}
    for row in calls:
        by_fusion.setdefault(row["fusion"], set()).add(row["tool"])
    out = []
    for fusion, tools in by_fusion.items():
        out.append({
            "fusion": fusion,
            "support_tools": ",".join(sorted(tools)),
            "n_tools": len(tools),
            "status": "MULTI_TOOL" if len(tools) >= 2 else "SINGLE_TOOL",
        })
    return sorted(out, key=lambda r: (-int(r["n_tools"]), r["fusion"]))


def suggestions(fastq1: str | None, fastq2: str | None, bam: str | None) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    if fastq1 and fastq2:
        out += [
            {"tool": "EasyFuse", "command": f"bash scripts/run_easyfuse_sample.sh --fastq1 {fastq1} --fastq2 {fastq2} --outdir results/fusion/easyfuse"},
            {"tool": "STAR-Fusion", "command": f"STAR-Fusion --left_fq {fastq1} --right_fq {fastq2} --genome_lib_dir $CTAT_GENOME_LIB --output_dir results/fusion/star-fusion"},
            {"tool": "FusionCatcher", "command": f"fusioncatcher -d $FUSIONCATCHER_DB -i {fastq1},{fastq2} -o results/fusion/fusioncatcher"},
        ]
    if bam:
        out.append({"tool": "Arriba", "command": f"bash scripts/run_arriba_sample.sh --bam {bam} --outdir results/fusion/arriba"})
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Run/review RNA fusion workflows across EasyFuse, STAR-Fusion, Arriba, and FusionCatcher")
    ap.add_argument("--result-dir", action="append", default=[])
    ap.add_argument("--file", action="append", default=[])
    ap.add_argument("--fastq1")
    ap.add_argument("--fastq2")
    ap.add_argument("--bam")
    ap.add_argument("--sample-id")
    ap.add_argument("--project-root", default=".")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)

    root = Path(args.project_root).resolve()
    outdir = ensure_dir(args.outdir)
    search_paths = [Path(x) for x in args.result_dir if x] + [Path(x) for x in args.file if x]
    if not search_paths:
        for rel in ["results/fusion", "results/easyfuse", "results/arriba", "results/star-fusion", "results/fusioncatcher", "results/llm_agent_web/neoag-sliding-run/run-full/parsed"]:
            p = root / rel
            if p.exists():
                search_paths.append(p)
    search_paths.append(outdir.parent)

    summaries, calls = collect_results(search_paths, sample_id=args.sample_id)
    cons = consensus(calls)
    run_suggestions = suggestions(args.fastq1, args.fastq2, args.bam)
    tool_counts = Counter(row["tool"] for row in calls)

    write_tsv(outdir / "fusion_tool_summary.tsv", summaries, ["tool", "source_file", "rows", "n_unique_fusions", "top_fusions"])
    write_tsv(outdir / "fusion_consensus.tsv", cons, ["fusion", "support_tools", "n_tools", "status"])
    write_tsv(outdir / "fusion_run_suggestions.tsv", run_suggestions, ["tool", "command"])
    write_json(outdir / "fusion_recommendation.json", {"sample_id": args.sample_id or "", "n_calls": len(calls), "tool_counts": dict(tool_counts), "consensus": cons[:100]})

    md = [
        "# RNA fusion cross-tool workflow review",
        "",
        f"Sample filter: `{args.sample_id}`" if args.sample_id else "Sample filter: not set",
        "",
        "## Tool summary",
        "| Tool | Rows | Unique fusions | Source |",
        "| --- | ---: | ---: | --- |",
    ]
    if summaries:
        for row in summaries[:80]:
            md.append(f"| {row['tool']} | {row['rows']} | {row['n_unique_fusions']} | `{row['source_file']}` |")
    else:
        md.append("| NA | 0 | 0 | no parsed fusion result found |")
    md += ["", "## Cross-tool fusion consensus", "| Fusion | Tools | N tools | Status |", "| --- | --- | ---: | --- |"]
    if cons:
        for row in cons[:50]:
            md.append(f"| {row['fusion']} | {row['support_tools']} | {row['n_tools']} | {row['status']} |")
    else:
        md.append("| NA |  | 0 | MISSING |")
    md += [
        "",
        "## Interpretation rules",
        "- Multi-tool supported fusions are stronger candidates than single-tool calls, but still require read-level review and optional orthogonal validation.",
        "- EasyFuse is a convenient Nextflow workflow; STAR-Fusion and FusionCatcher use CTAT/reference databases; Arriba is typically driven by STAR-aligned BAM/chimeric evidence.",
        "- Discordance is common because callers differ in filtering, reference annotations, and support thresholds.",
        "- For neoantigen use, prioritize fusions with in-frame ORF evidence, expression support, tumor specificity, and peptide-generation potential.",
        "",
        "## Suggested missing runs",
    ]
    if run_suggestions:
        for item in run_suggestions:
            md.append(f"- {item['tool']}: `{item['command']}`")
    else:
        md.append("- No RNA FASTQ pair or BAM was provided, so no run commands were generated.")
    md += [
        "",
        "## Outputs",
        f"- tool summary: `{outdir / 'fusion_tool_summary.tsv'}`",
        f"- consensus: `{outdir / 'fusion_consensus.tsv'}`",
        f"- recommendation: `{outdir / 'fusion_recommendation.json'}`",
        f"- run suggestions: `{outdir / 'fusion_run_suggestions.tsv'}`",
    ]
    (outdir / "fusion_rna_review.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
