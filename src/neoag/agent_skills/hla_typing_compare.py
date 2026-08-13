from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from .common import ensure_dir, write_json, write_tsv

HLA_RE = re.compile(r"(?:HLA[-_])?((?:A|B|C|DRB1|DQB1|DPB1|DQA1|DPA1)\*[0-9]{2,3}(?::[0-9A-Za-z]{2,3}){0,4})", re.I)
LOCUS_ORDER = ["A", "B", "C", "DRB1", "DQB1", "DPB1", "DQA1", "DPA1"]
DISPLAY = {"DRB1": "DR", "DQB1": "DQ", "DPB1": "DP", "DQA1": "DQA", "DPA1": "DPA"}


def norm_allele(raw: str) -> str:
    allele = raw.strip().upper().replace("HLA-", "")
    if "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    parts = [p for p in rest.split(":") if p]
    return gene + "*" + ":".join(parts)


def locus_of(allele: str) -> str:
    return allele.split("*", 1)[0]


def lowres(allele: str, fields: int = 2) -> str:
    if "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    parts = rest.split(":")[:fields]
    return gene + "*" + ":".join(parts)


def highres(allele: str) -> str:
    """Keep all reported fields while removing terminal expression/group suffixes."""
    return re.sub(r"(?<=[0-9])[NLSQGAP]+$", "", norm_allele(allele))


def infer_tool(path: Path) -> str:
    s = str(path).lower()
    if "optitype" in s:
        return "OptiType"
    if "hla-la" in s or "hlala" in s or "bestguess" in s:
        return "HLA-LA"
    if "spechla" in s or "spec_hla" in s:
        return "SpecHLA"
    return "Unknown"


def _looks_like_hla_result(path: Path) -> bool:
    name = path.name.lower()
    full = str(path).lower()
    return any(x in full for x in ["optitype", "spechla", "hla-la", "hlala", "hla_typing", "bestguess", "hla"]) and not any(x in name for x in ["ranked_peptides", "ranked_events", "raw_peptides", "raw_events", "agent.log"])


def candidate_files(paths: list[Path], sample_id: str | None = None, max_files: int = 500) -> list[Path]:
    pats = ["*.txt", "*.tsv", "*.csv", "*.json", "*.out", "*.hla", "*.HLA"]
    files: list[Path] = []
    for base in paths:
        if not base or not base.exists():
            continue
        if base.is_file():
            files.append(base)
        else:
            for pat in pats:
                files.extend(p for p in base.rglob(pat) if p.is_file() and p.stat().st_size < 20_000_000 and _looks_like_hla_result(p))
    seen: set[str] = set()
    out: list[Path] = []
    sid = sample_id.lower() if sample_id else ""
    for p in files:
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
    return sorted(out, key=lambda p: p.stat().st_mtime, reverse=True)[:max_files]


def parse_file(path: Path) -> list[str]:
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return []
    alleles = [norm_allele(m.group(1)) for m in HLA_RE.finditer(text)]
    # Keep order and unique alleles. Skip obvious reference/header-like duplicates by limiting per locus later.
    seen: set[str] = set()
    out: list[str] = []
    for allele in alleles:
        if allele not in seen:
            seen.add(allele)
            out.append(allele)
    return out


def collect_typing(paths: list[Path], sample_id: str | None = None) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in candidate_files(paths, sample_id=sample_id):
        alleles = parse_file(path)
        if not alleles:
            continue
        tool = infer_tool(path)
        by_locus: dict[str, list[str]] = defaultdict(list)
        for allele in alleles:
            loc = locus_of(allele)
            if loc in LOCUS_ORDER and allele not in by_locus[loc] and len(by_locus[loc]) < 2:
                by_locus[loc].append(allele)
        for loc, vals in by_locus.items():
            rows.append({
                "tool": tool,
                "locus": loc,
                "allele1": vals[0] if vals else "",
                "allele2": vals[1] if len(vals) > 1 else "",
                "lowres1": lowres(vals[0]) if vals else "",
                "lowres2": lowres(vals[1]) if len(vals) > 1 else "",
                "highres1": highres(vals[0]) if vals else "",
                "highres2": highres(vals[1]) if len(vals) > 1 else "",
                "source_file": str(path),
            })
    # Prefer known tools and recent files; keep first per tool/locus/source.
    unique: list[dict[str, Any]] = []
    seen_keys: set[tuple[str, str, str]] = set()
    for row in rows:
        key = (row["tool"], row["locus"], row["source_file"])
        if key not in seen_keys:
            seen_keys.add(key)
            unique.append(row)
    return unique


def consensus(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_locus: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_locus[row["locus"]].append(row)
    out: list[dict[str, Any]] = []
    for loc in LOCUS_ORDER:
        locus_rows = by_locus.get(loc, [])
        if not locus_rows:
            out.append({
                "locus": DISPLAY.get(loc, loc),
                "consensus_lowres": "",
                "consensus_2field": "",
                "consensus_high_resolution": "",
                "support": "0",
                "high_resolution_support": "0",
                "status": "MISSING",
                "high_resolution_status": "MISSING",
                "details": "no tool result",
            })
            continue
        pairs = []
        high_pairs = []
        for row in locus_rows:
            vals = sorted([x for x in [row.get("lowres1"), row.get("lowres2")] if x])
            pairs.append(" / ".join(vals))
            high_vals = sorted([x for x in [row.get("highres1"), row.get("highres2")] if x])
            high_pairs.append(" / ".join(high_vals))
        counts = Counter(pairs)
        best, support = counts.most_common(1)[0]
        high_counts = Counter(high_pairs)
        high_best, high_support = high_counts.most_common(1)[0]
        n_tools = len({r["tool"] for r in locus_rows})
        if support >= 2:
            status = "CONSENSUS"
        elif n_tools == 1:
            status = "SINGLE_TOOL"
        else:
            status = "DISCORDANT"
        out.append({
            "locus": DISPLAY.get(loc, loc),
            "consensus_lowres": best,
            "consensus_2field": best,
            "consensus_high_resolution": high_best,
            "support": f"{support}/{len(pairs)}",
            "high_resolution_support": f"{high_support}/{len(high_pairs)}",
            "status": status,
            "high_resolution_status": "CONSENSUS" if high_support >= 2 else "SINGLE_TOOL" if n_tools == 1 else "DISCORDANT",
            "details": "; ".join(f"{r['tool']}={r.get('lowres1','')}/{r.get('lowres2','')}" for r in locus_rows),
        })
    return out


def suggestions(bam: str | None, fastq1: str | None, fastq2: str | None) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    if bam:
        out.append({"tool": "SpecHLA", "command": f"run neoag-spechla-from-bam with --bam {bam}"})
        out.append({"tool": "HLA-LA", "command": f"bash scripts/run_hla_la_container.sh --bam {bam} --outdir results/hla_typing/hla-la"})
        out.append({"tool": "OptiType", "command": f"optitype can use BAM/FASTQ for HLA-I only; prefer FASTQ if available, otherwise configure BAM mode for {bam}"})
    if fastq1 and fastq2:
        out.append({"tool": "OptiType", "command": f"optitype -i {fastq1} {fastq2} --dna -o results/hla_typing/optitype"})
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Compare HLA typing calls across OptiType, SpecHLA, and HLA-LA")
    ap.add_argument("--result-dir", action="append", default=[])
    ap.add_argument("--file", action="append", default=[])
    ap.add_argument("--bam")
    ap.add_argument("--fastq1")
    ap.add_argument("--fastq2")
    ap.add_argument("--sample-id")
    ap.add_argument("--project-root", default=".")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)

    root = Path(args.project_root).resolve()
    outdir = ensure_dir(args.outdir)
    explicit_search_paths = bool(args.result_dir or args.file)
    search_paths = [Path(x) for x in args.result_dir if x] + [Path(x) for x in args.file if x]
    if not search_paths:
        for rel in ["results/hla_typing", "results/agent_hla_typing_compare_smoke", "results/llm_agent_web"]:
            p = root / rel
            if p.exists():
                search_paths.append(p)
    search_paths.append(outdir.parent)
    # Explicit result paths define the sample boundary. Tool outputs often use
    # fixed filenames (for example result.tsv or hla.result.txt) and need not
    # repeat the sample ID in either path or content.
    rows = collect_typing(search_paths, sample_id=None if explicit_search_paths else args.sample_id)
    cons = consensus(rows)
    run_suggestions = suggestions(args.bam, args.fastq1, args.fastq2)

    standardized_columns = ["tool", "locus", "allele1", "allele2", "lowres1", "lowres2", "highres1", "highres2", "source_file"]
    write_tsv(outdir / "hla_typing_tool_summary.tsv", rows, standardized_columns)
    write_tsv(outdir / "hla_typing.standardized.tsv", rows, standardized_columns)
    consensus_columns = ["locus", "consensus_lowres", "consensus_2field", "consensus_high_resolution", "support", "high_resolution_support", "status", "high_resolution_status", "details"]
    write_tsv(outdir / "hla_typing_consensus.tsv", cons, consensus_columns)
    write_tsv(outdir / "hla_typing_run_suggestions.tsv", run_suggestions, ["tool", "command"])
    write_json(outdir / "hla_typing_recommendation.json", {"sample_id": args.sample_id or "", "n_tool_locus_calls": len(rows), "consensus": cons})
    recommended: list[str] = []
    for row in cons:
        if row["locus"] in {"A", "B", "C"} and row["consensus_2field"]:
            recommended.extend(x.strip() for x in str(row["consensus_2field"]).split("/") if x.strip())
    (outdir / "recommended_hla.txt").write_text("\n".join(dict.fromkeys(recommended)) + ("\n" if recommended else ""), encoding="utf-8")

    md = [
        "# HLA typing cross-tool comparison",
        "",
        f"Sample filter: `{args.sample_id}`" if args.sample_id else "Sample filter: not set",
        "",
        "## Consensus by locus",
        "| Locus | 2-field consensus | High-resolution consensus | Support | Status | Details |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in cons:
        md.append(f"| {row['locus']} | {row['consensus_2field']} | {row['consensus_high_resolution']} | {row['support']} | {row['status']} | {row['details']} |")
    md += [
        "",
        "## Tool-level calls",
        "| Tool | Locus | Allele 1 | Allele 2 | Source |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in rows[:80]:
        md.append(f"| {row['tool']} | {DISPLAY.get(row['locus'], row['locus'])} | {row['allele1']} | {row['allele2']} | `{row['source_file']}` |")
    if len(rows) > 80:
        md.append(f"\nOnly first 80 calls shown; total {len(rows)} tool-locus calls.")
    md += [
        "",
        "## Interpretation rules",
        "- A/B/C can be compared across OptiType, SpecHLA, and HLA-LA; OptiType is HLA-I only.",
        "- DR/DQ/DP should be judged mainly from SpecHLA and HLA-LA; OptiType missing values are not discordance.",
        "- If tools agree at two-field resolution but differ at four-field resolution, report low-resolution concordance and high-resolution uncertainty.",
        "- Discordance should trigger review of input type, read length, tumor/normal source, reference build, graph/db version, and allele ambiguity.",
        "",
        "## Suggested missing runs",
    ]
    if run_suggestions:
        for item in run_suggestions:
            md.append(f"- {item['tool']}: `{item['command']}`")
    else:
        md.append("- No BAM/FASTQ was provided, so no run commands were generated.")
    md += [
        "",
        "## Outputs",
        f"- tool summary: `{outdir / 'hla_typing_tool_summary.tsv'}`",
        f"- consensus: `{outdir / 'hla_typing_consensus.tsv'}`",
        f"- recommendation: `{outdir / 'hla_typing_recommendation.json'}`",
        f"- run suggestions: `{outdir / 'hla_typing_run_suggestions.tsv'}`",
    ]
    (outdir / "hla_typing_compare.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    (outdir / "hla_typing_review.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
