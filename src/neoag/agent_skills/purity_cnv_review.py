from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
from pathlib import Path
from typing import Any

from .common import ensure_dir, safe_float, write_json, write_tsv

TOOL_PATTERNS = {
    "FACETS": ["*facets*", "*cncf*", "*purity*"],
    "PURPLE": ["*purple*", "*amber*", "*cobalt*"],
    "Sequenza": ["*sequenza*", "*seqz*", "*cellularity*"],
    "ASCAT": ["*ascat*"],
}

PURITY_KEYS = ["purity", "cellularity", "tumour_content", "tumor_content", "tumour_purity", "tumor_purity"]
PLOIDY_KEYS = ["ploidy", "psi"]


def _file_matches_sample(path: Path, sample_id: str | None) -> bool:
    if not sample_id:
        return True
    sid = sample_id.lower()
    if sid in str(path).lower():
        return True
    try:
        text = path.read_text(encoding="utf-8", errors="replace")[:200_000].lower()
        return sid in text
    except Exception:
        return False


def _candidate_files(paths: list[Path], sample_id: str | None = None) -> list[Path]:
    files: list[Path] = []
    for base in paths:
        if not base or not base.exists():
            continue
        if base.is_file():
            files.append(base)
            continue
        for pat in ["*.tsv", "*.csv", "*.txt", "*.out", "*.json", "*.yaml", "*.yml"]:
            files.extend(p for p in base.rglob(pat) if p.is_file() and p.stat().st_size < 50_000_000)
    seen: set[str] = set()
    out: list[Path] = []
    for p in files:
        key = str(p.resolve())
        if key not in seen and _file_matches_sample(p, sample_id):
            seen.add(key)
            out.append(p)
    return sorted(out, key=lambda p: p.stat().st_mtime, reverse=True)


def _tool_for_file(path: Path) -> str | None:
    name = str(path).lower()
    if "facets" in name or "cncf" in name:
        return "FACETS"
    if "purple" in name or "amber" in name or "cobalt" in name:
        return "PURPLE"
    if "sequenza" in name or "seqz" in name or "cellularity" in name:
        return "Sequenza"
    if "ascat" in name:
        return "ASCAT"
    return None


def _first_numeric_by_keys(row: dict[str, Any], keys: list[str]) -> float | None:
    lower = {str(k).strip().lower(): v for k, v in row.items()}
    for key in keys:
        for col, val in lower.items():
            if key == col or key in col:
                parsed = safe_float(val, None)
                if parsed is not None:
                    return parsed
    return None


def _parse_json(path: Path) -> tuple[float | None, float | None, str]:
    try:
        data = json.loads(path.read_text(encoding="utf-8", errors="replace"))
    except Exception:
        return None, None, "json parse failed"
    found: list[dict[str, Any]] = []
    def walk(x: Any):
        if isinstance(x, dict):
            found.append(x)
            for v in x.values():
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
    walk(data)
    for row in found:
        purity = _first_numeric_by_keys(row, PURITY_KEYS)
        ploidy = _first_numeric_by_keys(row, PLOIDY_KEYS)
        if purity is not None or ploidy is not None:
            return purity, ploidy, "json keys"
    return None, None, "no purity/ploidy keys"


def _parse_table(path: Path) -> tuple[float | None, float | None, str]:
    text = path.read_text(encoding="utf-8", errors="replace")[:200_000]
    lines = [x for x in text.splitlines() if x.strip()]
    if not lines:
        return None, None, "empty"
    if path.suffix.lower() == ".json":
        return _parse_json(path)

    # key=value or key: value summaries.
    purity = None
    ploidy = None
    for line in lines[:300]:
        m = re.search(r"(?i)\b(purity|cellularity|tumou?r[_ -]?content|tumou?r[_ -]?purity)\b\s*[:=]\s*([0-9.]+)", line)
        if m and purity is None:
            purity = safe_float(m.group(2), None)
        m = re.search(r"(?i)\b(ploidy|psi)\b\s*[:=]\s*([0-9.]+)", line)
        if m and ploidy is None:
            ploidy = safe_float(m.group(2), None)
        if purity is not None and ploidy is not None:
            return purity, ploidy, "key-value text"

    delimiter = "\t" if "\t" in lines[0] else ","
    try:
        reader = csv.DictReader(lines, delimiter=delimiter)
        for row in reader:
            purity = _first_numeric_by_keys(row, PURITY_KEYS)
            ploidy = _first_numeric_by_keys(row, PLOIDY_KEYS)
            if purity is not None or ploidy is not None:
                return purity, ploidy, "table columns"
    except Exception:
        pass
    return purity, ploidy, "not found"


def collect_tool_results(paths: list[Path], sample_id: str | None = None) -> list[dict[str, Any]]:
    best: dict[str, dict[str, Any]] = {}
    for path in _candidate_files(paths, sample_id=sample_id):
        tool = _tool_for_file(path)
        if not tool:
            continue
        purity, ploidy, method = _parse_table(path)
        if purity is None and ploidy is None:
            continue
        row = {
            "tool": tool,
            "status": "FOUND",
            "purity": purity if purity is not None else "",
            "ploidy": ploidy if ploidy is not None else "",
            "source_file": str(path),
            "parse_method": method,
            "notes": "",
        }
        current = best.get(tool)
        if current is None or (current.get("purity") == "" and purity is not None):
            best[tool] = row
    rows = [best[t] for t in ["FACETS", "PURPLE", "Sequenza", "ASCAT"] if t in best]
    for tool in ["FACETS", "PURPLE", "Sequenza", "ASCAT"]:
        if tool not in best:
            rows.append({"tool": tool, "status": "MISSING", "purity": "", "ploidy": "", "source_file": "", "parse_method": "", "notes": "no parsed result found"})
    return rows


def consensus(rows: list[dict[str, Any]]) -> dict[str, Any]:
    vals = []
    for row in rows:
        val = safe_float(row.get("purity"), None)
        if val is not None:
            vals.append((row["tool"], val))
    if not vals:
        return {"status": "NO_PURITY", "recommended_purity": "", "range": "", "interpretation": "No parsed purity value found."}
    nums = [v for _, v in vals]
    med = statistics.median(nums)
    spread = max(nums) - min(nums) if len(nums) > 1 else 0.0
    if len(nums) == 1:
        interp = f"Only {vals[0][0]} produced a parsed purity; use cautiously and cross-check with another tool."
        status = "SINGLE_TOOL"
    elif spread <= 0.08:
        interp = "Tools are broadly concordant; median purity is a reasonable working value."
        status = "CONCORDANT"
    elif spread <= 0.18:
        interp = "Tools differ moderately; use a range and inspect BAF/depth plots before choosing one value."
        status = "MODERATE_DISCORDANCE"
    else:
        interp = "Tools are strongly discordant; do not use a single purity without reviewing QC, SNP references, coverage, and CNV signal."
        status = "STRONG_DISCORDANCE"
    return {
        "status": status,
        "recommended_purity": round(med, 4),
        "range": f"{min(nums):.4f}-{max(nums):.4f}",
        "n_tools": len(nums),
        "tool_values": dict(vals),
        "interpretation": interp,
    }


def select_recommended_tool(rows: list[dict[str, Any]], recommended_purity: Any) -> dict[str, Any] | None:
    target = safe_float(recommended_purity, None)
    candidates = [row for row in rows if safe_float(row.get("purity"), None) is not None]
    if target is None or not candidates:
        return None
    with_ploidy = [row for row in candidates if safe_float(row.get("ploidy"), None) is not None]
    if with_ploidy:
        candidates = with_ploidy
    # Prefer a completed robust FACETS result on ties, then the estimate nearest the consensus median.
    return min(
        candidates,
        key=lambda row: (
            abs(float(row["purity"]) - target),
            0 if row.get("tool") == "FACETS" and "robust" in str(row.get("source_file", "")).lower() else 1,
        ),
    )


def find_segment_file(source_file: str) -> Path | None:
    source = Path(source_file)
    roots = [source.parent, source.parent.parent]
    patterns = ["cnv_segments.tsv", "*.purple.cnv.somatic.tsv", "*segments*.tsv", "*cncf*.tsv"]
    for root in roots:
        if not root.is_dir():
            continue
        for pattern in patterns:
            hit = next((p for p in sorted(root.glob(pattern)) if p.is_file() and p.stat().st_size > 0), None)
            if hit:
                return hit
    for pattern in patterns:
        hit = next((p for p in sorted(source.parent.rglob(pattern)) if p.is_file() and p.stat().st_size > 0), None)
        if hit:
            return hit
    return None


def normalize_segment_file(source: Path, output: Path) -> int:
    aliases = {
        "chromosome": ("chromosome", "chrom", "chr"),
        "start": ("start", "start.pos", "loc.start"),
        "end": ("end", "end.pos", "loc.end"),
        "total_cn": ("total_cn", "copyNumber", "CNt", "tcn.em"),
        "major_cn": ("major_cn", "majorAlleleCopyNumber", "A", "lcn.em"),
        "minor_cn": ("minor_cn", "minorAlleleCopyNumber", "B"),
        "loh_status": ("loh_status", "LOH", "status"),
    }

    def value(row: dict[str, str], field: str) -> str:
        return next((row[key] for key in aliases[field] if key in row and row[key] not in (None, "")), "")

    rows: list[dict[str, str]] = []
    with source.open(encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chromosome = value(row, "chromosome")
            start, end = value(row, "start"), value(row, "end")
            if not chromosome or not start or not end:
                continue
            minor = value(row, "minor_cn")
            loh = value(row, "loh_status")
            if not loh:
                minor_value = safe_float(minor)
                loh = "LOH" if minor_value is not None and minor_value <= 0.25 else "RETAINED" if minor_value is not None else "UNASSESSED"
            rows.append({
                "chromosome": chromosome,
                "start": start,
                "end": end,
                "total_cn": value(row, "total_cn"),
                "major_cn": value(row, "major_cn"),
                "minor_cn": minor,
                "loh_status": loh,
                "source_file": str(source),
            })
    write_tsv(output, rows, ["chromosome", "start", "end", "total_cn", "major_cn", "minor_cn", "loh_status", "source_file"])
    return len(rows)


def command_suggestions(tumor_bam: str | None, normal_bam: str | None) -> list[dict[str, str]]:
    if not tumor_bam or not normal_bam:
        return []
    return [
        {"tool": "FACETS", "command": f"bash scripts/run_facets_sample.sh --tumor-bam {tumor_bam} --normal-bam {normal_bam} --outdir results/purity_cnv/facets"},
        {"tool": "Sequenza", "command": f"bash scripts/run_sequenza_sample_by_chrom.sh --tumor-bam {tumor_bam} --normal-bam {normal_bam} --outdir results/purity_cnv/sequenza"},
        {"tool": "PURPLE", "command": f"bash scripts/run_purple_suite_container.sh --tumor-bam {tumor_bam} --normal-bam {normal_bam} --outdir results/purity_cnv/purple"},
    ]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Review FACETS/PURPLE/Sequenza/ASCAT purity and CNV outputs")
    ap.add_argument("--result-dir", action="append", default=[])
    ap.add_argument("--file", action="append", default=[])
    ap.add_argument("--tumor-bam")
    ap.add_argument("--normal-bam")
    ap.add_argument("--sample-id", help="Only parse result files whose path or content mentions this sample ID")
    ap.add_argument("--project-root", default=".")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)

    project_root = Path(args.project_root).resolve()
    outdir = ensure_dir(args.outdir)
    search_paths = [Path(x) for x in args.result_dir if x]
    search_paths += [Path(x) for x in args.file if x]
    search_paths += [outdir.parent, project_root / "results"]

    rows = collect_tool_results(search_paths, sample_id=args.sample_id)
    cons = consensus(rows)
    selected = select_recommended_tool(rows, cons.get("recommended_purity"))
    suggestions = command_suggestions(args.tumor_bam, args.normal_bam)

    write_tsv(outdir / "purity_cnv_tool_summary.tsv", rows, ["tool", "status", "purity", "ploidy", "source_file", "parse_method", "notes"])
    write_tsv(outdir / "purity_cnv_consensus.tsv", [{
        "sample_id": args.sample_id or "",
        "status": cons.get("status", ""),
        "recommended_purity": cons.get("recommended_purity", ""),
        "purity_range": cons.get("range", ""),
        "n_tools": cons.get("n_tools", 0),
        "interpretation": cons.get("interpretation", ""),
    }], ["sample_id", "status", "recommended_purity", "purity_range", "n_tools", "interpretation"])
    recommended_row = {
        "sample_id": args.sample_id or "",
        "purity": selected.get("purity", "") if selected else "",
        "ploidy": selected.get("ploidy", "") if selected else "",
        "evidence_tool": selected.get("tool", "UNASSESSED") if selected else "UNASSESSED",
        "evidence_file": selected.get("source_file", "") if selected else "",
        "consensus_status": cons.get("status", "NO_PURITY"),
    }
    write_tsv(outdir / "recommended_purity.tsv", [recommended_row], list(recommended_row))
    segment_source = find_segment_file(str(selected.get("source_file", ""))) if selected else None
    if segment_source:
        normalized_segment_count = normalize_segment_file(segment_source, outdir / "recommended_cnv_segments.tsv")
    else:
        normalized_segment_count = 0
        write_tsv(outdir / "recommended_cnv_segments.tsv", [], ["chromosome", "start", "end", "total_cn", "major_cn", "minor_cn", "loh_status", "source_file"])
    write_json(outdir / "recommended_cnv_segments.provenance.json", {
        "selected_tool": selected.get("tool", "UNASSESSED") if selected else "UNASSESSED",
        "source_file": str(segment_source or ""),
        "normalized_segment_count": normalized_segment_count,
        "selection_rule": "purity estimate nearest cross-tool median; robust FACETS preferred on ties",
    })
    write_tsv(outdir / "purity_cnv_run_suggestions.tsv", suggestions, ["tool", "command"])
    write_json(outdir / "purity_recommendation.json", cons)

    md = [
        "# Purity/CNV cross-tool review",
        "",
        f"Sample filter: `{args.sample_id}`" if args.sample_id else "Sample filter: not set",
        "",
        "## Consensus",
        f"- status: {cons.get('status')}",
        f"- recommended purity: {cons.get('recommended_purity') or 'NA'}",
        f"- purity range: {cons.get('range') or 'NA'}",
        f"- interpretation: {cons.get('interpretation')}",
        "",
        "## Tool summary",
        "| Tool | Status | Purity | Ploidy | Source |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        md.append(f"| {row['tool']} | {row['status']} | {row.get('purity','')} | {row.get('ploidy','')} | `{row.get('source_file','')}` |")
    md += [
        "",
        "## Recommended interpretation rules",
        "- Concordant values across FACETS/PURPLE/Sequenza/ASCAT are stronger than any single tool.",
        "- Strong disagreement should trigger review of SNP reference, GC/reference build, coverage, BAF separation, and low-purity signal.",
        "- FACETS can be sensitive to cval and SNP set; PURPLE depends on AMBER/COBALT references; Sequenza depends on FASTA/GC wiggle consistency; ASCAT depends on loci/alleles build and chr naming.",
        "- Use this consensus as a computational working value for CCF/HLA LOH/scoring, not as a clinical assertion.",
        "",
        "## Suggested missing runs",
    ]
    if suggestions:
        for item in suggestions:
            md.append(f"- {item['tool']}: `{item['command']}`")
    else:
        md.append("- No tumor/normal BAM pair was provided, so no run commands were generated.")
    md += [
        "",
        "## Outputs",
        f"- tool summary: `{outdir / 'purity_cnv_tool_summary.tsv'}`",
        f"- recommendation: `{outdir / 'purity_recommendation.json'}`",
        f"- run suggestions: `{outdir / 'purity_cnv_run_suggestions.tsv'}`",
    ]
    (outdir / "purity_cnv_review.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
