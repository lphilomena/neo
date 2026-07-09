from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
from pathlib import Path
from statistics import median
from typing import Any

from .common import ensure_dir, safe_float, write_json, write_tsv


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def looks_like_expression(path: Path) -> bool:
    name = path.name.lower()
    s = str(path).lower()
    if path.stat().st_size > 200_000_000:
        return False
    if any(x in name for x in ["ranked_peptides", "ranked_events", "raw_peptides", "agent.log"]):
        return False
    return any(x in name for x in ["tpm", "expression", "gene_abund", "isoforms.results", "genes.results", "quant.sf"]) or "/expression" in s


def is_fastq(path: Path) -> bool:
    n = path.name.lower()
    return n.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))


def find_expression_files(paths: list[Path], sample_id: str | None = None) -> list[Path]:
    files: list[Path] = []
    for base in paths:
        if not base or not base.exists():
            continue
        if base.is_file():
            if looks_like_expression(base):
                files.append(base)
        else:
            for pat in ["*.tsv", "*.csv", "*.txt", "*.sf", "*.results"]:
                files.extend(p for p in base.rglob(pat) if p.is_file() and looks_like_expression(p))
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
    return out[:200]


def parse_expression(path: Path, limit: int = 200000) -> tuple[list[dict[str, str]], str]:
    rows: list[dict[str, str]] = []
    with open_text(path) as fh:
        first = fh.readline()
        if not first:
            return [], "empty"
        delim = "\t" if "\t" in first else ","
        fh.seek(0)
        reader = csv.DictReader(fh, delimiter=delim)
        fields = reader.fieldnames or []
        gene_cols = [c for c in fields if c.lower() in {"gene", "gene_id", "gene_name", "target_id", "name"} or "gene" in c.lower()]
        tpm_cols = [c for c in fields if c.lower() == "tpm" or c.lower().endswith("_tpm") or "tpm" in c.lower()]
        if not tpm_cols:
            return [], "no TPM column"
        gene_col = gene_cols[0] if gene_cols else fields[0]
        tpm_col = tpm_cols[0]
        for i, row in enumerate(reader):
            if i >= limit:
                break
            gene = row.get(gene_col, "")
            tpm = row.get(tpm_col, "")
            if gene:
                rows.append({"gene_id": gene, "tpm": tpm, "source_file": str(path)})
    return rows, f"gene_col={gene_col};tpm_col={tpm_col}"


def summarize_tpm(rows: list[dict[str, str]]) -> dict[str, Any]:
    vals = [safe_float(r.get("tpm"), None) for r in rows]
    nums = [v for v in vals if v is not None]
    expressed = [v for v in nums if v >= 1.0]
    if not nums:
        return {"n_genes": 0, "n_expressed_tpm_ge_1": 0, "median_tpm": "", "max_tpm": ""}
    return {"n_genes": len(nums), "n_expressed_tpm_ge_1": len(expressed), "median_tpm": round(median(nums), 4), "max_tpm": round(max(nums), 4)}


def fastq_summary(paths: list[Path]) -> list[dict[str, Any]]:
    out = []
    for p in paths:
        if is_fastq(p):
            out.append({"file": str(p), "size_bytes": p.stat().st_size, "role": "RNA_FASTQ"})
        elif p.name.lower().endswith(".bam"):
            out.append({"file": str(p), "size_bytes": p.stat().st_size, "role": "RNA_BAM"})
    return out


def suggestions(fastq1: str | None, fastq2: str | None, bam: str | None, outprefix: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if fastq1 and fastq2:
        rows += [
            {"mode": "Salmon", "command": f"bash scripts/run_salmon_fastq_to_tpm.sh --fastq1 {fastq1} --fastq2 {fastq2} --outdir {outprefix}/salmon"},
            {"mode": "RSEM", "command": f"bash scripts/run_rsem_fastq_to_tpm.sh --fastq1 {fastq1} --fastq2 {fastq2} --outdir {outprefix}/rsem"},
            {"mode": "Kallisto", "command": f"kallisto quant -i $KALLISTO_INDEX -o {outprefix}/kallisto {fastq1} {fastq2}"},
        ]
    if bam:
        rows.append({"mode": "BAM-counts-to-TPM", "command": f"featureCounts/RSEM from BAM {bam}, then scripts/counts_to_tpm_gtf.py with matching GTF"})
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Prepare or review RNA expression TPM evidence from RNA FASTQ/BAM/existing TPM files")
    ap.add_argument("--result-dir", action="append", default=[])
    ap.add_argument("--file", action="append", default=[])
    ap.add_argument("--fastq1")
    ap.add_argument("--fastq2")
    ap.add_argument("--bam")
    ap.add_argument("--sample-id")
    ap.add_argument("--method", choices=["auto", "salmon", "rsem"], default="auto")
    ap.add_argument("--execute", action="store_true", help="Actually run Salmon or RSEM when FASTQ pair and references are available")
    ap.add_argument("--project-root", default=".")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args(argv)

    root = Path(args.project_root).resolve()
    outdir = ensure_dir(args.outdir)
    input_paths = [Path(x) for x in args.result_dir if x] + [Path(x) for x in args.file if x]
    if args.fastq1:
        input_paths.append(Path(args.fastq1))
    if args.fastq2:
        input_paths.append(Path(args.fastq2))
    if args.bam:
        input_paths.append(Path(args.bam))
    if not input_paths:
        for rel in ["results/expression", "results/rna", "results/llm_agent_web/neoag-sliding-run/run-full/parsed"]:
            p = root / rel
            if p.exists():
                input_paths.append(p)

    execution_record: dict[str, Any] = {"executed": False, "method": args.method, "returncode": "", "cmd": "", "stdout": "", "stderr": ""}
    if args.execute and args.fastq1 and args.fastq2 and args.method in {"auto", "salmon"}:
        salmon_out = outdir / "salmon"
        cmd = [
            "bash", str(root / "scripts/run_salmon_fastq_to_tpm.sh"),
            "--fastq1", args.fastq1,
            "--fastq2", args.fastq2,
            "--outdir", str(salmon_out),
        ]
        if args.sample_id:
            cmd += ["--sample-id", args.sample_id]
        proc = subprocess.run(cmd, text=True, capture_output=True, cwd=str(root))
        execution_record.update({
            "executed": True,
            "method": "salmon",
            "returncode": proc.returncode,
            "cmd": " ".join(cmd),
            "stdout": proc.stdout[-4000:],
            "stderr": proc.stderr[-4000:],
        })
        input_paths.append(salmon_out)
    elif args.execute and args.fastq1 and args.fastq2 and args.method == "rsem":
        rsem_out = outdir / "rsem"
        cmd = [
            "bash", str(root / "scripts/run_rsem_fastq_to_tpm.sh"),
            "--fastq1", args.fastq1,
            "--fastq2", args.fastq2,
            "--outdir", str(rsem_out),
        ]
        if args.sample_id:
            cmd += ["--sample-id", args.sample_id]
        proc = subprocess.run(cmd, text=True, capture_output=True, cwd=str(root))
        execution_record.update({
            "executed": True,
            "method": "rsem",
            "returncode": proc.returncode,
            "cmd": " ".join(cmd),
            "stdout": proc.stdout[-4000:],
            "stderr": proc.stderr[-4000:],
        })
        input_paths.append(rsem_out)

    expr_files = find_expression_files(input_paths, sample_id=args.sample_id)
    normalized_rows: list[dict[str, str]] = []
    expr_summary_rows = []
    for path in expr_files[:20]:
        rows, method = parse_expression(path)
        stats = summarize_tpm(rows)
        expr_summary_rows.append({"source_file": str(path), "parse_method": method, **stats})
        if not normalized_rows and rows:
            normalized_rows = rows

    fq_summary = fastq_summary(input_paths)
    run_suggestions = suggestions(args.fastq1, args.fastq2, args.bam, str(outdir / "run"))

    write_tsv(outdir / "expression_input_summary.tsv", expr_summary_rows, ["source_file", "parse_method", "n_genes", "n_expressed_tpm_ge_1", "median_tpm", "max_tpm"])
    write_tsv(outdir / "rna_input_files.tsv", fq_summary, ["file", "size_bytes", "role"])
    write_tsv(outdir / "rna_tpm_run_suggestions.tsv", run_suggestions, ["mode", "command"])
    if normalized_rows:
        write_tsv(outdir / "gene_tpm.tsv", normalized_rows, ["gene_id", "tpm", "source_file"])
    write_json(outdir / "rna_expression_recommendation.json", {"sample_id": args.sample_id or "", "n_expression_files": len(expr_files), "has_gene_tpm": bool(normalized_rows), "execution": execution_record, "suggestions": run_suggestions})

    md = [
        "# RNA FASTQ/BAM to TPM review",
        "",
        f"Sample filter: `{args.sample_id}`" if args.sample_id else "Sample filter: not set",
        "",
        "## Existing expression files",
        "| Source | Parse method | Genes | TPM>=1 | Median TPM | Max TPM |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    if expr_summary_rows:
        for row in expr_summary_rows[:20]:
            md.append(f"| `{row['source_file']}` | {row['parse_method']} | {row['n_genes']} | {row['n_expressed_tpm_ge_1']} | {row['median_tpm']} | {row['max_tpm']} |")
    else:
        md.append("| NA | no TPM/expression file parsed | 0 | 0 |  |  |")
    md += ["", "## RNA input files", "| File | Size | Role |", "| --- | ---: | --- |"]
    if fq_summary:
        for row in fq_summary:
            md.append(f"| `{row['file']}` | {row['size_bytes']} | {row['role']} |")
    else:
        md.append("| NA | 0 | no FASTQ/BAM provided |")
    md += [
        "",
        "## Execution",
        f"- executed: {execution_record.get('executed')}",
        f"- method: {execution_record.get('method')}",
        f"- returncode: {execution_record.get('returncode')}",
    ]
    if execution_record.get("stderr"):
        md.append(f"- stderr: `{str(execution_record.get('stderr'))[:500]}`")
    md += [
        "",
        "## Suggested TPM generation commands",
    ]
    if run_suggestions:
        for row in run_suggestions:
            md.append(f"- {row['mode']}: `{row['command']}`")
    else:
        md.append("- No RNA FASTQ pair or BAM was provided, so no run commands were generated.")
    md += [
        "",
        "## Interpretation rules",
        "- Existing gene TPM can be used as expression evidence for neoantigen scoring after gene ID/name consistency checks.",
        "- FASTQ-based quantification requires a matching transcriptome/reference index and annotation version.",
        "- BAM-based expression depends on whether the BAM is genome-aligned or transcriptome-aligned; use a matching GTF/reference.",
        "- TPM evidence should be treated as support, not proof of peptide presentation.",
        "",
        "## Outputs",
        f"- expression summary: `{outdir / 'expression_input_summary.tsv'}`",
        f"- normalized gene TPM: `{outdir / 'gene_tpm.tsv'}`" if normalized_rows else "- normalized gene TPM: not generated",
        f"- run suggestions: `{outdir / 'rna_tpm_run_suggestions.tsv'}`",
        f"- recommendation: `{outdir / 'rna_expression_recommendation.json'}`",
    ]
    (outdir / "rna_expression_summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
