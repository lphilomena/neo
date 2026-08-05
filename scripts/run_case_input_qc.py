#!/usr/bin/env python3
"""Inventory and validate tumor/normal DNA, somatic VCF, and optional RNA inputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


NORMAL_TOKENS = ("blood", "normal", "germline", "pbmc", "buffy")
TUMOR_TOKENS = ("tumor", "tumour", "liver", "lesion", "met", "cancer")
RNA_TOKENS = ("rna", "transcriptome")


def run(cmd: list[str], log: Path, check: bool = False) -> subprocess.CompletedProcess[str]:
    proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(proc.stdout or "", encoding="utf-8")
    if check and proc.returncode:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(cmd)}; see {log}")
    return proc


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def index_for(path: Path) -> Path | None:
    candidates = [Path(str(path) + ".bai"), path.with_suffix(".bai")]
    if path.suffix.lower() == ".cram":
        candidates.insert(0, Path(str(path) + ".crai"))
    return next((p for p in candidates if p.is_file() and p.stat().st_size > 0), None)


def discover(root: Path) -> tuple[list[Path], list[Path], list[Path]]:
    bams: list[Path] = []
    vcfs: list[Path] = []
    rna: list[Path] = []
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        low = str(path).lower()
        if any(part in low for part in ("/.git/", "/work/", "/results/")):
            continue
        if path.suffix.lower() in {".bam", ".cram"}:
            (rna if any(t in low for t in RNA_TOKENS) else bams).append(path)
        elif low.endswith((".vcf", ".vcf.gz", ".bcf")):
            vcfs.append(path)
        elif low.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")) and any(t in low for t in RNA_TOKENS):
            rna.append(path)
    return sorted(bams), sorted(vcfs), sorted(rna)


def choose_pair(bams: list[Path]) -> tuple[Path | None, Path | None, str]:
    normal = [p for p in bams if any(t in str(p).lower() for t in NORMAL_TOKENS)]
    tumor = [p for p in bams if p not in normal and any(t in str(p).lower() for t in TUMOR_TOKENS)]
    if len(normal) == 1 and len(tumor) == 1:
        return tumor[0], normal[0], "AUTO_NAME_TOKENS"
    if len(normal) == 1 and len(bams) == 2:
        return next(p for p in bams if p != normal[0]), normal[0], "AUTO_TWO_BAMS"
    return None, None, "EXPLICIT_SELECTION_REQUIRED"


def choose_vcf(vcfs: list[Path]) -> Path | None:
    ranked = sorted(
        vcfs,
        key=lambda p: (
            "somatic" not in str(p).lower(),
            "pass" not in str(p).lower(),
            len(str(p)),
        ),
    )
    return ranked[0] if ranked else None


def parse_header(text: str) -> dict[str, object]:
    lengths: dict[str, int] = {}
    samples: set[str] = set()
    read_groups: set[str] = set()
    platforms: set[str] = set()
    for line in text.splitlines():
        if line.startswith("@SQ"):
            fields = dict(x.split(":", 1) for x in line.split("\t")[1:] if ":" in x)
            if fields.get("SN") and str(fields.get("LN", "")).isdigit():
                lengths[fields["SN"]] = int(fields["LN"])
        elif line.startswith("@RG"):
            fields = dict(x.split(":", 1) for x in line.split("\t")[1:] if ":" in x)
            if fields.get("SM"):
                samples.add(fields["SM"])
            if fields.get("ID"):
                read_groups.add(fields["ID"])
            if fields.get("PL"):
                platforms.add(fields["PL"])
    chr1 = lengths.get("chr1", lengths.get("1"))
    chr2 = lengths.get("chr2", lengths.get("2"))
    build = "GRCh38" if chr1 == 248956422 and chr2 == 242193529 else "NON_GRCH38_OR_UNRESOLVED"
    return {
        "genome_build": build,
        "sample_names": ",".join(sorted(samples)),
        "read_groups": ",".join(sorted(read_groups)),
        "platforms": ",".join(sorted(platforms)),
        "contigs": len(lengths),
    }


def parse_flagstat(text: str) -> dict[str, str]:
    out = {"total_reads": "", "mapped_pct": "", "duplicate_pct": "", "properly_paired_pct": ""}
    for line in text.splitlines():
        if " in total " in line:
            out["total_reads"] = line.split()[0]
        elif " mapped (" in line:
            match = re.search(r"\(([0-9.]+)%", line)
            out["mapped_pct"] = match.group(1) if match else ""
        elif " duplicates" in line:
            match = re.search(r"\(([0-9.]+)%", line)
            out["duplicate_pct"] = match.group(1) if match else ""
        elif " properly paired (" in line:
            match = re.search(r"\(([0-9.]+)%", line)
            out["properly_paired_pct"] = match.group(1) if match else ""
    return out


def mean_coverage(path: Path) -> str:
    total_bases = 0
    weighted_depth = 0.0
    try:
        with path.open(encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                length = int(row.get("endpos", 0)) - int(row.get("startpos", 0)) + 1
                total_bases += max(length, 0)
                weighted_depth += max(length, 0) * float(row.get("meandepth", 0) or 0)
    except (OSError, TypeError, ValueError):
        return ""
    return f"{weighted_depth / total_bases:.4f}" if total_bases else ""


def inspect_bam(path: Path, role: str, outdir: Path, threads: int, quick: bool) -> dict[str, object]:
    logs = outdir / "logs"
    stem = f"{role}.{path.name}"
    quickcheck = run(["samtools", "quickcheck", "-v", str(path)], logs / f"{stem}.quickcheck.log")
    header_proc = run(["samtools", "view", "-H", str(path)], logs / f"{stem}.header.sam")
    header = parse_header(header_proc.stdout)
    row: dict[str, object] = {
        "role": role,
        "bam": str(path),
        "size_bytes": path.stat().st_size,
        "index": str(index_for(path) or ""),
        "quickcheck": "PASS" if quickcheck.returncode == 0 else "FAIL",
        **header,
        "total_reads": "UNASSESSED",
        "mapped_pct": "UNASSESSED",
        "duplicate_pct": "UNASSESSED",
        "properly_paired_pct": "UNASSESSED",
        "mean_depth": "UNASSESSED",
        "contamination": "UNASSESSED",
    }
    if not quick:
        flagstat = run(["samtools", "flagstat", "-@", str(threads), str(path)], logs / f"{stem}.flagstat.txt")
        row.update(parse_flagstat(flagstat.stdout))
        coverage_path = logs / f"{stem}.coverage.tsv"
        run(["samtools", "coverage", "-o", str(coverage_path), str(path)], logs / f"{stem}.coverage.log")
        row["mean_depth"] = mean_coverage(coverage_path) or "UNASSESSED"
    return row


def estimate_contamination(path: Path, role: str, outdir: Path, sites: Path) -> str:
    gatk = shutil.which("gatk")
    if not gatk or not sites.is_file():
        return "UNASSESSED"
    metrics = outdir / "contamination"
    metrics.mkdir(parents=True, exist_ok=True)
    pileups = metrics / f"{role}.pileups.table"
    result = metrics / f"{role}.contamination.table"
    pileup_proc = run(
        [gatk, "--java-options", f"-Xmx4g -Djava.io.tmpdir={metrics}", "GetPileupSummaries", "-I", str(path), "-V", str(sites), "-L", str(sites), "-O", str(pileups)],
        outdir / f"logs/{role}.contamination.getpileup.log",
    )
    if pileup_proc.returncode:
        return "UNASSESSED"
    calc_proc = run(
        [gatk, "--java-options", "-Xmx2g", "CalculateContamination", "-I", str(pileups), "-O", str(result)],
        outdir / f"logs/{role}.contamination.calculate.log",
    )
    if calc_proc.returncode or not result.is_file():
        return "UNASSESSED"
    with result.open(encoding="utf-8") as handle:
        rows = [line.rstrip().split("\t") for line in handle if line.strip() and not line.startswith("#")]
    if len(rows) < 2 or "contamination" not in rows[0]:
        return "UNASSESSED"
    return rows[1][rows[0].index("contamination")]


def vcf_samples(path: Path) -> tuple[list[str], str]:
    if shutil.which("bcftools"):
        proc = subprocess.run(["bcftools", "query", "-l", str(path)], text=True, capture_output=True)
        if proc.returncode == 0:
            return [x for x in proc.stdout.splitlines() if x], "PASS"
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("#CHROM"):
                    return line.rstrip().split("\t")[9:], "PASS"
    except OSError:
        return [], "FAIL"
    return [], "FAIL"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tumor-bam")
    ap.add_argument("--normal-bam")
    ap.add_argument("--somatic-vcf")
    ap.add_argument("--reference", default=os.environ.get("BAM_MATCHER_REFERENCE") or os.environ.get("NEOAG_REFERENCE_FASTA"))
    ap.add_argument("--bam-matcher-loci")
    ap.add_argument("--contamination-sites", default=os.environ.get("NEOAG_CONTAMINATION_SITES"))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--quick", action="store_true", help="Skip full flagstat/coverage; mark those fields UNASSESSED")
    ap.add_argument("--skip-bam-matcher", action="store_true")
    args = ap.parse_args()

    root = Path(args.input_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    if not root.is_dir():
        raise SystemExit(f"ERROR: input directory missing: {root}")
    if not shutil.which("samtools"):
        raise SystemExit("ERROR: samtools is required")

    bams, vcfs, rna_files = discover(root)
    auto_tumor, auto_normal, pairing_method = choose_pair(bams)
    tumor = Path(args.tumor_bam).resolve() if args.tumor_bam else auto_tumor
    normal = Path(args.normal_bam).resolve() if args.normal_bam else auto_normal
    if args.tumor_bam and args.normal_bam:
        pairing_method = "EXPLICIT"
    vcf = Path(args.somatic_vcf).resolve() if args.somatic_vcf else choose_vcf(vcfs)

    inventory: list[dict[str, object]] = []
    for path in [*bams, *vcfs, *rna_files]:
        kind = "DNA_BAM" if path in bams else "SOMATIC_VCF_CANDIDATE" if path in vcfs else "RNA"
        inventory.append({"type": kind, "path": str(path), "size_bytes": path.stat().st_size, "index": str(index_for(path) or "")})
    write_tsv(outdir / "input_inventory.tsv", inventory, ["type", "path", "size_bytes", "index"])

    bam_rows: list[dict[str, object]] = []
    critical: list[str] = []
    for role, path in (("tumor", tumor), ("normal", normal)):
        if path is None or not path.is_file():
            critical.append(f"{role}_bam_missing_or_ambiguous")
            continue
        bam_rows.append(inspect_bam(path, role, outdir, args.threads, args.quick))
    write_tsv(
        outdir / "bam_qc.tsv",
        bam_rows,
        ["role", "bam", "size_bytes", "index", "quickcheck", "genome_build", "sample_names", "read_groups", "platforms", "contigs", "total_reads", "mapped_pct", "duplicate_pct", "properly_paired_pct", "mean_depth", "contamination"],
    )
    for row in bam_rows:
        if not row["index"]:
            critical.append(f"{row['role']}_bam_index_missing")
        if row["quickcheck"] != "PASS":
            critical.append(f"{row['role']}_bam_quickcheck_failed")
        if row["genome_build"] != "GRCh38":
            critical.append(f"{row['role']}_bam_not_grch38")
    if not args.quick:
        contamination_sites = Path(args.contamination_sites).resolve() if args.contamination_sites else None
        for row in bam_rows:
            role = str(row["role"])
            value = estimate_contamination(Path(str(row["bam"])), role, outdir, contamination_sites) if contamination_sites else "UNASSESSED"
            row["contamination"] = value
            if value == "UNASSESSED":
                critical.append(f"{role}_contamination_unassessed")
        write_tsv(
            outdir / "bam_qc.tsv",
            bam_rows,
            ["role", "bam", "size_bytes", "index", "quickcheck", "genome_build", "sample_names", "read_groups", "platforms", "contigs", "total_reads", "mapped_pct", "duplicate_pct", "properly_paired_pct", "mean_depth", "contamination"],
        )

    samples: list[str] = []
    vcf_status = "MISSING"
    if vcf and vcf.is_file():
        samples, vcf_status = vcf_samples(vcf)
    else:
        critical.append("somatic_vcf_missing_or_ambiguous")
    write_tsv(outdir / "vcf_qc.tsv", [{"vcf": str(vcf or ""), "status": vcf_status, "samples": ",".join(samples)}], ["vcf", "status", "samples"])
    if vcf_status != "PASS":
        critical.append("somatic_vcf_header_unreadable")

    identity_status = "UNASSESSED"
    identity_detail = "bam-matcher not requested or unavailable"
    reference = Path(args.reference).resolve() if args.reference else None
    loci = Path(args.bam_matcher_loci).resolve() if args.bam_matcher_loci else None
    if not args.skip_bam_matcher and tumor and normal and reference and loci and shutil.which("bam-matcher"):
        identity_dir = outdir / "bam_matcher"
        proc = run(
            [
                str(Path(__file__).resolve().parent / "run_bam_matcher_pair.sh"),
                "--bam1", str(normal), "--bam2", str(tumor), "--reference", str(reference),
                "--loci", str(loci), "--outdir", str(identity_dir),
            ],
            outdir / "logs/bam_matcher.driver.log",
        )
        identity_result = identity_dir / "sample_identity.tsv"
        identity_detail = str(identity_result) if identity_result.is_file() else str(outdir / "logs/bam_matcher.driver.log")
        if proc.returncode == 0:
            identity_status = "PASS"
        elif identity_result.is_file() and identity_result.stat().st_size:
            identity_status = "FAIL"
            critical.append("tumor_normal_identity_mismatch")
        else:
            identity_status = "ERROR"
            critical.append("tumor_normal_identity_tool_failed")
    elif not args.skip_bam_matcher:
        critical.append("tumor_normal_identity_unassessed")

    pairing = [{
        "tumor_bam": str(tumor or ""),
        "normal_bam": str(normal or ""),
        "selection_method": pairing_method,
        "identity_status": identity_status,
        "identity_detail": identity_detail,
    }]
    write_tsv(outdir / "sample_pairing.tsv", pairing, list(pairing[0]))

    status = "PASS" if not critical else "INCOMPLETE" if all("unassessed" in x or "ambiguous" in x for x in critical) else "FAIL"
    summary = {"status": status, "input_dir": str(root), "critical_findings": sorted(set(critical)), "quick_mode": args.quick}
    (outdir / "input_status.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    review = [
        "# Input and pairing QC",
        "",
        f"- status: **{status}**",
        f"- tumor BAM: `{tumor or 'UNRESOLVED'}`",
        f"- normal BAM: `{normal or 'UNRESOLVED'}`",
        f"- somatic VCF: `{vcf or 'UNRESOLVED'}`",
        f"- tumor-normal identity: **{identity_status}**",
        f"- full metrics: **{'no (quick mode)' if args.quick else 'yes'}**",
        "",
        "## Findings",
    ]
    review.extend(f"- {item}" for item in sorted(set(critical)))
    if not critical:
        review.append("- No blocking input or pairing finding.")
    (outdir / "input_qc_review.md").write_text("\n".join(review) + "\n", encoding="utf-8")
    if status != "FAIL":
        (outdir / ".complete").write_text("input_qc_complete\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))
    return 0 if status != "FAIL" else 4


if __name__ == "__main__":
    sys.exit(main())
