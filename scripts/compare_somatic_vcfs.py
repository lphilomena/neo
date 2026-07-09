#!/usr/bin/env python3
"""Compare two somatic pass VCFs; export overlap/unique TSV + gene/VAF summaries."""

from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

VAF_BINS = [
    ("0.00-0.05", 0.0, 0.05),
    ("0.05-0.10", 0.05, 0.10),
    ("0.10-0.20", 0.10, 0.20),
    ("0.20-0.50", 0.20, 0.50),
    (">0.50", 0.50, 1.01),
]

OUT_FIELDS = [
    "set",
    "chrom",
    "pos",
    "ref",
    "alt",
    "variant_key",
    "variant_type",
    "gene",
    "consequence",
    "hgvsp",
    "transcript_id",
    "filter_946",
    "filter_383",
    "tumor_sample_946",
    "tumor_sample_383",
    "normal_sample",
    "tumor_vaf_946",
    "tumor_vaf_383",
    "tumor_dp_946",
    "tumor_dp_383",
    "normal_vaf_946",
    "normal_vaf_383",
    "vaf_delta_tumor",
]


def open_maybe_gz(path: Path):
    if path.read_bytes()[:2] == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open(encoding="utf-8", errors="replace")


def variant_type(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    return "indel"


def parse_format_value(fmt: str, sample: str, key: str) -> str:
    keys = fmt.split(":")
    if key not in keys:
        return ""
    idx = keys.index(key)
    parts = sample.split(":")
    if idx >= len(parts):
        return ""
    return parts[idx]


def parse_vaf(fmt: str, sample: str) -> tuple[str, str]:
    af = parse_format_value(fmt, sample, "AF")
    if af:
        vals = [v for v in af.split(",") if v not in {"", "."}]
        if vals:
            return vals[0], ""
    ad = parse_format_value(fmt, sample, "AD")
    if not ad or "," not in ad:
        return "", parse_format_value(fmt, sample, "DP")
    ref_ad, alt_ad, *rest = ad.split(",")
    try:
        ref_n, alt_n = int(ref_ad), int(alt_ad)
        total = ref_n + alt_n
        if total > 0:
            return f"{alt_n / total:.4f}", str(total)
    except ValueError:
        pass
    return "", parse_format_value(fmt, sample, "DP")


def parse_csq_gene(csq_field: str) -> dict[str, str]:
    if not csq_field:
        return {}
    parts = csq_field.split("|")
    # VEP CSQ: Allele|Consequence|...|SYMBOL|...|Feature|...|HGVSp at ~10
    gene = parts[3] if len(parts) > 3 else ""
    transcript = parts[6] if len(parts) > 6 else ""
    hgvsp = parts[10] if len(parts) > 10 else ""
    consequence = parts[1] if len(parts) > 1 else ""
    if not gene and len(parts) > 4:
        gene = parts[4]
    return {
        "gene": gene,
        "consequence": consequence.split(",")[0] if consequence else "",
        "hgvsp": hgvsp,
        "transcript_id": transcript,
    }


def load_vep_gene_index(path: Path) -> dict[tuple[str, str, str, str], dict[str, str]]:
    index: dict[tuple[str, str, str, str], dict[str, str]] = {}
    if not path.is_file():
        return index
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _id, ref, alt, _qual, _filt, info = parts[:8]
            csq = ""
            for item in info.split(";"):
                if item.startswith("CSQ="):
                    csq = item[4:]
                    break
            if not csq:
                continue
            key = (chrom, pos, ref, alt)
            if key in index:
                continue
            for entry in csq.split(","):
                ann = parse_csq_gene(entry)
                if ann.get("gene"):
                    index[key] = ann
                    break
    return index


def load_somatic_vcf(
    path: Path,
    *,
    tumor_sample: str,
    normal_sample: str,
    label: str,
) -> dict[tuple[str, str, str, str], dict[str, str]]:
    rows: dict[tuple[str, str, str, str], dict[str, str]] = {}
    with open_maybe_gz(path) as fh:
        header: list[str] = []
        header_fmt = ""
        for line in fh:
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                header_fmt = header[8] if len(header) > 8 else ""
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos, _vid, ref, alt = parts[:5]
            filt = parts[6]
            key = (chrom, pos, ref, alt)
            fmt = parts[8] if header_fmt == "FORMAT" else header_fmt
            sample_cols = {header[i]: parts[i] for i in range(9, len(header))}
            tumor = sample_cols.get(tumor_sample, "")
            normal = sample_cols.get(normal_sample, "")
            tvaf, tdp = parse_vaf(fmt, tumor)
            nvaf, _ndp = parse_vaf(fmt, normal)
            rows[key] = {
                f"filter_{label}": filt,
                f"tumor_vaf_{label}": tvaf,
                f"tumor_dp_{label}": tdp,
                f"normal_vaf_{label}": nvaf,
            }
    return rows


def vaf_bin(vaf: str) -> str:
    if not vaf:
        return "missing"
    try:
        x = float(vaf)
    except ValueError:
        return "missing"
    for name, lo, hi in VAF_BINS:
        if lo <= x < hi:
            return name
    return "missing"


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf-946", required=True)
    ap.add_argument("--vcf-383", required=True)
    ap.add_argument("--vep-946")
    ap.add_argument("--vep-383")
    ap.add_argument("--tumor-946", default="ML150006946_L01_137")
    ap.add_argument("--tumor-383", default="M1ML150017383_L01_438")
    ap.add_argument("--normal", default="ML150006927_L01_470")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    path946 = Path(args.vcf_946)
    path383 = Path(args.vcf_383)
    outdir = Path(args.outdir)
    vep946 = Path(args.vep_946) if args.vep_946 else None
    vep383 = Path(args.vep_383) if args.vep_383 else None

    v946 = load_somatic_vcf(path946, tumor_sample=args.tumor_946, normal_sample=args.normal, label="946")
    v383 = load_somatic_vcf(path383, tumor_sample=args.tumor_383, normal_sample=args.normal, label="383")
    gene946 = load_vep_gene_index(vep946) if vep946 else {}
    gene383 = load_vep_gene_index(vep383) if vep383 else {}
    gene_index = dict(gene946)
    for k, v in gene383.items():
        gene_index.setdefault(k, v)

    keys946 = set(v946)
    keys383 = set(v383)
    overlap = keys946 & keys383
    only946 = keys946 - keys383
    only383 = keys383 - keys946

    def build_row(key: tuple[str, str, str, str], set_name: str) -> dict[str, str]:
        chrom, pos, ref, alt = key
        ann = gene_index.get(key) or {}
        r946 = v946.get(key, {})
        r383 = v383.get(key, {})
        tv946 = r946.get("tumor_vaf_946", "")
        tv383 = r383.get("tumor_vaf_383", "")
        delta = ""
        if tv946 and tv383:
            try:
                delta = f"{float(tv383) - float(tv946):.4f}"
            except ValueError:
                delta = ""
        return {
            "set": set_name,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "variant_key": f"{chrom}:{pos}:{ref}>{alt}",
            "variant_type": variant_type(ref, alt),
            "gene": ann.get("gene", ""),
            "consequence": ann.get("consequence", ""),
            "hgvsp": ann.get("hgvsp", ""),
            "transcript_id": ann.get("transcript_id", ""),
            "filter_946": r946.get("filter_946", ""),
            "filter_383": r383.get("filter_383", ""),
            "tumor_sample_946": args.tumor_946,
            "tumor_sample_383": args.tumor_383,
            "normal_sample": args.normal,
            "tumor_vaf_946": tv946,
            "tumor_vaf_383": tv383,
            "tumor_dp_946": r946.get("tumor_dp_946", ""),
            "tumor_dp_383": r383.get("tumor_dp_383", ""),
            "normal_vaf_946": r946.get("normal_vaf_946", ""),
            "normal_vaf_383": r383.get("normal_vaf_383", ""),
            "vaf_delta_tumor": delta,
        }

    all_rows: list[dict[str, str]] = []
    for key in sorted(overlap):
        all_rows.append(build_row(key, "overlap"))
    for key in sorted(only946):
        all_rows.append(build_row(key, "only_ML150006946"))
    for key in sorted(only383):
        all_rows.append(build_row(key, "only_M1ML150017383"))

    write_tsv(outdir / "somatic_variants_compare.all.tsv", all_rows, OUT_FIELDS)
    write_tsv(outdir / "somatic_variants_compare.overlap.tsv", [r for r in all_rows if r["set"] == "overlap"], OUT_FIELDS)
    write_tsv(outdir / "somatic_variants_compare.only_ML150006946.tsv", [r for r in all_rows if r["set"] == "only_ML150006946"], OUT_FIELDS)
    write_tsv(outdir / "somatic_variants_compare.only_M1ML150017383.tsv", [r for r in all_rows if r["set"] == "only_M1ML150017383"], OUT_FIELDS)

    # VAF stratification
    vaf_summary: list[dict[str, str]] = []
    for set_name, keys, vaf_key in [
        ("overlap", overlap, None),
        ("only_ML150006946", only946, "tumor_vaf_946"),
        ("only_M1ML150017383", only383, "tumor_vaf_383"),
    ]:
        type_ctr: Counter[str] = Counter()
        bin_ctr: Counter[str] = Counter()
        for key in keys:
            ref, alt = key[2], key[3]
            type_ctr[variant_type(ref, alt)] += 1
            if set_name == "overlap":
                vaf = v946.get(key, {}).get("tumor_vaf_946") or v383.get(key, {}).get("tumor_vaf_383", "")
            elif vaf_key == "tumor_vaf_946":
                vaf = v946.get(key, {}).get("tumor_vaf_946", "")
            else:
                vaf = v383.get(key, {}).get("tumor_vaf_383", "")
            bin_ctr[vaf_bin(vaf)] += 1
        for vt, n in sorted(type_ctr.items()):
            vaf_summary.append({"set": set_name, "stratum_type": "variant_type", "stratum": vt, "count": str(n)})
        for b, n in sorted(bin_ctr.items()):
            vaf_summary.append({"set": set_name, "stratum_type": "tumor_vaf_bin", "stratum": b, "count": str(n)})

    write_tsv(outdir / "somatic_variants_compare.vaf_type_summary.tsv", vaf_summary, ["set", "stratum_type", "stratum", "count"])

    # Gene stratification (annotated only)
    gene_rows: list[dict[str, str]] = []
    gene_ctr: dict[str, Counter[str]] = defaultdict(Counter)
    for row in all_rows:
        gene = row.get("gene") or "UNKNOWN"
        gene_ctr[row["set"]][gene] += 1
    for set_name in ("overlap", "only_ML150006946", "only_M1ML150017383"):
        for gene, n in gene_ctr[set_name].most_common():
            gene_rows.append({"set": set_name, "gene": gene, "count": str(n)})
    write_tsv(outdir / "somatic_variants_compare.gene_summary.tsv", gene_rows, ["set", "gene", "count"])

    # overlap VAF delta summary
    delta_rows: list[dict[str, str]] = []
    for row in all_rows:
        if row["set"] != "overlap" or not row["vaf_delta_tumor"]:
            continue
        delta_rows.append({
            "variant_key": row["variant_key"],
            "gene": row["gene"],
            "tumor_vaf_946": row["tumor_vaf_946"],
            "tumor_vaf_383": row["tumor_vaf_383"],
            "vaf_delta_tumor": row["vaf_delta_tumor"],
            "consequence": row["consequence"],
        })
    write_tsv(outdir / "somatic_variants_compare.overlap_vaf_delta.tsv", delta_rows, list(delta_rows[0].keys()) if delta_rows else ["variant_key"])

    counts = {
        "ML150006946_total": len(keys946),
        "M1ML150017383_total": len(keys383),
        "overlap": len(overlap),
        "only_ML150006946": len(only946),
        "only_M1ML150017383": len(only383),
        "gene_annotated_overlap": sum(1 for k in overlap if gene_index.get(k)),
    }
    summary_path = outdir / "somatic_variants_compare.counts.tsv"
    write_tsv(summary_path, [{"metric": k, "value": str(v)} for k, v in counts.items()], ["metric", "value"])

    print("Wrote:", outdir)
    for k, v in counts.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
