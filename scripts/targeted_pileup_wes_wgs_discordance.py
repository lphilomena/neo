#!/usr/bin/env python3
"""Targeted BAM pileup for WES/WGS coding call-set discordance."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import pysam


def classify_variant(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == 1 and alt.startswith(ref):
        return "INSERTION"
    if len(alt) == 1 and ref.startswith(alt):
        return "DELETION"
    if len(ref) == len(alt):
        return "MNV"
    return "COMPLEX"


def bam_contig(bam: pysam.AlignmentFile, chrom: str) -> str | None:
    candidates = [chrom, f"chr{chrom}"]
    if chrom in {"MT", "M"}:
        candidates.extend(["MT", "M", "chrM", "chrMT"])
    refs = set(bam.references)
    return next((candidate for candidate in candidates if candidate in refs), None)


def count_site(
    bam: pysam.AlignmentFile,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    min_mapq: int,
    min_baseq: int,
) -> dict:
    contig = bam_contig(bam, chrom)
    variant_type = classify_variant(ref, alt)
    result = {
        "contig_in_bam": contig or "",
        "variant_type": variant_type,
        "depth": 0,
        "ref_count": 0,
        "alt_count": 0,
        "other_count": 0,
        "alt_vaf": "",
        "count_status": "PASS",
    }
    if contig is None:
        result["count_status"] = "CONTIG_MISSING"
        return result
    if variant_type in {"MNV", "COMPLEX"}:
        result["count_status"] = "COMPLEX_REQUIRES_HAPLOTYPE_REVIEW"

    found = False
    for column in bam.pileup(
        contig,
        pos - 1,
        pos,
        truncate=True,
        stepper="all",
        min_base_quality=0,
        min_mapping_quality=0,
        max_depth=100000,
        ignore_overlaps=True,
        ignore_orphans=False,
    ):
        if column.reference_pos != pos - 1:
            continue
        found = True
        for pileup_read in column.pileups:
            read = pileup_read.alignment
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
                or read.is_duplicate
                or read.is_qcfail
                or read.mapping_quality < min_mapq
                or pileup_read.is_refskip
                or pileup_read.query_position is None
            ):
                continue
            query_pos = pileup_read.query_position
            qualities = read.query_qualities
            if qualities is not None and qualities[query_pos] < min_baseq:
                continue
            base = read.query_sequence[query_pos].upper()
            result["depth"] += 1

            if variant_type == "SNV":
                if base == alt:
                    result["alt_count"] += 1
                elif base == ref:
                    result["ref_count"] += 1
                else:
                    result["other_count"] += 1
            elif variant_type == "INSERTION":
                inserted = alt[1:]
                observed = ""
                if pileup_read.indel == len(inserted):
                    observed = read.query_sequence[query_pos + 1 : query_pos + 1 + len(inserted)].upper()
                if observed == inserted:
                    result["alt_count"] += 1
                elif pileup_read.indel == 0 and base == ref[0]:
                    result["ref_count"] += 1
                else:
                    result["other_count"] += 1
            elif variant_type == "DELETION":
                deletion_length = len(ref) - 1
                if pileup_read.indel == -deletion_length:
                    result["alt_count"] += 1
                elif pileup_read.indel == 0 and base == ref[0]:
                    result["ref_count"] += 1
                else:
                    result["other_count"] += 1
            else:
                result["other_count"] += 1
    if not found:
        result["count_status"] = "NO_PILEUP" if result["count_status"] == "PASS" else result["count_status"]
    if result["depth"]:
        result["alt_vaf"] = f"{result['alt_count'] / result['depth']:.6f}"
    return result


def evidence_class(source_status: str, wes: dict, wgs: dict, normal: dict) -> str:
    source = wes if source_status == "WES_ONLY" else wgs
    other = wgs if source_status == "WES_ONLY" else wes
    if source["variant_type"] in {"INSERTION", "DELETION", "MNV", "COMPLEX"} and source["alt_count"] == 0:
        return "SOURCE_INDEL_NOT_REPRODUCED_REASSEMBLY_REQUIRED"
    if source["variant_type"] == "SNV" and source["alt_count"] == 0:
        return "SOURCE_PASS_NOT_REPRODUCED_BY_PILEUP"
    if source["alt_count"] < 3:
        return "SOURCE_WEAK_EXACT_PILEUP_SUPPORT"
    if other["count_status"] == "CONTIG_MISSING":
        return "OTHER_CONTIG_MISSING"
    if other["depth"] < 10:
        return "OTHER_LOW_OR_NO_COVERAGE"
    normal_vaf = float(normal["alt_vaf"] or 0)
    if normal["alt_count"] >= 3 and normal_vaf >= 0.02:
        return "NORMAL_SUPPORT_REVIEW"
    other_vaf = float(other["alt_vaf"] or 0)
    if other["alt_count"] >= 3 and other_vaf >= 0.02:
        return "ALT_PRESENT_BELOW_PASS_OR_CALLER_DIFFERENCE"
    if other["count_status"] == "COMPLEX_REQUIRES_HAPLOTYPE_REVIEW":
        return "COMPLEX_REQUIRES_HAPLOTYPE_REVIEW"
    if other["depth"] >= 10 and other["alt_count"] == 0:
        source_vaf = source["alt_count"] / source["depth"] if source["depth"] else 0
        zero_probability = (1 - source_vaf) ** other["depth"]
        if zero_probability > 0.05:
            return "OTHER_COVERED_BUT_LIMITED_POWER_AT_SOURCE_VAF"
        return "COVERED_NO_ALT_SAMPLE_OR_ASSAY_DIFFERENCE"
    return "WEAK_OR_ABSENT_ALT_REVIEW"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--comparison-tsv", type=Path, required=True)
    parser.add_argument("--wes-tumor-bam", type=Path, required=True)
    parser.add_argument("--wgs-tumor-bam", type=Path, required=True)
    parser.add_argument("--normal-bam", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-baseq", type=int, default=20)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    with args.comparison_tsv.open() as handle:
        comparison_rows = list(csv.DictReader(handle, delimiter="\t"))
        input_rows = [
            row for row in comparison_rows
            if row["comparison_status"] in {"WES_ONLY", "WGS_ONLY"}
        ]

    bams = {
        "wes_tumor": pysam.AlignmentFile(str(args.wes_tumor_bam), "rb"),
        "wgs_tumor": pysam.AlignmentFile(str(args.wgs_tumor_bam), "rb"),
        "normal": pysam.AlignmentFile(str(args.normal_bam), "rb"),
    }
    long_rows = []
    wide_rows = []
    try:
        for index, row in enumerate(input_rows, 1):
            counts = {}
            for label, bam in bams.items():
                count = count_site(
                    bam,
                    row["chrom"],
                    int(row["pos"]),
                    row["ref"],
                    row["alt"],
                    args.min_mapq,
                    args.min_baseq,
                )
                counts[label] = count
                long_rows.append({
                    "comparison_status": row["comparison_status"],
                    "variant_key": row["variant_key"],
                    "genes": row["genes"],
                    "coding_class": row["coding_class"],
                    "sample": label,
                    **count,
                })
            evidence = evidence_class(row["comparison_status"], counts["wes_tumor"], counts["wgs_tumor"], counts["normal"])
            source_label = "wes_tumor" if row["comparison_status"] == "WES_ONLY" else "wgs_tumor"
            other_label = "wgs_tumor" if row["comparison_status"] == "WES_ONLY" else "wes_tumor"
            source_count = counts[source_label]
            other_count = counts[other_label]
            source_pileup_vaf = source_count["alt_count"] / source_count["depth"] if source_count["depth"] else 0
            zero_probability = (
                (1 - source_pileup_vaf) ** other_count["depth"]
                if other_count["depth"] and source_pileup_vaf else ""
            )
            wide = {
                "comparison_status": row["comparison_status"],
                "variant_key": row["variant_key"],
                "genes": row["genes"],
                "coding_class": row["coding_class"],
                "hgvsp": row["hgvsp"],
                "source_vcf_tumor_ad": row["tumor_ad"],
                "source_vcf_tumor_af": row["tumor_af"],
                "other_zero_alt_probability_at_source_pileup_vaf": (
                    f"{zero_probability:.6g}" if isinstance(zero_probability, float) else ""
                ),
                "pileup_interpretation": evidence,
            }
            for label, count in counts.items():
                for field in ("depth", "ref_count", "alt_count", "alt_vaf", "count_status"):
                    wide[f"{label}_{field}"] = count[field]
            wide_rows.append(wide)
            if index % 25 == 0:
                print(f"processed {index}/{len(input_rows)}", flush=True)
    finally:
        for bam in bams.values():
            bam.close()

    long_fields = [
        "comparison_status", "variant_key", "genes", "coding_class", "sample",
        "contig_in_bam", "variant_type", "depth", "ref_count", "alt_count",
        "other_count", "alt_vaf", "count_status",
    ]
    with (args.outdir / "discordant_targeted_pileup_long.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=long_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(long_rows)

    wide_fields = list(wide_rows[0]) if wide_rows else []
    with (args.outdir / "discordant_targeted_pileup.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=wide_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(wide_rows)

    by_variant = {row["variant_key"]: row for row in wide_rows}
    evidence_rows = []
    evidence_fields = [
        "variant_key", "chrom", "pos", "ref", "alt", "genes", "coding_class",
        "comparison_status", "cross_platform_status", "source_vcf_tumor_ad",
        "source_vcf_tumor_af", "wes_tumor_depth", "wes_tumor_alt_count",
        "wes_tumor_alt_vaf", "wgs_tumor_depth", "wgs_tumor_alt_count",
        "wgs_tumor_alt_vaf", "normal_depth", "normal_alt_count", "normal_alt_vaf",
        "other_zero_alt_probability_at_source_pileup_vaf",
    ]
    for row in comparison_rows:
        if row["comparison_status"] == "COMMON":
            evidence_rows.append({
                "variant_key": row["variant_key"],
                "chrom": row["chrom"],
                "pos": row["pos"],
                "ref": row["ref"],
                "alt": row["alt"],
                "genes": row["genes"],
                "coding_class": row["coding_class"],
                "comparison_status": "COMMON",
                "cross_platform_status": "CROSS_PLATFORM_PASS_CONCORDANT",
                "source_vcf_tumor_ad": row["tumor_ad"],
                "source_vcf_tumor_af": row["tumor_af"],
            })
            continue
        pileup = by_variant[row["variant_key"]]
        evidence_rows.append({
            "variant_key": row["variant_key"],
            "chrom": row["chrom"],
            "pos": row["pos"],
            "ref": row["ref"],
            "alt": row["alt"],
            "genes": row["genes"],
            "coding_class": row["coding_class"],
            "comparison_status": row["comparison_status"],
            "cross_platform_status": pileup["pileup_interpretation"],
            **{field: pileup.get(field, "") for field in evidence_fields if field not in {
                "variant_key", "chrom", "pos", "ref", "alt", "genes", "coding_class",
                "comparison_status", "cross_platform_status",
            }},
        })
    with (args.outdir / "cross_platform_evidence.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=evidence_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(evidence_rows)

    category_counts = dict(sorted(Counter(row["pileup_interpretation"] for row in wide_rows).items()))
    by_source = {}
    for source in ("WES_ONLY", "WGS_ONLY"):
        subset = [row for row in wide_rows if row["comparison_status"] == source]
        by_source[source] = {
            "n": len(subset),
            "categories": dict(sorted(Counter(row["pileup_interpretation"] for row in subset).items())),
        }
    summary = {
        "n_discordant": len(wide_rows),
        "min_mapq": args.min_mapq,
        "min_baseq": args.min_baseq,
        "category_counts": category_counts,
        "by_source": by_source,
        "bam_inputs": {key: str(value) for key, value in {
            "wes_tumor": args.wes_tumor_bam,
            "wgs_tumor": args.wgs_tumor_bam,
            "normal": args.normal_bam,
        }.items()},
    }
    (args.outdir / "discordant_targeted_pileup_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (args.outdir / "discordant_targeted_pileup_summary.tsv").open("w") as handle:
        handle.write("source\tcategory\tcount\n")
        for source, source_summary in by_source.items():
            for category, count in source_summary["categories"].items():
                handle.write(f"{source}\t{category}\t{count}\n")

    md = [
        "# WES/WGS discordant coding-variant targeted pileup",
        "",
        f"Discordant variants reviewed: {len(wide_rows)}",
        f"Read filters: mapping quality >= {args.min_mapq}; base quality >= {args.min_baseq}; duplicates, QC-fail, secondary and supplementary alignments excluded.",
        "",
        "| Source | Interpretation | Count |",
        "|---|---|---:|",
    ]
    for source, source_summary in by_source.items():
        for category, count in source_summary["categories"].items():
            md.append(f"| {source} | {category} | {count} |")
    md.extend([
        "",
        "## Interpretation boundary",
        "",
        "Ordinary pileup does not reproduce every local-assembly indel. SOURCE_INDEL_NOT_REPRODUCED_REASSEMBLY_REQUIRED events require Mutect2 assembly-region or IGV haplotype review.",
        "Only source-reproduced variants with adequate coverage and no ALT support in the other tumor are evidence for sample, time-point or assay discordance.",
        "A PASS-set difference is not, by itself, evidence of biological gain or loss.",
    ])
    (args.outdir / "discordant_targeted_pileup_summary.md").write_text("\n".join(md) + "\n")


if __name__ == "__main__":
    main()
