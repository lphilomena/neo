#!/usr/bin/env python3
"""Compare coding variants from matched WES and WGS VEP-annotated VCFs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path


PROTEIN_ALTERING = {
    "frameshift_variant",
    "inframe_deletion",
    "inframe_insertion",
    "missense_variant",
    "protein_altering_variant",
    "start_lost",
    "stop_gained",
    "stop_lost",
}
CODING_TERMS = PROTEIN_ALTERING | {
    "coding_sequence_variant",
    "incomplete_terminal_codon_variant",
    "synonymous_variant",
    "stop_retained_variant",
}
SPLICE_TERMS = {
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_region_variant",
}
IMPACT_ORDER = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1, "": 0}


def open_text(path: Path):
    with path.open("rb") as handle:
        is_gzip = handle.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if is_gzip else path.open()


def canonical_key(chrom: str, pos: int, ref: str, alt: str) -> tuple[str, int, str, str]:
    chrom = chrom[3:] if chrom.startswith("chr") else chrom
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return chrom, pos, ref, alt


def variant_type(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) == len(alt):
        return "MNV"
    return "InDel"


def parse_info(value: str) -> dict[str, str]:
    result = {}
    for item in value.split(";"):
        key, sep, val = item.partition("=")
        result[key] = val if sep else ""
    return result


def sample_metrics(fmt: str, sample: str) -> tuple[str, str, str]:
    values = dict(zip(fmt.split(":"), sample.split(":")))
    dp = values.get("DP", "")
    ad = values.get("AD", "")
    af = values.get("AF", "")
    if not af and ad and dp not in {"", "0", "."}:
        try:
            alt_depth = float(ad.split(",")[1])
            af = f"{alt_depth / float(dp):.6g}"
        except (IndexError, ValueError, ZeroDivisionError):
            pass
    return dp, ad, af


def split_terms(value: str) -> set[str]:
    return {term for term in value.split("&") if term}


def parse_vep_vcf(path: Path, source: str, tumor_sample: str) -> tuple[dict, dict]:
    csq_fields: list[str] = []
    samples: list[str] = []
    variants: dict[tuple[str, int, str, str], dict] = {}

    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                match = re.search(r"Format: ([^\"]+)", line)
                if match:
                    csq_fields = match.group(1).split("|")
            elif line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            elif line.startswith("#"):
                continue
            else:
                fields = line.rstrip("\n").split("\t")
                chrom, pos_s, _vid, ref, alt_text, _qual, filt, info_text = fields[:8]
                if filt not in {"PASS", "."}:
                    continue
                if tumor_sample not in samples:
                    raise ValueError(f"Tumor sample {tumor_sample!r} not found in {path}; samples={samples}")
                sample = fields[9 + samples.index(tumor_sample)]
                fmt = fields[8]
                dp, ad, af = sample_metrics(fmt, sample)
                info = parse_info(info_text)
                annotations = []
                for raw_csq in info.get("CSQ", "").split(","):
                    if not raw_csq:
                        continue
                    values = raw_csq.split("|")
                    annotations.append(dict(zip(csq_fields, values)))

                alts = alt_text.split(",")
                for alt in alts:
                    key = canonical_key(chrom, int(pos_s), ref, alt)
                    relevant = []
                    for ann in annotations:
                        terms = split_terms(ann.get("Consequence", ""))
                        if ann.get("BIOTYPE") != "protein_coding":
                            continue
                        if not (ann.get("CDS_position") or terms & (CODING_TERMS | SPLICE_TERMS)):
                            continue
                        relevant.append(ann)
                    if not relevant:
                        continue

                    all_terms = set().union(*(split_terms(a.get("Consequence", "")) for a in relevant))
                    impacts = {a.get("IMPACT", "") for a in relevant}
                    top_impact = max(impacts, key=lambda x: IMPACT_ORDER.get(x, 0), default="")
                    if all_terms & PROTEIN_ALTERING:
                        coding_class = "PROTEIN_ALTERING"
                    elif all_terms & {"splice_acceptor_variant", "splice_donor_variant"}:
                        coding_class = "ESSENTIAL_SPLICE"
                    elif "synonymous_variant" in all_terms:
                        coding_class = "SYNONYMOUS"
                    else:
                        coding_class = "OTHER_CODING_OR_SPLICE"

                    row = {
                        "variant_key": f"{key[0]}:{key[1]}:{key[2]}:{key[3]}",
                        "chrom": key[0],
                        "pos": key[1],
                        "ref": key[2],
                        "alt": key[3],
                        "variant_type": variant_type(key[2], key[3]),
                        "source": source,
                        "filter": filt,
                        "tumor_dp": dp,
                        "tumor_ad": ad,
                        "tumor_af": af,
                        "genes": ";".join(sorted({a.get("SYMBOL", "") for a in relevant if a.get("SYMBOL")})),
                        "gene_ids": ";".join(sorted({a.get("Gene", "") for a in relevant if a.get("Gene")})),
                        "transcripts": ";".join(sorted({a.get("Feature", "") for a in relevant if a.get("Feature")})),
                        "consequences": ";".join(sorted(all_terms)),
                        "highest_impact": top_impact,
                        "coding_class": coding_class,
                        "hgvsc": ";".join(sorted({a.get("HGVSc", "") for a in relevant if a.get("HGVSc")})),
                        "hgvsp": ";".join(sorted({a.get("HGVSp", "") for a in relevant if a.get("HGVSp")})),
                    }
                    variants[key] = row

    metadata = {"path": str(path), "tumor_sample": tumor_sample, "vcf_samples": samples}
    return variants, metadata


FIELDS = [
    "variant_key", "chrom", "pos", "ref", "alt", "variant_type", "source", "filter", "tumor_dp",
    "tumor_ad", "tumor_af", "genes", "gene_ids", "transcripts", "consequences",
    "highest_impact", "coding_class", "hgvsc", "hgvsp",
]


def write_tsv(path: Path, rows: list[dict], fields: list[str] = FIELDS) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--wes-vep-vcf", type=Path, required=True)
    parser.add_argument("--wgs-vep-vcf", type=Path, required=True)
    parser.add_argument("--wes-tumor-sample", required=True)
    parser.add_argument("--wgs-tumor-sample", required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    wes, wes_meta = parse_vep_vcf(args.wes_vep_vcf, "WES", args.wes_tumor_sample)
    wgs, wgs_meta = parse_vep_vcf(args.wgs_vep_vcf, "WGS", args.wgs_tumor_sample)
    common = set(wes) & set(wgs)
    wes_only = set(wes) - set(wgs)
    wgs_only = set(wgs) - set(wes)
    wes_protein = {
        key for key, row in wes.items()
        if row["coding_class"] == "PROTEIN_ALTERING" and row["variant_type"] in {"SNV", "InDel"}
    }
    wgs_protein = {
        key for key, row in wgs.items()
        if row["coding_class"] == "PROTEIN_ALTERING" and row["variant_type"] in {"SNV", "InDel"}
    }
    wes_protein_mnv = {
        key for key, row in wes.items()
        if row["coding_class"] == "PROTEIN_ALTERING" and row["variant_type"] == "MNV"
    }
    wgs_protein_mnv = {
        key for key, row in wgs.items()
        if row["coding_class"] == "PROTEIN_ALTERING" and row["variant_type"] == "MNV"
    }
    protein_common = wes_protein & wgs_protein
    protein_wes_only = wes_protein - wgs_protein
    protein_wgs_only = wgs_protein - wes_protein

    write_tsv(args.outdir / "wes_coding_variants.tsv", [wes[k] for k in sorted(wes)])
    write_tsv(args.outdir / "wgs_coding_variants.tsv", [wgs[k] for k in sorted(wgs)])
    write_tsv(args.outdir / "common_coding_variants.tsv", [wes[k] for k in sorted(common)])
    write_tsv(args.outdir / "wes_only_coding_variants.tsv", [wes[k] for k in sorted(wes_only)])
    write_tsv(args.outdir / "wgs_only_coding_variants.tsv", [wgs[k] for k in sorted(wgs_only)])
    write_tsv(args.outdir / "wes_protein_altering_snv_indel.tsv", [wes[k] for k in sorted(wes_protein)])
    write_tsv(args.outdir / "wgs_protein_altering_snv_indel.tsv", [wgs[k] for k in sorted(wgs_protein)])
    write_tsv(args.outdir / "common_protein_altering_snv_indel.tsv", [wes[k] for k in sorted(protein_common)])
    write_tsv(args.outdir / "wes_only_protein_altering_snv_indel.tsv", [wes[k] for k in sorted(protein_wes_only)])
    write_tsv(args.outdir / "wgs_only_protein_altering_snv_indel.tsv", [wgs[k] for k in sorted(protein_wgs_only)])

    comparison_fields = ["comparison_status"] + FIELDS
    comparison = []
    for key in sorted(set(wes) | set(wgs)):
        if key in common:
            status = "COMMON"
            base = dict(wes[key])
            base["source"] = "WES;WGS"
            base["tumor_dp"] = f"WES={wes[key]['tumor_dp']};WGS={wgs[key]['tumor_dp']}"
            base["tumor_ad"] = f"WES={wes[key]['tumor_ad']};WGS={wgs[key]['tumor_ad']}"
            base["tumor_af"] = f"WES={wes[key]['tumor_af']};WGS={wgs[key]['tumor_af']}"
        elif key in wes:
            status, base = "WES_ONLY", dict(wes[key])
        else:
            status, base = "WGS_ONLY", dict(wgs[key])
        base["comparison_status"] = status
        comparison.append(base)
    write_tsv(args.outdir / "wes_wgs_coding_variant_comparison.tsv", comparison, comparison_fields)
    protein_union = wes_protein | wgs_protein
    protein_variant_keys = {
        f"{chrom}:{pos}:{ref}:{alt}" for chrom, pos, ref, alt in protein_union
    }
    write_tsv(
        args.outdir / "protein_altering_snv_indel_comparison.tsv",
        [row for row in comparison if row["variant_key"] in protein_variant_keys],
        comparison_fields,
    )

    def class_counts(rows: dict) -> dict:
        return dict(sorted(Counter(row["coding_class"] for row in rows.values()).items()))

    af_pairs = []
    for key in common:
        try:
            af_pairs.append((float(wes[key]["tumor_af"].split(",")[0]), float(wgs[key]["tumor_af"].split(",")[0])))
        except (ValueError, IndexError):
            pass
    af_correlation = None
    if len(af_pairs) > 1:
        wes_mean = sum(a for a, _ in af_pairs) / len(af_pairs)
        wgs_mean = sum(b for _, b in af_pairs) / len(af_pairs)
        numerator = sum((a - wes_mean) * (b - wgs_mean) for a, b in af_pairs)
        denominator = (
            sum((a - wes_mean) ** 2 for a, _ in af_pairs)
            * sum((b - wgs_mean) ** 2 for _, b in af_pairs)
        ) ** 0.5
        af_correlation = round(numerator / denominator, 6) if denominator else None
    else:
        wes_mean = wgs_mean = None

    summary = {
        "wes": {"coding_variants": len(wes), "class_counts": class_counts(wes), **wes_meta},
        "wgs": {"coding_variants": len(wgs), "class_counts": class_counts(wgs), **wgs_meta},
        "common": len(common),
        "common_class_counts": class_counts({key: wes[key] for key in common}),
        "wes_only": len(wes_only),
        "wes_only_class_counts": class_counts({key: wes[key] for key in wes_only}),
        "wgs_only": len(wgs_only),
        "wgs_only_class_counts": class_counts({key: wgs[key] for key in wgs_only}),
        "protein_altering_snv_indel": {
            "wes": len(wes_protein),
            "wgs": len(wgs_protein),
            "common": len(protein_common),
            "wes_only": len(protein_wes_only),
            "wgs_only": len(protein_wgs_only),
            "definition": "SNV/InDel with a VEP protein-coding consequence in PROTEIN_ALTERING",
        },
        "protein_altering_mnv_excluded": {
            "wes": len(wes_protein_mnv),
            "wgs": len(wgs_protein_mnv),
            "definition": "Equal-length multi-nucleotide substitutions are MNVs, not InDels",
        },
        "wes_recall_of_wgs_coding_calls": round(len(common) / len(wgs), 6) if wgs else None,
        "wgs_recall_of_wes_coding_calls": round(len(common) / len(wes), 6) if wes else None,
        "common_af_pairs": len(af_pairs),
        "common_af_pearson": af_correlation,
        "common_wes_mean_af": round(wes_mean, 6) if wes_mean is not None else None,
        "common_wgs_mean_af": round(wgs_mean, 6) if wgs_mean is not None else None,
        "interpretation_boundary": (
            "Call-set concordance only. A platform-specific call is not proof that the other platform "
            "is reference at that locus; callable coverage and targeted pileup are required."
        ),
    }
    (args.outdir / "wes_wgs_coding_comparison.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n"
    )
    with (args.outdir / "wes_wgs_coding_summary.tsv").open("w") as handle:
        handle.write("metric\tvalue\n")
        for key in (
            "common", "wes_only", "wgs_only", "wes_recall_of_wgs_coding_calls",
            "wgs_recall_of_wes_coding_calls", "common_af_pairs", "common_af_pearson",
            "common_wes_mean_af", "common_wgs_mean_af",
        ):
            handle.write(f"{key}\t{summary[key]}\n")
        for source in ("wes", "wgs"):
            handle.write(f"{source}_coding_variants\t{summary[source]['coding_variants']}\n")
            for klass, count in summary[source]["class_counts"].items():
                handle.write(f"{source}_{klass.lower()}\t{count}\n")
        for key, value in summary["protein_altering_snv_indel"].items():
            if key != "definition":
                handle.write(f"protein_altering_{key}\t{value}\n")
        for source, value in summary["protein_altering_mnv_excluded"].items():
            if source != "definition":
                handle.write(f"protein_altering_mnv_excluded_{source}\t{value}\n")

    md = [
        "# WES/WGS coding variant comparison",
        "",
        "| Metric | Count |",
        "|---|---:|",
        f"| WES coding/splice variants | {len(wes)} |",
        f"| WGS coding/splice variants | {len(wgs)} |",
        f"| Common exact normalized variants | {len(common)} |",
        f"| WES only | {len(wes_only)} |",
        f"| WGS only | {len(wgs_only)} |",
        f"| Common protein-altering variants | {summary['common_class_counts'].get('PROTEIN_ALTERING', 0)} |",
        "",
        "## Protein-altering SNV/InDel subset used for biological comparison",
        "",
        "| Metric | Count |",
        "|---|---:|",
        f"| WES protein-altering SNV/InDel | {len(wes_protein)} |",
        f"| WGS protein-altering SNV/InDel | {len(wgs_protein)} |",
        f"| Common protein-altering SNV/InDel | {len(protein_common)} |",
        f"| WES-only protein-altering SNV/InDel | {len(protein_wes_only)} |",
        f"| WGS-only protein-altering SNV/InDel | {len(protein_wgs_only)} |",
        f"| Protein-altering MNV excluded from SNV/InDel subset, WES / WGS | {len(wes_protein_mnv)} / {len(wgs_protein_mnv)} |",
        f"| Common-call VAF Pearson correlation | {af_correlation if af_correlation is not None else 'NA'} |",
        f"| Common-call mean VAF, WES / WGS | {summary['common_wes_mean_af']} / {summary['common_wgs_mean_af']} |",
        "",
        "Both inputs are existing VEP-annotated PASS VCFs. The comparison uses normalized chromosome/position/ref/alt keys and protein-coding transcript annotations.",
        "",
        "## Interpretation boundary",
        "",
        "A platform-specific call means it is absent from the other PASS call set. It does not prove a reference genotype there. Review callable coverage and tumor/normal pileups before assigning technical or biological discordance.",
    ]
    (args.outdir / "wes_wgs_coding_summary.md").write_text("\n".join(md) + "\n")


if __name__ == "__main__":
    main()
