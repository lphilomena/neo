#!/usr/bin/env python3
"""Count RNA ref/alt support at somatic VCF sites (pysam pileup)."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--somatic-vcf", type=Path, required=True)
    p.add_argument("--rna-bam", type=Path, required=True)
    p.add_argument("--output-tsv", type=Path, required=True)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--min-baseq", type=int, default=13)
    return p.parse_args()


def iter_somatic_sites(vcf: Path):
    with open_text(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, pos, _vid, ref, alt = parts[:5]
            if not alt or alt == ".":
                continue
            alt_allele = alt.split(",")[0]
            yield chrom, int(pos), ref, alt_allele


def normalize_chrom(chrom: str, references: set[str]) -> str:
    if chrom in references:
        return chrom
    if chrom.startswith("chr") and chrom[3:] in references:
        return chrom[3:]
    if f"chr{chrom}" in references:
        return f"chr{chrom}"
    return chrom


def variant_type(ref: str, alt: str) -> str:
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(alt) > len(ref) and alt.startswith(ref):
        return "INS"
    if len(ref) > len(alt) and ref.startswith(alt):
        return "DEL"
    if len(ref) == len(alt):
        return "MNV"
    return "COMPLEX"


def classify_observation(*, query_sequence: str, query_position: int, indel: int, ref: str, alt: str) -> str:
    """Classify one pileup observation at the VCF anchor as ref, alt or other."""
    kind = variant_type(ref, alt)
    ref_u = ref.upper()
    alt_u = alt.upper()
    if kind == "SNV":
        base = query_sequence[query_position].upper()
        return "ref" if base == ref_u else ("alt" if base == alt_u else "other")
    if kind == "INS":
        expected = len(alt_u) - len(ref_u)
        if indel == expected:
            inserted = query_sequence[query_position + 1 : query_position + 1 + expected].upper()
            return "alt" if inserted == alt_u[len(ref_u) :] else "other"
        return "ref" if indel == 0 and query_sequence[query_position].upper() == ref_u[0] else "other"
    if kind == "DEL":
        expected = -(len(ref_u) - len(alt_u))
        if indel == expected:
            return "alt"
        return "ref" if indel == 0 and query_sequence[query_position].upper() == ref_u[0] else "other"
    if kind == "MNV" and indel == 0:
        observed = query_sequence[query_position : query_position + len(ref_u)].upper()
        return "ref" if observed == ref_u else ("alt" if observed == alt_u else "other")
    return "other"


def count_alleles(
    bam,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    *,
    references: set[str],
    min_mapq: int,
    min_baseq: int,
) -> tuple[int, int, int]:
  import pysam

  chrom = normalize_chrom(chrom, references)
  ref_count = 0
  alt_count = 0
  try:
      for col in bam.pileup(
          chrom,
          pos - 1,
          pos,
          truncate=True,
          max_depth=50_000,
          min_base_quality=min_baseq,
          min_mapping_quality=min_mapq,
      ):
          if col.reference_pos != pos - 1:
              continue
          for pr in col.pileups:
              if pr.is_refskip or pr.is_del or pr.query_position is None:
                  continue
              observation = classify_observation(
                  query_sequence=pr.alignment.query_sequence,
                  query_position=pr.query_position,
                  indel=pr.indel,
                  ref=ref,
                  alt=alt,
              )
              if observation == "ref":
                  ref_count += 1
              elif observation == "alt":
                  alt_count += 1
  except (ValueError, OSError):
      return 0, 0, 0
  depth = ref_count + alt_count
  return ref_count, alt_count, depth


def main() -> None:
    args = parse_args()
    import pysam

    bam = pysam.AlignmentFile(str(args.rna_bam), "rb")
    references = set(bam.references)
    rows: list[dict[str, str]] = []
    seen = 0
    with_support = 0
    for chrom, pos, ref, alt in iter_somatic_sites(args.somatic_vcf):
        seen += 1
        ref_n, alt_n, depth = count_alleles(
            bam,
            chrom,
            pos,
            ref,
            alt,
            references=references,
            min_mapq=args.min_mapq,
            min_baseq=args.min_baseq,
        )
        vaf = (alt_n / depth) if depth > 0 else 0.0
        if alt_n > 0:
            with_support += 1
        rows.append(
            {
                "chrom": chrom,
                "pos": str(pos),
                "ref": ref,
                "alt": alt,
                "rna_ref_reads": str(ref_n),
                "rna_alt_reads": str(alt_n),
                "rna_depth": str(depth),
                "rna_vaf": f"{vaf:.4f}" if depth > 0 else "",
                "variant_key": f"{chrom}:{pos}{ref}>{alt}",
                "variant_type": variant_type(ref, alt),
                "count_status": "ASSESSED" if variant_type(ref, alt) != "COMPLEX" else "UNASSESSED_COMPLEX",
            }
        )
        if seen % 500 == 0:
            print(f"... {seen} sites", flush=True)
    bam.close()

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys()) if rows else []
    with args.output_tsv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {args.output_tsv} sites={seen} with_alt_support={with_support}")


if __name__ == "__main__":
    main()
