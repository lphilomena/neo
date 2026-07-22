#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

from neoag.haplotype import parse_variant, phase_bam_region
from neoag.utils import write_tsv


def main() -> None:
    parser = argparse.ArgumentParser(description="Read-backed phasing for nearby variants")
    parser.add_argument("--bam", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--variant", action="append", required=True, help="POS:REF>ALT; repeat for each variant")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    variants = [parse_variant(token) for token in args.variant]
    reads, summary = phase_bam_region(
        args.bam, args.chrom, variants, samtools=args.samtools,
        min_mapq=args.min_mapq, min_baseq=args.min_baseq,
    )
    write_tsv(outdir / "read_haplotypes.tsv", reads, ["read_name", "haplotype", "n_alignments", "max_mapq", "flags"])
    count_rows = [{"haplotype": hap, "fragment_count": n} for hap, n in summary["haplotype_counts"].items()]
    write_tsv(outdir / "haplotype_counts.tsv", count_rows, ["haplotype", "fragment_count"])
    (outdir / "phasing_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
