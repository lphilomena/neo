#!/usr/bin/env python3
"""Normalize RNA splice sources with v0.4.4 canonical junction provenance."""

from __future__ import annotations

import argparse
from pathlib import Path

from neoag.splice.normalization import normalize_splice_sources


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Canonicalize RegTools/SNAF/SpliceMutr junctions, prevent cross-junction "
            "read leakage, preserve every source record, and emit raw event/peptide tables."
        )
    )
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--profile", default="default")
    primary = parser.add_mutually_exclusive_group(required=True)
    primary.add_argument("--junctions", type=Path, help="RegTools BED/TSV junction evidence")
    primary.add_argument("--star-sj", type=Path, help="STAR SJ.out.tab from the matched sample")
    parser.add_argument("--snaf", type=Path)
    parser.add_argument("--splicemutr", type=Path)
    parser.add_argument("--normal-junctions", type=Path)
    parser.add_argument(
        "--annotation-gtf", type=Path,
        help="GTF used to resolve strand only from exact same-transcript exon boundaries",
    )
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--genome-build", default="GRCh38")
    parser.add_argument("--junction-coordinate-system", default="auto")
    parser.add_argument("--snaf-coordinate-system", default="auto")
    parser.add_argument("--splicemutr-coordinate-system", default="auto")
    parser.add_argument("--normal-coordinate-system", default="auto")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="fail when a primary junction row cannot be normalized",
    )
    parser.add_argument(
        "--candidate-only",
        action="store_true",
        help="emit only neoantigen candidates and their exact linked primary junction evidence",
    )
    args = parser.parse_args()

    junctions = args.star_sj or args.junctions
    junction_tool = "STAR" if args.star_sj else "RegTools"
    coordinate_system = "star_sj" if args.star_sj else args.junction_coordinate_system

    normalize_splice_sources(
        sample_id=args.sample_id,
        profile_name=args.profile,
        junctions=junctions,
        junction_tool=junction_tool,
        snaf=args.snaf,
        splicemutr=args.splicemutr,
        normal_junctions=args.normal_junctions,
        annotation_gtf=args.annotation_gtf,
        outdir=args.outdir,
        genome_build=args.genome_build,
        junction_coordinate_system=coordinate_system,
        snaf_coordinate_system=args.snaf_coordinate_system,
        splicemutr_coordinate_system=args.splicemutr_coordinate_system,
        normal_coordinate_system=args.normal_coordinate_system,
        strict=args.strict,
        candidate_only=args.candidate_only,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
