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
    parser.add_argument("--junctions", required=True, type=Path)
    parser.add_argument("--snaf", type=Path)
    parser.add_argument("--splicemutr", type=Path)
    parser.add_argument("--normal-junctions", type=Path)
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
    args = parser.parse_args()

    normalize_splice_sources(
        sample_id=args.sample_id,
        profile_name=args.profile,
        junctions=args.junctions,
        snaf=args.snaf,
        splicemutr=args.splicemutr,
        normal_junctions=args.normal_junctions,
        outdir=args.outdir,
        genome_build=args.genome_build,
        junction_coordinate_system=args.junction_coordinate_system,
        snaf_coordinate_system=args.snaf_coordinate_system,
        splicemutr_coordinate_system=args.splicemutr_coordinate_system,
        normal_coordinate_system=args.normal_coordinate_system,
        strict=args.strict,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
