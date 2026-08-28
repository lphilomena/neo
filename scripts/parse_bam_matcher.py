#!/usr/bin/env python3
from __future__ import annotations

import argparse

from neoag.sample_identity.bam_matcher import parse_bam_matcher_short, write_identity_tsv


def main() -> int:
    parser = argparse.ArgumentParser(description="Normalize BAM-matcher short output")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--fail-on-mismatch", action="store_true")
    args = parser.parse_args()
    parsed = parse_bam_matcher_short(args.input)
    write_identity_tsv(parsed, args.output)
    if args.fail_on_mismatch and parsed["sample_identity_status"] == "MISMATCH":
        return 42
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
