#!/usr/bin/env python3
"""Refresh variant_peptides.annotated.tsv with predictor/immunogenicity outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from neoag_v03.adapters.peptide_netmhcpan import annotate_variant_peptide_tsv


def existing(path: str | None) -> str | None:
    if not path:
        return None
    p = Path(path)
    return str(p) if p.is_file() else None


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Update variant_peptides.annotated.tsv after binding/immunogenicity scoring."
    )
    parser.add_argument("--variant-peptides", required=True, help="Input variant_peptides.tsv")
    parser.add_argument("--output", required=True, help="Output variant_peptides.annotated.tsv")
    parser.add_argument("--hla-alleles", required=True, help="Comma-separated HLA alleles")
    parser.add_argument("--netmhcpan-xls", help="NetMHCpan xls output")
    parser.add_argument("--mhcflurry-csv", help="MHCflurry csv output")
    parser.add_argument("--netmhcstabpan-tsv", help="NetMHCstabpan evidence TSV")
    parser.add_argument("--prime-tsv", help="PRIME output TSV")
    parser.add_argument("--bigmhc-im-tsv", help="BigMHC_IM output TSV")
    parser.add_argument("--iedb-immunogenicity-tsv", help="IEDB immunogenicity evidence TSV")
    args = parser.parse_args()

    hla_alleles = [a.strip() for a in args.hla_alleles.split(",") if a.strip()]
    result = annotate_variant_peptide_tsv(
        args.variant_peptides,
        hla_alleles,
        netmhcpan_xls=existing(args.netmhcpan_xls),
        mhcflurry_csv=existing(args.mhcflurry_csv),
        netmhcstabpan_tsv=existing(args.netmhcstabpan_tsv),
        prime_tsv=existing(args.prime_tsv),
        bigmhc_im_tsv=existing(args.bigmhc_im_tsv),
        iedb_immunogenicity_tsv=existing(args.iedb_immunogenicity_tsv),
        output_tsv=args.output,
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
