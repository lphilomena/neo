#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

from neoag.input_router import build_raw_intermediates
from neoag.tools.upstream import run_upstream


HLA_RE = re.compile(r"(?:HLA-)?(?:A|B|C)\*[0-9]{2,3}(?::[0-9A-Z]{2,3}){1,4}", re.I)


def read_hla(path: Path) -> list[str]:
    values = HLA_RE.findall(path.read_text(encoding="utf-8", errors="replace"))
    return list(dict.fromkeys(value.upper() if value.upper().startswith("HLA-") else "HLA-" + value.upper() for value in values))


def main() -> int:
    parser = argparse.ArgumentParser(description="Normalize one candidate source using a runtime HLA consensus")
    parser.add_argument("--mode", choices=["snv", "fusion", "splice"], required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--hla-file", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--reference-fasta", default="")
    parser.add_argument("--vep-cache", default="")
    parser.add_argument("--vep-plugins", default=os.environ.get("NEOAG_VEP_PLUGINS", ""))
    parser.add_argument("--normal-proteome", default="")
    args = parser.parse_args()

    source = Path(args.input)
    hla_file = Path(args.hla_file)
    if not source.is_file() or not hla_file.is_file():
        raise SystemExit("input and HLA consensus must exist")
    alleles = read_hla(hla_file)
    if not alleles:
        raise SystemExit("HLA consensus contains no class-I alleles")
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    inputs: dict[str, object] = {"hla_alleles": alleles}
    if args.mode == "snv":
        inputs.update({
            "entry_mode": "snv_indel",
            "variants_vcf": str(source.resolve()),
            "variant_peptide_extraction": True,
            "auto_vep_annotate": True,
            "reference_fasta": args.reference_fasta,
            "vep_cache": args.vep_cache,
            "vep_plugins": args.vep_plugins,
            "normal_proteome_fasta": args.normal_proteome,
        })
        cfg = {"sample": {"id": args.sample_id, "profile": "default"}, "tools": {"enabled": []}, "inputs": inputs}
        config = outdir / "candidate_upstream.runtime.toml"
        config.write_text(
            "[sample]\n"
            f'id = "{args.sample_id}"\nprofile = "default"\n\n'
            "[tools]\nenabled = []\n\n[inputs]\nentry_mode = \"snv_indel\"\n"
            f'variants_vcf = "{source.resolve()}"\nvariant_peptide_extraction = true\nauto_vep_annotate = true\n'
            f'hla_alleles = [{", ".join(repr(value) for value in alleles)}]\n'
            f'reference_fasta = "{args.reference_fasta}"\nvep_cache = "{args.vep_cache}"\n'
            f'vep_plugins = "{args.vep_plugins}"\n'
            f'normal_proteome_fasta = "{args.normal_proteome}"\n',
            encoding="utf-8",
        )
        run_upstream(config, outdir)
    else:
        key = "fusion_tsv" if args.mode == "fusion" else "splice_junction_tsv"
        mode = "fusion" if args.mode == "fusion" else "splice_junction"
        inputs.update({"entry_mode": mode, key: str(source.resolve())})
        cfg = {"sample": {"id": args.sample_id, "profile": "default"}, "inputs": inputs}
        build_raw_intermediates(cfg, outdir, root=Path.cwd())
    for name in ("raw_events.tsv", "raw_peptides.tsv"):
        path = outdir / "parsed" / name
        if not path.is_file():
            raise SystemExit(f"candidate normalization did not create {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
