#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
from pathlib import Path

from neoag.input_router import build_raw_intermediates
from neoag.tools.upstream import run_upstream


HLA_RE = re.compile(r"(?:HLA-)?(?:A|B|C)\*[0-9]{2,3}(?::[0-9A-Z]{2,3}){1,4}", re.I)

NORMAL_SAMPLE_TOKENS = ("normal", "blood", "germline", "control")
TUMOR_SAMPLE_TOKENS = ("tumor", "tumour", "cancer", "lesion")


def read_hla(path: Path) -> list[str]:
    values = HLA_RE.findall(path.read_text(encoding="utf-8", errors="replace"))
    return list(dict.fromkeys(value.upper() if value.upper().startswith("HLA-") else "HLA-" + value.upper() for value in values))


def discover_vep_cache(explicit: str, reference_fasta: str) -> str:
    reference_fasta = reference_fasta or os.environ.get("NEOAG_REFERENCE_FASTA", "")
    candidates = [
        explicit,
        os.environ.get("NEOAG_VEP_CACHE", ""),
        str(Path(os.environ["OPEN_NEO_ASSET_ROOT"]) / "data/vep") if os.environ.get("OPEN_NEO_ASSET_ROOT") else "",
        str(Path(os.environ["OPEN_NEO_REFERENCE_ROOT"]) / "data/vep") if os.environ.get("OPEN_NEO_REFERENCE_ROOT") else "",
    ]
    if reference_fasta:
        candidates.extend(str(parent / "vep") for parent in Path(reference_fasta).resolve().parents)
    for value in dict.fromkeys(item for item in candidates if item):
        path = Path(value).expanduser()
        if (path / "homo_sapiens").is_dir():
            return str(path.resolve())
        if path.name.endswith("_GRCh38") and path.parent.name == "homo_sapiens":
            return str(path.parent.parent.resolve())
    return ""


def _vcf_samples(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    return []


def resolve_vcf_sample_roles(
    path: Path,
    *,
    tumor_sample_name: str = "",
    normal_sample_name: str = "",
) -> tuple[str, str]:
    """Resolve paired VCF roles conservatively; never guess from column order alone."""
    samples = _vcf_samples(path)
    tumor = tumor_sample_name.strip()
    normal = normal_sample_name.strip()
    for role, value in (("tumor", tumor), ("normal", normal)):
        if value and value not in samples:
            raise SystemExit(f"{role} VCF sample {value!r} is not present; samples={samples}")
    if len(samples) == 2 and tumor and not normal:
        normal = next((sample for sample in samples if sample != tumor), "")
    if len(samples) == 2 and normal and not tumor:
        tumor = next((sample for sample in samples if sample != normal), "")
    if not normal:
        matches = [sample for sample in samples if any(token in sample.lower() for token in NORMAL_SAMPLE_TOKENS)]
        if len(matches) == 1:
            normal = matches[0]
    if not tumor:
        matches = [sample for sample in samples if any(token in sample.lower() for token in TUMOR_SAMPLE_TOKENS)]
        if len(matches) == 1:
            tumor = matches[0]
    if len(samples) == 2 and tumor and not normal:
        normal = next((sample for sample in samples if sample != tumor), "")
    if len(samples) == 2 and normal and not tumor:
        tumor = next((sample for sample in samples if sample != normal), "")
    if tumor and normal and tumor == normal:
        raise SystemExit("tumor and normal VCF sample names must differ")
    return tumor, normal


def main() -> int:
    parser = argparse.ArgumentParser(description="Normalize one candidate source using a runtime HLA consensus")
    parser.add_argument("--mode", choices=["snv", "fusion", "splice"], required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--hla-file", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--tumor-sample-name", default="")
    parser.add_argument("--normal-sample-name", default="")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--reference-fasta", default="")
    parser.add_argument("--vep-cache", default="")
    parser.add_argument("--vep-plugins", default=os.environ.get("NEOAG_VEP_PLUGINS", ""))
    parser.add_argument("--vep-bin", default=os.environ.get("NEOAG_VEP_BIN", ""))
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
        reference_fasta = args.reference_fasta or os.environ.get("NEOAG_REFERENCE_FASTA", "")
        vep_cache = discover_vep_cache(args.vep_cache, reference_fasta)
        if not vep_cache:
            raise SystemExit(
                "VEP cache not found; provide --vep-cache or set NEOAG_VEP_CACHE/OPEN_NEO_ASSET_ROOT"
            )
        tumor_sample_name, normal_sample_name = resolve_vcf_sample_roles(
            source,
            tumor_sample_name=args.tumor_sample_name,
            normal_sample_name=args.normal_sample_name,
        )
        inputs.update({
            "entry_mode": "snv_indel",
            "variants_vcf": str(source.resolve()),
            "variant_peptide_extraction": True,
            "auto_vep_annotate": True,
            "reference_fasta": reference_fasta,
            "vep_cache": vep_cache,
            "vep_plugins": args.vep_plugins,
            "vep_bin": args.vep_bin,
            "normal_proteome_fasta": args.normal_proteome,
            "tumor_sample_name": tumor_sample_name,
            "normal_sample_name": normal_sample_name,
        })
        cfg = {"sample": {"id": args.sample_id, "profile": "default"}, "tools": {"enabled": []}, "inputs": inputs}
        config = outdir / "candidate_upstream.runtime.toml"
        config.write_text(
            "[sample]\n"
            f'id = "{args.sample_id}"\nprofile = "default"\n\n'
            "[tools]\nenabled = []\n\n[inputs]\nentry_mode = \"snv_indel\"\n"
            f'variants_vcf = "{source.resolve()}"\nvariant_peptide_extraction = true\nauto_vep_annotate = true\n'
            f'hla_alleles = [{", ".join(repr(value) for value in alleles)}]\n'
            f'reference_fasta = "{reference_fasta}"\nvep_cache = "{vep_cache}"\n'
            f'vep_plugins = "{args.vep_plugins}"\n'
            f'vep_bin = "{args.vep_bin}"\n'
            f'normal_proteome_fasta = "{args.normal_proteome}"\n'
            f'tumor_sample_name = {json.dumps(tumor_sample_name)}\n'
            f'normal_sample_name = {json.dumps(normal_sample_name)}\n',
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
