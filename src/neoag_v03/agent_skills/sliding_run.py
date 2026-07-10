from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import subprocess
from pathlib import Path
from typing import Iterable


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _read_hla_alleles(path: Path) -> list[str]:
    alleles: list[str] = []
    text = path.read_text(encoding="utf-8", errors="replace")
    for token in re.split(r"[\s,;]+", text):
        token = token.strip().strip('"').strip("'")
        if not token or token.startswith("#"):
            continue
        if token.upper().startswith("HLA-") and "*" in token:
            alleles.append(token)
    seen: set[str] = set()
    out: list[str] = []
    for allele in alleles:
        if allele not in seen:
            seen.add(allele)
            out.append(allele)
    return out


def _sample_id_from_vcf(path: Path) -> str:
    name = path.name
    for suffix in (".vcf.gz", ".vcf", ".gz"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    name = re.sub(r"\.somatic\..*$", "", name)
    name = re.sub(r"\.align$", "", name)
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name) or "SAMPLE001"


def _tumor_sample_from_vcf(path: Path) -> str | None:
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                samples = parts[9:]
                return samples[0] if samples else None
    return None


def _toml_array(values: Iterable[str]) -> str:
    return "[" + ", ".join(json.dumps(v) for v in values) + "]"


def _default_reference_fasta(root: Path) -> Path | None:
    candidates: list[Path] = []
    for key in ("NEOAG_REFERENCE_FASTA", "REFERENCE_FASTA"):
        value = os.environ.get(key)
        if value:
            candidates.append(Path(value))
    candidates.append(root / "data/ref/hg38/Homo_sapiens_assembly38.fasta")
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def _write_config(*, config_path: Path, sample_id: str, profile: str, variants_vcf: Path, tumor_sample_name: str, hla_alleles: list[str], tools_stub: bool, immunogenicity_stub: bool, reference_fasta: Path | None = None) -> None:
    config_path.parent.mkdir(parents=True, exist_ok=True)
    enabled = ["netmhcpan", "mhcflurry"]
    reference_line = f'reference_fasta = {json.dumps(str(reference_fasta))}\n' if reference_fasta else ""
    text = f'''[sample]
id = {json.dumps(sample_id)}
profile = {json.dumps(profile)}

[tools]
stub = {str(bool(tools_stub)).lower()}
enabled = {_toml_array(enabled)}
immunogenicity_stub = {str(bool(immunogenicity_stub)).lower()}

[tools.executables]
netmhcpan = "{Path('scripts/run_netmhcpan_container.sh').resolve()}"
mhcflurry = "mhcflurry-predict"

[inputs]
entry_mode = "snv_indel"
variant_peptide_extraction = true
variants_vcf = {json.dumps(str(variants_vcf))}
tumor_sample_name = {json.dumps(tumor_sample_name)}
hla_alleles = {_toml_array(hla_alleles)}
extract_appm_from_vcf = true
{reference_line}normal_expression = "resources/normal_expression.example.tsv"
normal_hla_ligands = "resources/normal_hla_ligands.example.tsv"
'''
    config_path.write_text(text, encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Run NeoAg sliding-window SNV/InDel workflow from VCF + HLA file")
    ap.add_argument("--variants-vcf", required=True)
    ap.add_argument("--hla", required=True, help="HLA allele text file")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--project-root", default=".")
    ap.add_argument("--sample-id")
    ap.add_argument("--tumor-sample-name")
    ap.add_argument("--profile", default="default")
    ap.add_argument("--tools-stub", action="store_true")
    ap.add_argument("--immunogenicity-stub", action="store_true")
    args = ap.parse_args(argv)

    root = Path(args.project_root).resolve()
    variants_vcf = Path(args.variants_vcf)
    hla_path = Path(args.hla)
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if not variants_vcf.is_file():
        raise SystemExit(f"variants VCF not found: {variants_vcf}")
    if not hla_path.is_file():
        raise SystemExit(f"HLA file not found: {hla_path}")
    hla_alleles = _read_hla_alleles(hla_path)
    if not hla_alleles:
        raise SystemExit(f"No HLA alleles parsed from: {hla_path}")

    sample_id = args.sample_id or _sample_id_from_vcf(variants_vcf)
    tumor_sample_name = args.tumor_sample_name or _tumor_sample_from_vcf(variants_vcf) or sample_id
    reference_fasta = _default_reference_fasta(root)
    config_path = outdir / "run.sliding.private.toml"
    _write_config(
        config_path=config_path,
        sample_id=sample_id,
        profile=args.profile,
        variants_vcf=variants_vcf,
        tumor_sample_name=tumor_sample_name,
        hla_alleles=hla_alleles,
        tools_stub=args.tools_stub,
        immunogenicity_stub=args.immunogenicity_stub,
        reference_fasta=reference_fasta,
    )

    run_outdir = outdir / "run-full"
    log_path = outdir / "run-full.log"
    cmd = f"set -euo pipefail; cd {root}; source conf/tools.env.sh; bin/neoag-v03 run-full --config {config_path} --outdir {run_outdir}"
    proc = subprocess.run(["bash", "-lc", cmd], text=True, capture_output=True)
    log_path.write_text(proc.stdout + ("\n--- STDERR ---\n" if proc.stderr else "") + proc.stderr, encoding="utf-8")

    summary = {
        "config": str(config_path),
        "run_outdir": str(run_outdir),
        "log": str(log_path),
        "sample_id": sample_id,
        "tumor_sample_name": tumor_sample_name,
        "hla_alleles": hla_alleles,
        "reference_fasta": str(reference_fasta) if reference_fasta else None,
        "returncode": proc.returncode,
    }
    (outdir / "sliding_run_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return proc.returncode


if __name__ == "__main__":
    raise SystemExit(main())
