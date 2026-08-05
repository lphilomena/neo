#!/usr/bin/env python3
"""Generate a production run-full TOML from recommended upstream evidence."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


def q(value: str | Path) -> str:
    return json.dumps(str(value))


def required(path: str | None, label: str) -> Path:
    if not path or not Path(path).is_file() or Path(path).stat().st_size == 0:
        raise SystemExit(f"ERROR: {label} missing: {path}")
    return Path(path).resolve()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--tumor-sample-name", required=True)
    ap.add_argument("--normal-sample-name", required=True)
    ap.add_argument("--hla-file", required=True)
    ap.add_argument("--purity", required=True)
    ap.add_argument("--cnv", required=True)
    ap.add_argument("--hla-loh", required=True)
    ap.add_argument("--expression")
    ap.add_argument("--transcript-expression")
    ap.add_argument("--rna-vaf")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    vcf = required(args.vcf, "somatic VCF")
    hla_file = required(args.hla_file, "recommended HLA")
    purity = required(args.purity, "recommended purity")
    cnv = required(args.cnv, "recommended CNV")
    hla_loh = required(args.hla_loh, "recommended HLA LOH")
    hla = [line.strip() for line in hla_file.read_text().splitlines() if line.strip() and not line.startswith("#")]
    if not hla:
        raise SystemExit("ERROR: recommended HLA file is empty")

    ref_bundle = Path(os.environ.get("NEOAG_REF_BUNDLE", "/root/neo/refs/neodata4git"))
    reference = Path(os.environ.get("NEOAG_REFERENCE_FASTA", str(ref_bundle / "data/ref/hg38/Homo_sapiens_assembly38.fasta")))
    vep_cache = Path(os.environ.get("NEOAG_VEP_CACHE", str(ref_bundle / "data/vep")))
    vep_plugins = Path(os.environ.get("NEOAG_VEP_PLUGINS", str(ref_bundle / "work/vep_plugins")))
    proteome = Path(os.environ.get("NEOAG_NORMAL_PROTEOME_FASTA", str(ref_bundle / "data/normal/proteome/Homo_sapiens.GRCh38.pep.all.fa")))
    normal_expression = ref_bundle / "data/normal/expression/normal_expression.gtex_v11_hpa_hspc.tsv"
    normal_ligands = ref_bundle / "data/normal/ligandome/normal_ms_ligands.tsv"
    normal_junctions = ref_bundle / "data/normal/junctions/normal_junctions.GRCh38.tsv.gz"
    project = Path(os.environ.get("NEOAG_PROJECT_ROOT", Path(__file__).resolve().parents[1]))
    outdir = Path(args.outdir).resolve()

    for path, label in ((reference, "GRCh38 FASTA"), (vep_cache, "VEP cache"), (proteome, "normal proteome"), (normal_expression, "normal expression"), (normal_ligands, "normal HLA ligandome")):
        if not path.exists():
            raise SystemExit(f"ERROR: {label} missing: {path}")

    netmhcpan = os.environ.get("NEOAG_NETMHCPAN_BIN", str(project / "scripts/run_netmhcpan_container.sh"))
    mhcflurry = os.environ.get("MHCFLURRY_BIN", "mhcflurry-predict")
    prime = os.environ.get("NEOAG_PRIME_BIN", os.environ.get("PRIME_HOME", "") + "/PRIME")
    mixmhcpred = os.environ.get("MIXMHCPRED_BIN", os.environ.get("MIXMHCPRED_HOME", "") + "/MixMHCpred")
    rules = project / "configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml"

    lines = [
        "[sample]", f"id = {q(args.sample_id)}", 'profile = "sarcoma_rna_supported_v2_provisional"' if args.expression else 'profile = "sarcoma"', f"outdir = {q(outdir)}", "",
        "[tools]", "stub = false", "enabled = []", "immunogenicity_stub = false", "",
        "[tools.executables]", f"netmhcpan = {q(netmhcpan)}", f"mhcflurry = {q(mhcflurry)}", f"prime = {q(prime)}", f"mixmhcpred = {q(mixmhcpred)}", "",
        "[inputs]", 'entry_mode = "snv_indel"', "variant_peptide_extraction = true", "auto_vep_annotate = true", "vep_online = false", "vep_fork = 6", 'vep_cache_version = "105"',
        f"variants_vcf = {q(vcf)}", f"tumor_sample_name = {q(args.tumor_sample_name)}", f"normal_sample_name = {q(args.normal_sample_name)}",
        f"reference_fasta = {q(reference)}", f"gencode_gtf = {q(ref_bundle / 'data/ref/hg38/gencode.gtf')}", f"vep_cache = {q(vep_cache)}", f"vep_plugins = {q(vep_plugins)}", f"vep_bin = {q(os.environ.get('NEOAG_VEP_BIN', project / 'bin/vep-neoag'))}",
        "hla_alleles = [" + ", ".join(q(x) for x in hla) + "]", "variant_peptide_length_min = 8", "variant_peptide_length_max = 11",
        f"normal_proteome_fasta = {q(proteome)}", f"reference_proteome = {q(proteome)}", "variant_peptide_annotate_normal_only = true", "extract_appm_from_vcf = true",
        f"purity_tsv = {q(purity)}", f"cnv_tsv = {q(cnv)}", f"hla_loh_tsv = {q(hla_loh)}",
        f"normal_expression = {q(normal_expression)}", f"normal_hla_ligands = {q(normal_ligands)}", f"normal_junctions = {q(normal_junctions)}", f"evidence_consensus_rules = {q(rules)}",
    ]
    optional = (("expression", args.expression), ("transcript_expression_tsv", args.transcript_expression), ("rna_vaf_tsv", args.rna_vaf))
    for key, value in optional:
        if value:
            lines.append(f"{key} = {q(required(value, key))}")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
