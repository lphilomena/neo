#!/usr/bin/env python3
"""Generate a production manifest that reuses completed upstream tool results."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def q(value) -> str:
    return json.dumps(str(value))


def require(path: str | None, label: str) -> str:
    if not path or not Path(path).exists():
        raise SystemExit(f"{label} missing: {path}")
    return str(Path(path).resolve())


def require_path_list(value: str | None, label: str) -> str:
    """Validate a comma-separated path list while preserving lane order."""
    paths = [item.strip() for item in str(value or "").split(",") if item.strip()]
    if not paths:
        raise SystemExit(f"{label} missing: {value}")
    resolved = [require(path, label) for path in paths]
    return ",".join(resolved)


def resolve_result_file(path: str | None, label: str, names: tuple[str, ...]) -> str:
    resolved = Path(require(path, label))
    if resolved.is_file():
        return str(resolved)
    for name in names:
        matches = sorted(candidate for candidate in resolved.rglob(name) if candidate.is_file() and candidate.stat().st_size > 0)
        if matches:
            return str(matches[0].resolve())
    raise SystemExit(f"{label} has no recognized result file below {resolved}: {', '.join(names)}")


def stage(lines, name, *, outputs, source="", command="", required=True, depends=None):
    lines += ["", f"[stages.{name}]", f"required = {str(required).lower()}"]
    if source: lines.append(f"source = {q(source)}")
    if depends: lines.append("depends_on = [" + ", ".join(q(value) for value in depends) + "]")
    if command: lines.append(f"command = {q(command)}")
    lines.append(f"[stages.{name}.outputs]")
    lines.extend(f"{key} = {q(value)}" for key, value in outputs.items())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-root", default=Path(__file__).resolve().parents[1])
    ap.add_argument("--sample-id", required=True); ap.add_argument("--outdir", required=True); ap.add_argument("--output", required=True)
    ap.add_argument("--profile", default="profiles/sarcoma_rna_supported_v2_provisional.toml")
    ap.add_argument(
        "--evidence-consensus-rules",
        default="configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml",
        help="Independent Evidence-consensus R1-R4 rules; does not replace the weighted profile.",
    )
    ap.add_argument("--reference-fasta")
    ap.add_argument("--hla-file"); ap.add_argument("--optitype"); ap.add_argument("--spechla-typing"); ap.add_argument("--hla-la"); ap.add_argument("--somatic-vcf")
    ap.add_argument("--facets"); ap.add_argument("--sequenza"); ap.add_argument("--purple"); ap.add_argument("--ascat")
    ap.add_argument("--purity"); ap.add_argument("--cnv")
    ap.add_argument("--lohhla", required=True); ap.add_argument("--spechla-loh", required=True); ap.add_argument("--hla-loh")
    ap.add_argument("--expression"); ap.add_argument("--transcript-expression")
    ap.add_argument("--rna-fastq1", help="Tumor RNA FASTQ R1; comma-separate multiple lanes")
    ap.add_argument("--rna-fastq2", help="Tumor RNA FASTQ R2; comma-separate multiple lanes")
    ap.add_argument("--rna-bam", help="Existing coordinate-sorted tumor RNA BAM")
    ap.add_argument("--rna-vaf", help="Existing RNA ref/alt/depth/VAF table to reuse")
    ap.add_argument("--star-index", help="GRCh38 STAR index used when RNA FASTQ is supplied")
    ap.add_argument("--gencode-gtf", help="GENCODE GTF matching the STAR index and somatic VCF build")
    ap.add_argument("--star-executable", default="", help="Optional explicit STAR executable")
    ap.add_argument("--samtools-executable", default="samtools", help="samtools used to index RNA BAM")
    ap.add_argument("--rna-threads", type=int, default=16)
    ap.add_argument("--easyfuse"); ap.add_argument("--star-fusion"); ap.add_argument("--arriba")
    ap.add_argument("--junctions"); ap.add_argument("--snaf"); ap.add_argument("--splicemutr"); ap.add_argument("--normal-junctions")
    ap.add_argument("--normal-expression"); ap.add_argument("--normal-hla-ligands"); ap.add_argument("--reference-proteome")
    ap.add_argument("--netchop-executable", default="netChop"); ap.add_argument("--netchop-home", default="")
    ap.add_argument(
        "--skip-netmhcstabpan",
        action="store_true",
        help="Omit NetMHCstabpan when no approved local installation is available.",
    )
    args = ap.parse_args()
    if bool(args.rna_fastq1) != bool(args.rna_fastq2):
        raise SystemExit("--rna-fastq1 and --rna-fastq2 must be supplied together")
    if args.rna_threads < 1:
        raise SystemExit("--rna-threads must be a positive integer")
    root = Path(args.project_root).resolve()
    profile_path = Path(args.profile)
    if not profile_path.is_absolute():
        profile_path = root / profile_path
    profile = require(str(profile_path), "ranking profile")
    consensus_rules_path = Path(args.evidence_consensus_rules)
    if not consensus_rules_path.is_absolute():
        consensus_rules_path = root / consensus_rules_path
    evidence_consensus_rules = require(str(consensus_rules_path), "evidence-consensus rules")
    reference_fasta = require(args.reference_fasta, "reference FASTA") if args.reference_fasta else ""
    generated_hla = not bool(args.hla_file)
    if generated_hla:
        optitype = require(args.optitype, "OptiType result")
        spechla_typing = require(args.spechla_typing, "SpecHLA typing result")
        hla = "{outdir}/evidence/hla_typing/recommended_hla.txt"
    else:
        hla = require(args.hla_file, "HLA consensus")
        optitype = spechla_typing = ""
    if bool(args.purity) != bool(args.cnv):
        raise SystemExit("--purity and --cnv must be supplied together")
    generated_purity = not bool(args.purity and args.cnv)
    purity = "{outdir}/evidence/purity_cnv/recommended_purity.tsv" if generated_purity else require(args.purity, "purity consensus")
    cnv = "{outdir}/evidence/purity_cnv/recommended_cnv_segments.tsv" if generated_purity else require(args.cnv, "CNV consensus")
    purity_tools = {
        tool: require(getattr(args, tool), f"{tool} result")
        for tool in ("facets", "sequenza", "purple", "ascat")
        if getattr(args, tool)
    }
    if generated_purity and len(purity_tools) < 2:
        raise SystemExit(
            "At least two valid purity/CNV tool results are required from "
            "--facets, --sequenza, --purple, and --ascat"
        )
    lohhla = resolve_result_file(args.lohhla, "LOHHLA", ("hla_loh.tsv", "*HLAlossPrediction_CI*.txt", "*HLAlossPrediction_CI*"))
    spechla_loh = resolve_result_file(args.spechla_loh, "SpecHLA LOH", ("hla_loh.tsv", "spechla_hla_loh.tsv"))
    generated_hla_loh = not bool(args.hla_loh)
    hla_loh = "{outdir}/evidence/hla_loh/hla_loh_consensus.tsv" if generated_hla_loh else require(args.hla_loh, "HLA LOH consensus")
    presentation_predictors = ["netmhcpan", "mhcflurry", "netchop"]
    production_limitations = []
    if not args.skip_netmhcstabpan:
        presentation_predictors.insert(2, "netmhcstabpan")
    else:
        production_limitations.append("NETMHCSTABPAN_LOCAL_UNAVAILABLE")
    predictor_toml = "[" + ", ".join(q(tool) for tool in presentation_predictors) + "]"
    lines = ["# Generated from completed upstream tool results.", "[run]", f"sample_id = {q(args.sample_id)}", f"profile = {q(profile)}", f"outdir = {q(Path(args.outdir).resolve())}", f"hla_file = {q(hla)}", "tools_stub = false", "immunogenicity_stub = false", f"presentation_predictors = {predictor_toml}", f"required_presentation_predictors = {predictor_toml}", 'reports = "patient,technical"', f"netchop_executable = {q(args.netchop_executable)}"]
    if production_limitations:
        limitations_toml = "[" + ", ".join(q(item) for item in production_limitations) + "]"
        lines.append(f"production_limitations = {limitations_toml}")
    if args.netchop_home: lines.append(f"netchop_home = {q(args.netchop_home)}")
    if generated_hla:
        lines += ["", "[run.required_tool_groups.hla_typing]", 'tools = ["optitype", "spechla"]', "min_successful = 2", "require_all_declared = true"]
    if generated_purity:
        purity_tool_names = ", ".join(q(tool) for tool in purity_tools)
        lines += ["", "[run.required_tool_groups.purity_cnv]", f"tools = [{purity_tool_names}]", "min_successful = 2", "require_all_declared = true"]
    lines += ["", "[run.required_tool_groups.hla_loh]", 'tools = ["lohhla", "spechla"]', "min_successful = 2", "require_all_declared = true"]
    hla_dependency = []
    if generated_hla:
        stage(lines, "hla_optitype", outputs={"optitype_result": optitype}, required=True)
        stage(lines, "hla_spechla", outputs={"spechla_result": spechla_typing}, required=True)
        hla_args = f"--result-dir {q(optitype)} --result-dir {q(spechla_typing)}"
        if args.hla_la:
            hla_args += f" --result-dir {q(require(args.hla_la, 'HLA-LA result'))}"
        command = f"PYTHONPATH={q(root / 'src')} {q(sys.executable)} -m neoag.agent_skills.hla_typing_compare {hla_args} --sample-id {q(args.sample_id)} --outdir {{outdir}}/evidence/hla_typing && PYTHONPATH={q(root / 'src')} {q(sys.executable)} {q(root / 'scripts/hla_consensus_to_file.py')} --consensus {{outdir}}/evidence/hla_typing/hla_typing_consensus.tsv --output {q(hla)}"
        stage(lines, "hla_typing_consensus", command=command, outputs={"hla_file": hla, "hla_consensus": "{outdir}/evidence/hla_typing/hla_typing_consensus.tsv"}, depends=["hla_optitype", "hla_spechla"])
        hla_dependency = ["hla_typing_consensus"]
    for tool, result_path in purity_tools.items():
        stage(lines, f"purity_{tool}", outputs={f"{tool}_result": result_path}, required=False)
    if generated_purity:
        result_args = " ".join(f"--result-dir {q(path)}" for path in purity_tools.values())
        purity_command = f"PYTHONPATH={q(root / 'src')} {q(sys.executable)} -m neoag.agent_skills.purity_cnv_review {result_args} --sample-id {q(args.sample_id)} --outdir {{outdir}}/evidence/purity_cnv"
        stage(lines, "purity_cnv_consensus", command=purity_command, outputs={"purity_consensus": purity, "cnv_consensus": cnv, "purity_tool_summary": "{outdir}/evidence/purity_cnv/purity_cnv_tool_summary.tsv"}, depends=[f"purity_{tool}" for tool in purity_tools])
    stage(lines, "hla_loh_lohhla", outputs={"lohhla_hla_loh": lohhla}, required=True)
    stage(lines, "hla_loh_spechla", outputs={"spechla_hla_loh": spechla_loh}, required=True)
    if generated_hla_loh:
        loh_command = f"PYTHONPATH={q(root / 'src')} {q(sys.executable)} {q(root / 'scripts/build_hla_loh_consensus.py')} --sample-id {q(args.sample_id)} --lohhla {q(lohhla)} --spechla {q(spechla_loh)} --outdir {{outdir}}/evidence/hla_loh"
        stage(lines, "hla_loh_consensus", command=loh_command, outputs={"hla_loh_consensus": hla_loh, "hla_loh_summary": "{outdir}/evidence/hla_loh/hla_loh_summary.json"}, depends=["hla_loh_lohhla", "hla_loh_spechla"])

    rna_vaf = ""
    if args.rna_vaf:
        rna_vaf = require(args.rna_vaf, "RNA VAF evidence")
        stage(
            lines,
            "rna_alt_vaf_input",
            source="RNA_ALLELE_COUNTS",
            outputs={"rna_vaf": rna_vaf},
            required=False,
        )
    elif args.rna_bam:
        if not args.somatic_vcf:
            raise SystemExit("--rna-bam requires --somatic-vcf for RNA allele counting")
        rna_bam = require(args.rna_bam, "RNA BAM")
        bam_path = Path(rna_bam)
        existing_bai = next((candidate for candidate in (
            Path(str(bam_path) + ".bai"), bam_path.with_suffix(".bai"),
        ) if candidate.is_file() and candidate.stat().st_size > 0), None)
        rna_bai = str(existing_bai or Path(str(bam_path) + ".bai"))
        index_command = "" if existing_bai else f"{q(args.samtools_executable)} index -@ {args.rna_threads} {q(rna_bam)}"
        stage(
            lines,
            "rna_bam_input",
            source="RNA_BAM",
            command=index_command,
            outputs={"rna_bam": rna_bam, "rna_bai": rna_bai},
            required=True,
        )
        rna_vaf = "{outdir}/rna/rna_alt_vaf.tsv"
        allele_command = (
            f"PYTHONPATH={q(root / 'src')} {q(sys.executable)} "
            f"{q(root / 'scripts/rna_allele_counts_pysam.py')} "
            f"--somatic-vcf {q(require(args.somatic_vcf, 'somatic VCF'))} "
            f"--rna-bam {q(rna_bam)} --output-tsv {q(rna_vaf)}"
        )
        stage(
            lines,
            "rna_alt_vaf",
            source="RNA_BAM_PILEUP",
            command=allele_command,
            outputs={"rna_vaf": rna_vaf},
            depends=["rna_bam_input"],
        )
    elif args.rna_fastq1 and args.rna_fastq2:
        if not args.somatic_vcf:
            raise SystemExit("RNA FASTQ allele counting requires --somatic-vcf")
        fastq1 = require_path_list(args.rna_fastq1, "RNA FASTQ R1")
        fastq2 = require_path_list(args.rna_fastq2, "RNA FASTQ R2")
        star_index = require(args.star_index, "STAR index")
        if not Path(star_index).is_dir():
            raise SystemExit(f"STAR index is not a directory: {star_index}")
        gencode_gtf = require(args.gencode_gtf, "GENCODE GTF")
        star_dir = "{outdir}/rna/star"
        rna_bam = f"{star_dir}/Aligned.sortedByCoord.out.bam"
        rna_bai = f"{rna_bam}.bai"
        star_env = f"NEOAG_STAR_BIN={q(require(args.star_executable, 'STAR executable'))} " if args.star_executable else ""
        star_command = (
            f"{star_env}bash {q(root / 'scripts/run_star_rna_fastq.sh')} "
            f"--fastq1 {q(fastq1)} --fastq2 {q(fastq2)} "
            f"--star-index {q(star_index)} --gtf {q(gencode_gtf)} "
            f"--sample-id {q(args.sample_id)} --outdir {q(star_dir)} --threads {args.rna_threads} "
            f"&& test -s {q(rna_bam)} && test -s {q(rna_bai)}"
        )
        stage(
            lines,
            "rna_star_alignment",
            source="RNA_FASTQ",
            command=star_command,
            outputs={"rna_bam": rna_bam, "rna_bai": rna_bai},
            required=True,
        )
        rna_vaf = "{outdir}/rna/rna_alt_vaf.tsv"
        allele_command = (
            f"PYTHONPATH={q(root / 'src')} {q(sys.executable)} "
            f"{q(root / 'scripts/rna_allele_counts_pysam.py')} "
            f"--somatic-vcf {q(require(args.somatic_vcf, 'somatic VCF'))} "
            f"--rna-bam {q(rna_bam)} --output-tsv {q(rna_vaf)}"
        )
        stage(
            lines,
            "rna_alt_vaf",
            source="RNA_BAM_PILEUP",
            command=allele_command,
            outputs={"rna_vaf": rna_vaf},
            depends=["rna_star_alignment"],
        )
    candidate_stages = []
    if args.somatic_vcf:
        reference_env = f"NEOAG_REFERENCE_FASTA={q(reference_fasta)} " if reference_fasta else ""
        command = f"{reference_env}PYTHONPATH={q(root / 'src')} {q(sys.executable)} {q(root / 'scripts/run_candidate_upstream.py')} --mode snv --input {q(require(args.somatic_vcf, 'somatic VCF'))} --hla-file {q(hla)} --sample-id {q(args.sample_id)} --outdir {{outdir}}/branches/snv"
        stage(lines, "snv_indel_candidates", source="SNV_INDEL", command=command, outputs={"raw_events": "{outdir}/branches/snv/parsed/raw_events.tsv", "raw_peptides": "{outdir}/branches/snv/parsed/raw_peptides.tsv"}, depends=hla_dependency)
        candidate_stages.append("snv_indel_candidates")
    if args.easyfuse or args.star_fusion or args.arriba:
        union_args = []
        easyfuse = require(args.easyfuse, "EasyFuse") if args.easyfuse else ""
        if easyfuse: union_args += ["--easyfuse", q(easyfuse)]
        if args.star_fusion: union_args += ["--star-fusion", q(require(args.star_fusion, "STAR-Fusion"))]
        if args.arriba: union_args += ["--arriba", q(require(args.arriba, "Arriba"))]
        command = f"PYTHONPATH={q(root / 'src')} {q(sys.executable)} {q(root / 'scripts/build_fusion_caller_union.py')} --sample-id {q(args.sample_id)} --profile {q(profile)} --hla-file {q(hla)} {' '.join(union_args)} --outdir {{outdir}}/branches/fusion/intermediates"
        stage(lines, "fusion_candidates", source="FusionCallerUnion", command=command, outputs={"raw_events": "{outdir}/branches/fusion/intermediates/raw_events.tsv", "raw_peptides": "{outdir}/branches/fusion/intermediates/raw_peptides.tsv", "fusion_union": "{outdir}/branches/fusion/intermediates/fusion_caller_union.tsv"}, depends=hla_dependency)
        candidate_stages.append("fusion_candidates")
        review = f"{q(sys.executable)} {q(root / 'scripts/review_rna_fusions.py')}"
        if easyfuse: review += f" --easyfuse {easyfuse}"
        if args.star_fusion: review += f" --star-fusion {require(args.star_fusion, 'STAR-Fusion')}"
        if args.arriba: review += f" --arriba {require(args.arriba, 'Arriba')}"
        review += " --outdir {outdir}/branches/fusion/consensus"
        stage(lines, "fusion_cross_validation", command=review, outputs={"fusion_consensus": "{outdir}/branches/fusion/consensus/fusion_consensus.tsv"}, required=True, depends=["fusion_candidates"])
    if args.junctions and (args.snaf or args.splicemutr):
        command = f"PYTHONPATH={q(root / 'src')} {q(sys.executable)} {q(root / 'scripts/normalize_rna_fusion_splice.py')} --sample-id {q(args.sample_id)} --profile {q(profile)} --junctions {q(require(args.junctions, 'junctions'))} --candidate-only"
        if args.snaf: command += f" --snaf {require(args.snaf, 'SNAF')}"
        if args.splicemutr: command += f" --splicemutr {require(args.splicemutr, 'SpliceMutr')}"
        if args.normal_junctions: command += f" --normal-junctions {require(args.normal_junctions, 'normal junctions')}"
        command += " --outdir {outdir}/branches/splice/intermediates"
        command += (
            f" && PYTHONPATH={q(root / 'src')} {q(sys.executable)} {q(root / 'scripts/filter_splice_production_candidates.py')}"
            " --events {outdir}/branches/splice/intermediates/raw_events.tsv"
            " --peptides {outdir}/branches/splice/intermediates/raw_peptides.tsv"
            " --consensus {outdir}/branches/splice/intermediates/splice_consensus.tsv"
            " --outdir {outdir}/branches/splice/production_selected"
            " --min-length 8 --max-length 12 --max-source-binding-rank 2.0"
        )
        stage(lines, "splice_candidates", source="SpliceConsensus", command=command, outputs={"raw_events": "{outdir}/branches/splice/production_selected/raw_events.tsv", "raw_peptides": "{outdir}/branches/splice/production_selected/raw_peptides.tsv", "production_filter_summary": "{outdir}/branches/splice/production_selected/production_filter_summary.json"}, depends=hla_dependency)
        candidate_stages.append("splice_candidates")
    if not candidate_stages: raise SystemExit("At least one SNV/Fusion/Splice candidate source is required")
    lines += [
        "", "[evidence]", f"purity = {q(purity)}", f"cnv = {q(cnv)}",
        f"hla_loh = {q(hla_loh)}",
        f"evidence_consensus_rules = {q(evidence_consensus_rules)}",
    ]
    if rna_vaf:
        lines.append(f"rna_vaf = {q(rna_vaf)}")
    for key in ("expression", "transcript_expression", "normal_junctions", "normal_expression", "normal_hla_ligands", "reference_proteome"):
        value = getattr(args, key)
        if value: lines.append(f"{key} = {q(require(value, key))}")
    Path(args.output).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__": raise SystemExit(main())
