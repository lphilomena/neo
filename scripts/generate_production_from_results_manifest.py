#!/usr/bin/env python3
"""Generate a production manifest that reuses completed upstream tool results."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


STAR_INDEX_REQUIRED_FILES = ("Genome", "SA", "SAindex", "genomeParameters.txt")


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


def inspect_star_index(path: str | Path | None, reference_fasta: str = "") -> dict[str, object]:
    """Validate a reusable STAR index without modifying it."""
    index = Path(str(path or "")).expanduser()
    missing = [name for name in STAR_INDEX_REQUIRED_FILES if not (index / name).is_file() or (index / name).stat().st_size == 0]
    result: dict[str, object] = {
        "path": str(index),
        "status": "VALID" if index.is_dir() and not missing else "INVALID",
        "missing_files": missing,
        "reference_check": "UNASSESSED",
        "parameters": {},
    }
    parameters = index / "genomeParameters.txt"
    if parameters.is_file():
        parsed: dict[str, str] = {}
        for line in parameters.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split(None, 1)
            if len(fields) == 2:
                parsed[fields[0]] = fields[1].strip()
        result["parameters"] = parsed
    fasta = Path(reference_fasta) if reference_fasta else None
    fai = Path(str(fasta) + ".fai") if fasta else None
    chr_lengths = index / "chrNameLength.txt"
    if result["status"] == "VALID" and fai and fai.is_file() and chr_lengths.is_file():
        expected = {
            row.split("\t", 2)[0]: row.split("\t", 2)[1]
            for row in fai.read_text(encoding="utf-8").splitlines()
            if row.count("\t") >= 1
        }
        observed = {
            row.split("\t", 2)[0]: row.split("\t", 2)[1]
            for row in chr_lengths.read_text(encoding="utf-8").splitlines()
            if row.count("\t") >= 1
        }
        if expected == observed:
            result["reference_check"] = "CONTIG_LENGTHS_MATCH"
        else:
            result["reference_check"] = "CONTIG_LENGTHS_MISMATCH"
            result["status"] = "INVALID"
    elif result["status"] == "VALID":
        result["reference_check"] = "CORE_FILES_ONLY"
    return result


def easyfuse_star_candidates(explicit: str | None) -> list[str]:
    candidates = [explicit or "", os.environ.get("EASYFUSE_STAR_INDEX", "")]
    for root in (os.environ.get("NEOAG_EASYFUSE_REF", ""), os.environ.get("EASYFUSE_REF", "")):
        if root:
            candidates.extend([
                str(Path(root) / "starfusion_index/ref_genome.fa.star.idx"),
                str(Path(root) / "star_index"),
            ])
    return list(dict.fromkeys(value for value in candidates if value))


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
    ap.add_argument("--secondary-rna-events", help="Standardized raw_events/all_tool_results from another tumor site")
    ap.add_argument("--secondary-sample-id", default="SECONDARY_RNA")
    ap.add_argument("--secondary-identity-status", default="UNASSESSED", help="CONFIRMED only after RNA-DNA fingerprint review")
    ap.add_argument("--rna-fastq1", help="Tumor RNA FASTQ R1; comma-separate multiple lanes")
    ap.add_argument("--rna-fastq2", help="Tumor RNA FASTQ R2; comma-separate multiple lanes")
    ap.add_argument("--rna-bam", help="Existing coordinate-sorted tumor RNA BAM")
    ap.add_argument("--rna-vaf", help="Existing RNA ref/alt/depth/VAF table to reuse")
    ap.add_argument("--star-index", help="GRCh38 STAR index used when RNA FASTQ is supplied")
    ap.add_argument("--easyfuse-star-index", help="Reusable EasyFuse/STAR-Fusion STAR index; used only after validation")
    ap.add_argument("--star-index-build-dir", help="Destination for a newly built STAR index when reusable indexes fail validation")
    ap.add_argument("--star-sjdb-overhang", type=int, default=149)
    ap.add_argument("--gencode-gtf", help="GENCODE GTF matching the STAR index and somatic VCF build")
    ap.add_argument("--star-executable", default="", help="Optional explicit STAR executable")
    ap.add_argument("--samtools-executable", default="samtools", help="samtools used to index RNA BAM")
    ap.add_argument("--rna-threads", type=int, default=16)
    ap.add_argument("--easyfuse"); ap.add_argument("--easyfuse-unfiltered")
    ap.add_argument("--diagnostic-fusion-whitelist", action="append", default=[], help="Exact disease-defining fusion label eligible for audited rescue; repeatable")
    ap.add_argument("--disable-diagnostic-fusion-rescue", action="store_true")
    ap.add_argument("--star-fusion"); ap.add_argument("--arriba")
    ap.add_argument("--fusioncatcher"); ap.add_argument("--jaffal")
    ap.add_argument("--fusion-caller-root", action="append", default=[], help="Directory containing completed fusion caller outputs; repeatable")
    ap.add_argument("--normal-readthrough", help="Normal/read-through fusion background table for review")
    ap.add_argument("--junctions"); ap.add_argument("--star-sj"); ap.add_argument("--snaf"); ap.add_argument("--splicemutr"); ap.add_argument("--normal-junctions")
    ap.add_argument("--normal-expression"); ap.add_argument("--normal-hla-ligands"); ap.add_argument("--reference-proteome")
    ap.add_argument("--prime-evidence", help="Existing normalized PRIME evidence TSV to reuse")
    ap.add_argument("--bigmhc-evidence", help="Existing normalized BigMHC_IM evidence TSV to reuse")
    ap.add_argument("--deepimmuno-evidence", help="Existing normalized DeepImmuno evidence TSV to reuse")
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
    if args.star_sjdb_overhang < 1:
        raise SystemExit("--star-sjdb-overhang must be a positive integer")
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
        gencode_gtf = require(args.gencode_gtf, "GENCODE GTF")
        inspected: list[dict[str, object]] = []
        star_index = ""
        star_index_source = ""
        candidates = []
        if args.star_index:
            candidates.append(("EXPLICIT", args.star_index))
        candidates.extend(("EASYFUSE", path) for path in easyfuse_star_candidates(args.easyfuse_star_index))
        for source, candidate in candidates:
            check = inspect_star_index(candidate, reference_fasta)
            check["source"] = source
            inspected.append(check)
            if check["status"] == "VALID":
                star_index = str(Path(candidate).resolve())
                star_index_source = source
                break

        validation_path = Path(args.output).resolve().parent / "star_index_validation.json"
        star_index_dependency: list[str] = []
        if not star_index:
            if not reference_fasta:
                detail = "; ".join(f"{row['source']}={row['path']}:{row['status']}" for row in inspected) or "no candidates"
                raise SystemExit(
                    "No valid STAR index was found and --reference-fasta was not supplied for rebuilding: " + detail
                )
            build_dir = args.star_index_build_dir or "{outdir}/rna/star_index"
            star_index = build_dir
            star_index_source = "BUILT_FOR_RUN"
            star_env = f"NEOAG_STAR_BIN={q(require(args.star_executable, 'STAR executable'))} " if args.star_executable else ""
            build_command = (
                f"{star_env}bash {q(root / 'scripts/build_star_index.sh')} "
                f"--reference-fasta {q(reference_fasta)} --gtf {q(gencode_gtf)} "
                f"--star-index {q(star_index)} --threads {args.rna_threads} "
                f"--sjdb-overhang {args.star_sjdb_overhang}"
            )
            stage(
                lines,
                "rna_star_index",
                source="STAR_INDEX_BUILD",
                command=build_command,
                outputs={
                    "star_index": star_index,
                    "genome_parameters": f"{star_index}/genomeParameters.txt",
                    "validation": str(validation_path),
                },
                required=True,
            )
            star_index_dependency = ["rna_star_index"]

        validation_path.parent.mkdir(parents=True, exist_ok=True)
        validation_payload = {
            "selected_path": star_index,
            "selected_source": star_index_source,
            "status": "PLANNED_BUILD" if star_index_source == "BUILT_FOR_RUN" else "VALIDATED_REUSE",
            "reference_fasta": reference_fasta,
            "gencode_gtf": gencode_gtf,
            "sjdb_overhang": args.star_sjdb_overhang,
            "candidates": inspected,
        }
        validation_path.write_text(json.dumps(validation_payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        if star_index_source != "BUILT_FOR_RUN":
            stage(
                lines,
                "rna_star_index",
                source=f"STAR_INDEX_REUSE_{star_index_source}",
                outputs={"star_index": star_index, "validation": str(validation_path)},
                required=True,
            )
            star_index_dependency = ["rna_star_index"]
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
            depends=star_index_dependency,
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
        vep_cache = os.environ.get("NEOAG_VEP_CACHE", "")
        if not vep_cache and args.normal_junctions:
            candidate_cache = Path(args.normal_junctions).resolve().parents[2] / "vep"
            if candidate_cache.is_dir():
                vep_cache = str(candidate_cache)
        vep_cache_arg = f" --vep-cache {q(vep_cache)}" if vep_cache and Path(vep_cache).is_dir() else ""
        command = f"{reference_env}PYTHONPATH={q(root / 'src')} {q(sys.executable)} {q(root / 'scripts/run_candidate_upstream.py')} --mode snv --input {q(require(args.somatic_vcf, 'somatic VCF'))} --hla-file {q(hla)} --sample-id {q(args.sample_id)} --outdir {{outdir}}/branches/snv{vep_cache_arg}"
        stage(lines, "snv_indel_candidates", source="SNV_INDEL", command=command, outputs={"raw_events": "{outdir}/branches/snv/parsed/raw_events.tsv", "raw_peptides": "{outdir}/branches/snv/parsed/raw_peptides.tsv"}, depends=hla_dependency)
        candidate_stages.append("snv_indel_candidates")
    if args.easyfuse or args.easyfuse_unfiltered or args.star_fusion or args.arriba or args.fusioncatcher or args.jaffal or args.fusion_caller_root:
        union_args = []
        easyfuse = require(args.easyfuse, "EasyFuse") if args.easyfuse else ""
        if easyfuse: union_args += ["--easyfuse", q(easyfuse)]
        if args.easyfuse_unfiltered: union_args += ["--easyfuse-unfiltered", q(require(args.easyfuse_unfiltered, "unfiltered EasyFuse"))]
        for fusion_label in args.diagnostic_fusion_whitelist:
            union_args += ["--diagnostic-fusion-whitelist", q(fusion_label)]
        if args.disable_diagnostic_fusion_rescue:
            union_args += ["--disable-diagnostic-fusion-rescue"]
        if args.star_fusion: union_args += ["--star-fusion", q(require(args.star_fusion, "STAR-Fusion"))]
        if args.arriba: union_args += ["--arriba", q(require(args.arriba, "Arriba"))]
        star_junction_source = args.star_sj or args.junctions
        if star_junction_source:
            chimeric = Path(star_junction_source).with_name("Chimeric.out.junction")
            if chimeric.is_file() and chimeric.stat().st_size > 0:
                union_args += ["--star-chimeric", q(str(chimeric))]
        if args.fusioncatcher: union_args += ["--fusioncatcher", q(require(args.fusioncatcher, "FusionCatcher"))]
        if args.jaffal: union_args += ["--jaffal", q(require(args.jaffal, "JAFFAL"))]
        for caller_root in args.fusion_caller_root:
            union_args += ["--caller-root", q(require(caller_root, "fusion caller result root"))]
        command = f"PYTHONPATH={q(root / 'src')} {q(sys.executable)} {q(root / 'scripts/build_fusion_caller_union.py')} --sample-id {q(args.sample_id)} --profile {q(profile)} --hla-file {q(hla)} {' '.join(union_args)} --outdir {{outdir}}/branches/fusion/intermediates"
        stage(lines, "fusion_candidates", source="FusionCallerUnion", command=command, outputs={"raw_events": "{outdir}/branches/fusion/intermediates/raw_events.tsv", "raw_peptides": "{outdir}/branches/fusion/intermediates/raw_peptides.tsv", "fusion_union": "{outdir}/branches/fusion/intermediates/fusion_caller_union.tsv", "fusion_consensus": "{outdir}/branches/fusion/intermediates/fusion_consensus.tsv", "diagnostic_fusion_rescue": "{outdir}/branches/fusion/intermediates/diagnostic_fusion_rescue.tsv"}, depends=hla_dependency)
        candidate_stages.append("fusion_candidates")
        review = f"{q(sys.executable)} {q(root / 'scripts/review_rna_fusions.py')}"
        if easyfuse: review += f" --easyfuse {q(easyfuse)}"
        if args.star_fusion: review += f" --star-fusion {q(require(args.star_fusion, 'STAR-Fusion'))}"
        if args.arriba: review += f" --arriba {q(require(args.arriba, 'Arriba'))}"
        if args.fusioncatcher: review += f" --fusioncatcher {q(require(args.fusioncatcher, 'FusionCatcher'))}"
        if args.jaffal: review += f" --jaffal {q(require(args.jaffal, 'JAFFAL'))}"
        for caller_root in args.fusion_caller_root:
            review += f" --caller-root {q(require(caller_root, 'fusion caller result root'))}"
        if args.normal_readthrough:
            review += f" --normal-readthrough {q(require(args.normal_readthrough, 'normal read-through background'))}"
        review += " --outdir {outdir}/branches/fusion/consensus"
        stage(lines, "fusion_cross_validation", command=review, outputs={"fusion_consensus": "{outdir}/branches/fusion/consensus/fusion_consensus.tsv"}, required=True, depends=["fusion_candidates"])
    if args.junctions and args.star_sj:
        raise SystemExit("Use only one primary splice evidence input: --junctions or --star-sj")
    if (args.junctions or args.star_sj) and (args.snaf or args.splicemutr):
        if args.star_sj:
            primary_arg = f"--star-sj {q(require(args.star_sj, 'STAR SJ.out.tab'))}"
        else:
            primary_arg = f"--junctions {q(require(args.junctions, 'junctions'))}"
        command = f"PYTHONPATH={q(root / 'src')} {q(sys.executable)} {q(root / 'scripts/normalize_rna_fusion_splice.py')} --sample-id {q(args.sample_id)} --profile {q(profile)} {primary_arg} --candidate-only"
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
    for key in ("prime_evidence", "bigmhc_evidence", "deepimmuno_evidence"):
        value = getattr(args, key)
        if value:
            lines.append(f"{key} = {q(require(value, key))}")
    if args.secondary_rna_events:
        lines.append(f"secondary_rna_events = {q(require(args.secondary_rna_events, 'secondary RNA events'))}")
        lines.append(f"secondary_sample_id = {q(args.secondary_sample_id)}")
        lines.append(f"secondary_identity_status = {q(args.secondary_identity_status.upper())}")
    for key in ("expression", "transcript_expression", "normal_junctions", "normal_expression", "normal_hla_ligands", "reference_proteome"):
        value = getattr(args, key)
        if value: lines.append(f"{key} = {q(require(value, key))}")
    Path(args.output).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__": raise SystemExit(main())
