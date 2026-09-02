#!/usr/bin/env python3
"""Generate a production manifest that reuses completed upstream tool results."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from neoag.cohort_rules import load_cohort_rule_contract, validate_cohort_rule_pair
from neoag.controlled_execution.io_utils import load_limited_yaml


STAR_INDEX_REQUIRED_FILES = ("Genome", "SA", "SAindex", "genomeParameters.txt")


def q(value) -> str:
    return json.dumps(str(value))


def require(path: str | None, label: str) -> str:
    if not path or not Path(path).exists():
        raise SystemExit(f"{label} missing: {path}")
    return str(Path(path).resolve())


def load_clinical_context(path: str | None) -> tuple[dict[str, str], str]:
    """Load explicit clinical metadata without inferring it from an analysis profile."""
    if not path:
        return {}, ""
    resolved = Path(require(path, "clinical context"))
    try:
        payload = load_limited_yaml(resolved)
    except Exception as exc:
        raise SystemExit(f"clinical context could not be parsed: {resolved}: {exc}") from exc
    if not isinstance(payload, dict):
        raise SystemExit(f"clinical context must contain a top-level mapping: {resolved}")
    context = {
        str(key): str(value).strip()
        for key, value in payload.items()
        if isinstance(value, (str, int, float, bool)) and str(value).strip()
    }
    clinical_keys = {"disease", "diagnosis", "disease_name", "cancer_type", "tumor_type"}
    if not clinical_keys.intersection(context):
        raise SystemExit(f"clinical context has no diagnosis field: {resolved}")
    return context, str(resolved)


def require_path_list(value: str | None, label: str) -> str:
    """Validate a comma-separated path list while preserving lane order."""
    paths = [item.strip() for item in str(value or "").split(",") if item.strip()]
    if not paths:
        raise SystemExit(f"{label} missing: {value}")
    resolved = [require(path, label) for path in paths]
    return ",".join(resolved)


def discover_vep_cache(reference_fasta: str, normal_junctions: str | None = None, explicit: str | None = None) -> str:
    candidates = [
        explicit or "",
        os.environ.get("NEOAG_VEP_CACHE", ""),
        str(Path(os.environ["OPEN_NEO_ASSET_ROOT"]) / "data/vep") if os.environ.get("OPEN_NEO_ASSET_ROOT") else "",
        str(Path(os.environ["OPEN_NEO_REFERENCE_ROOT"]) / "data/vep") if os.environ.get("OPEN_NEO_REFERENCE_ROOT") else "",
    ]
    for source in (reference_fasta, normal_junctions or ""):
        if source:
            candidates.extend(str(parent / "vep") for parent in Path(source).resolve().parents)
    for value in dict.fromkeys(item for item in candidates if item):
        path = Path(value).expanduser()
        if (path / "homo_sapiens").is_dir():
            return str(path.resolve())
        if path.name.endswith("_GRCh38") and path.parent.name == "homo_sapiens":
            return str(path.parent.parent.resolve())
    return ""


def discover_vep_plugins(
    reference_fasta: str,
    vep_cache: str,
    explicit: str | None = None,
) -> str:
    candidates = [
        explicit or "",
        os.environ.get("NEOAG_VEP_PLUGINS", ""),
        str(Path(os.environ["OPEN_NEO_ASSET_ROOT"]) / "work/vep_plugins")
        if os.environ.get("OPEN_NEO_ASSET_ROOT") else "",
        str(Path(os.environ["OPEN_NEO_REFERENCE_ROOT"]) / "work/vep_plugins")
        if os.environ.get("OPEN_NEO_REFERENCE_ROOT") else "",
    ]
    for source in (reference_fasta, vep_cache):
        if source:
            candidates.extend(str(parent / "work/vep_plugins") for parent in Path(source).resolve().parents)
    for value in dict.fromkeys(item for item in candidates if item):
        path = Path(value).expanduser()
        if all((path / filename).is_file() for filename in ("Wildtype.pm", "Frameshift.pm")):
            return str(path.resolve())
    return ""


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


def stage(lines, name, *, outputs, source="", command="", required=True, depends=None, cpus=None, memory_gb=None):
    lines += ["", f"[stages.{name}]", f"required = {str(required).lower()}"]
    if source: lines.append(f"source = {q(source)}")
    if depends: lines.append("depends_on = [" + ", ".join(q(value) for value in depends) + "]")
    if cpus is not None: lines.append(f"cpus = {int(cpus)}")
    if memory_gb is not None: lines.append(f"memory_gb = {float(memory_gb)}")
    if command: lines.append(f"command = {q(command)}")
    lines.append(f"[stages.{name}.outputs]")
    lines.extend(f"{key} = {q(value)}" for key, value in outputs.items())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-root", default=Path(__file__).resolve().parents[1])
    ap.add_argument("--sample-id", required=True); ap.add_argument("--outdir", required=True); ap.add_argument("--output", required=True)
    ap.add_argument("--profile", default="profiles/sarcoma_rna_supported_v2_provisional.toml")
    ap.add_argument(
        "--cohort-rule-set",
        default="configs/cohorts/dsrct_v1.toml",
        help="Versioned cohort contract that locks the ranking profile, evidence rules and report policy.",
    )
    ap.add_argument(
        "--clinical-context",
        help="Explicit YAML/JSON clinical context; its diagnosis is recorded separately from the cohort contract.",
    )
    ap.add_argument("--event-top-n", type=int, default=20)
    ap.add_argument("--candidate-top-n", type=int, default=100)
    ap.add_argument(
        "--evidence-consensus-rules",
        default="configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml",
        help="Independent Evidence-consensus R1-R4 rules; does not replace the weighted profile.",
    )
    ap.add_argument("--reference-fasta")
    ap.add_argument("--tumor-dna-bam", help="Explicit tumor DNA BAM; never inferred by sample order")
    ap.add_argument("--normal-dna-bam", help="Explicit matched-normal DNA BAM")
    ap.add_argument("--assay-type", choices=("WGS", "WES", "PANEL", "UNKNOWN"), default="UNKNOWN")
    ap.add_argument("--capture-bed", help="Required for formal WES/PANEL DNA-SV assessment")
    ap.add_argument("--genome-build", default="GRCh38")
    ap.add_argument("--sv-vcf", action="append", default=[], help="Existing DNA-SV VCF; repeatable")
    ap.add_argument("--sv-caller", action="append", default=[], help="Caller corresponding to each --sv-vcf")
    ap.add_argument("--skip-dna-sv", action="store_true", help="Explicitly leave DNA-SV unassessed")
    ap.add_argument("--sv-threads", type=int, default=8)
    ap.add_argument("--sv-memory-gb", type=float, default=48.0)
    ap.add_argument("--nextflow-executable", default="nextflow")
    ap.add_argument("--sv-nextflow-config", default="")
    ap.add_argument("--sv-nextflow-profile", default="")
    ap.add_argument("--bam-matcher-loci")
    ap.add_argument("--bam-matcher-reference")
    ap.add_argument("--skip-bam-matcher", action="store_true")
    ap.add_argument("--vep-cache", help="VEP cache root containing homo_sapiens/<version>_GRCh38")
    ap.add_argument("--vep-plugins", help="Directory containing Wildtype.pm and Frameshift.pm")
    ap.add_argument("--vep-bin", default=os.environ.get("NEOAG_VEP_BIN", ""), help="Installed VEP executable or wrapper")
    ap.add_argument("--hla-file"); ap.add_argument("--optitype"); ap.add_argument("--spechla-typing"); ap.add_argument("--hla-la"); ap.add_argument("--somatic-vcf")
    ap.add_argument("--tumor-sample-name", default="", help="Tumor sample column in a multi-sample somatic VCF")
    ap.add_argument("--normal-sample-name", default="", help="Matched-normal sample column in a multi-sample somatic VCF")
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
    ap.add_argument("--splice-rna-bam", help="Coordinate-sorted tumor RNA BAM used independently for splice read QC")
    ap.add_argument("--splice-star-sj", help="Matching tumor STAR SJ.out.tab used independently for splice read QC and PSI")
    ap.add_argument("--matched-normal-rna-bam", help="Optional matched-normal RNA BAM for coverage-aware junction exclusion")
    ap.add_argument("--matched-normal-star-sj", help="Optional matched-normal STAR SJ.out.tab matching --matched-normal-rna-bam")
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
    ap.add_argument(
        "--fusion-expressed-products",
        help="Exact-adjacency confirmed fusion transcript/ORF table used for formal peptide source-chain closure",
    )
    ap.add_argument(
        "--star-chimeric",
        action="append",
        default=[],
        help="STAR Chimeric.out.junction used for exact fusion breakpoint/read-name verification; repeatable",
    )
    ap.add_argument("--fusion-caller-root", action="append", default=[], help="Directory containing completed fusion caller outputs; repeatable")
    ap.add_argument("--normal-readthrough", help="Normal/read-through fusion background table for review")
    ap.add_argument("--junctions"); ap.add_argument("--star-sj"); ap.add_argument("--snaf"); ap.add_argument("--splicemutr"); ap.add_argument("--normal-junctions")
    ap.add_argument("--normal-junction-sqlite", help="Membership index built from --normal-junctions")
    ap.add_argument("--normal-expression"); ap.add_argument("--normal-hla-ligands"); ap.add_argument("--reference-proteome")
    ap.add_argument("--prime-evidence", help="Existing normalized PRIME evidence TSV to reuse")
    ap.add_argument("--bigmhc-evidence", help="Existing normalized BigMHC_IM evidence TSV to reuse")
    ap.add_argument("--deepimmuno-evidence", help="Existing normalized DeepImmuno evidence TSV to reuse")
    ap.add_argument("--netmhcstabpan-evidence", help="Existing normalized NetMHCstabpan evidence TSV to reuse")
    ap.add_argument("--netchop-executable", default="netChop"); ap.add_argument("--netchop-home", default="")
    ap.add_argument(
        "--skip-netmhcstabpan",
        action="store_true",
        help="Omit NetMHCstabpan when no approved local installation is available.",
    )
    args = ap.parse_args()
    if bool(args.rna_fastq1) != bool(args.rna_fastq2):
        raise SystemExit("--rna-fastq1 and --rna-fastq2 must be supplied together")
    if bool(args.matched_normal_rna_bam) != bool(args.matched_normal_star_sj):
        raise SystemExit("--matched-normal-rna-bam and --matched-normal-star-sj must be supplied together")
    if args.rna_threads < 1:
        raise SystemExit("--rna-threads must be a positive integer")
    if args.sv_threads < 1 or args.sv_memory_gb <= 0:
        raise SystemExit("--sv-threads and --sv-memory-gb must be positive")
    if bool(args.tumor_dna_bam) != bool(args.normal_dna_bam):
        raise SystemExit("--tumor-dna-bam and --normal-dna-bam must be supplied together")
    if bool(args.tumor_sample_name) != bool(args.normal_sample_name):
        raise SystemExit("--tumor-sample-name and --normal-sample-name must be supplied together")
    if args.sv_caller and len(args.sv_caller) != len(args.sv_vcf):
        raise SystemExit("--sv-caller must be repeated once for every --sv-vcf")
    if args.assay_type in {"WES", "PANEL"} and not args.skip_dna_sv and not args.capture_bed:
        raise SystemExit("WES/PANEL DNA-SV production requires --capture-bed; use --skip-dna-sv only for an explicit UNASSESSED result")
    if args.star_sjdb_overhang < 1:
        raise SystemExit("--star-sjdb-overhang must be a positive integer")
    if args.event_top_n < 1 or args.candidate_top_n < 1:
        raise SystemExit("--event-top-n and --candidate-top-n must be positive integers")
    splice_status_path = Path(args.output).resolve().parent / "splice_source_status.tsv"
    splice_status_path.parent.mkdir(parents=True, exist_ok=True)
    primary_splice = args.junctions or args.star_sj or ""
    splice_rows = [
        {
            "source": "junctions",
            "status": "AVAILABLE" if primary_splice and Path(primary_splice).exists() else "MISSING",
            "path": str(primary_splice),
            "reason": "primary tumor junction evidence" if primary_splice else "no --junctions/--star-sj input",
        },
        {
            "source": "snaf",
            "status": "AVAILABLE" if args.snaf and Path(args.snaf).exists() else "MISSING",
            "path": str(args.snaf or ""),
            "reason": "configured result" if args.snaf else "no completed sNAF result supplied",
        },
        {
            "source": "splicemutr",
            "status": "AVAILABLE" if args.splicemutr and Path(args.splicemutr).exists() else "MISSING",
            "path": str(args.splicemutr or ""),
            "reason": "configured result" if args.splicemutr else "no completed SpliceMutr result supplied",
        },
    ]
    with splice_status_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["source", "status", "path", "reason"], delimiter="\t")
        writer.writeheader()
        writer.writerows(splice_rows)
    if primary_splice and not (args.snaf or args.splicemutr):
        print(
            f"WARNING: splice junction evidence is present but neither sNAF nor SpliceMutr output was supplied; see {splice_status_path}",
            file=sys.stderr,
        )
    root = Path(args.project_root).resolve()
    profile_path = Path(args.profile)
    if not profile_path.is_absolute():
        profile_path = root / profile_path
    profile = require(str(profile_path), "ranking profile")
    consensus_rules_path = Path(args.evidence_consensus_rules)
    if not consensus_rules_path.is_absolute():
        consensus_rules_path = root / consensus_rules_path
    evidence_consensus_rules = require(str(consensus_rules_path), "evidence-consensus rules")
    contract_path = Path(args.cohort_rule_set)
    if not contract_path.is_absolute():
        contract_path = root / contract_path
    cohort_contract = load_cohort_rule_contract(require(str(contract_path), "cohort rule contract"))
    contract_mismatches = validate_cohort_rule_pair(
        cohort_contract,
        ranking_profile=profile,
        evidence_consensus_rules=evidence_consensus_rules,
    )
    if contract_mismatches:
        raise SystemExit(
            "Cohort rule contract mismatch; this run is not cohort-comparable: "
            + "; ".join(contract_mismatches)
        )
    clinical_context, clinical_context_source = load_clinical_context(args.clinical_context)
    reference_fasta = require(args.reference_fasta, "reference FASTA") if args.reference_fasta else ""
    tumor_dna_bam = require(args.tumor_dna_bam, "tumor DNA BAM") if args.tumor_dna_bam else ""
    normal_dna_bam = require(args.normal_dna_bam, "normal DNA BAM") if args.normal_dna_bam else ""
    capture_bed = require(args.capture_bed, "capture BED") if args.capture_bed else ""
    sv_vcfs = [require(value, "DNA-SV VCF") for value in args.sv_vcf]
    if args.assay_type in {"WGS", "WES", "PANEL"} and not args.skip_dna_sv and not (sv_vcfs or tumor_dna_bam):
        raise SystemExit("formal DNA-SV production requires existing --sv-vcf inputs or an explicit tumor/normal BAM pair")
    if tumor_dna_bam and not (args.tumor_sample_name and args.normal_sample_name):
        raise SystemExit("paired BAM DNA-SV calling requires explicit --tumor-sample-name and --normal-sample-name")
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
    if args.netmhcstabpan_evidence:
        netmhcstabpan_evidence = require(args.netmhcstabpan_evidence, "NetMHCstabpan evidence")
        presentation_predictors.insert(2, "netmhcstabpan")
    elif not args.skip_netmhcstabpan:
        netmhcstabpan_evidence = ""
        presentation_predictors.insert(2, "netmhcstabpan")
    else:
        netmhcstabpan_evidence = ""
        production_limitations.append("NETMHCSTABPAN_LOCAL_UNAVAILABLE")
    bam_pair = bool(tumor_dna_bam and normal_dna_bam)
    matcher_reference = require(args.bam_matcher_reference, "BAM-matcher reference") if args.bam_matcher_reference else reference_fasta
    matcher_loci = require(args.bam_matcher_loci, "BAM-matcher loci") if args.bam_matcher_loci else ""
    matcher_executable = shutil.which("bam-matcher")
    matcher_enabled = bam_pair and not args.skip_bam_matcher and bool(matcher_reference and matcher_loci and matcher_executable)
    if bam_pair and not matcher_enabled:
        production_limitations.append("BAM_MATCHER_UNASSESSED")
    if args.skip_dna_sv or not (sv_vcfs or bam_pair):
        production_limitations.append("DNA_SV_UNASSESSED")
    predictor_toml = "[" + ", ".join(q(tool) for tool in presentation_predictors) + "]"
    lines = [
        "# Generated from completed upstream tool results.",
        "[run]",
        f"sample_id = {q(args.sample_id)}",
        f"profile = {q(profile)}",
        f"outdir = {q(Path(args.outdir).resolve())}",
        f"hla_file = {q(hla)}",
        f"cohort_rule_set = {q(cohort_contract['path'])}",
        f"cohort_rule_set_id = {q(cohort_contract['id'])}",
        f"cohort_rule_set_version = {q(cohort_contract['version'])}",
        f"cohort_rule_set_sha256 = {q(cohort_contract['contract_sha256'])}",
        f"ranking_profile_sha256 = {q(cohort_contract['ranking_profile_sha256'])}",
        f"evidence_consensus_rules_sha256 = {q(cohort_contract['evidence_consensus_rules_sha256'])}",
        f"report_contract_version = {q(cohort_contract['report_contract_version'])}",
        f"release_audit_policy = {q(cohort_contract['release_audit_policy'])}",
        "cohort_comparability_required = true",
        "tools_stub = false",
        "immunogenicity_stub = false",
        f"presentation_predictors = {predictor_toml}",
        f"required_presentation_predictors = {predictor_toml}",
        'reports = "patient,technical"',
        f"event_top_n = {args.event_top_n}",
        f"candidate_top_n = {args.candidate_top_n}",
        f"netchop_executable = {q(args.netchop_executable)}",
        *( [f"clinical_context_source = {q(clinical_context_source)}"] if clinical_context_source else [] ),
        *( [f"clinical_context_schema = {q('open-neo-clinical-context-v1')}"] if clinical_context else [] ),
    ]
    if production_limitations:
        limitations_toml = "[" + ", ".join(q(item) for item in production_limitations) + "]"
        lines.append(f"production_limitations = {limitations_toml}")
    if args.netchop_home: lines.append(f"netchop_home = {q(args.netchop_home)}")
    if clinical_context:
        lines += ["", "[run.clinical_context]"]
        lines.extend(f"{key} = {q(value)}" for key, value in sorted(clinical_context.items()))
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

    pairing_dependency: list[str] = []
    sample_identity = ""
    if matcher_enabled:
        sample_identity = "{outdir}/sample_identity/bam_matcher/sample_identity.tsv"
        matcher_command = (
            f"bash {q(root / 'scripts/run_bam_matcher_pair.sh')} --bam1 {q(normal_dna_bam)} "
            f"--bam2 {q(tumor_dna_bam)} --reference {q(matcher_reference)} "
            f"--loci {q(matcher_loci)} --outdir {{outdir}}/sample_identity/bam_matcher"
        )
        stage(
            lines, "sample_identity_bam_matcher", command=matcher_command,
            outputs={"sample_identity": sample_identity, "raw_report": "{outdir}/sample_identity/bam_matcher/bam_matcher.short.tsv"},
            required=True, cpus=2, memory_gb=8,
        )
        pairing_dependency = ["sample_identity_bam_matcher"]

    rna_vaf = ""
    rna_bam = ""
    rna_bam_dependency: list[str] = []
    if args.rna_bam:
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
        rna_bam_dependency = ["rna_bam_input"]
    if args.rna_vaf:
        rna_vaf = require(args.rna_vaf, "RNA VAF evidence")
        stage(
            lines,
            "rna_alt_vaf_input",
            source="RNA_ALLELE_COUNTS",
            outputs={"rna_vaf": rna_vaf},
            required=False,
        )
    elif rna_bam:
        if not args.somatic_vcf:
            raise SystemExit("--rna-bam requires --somatic-vcf for RNA allele counting")
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
        rna_bam_dependency = ["rna_star_alignment"]
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
    dna_sv_stage = ""
    dna_sv_events = ""
    if not args.skip_dna_sv and (sv_vcfs or bam_pair):
        if not reference_fasta or not args.gencode_gtf:
            raise SystemExit("DNA-SV production requires --reference-fasta and --gencode-gtf")
        gencode_gtf = require(args.gencode_gtf, "GENCODE GTF")
        dna_sv_root = "{outdir}/branches/dna_sv"
        dna_sv_events = f"{dna_sv_root}/sv/sv_events.full.tsv"
        if sv_vcfs:
            callers = f" --callers {' '.join(q(value) for value in args.sv_caller)}" if args.sv_caller else ""
            sample_names = ""
            if args.tumor_sample_name and args.normal_sample_name:
                sample_names = f" --tumor-sample-name {q(args.tumor_sample_name)} --normal-sample-name {q(args.normal_sample_name)}"
            elif args.tumor_sample_name or args.normal_sample_name:
                raise SystemExit("DNA-SV VCF parsing requires both --tumor-sample-name and --normal-sample-name")
            mode = "sv-build-raw-wes" if args.assay_type in {"WES", "PANEL"} else "sv-build-raw"
            capture = f" --capture-bed {q(capture_bed)}" if capture_bed else ""
            command = (
                f"PYTHONPATH={q(root / 'src')} {q(sys.executable)} -m neoag.cli {mode} "
                f"--sample-id {q(args.sample_id)} --profile {q('sv_wes_phase1_5' if capture else 'sv_wgs_phase1')} "
                f"--sv-vcf {' '.join(q(value) for value in sv_vcfs)}{callers} "
                f"--reference-fasta {q(reference_fasta)} --gencode-gtf {q(gencode_gtf)} "
                f"--hla {q(hla)} --outdir {dna_sv_root} --genome-build {q(args.genome_build)}"
                f"{sample_names}{capture}"
            )
        else:
            config_arg = f" -c {q(require(args.sv_nextflow_config, 'SV Nextflow config'))}" if args.sv_nextflow_config else ""
            profile_arg = f" -profile {q(args.sv_nextflow_profile)}" if args.sv_nextflow_profile else ""
            capture_args = f" --wes_mode true --capture_bed {q(capture_bed)}" if capture_bed else ""
            sample_args = ""
            if args.tumor_sample_name and args.normal_sample_name:
                sample_args = f" --tumor_sample_name {q(args.tumor_sample_name)} --normal_sample_name {q(args.normal_sample_name)}"
            command = (
                f"{q(args.nextflow_executable)} run {q(root / 'workflows/sv_phase1_wgs.nf')} -resume{config_arg}{profile_arg} "
                f"--sample_id {q(args.sample_id)} --tumor_bam {q(tumor_dna_bam)} --normal_bam {q(normal_dna_bam)} "
                f"--reference_fasta {q(reference_fasta)} --gencode_gtf {q(gencode_gtf)} --hla {q(hla)} "
                f"--outdir {dna_sv_root} --run_scoring false --threads {args.sv_threads} "
                f"--genome_build {q(args.genome_build)}{sample_args}{capture_args}"
            )
        dna_sv_stage = "dna_sv_discovery"
        stage(
            lines, dna_sv_stage, source="DNA_SV", command=command,
            outputs={
                "raw_events": f"{dna_sv_root}/parsed/raw_events.tsv",
                "raw_peptides": f"{dna_sv_root}/parsed/raw_peptides.tsv",
                "sv_events_full": dna_sv_events,
                "sv_event_to_peptide": f"{dna_sv_root}/sv/sv_event_to_peptide.tsv",
            },
            depends=pairing_dependency + hla_dependency,
            cpus=args.sv_threads * (3 if not sv_vcfs else 1), memory_gb=args.sv_memory_gb,
        )
        candidate_stages.append(dna_sv_stage)
    if args.somatic_vcf:
        reference_env = f"NEOAG_REFERENCE_FASTA={q(reference_fasta)} " if reference_fasta else ""
        vep_cache = discover_vep_cache(reference_fasta, args.normal_junctions, args.vep_cache)
        vep_cache_arg = f" --vep-cache {q(vep_cache)}" if vep_cache and Path(vep_cache).is_dir() else ""
        vep_plugins = discover_vep_plugins(reference_fasta, vep_cache, args.vep_plugins)
        vep_plugins_arg = f" --vep-plugins {q(vep_plugins)}" if vep_plugins else ""
        vep_bin_arg = f" --vep-bin {q(args.vep_bin)}" if args.vep_bin else ""
        sample_role_args = ""
        if args.tumor_sample_name:
            sample_role_args += f" --tumor-sample-name {q(args.tumor_sample_name)}"
        if args.normal_sample_name:
            sample_role_args += f" --normal-sample-name {q(args.normal_sample_name)}"
        command = f"{reference_env}PYTHONPATH={q(root / 'src')} {q(sys.executable)} {q(root / 'scripts/run_candidate_upstream.py')} --mode snv --input {q(require(args.somatic_vcf, 'somatic VCF'))} --hla-file {q(hla)} --sample-id {q(args.sample_id)}{sample_role_args} --outdir {{outdir}}/branches/snv{vep_cache_arg}{vep_plugins_arg}{vep_bin_arg}"
        stage(lines, "snv_indel_candidates", source="SNV_INDEL", command=command, outputs={"raw_events": "{outdir}/branches/snv/parsed/raw_events.tsv", "raw_peptides": "{outdir}/branches/snv/parsed/raw_peptides.tsv"}, depends=list(dict.fromkeys(pairing_dependency + hla_dependency)))
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
        for chimeric_path in args.star_chimeric:
            union_args += ["--star-chimeric", q(require(chimeric_path, "STAR Chimeric.out.junction"))]
        star_junction_source = args.star_sj or args.junctions
        if star_junction_source:
            chimeric = Path(star_junction_source).with_name("Chimeric.out.junction")
            if chimeric.is_file() and chimeric.stat().st_size > 0 and str(chimeric) not in args.star_chimeric:
                union_args += ["--star-chimeric", q(str(chimeric))]
        if rna_bam:
            union_args += ["--rna-bam", q(rna_bam)]
        union_args += ["--samtools", q(args.samtools_executable)]
        if args.fusioncatcher: union_args += ["--fusioncatcher", q(require(args.fusioncatcher, "FusionCatcher"))]
        if args.jaffal: union_args += ["--jaffal", q(require(args.jaffal, "JAFFAL"))]
        if args.fusion_expressed_products:
            union_args += ["--fusion-expressed-products", q(require(args.fusion_expressed_products, "confirmed fusion expressed products"))]
        union_args += ["--genome-build", q(args.genome_build)]
        for caller_root in args.fusion_caller_root:
            union_args += ["--caller-root", q(require(caller_root, "fusion caller result root"))]
        command = f"PYTHONPATH={q(root / 'src')} {q(sys.executable)} {q(root / 'scripts/build_fusion_caller_union.py')} --sample-id {q(args.sample_id)} --profile {q(profile)} --hla-file {q(hla)} {' '.join(union_args)} --outdir {{outdir}}/branches/fusion/intermediates"
        stage(lines, "fusion_candidates", command=command, outputs={"raw_events": "{outdir}/branches/fusion/intermediates/raw_events.tsv", "raw_peptides": "{outdir}/branches/fusion/intermediates/raw_peptides.tsv", "fusion_union": "{outdir}/branches/fusion/intermediates/fusion_caller_union.tsv", "fusion_caller_availability": "{outdir}/branches/fusion/intermediates/fusion_caller_availability.tsv", "fusion_consensus": "{outdir}/branches/fusion/intermediates/fusion_consensus.tsv", "junction_verification": "{outdir}/branches/fusion/intermediates/junction_read_verification.tsv", "fusion_peptide_origin_chain": "{outdir}/branches/fusion/intermediates/fusion_peptide_origin_chain.tsv", "fusion_orf_completion_queue": "{outdir}/branches/fusion/intermediates/fusion_orf_completion_queue.tsv", "diagnostic_fusion_rescue": "{outdir}/branches/fusion/intermediates/diagnostic_fusion_rescue.tsv"}, depends=list(dict.fromkeys(pairing_dependency + hla_dependency + rna_bam_dependency)))
        link_command = (
            f"PYTHONPATH={q(root / 'src')} {q(sys.executable)} {q(root / 'scripts/link_dna_sv_rna_fusions.py')} "
            f"--fusion-events {{outdir}}/branches/fusion/intermediates/raw_events.tsv "
            f"--fusion-peptides {{outdir}}/branches/fusion/intermediates/raw_peptides.tsv "
            f"--fusion-union {{outdir}}/branches/fusion/intermediates/fusion_caller_union.tsv "
            f"--fusion-consensus {{outdir}}/branches/fusion/intermediates/fusion_consensus.tsv "
            + (f"--sv-events {dna_sv_events} " if dna_sv_events else "")
            + "--outdir {outdir}/branches/fusion/dna_sv_linked"
        )
        link_dependencies = ["fusion_candidates"] + ([dna_sv_stage] if dna_sv_stage else [])
        stage(
            lines, "fusion_dna_sv_link", source="FusionCallerUnion", command=link_command,
            outputs={
                "raw_events": "{outdir}/branches/fusion/dna_sv_linked/raw_events.tsv",
                "raw_peptides": "{outdir}/branches/fusion/dna_sv_linked/raw_peptides.tsv",
                "dna_sv_rna_links": "{outdir}/branches/fusion/dna_sv_linked/dna_sv_rna_fusion_links.tsv",
                "fusion_consensus": "{outdir}/branches/fusion/dna_sv_linked/fusion_consensus.tsv",
            }, depends=link_dependencies,
        )
        candidate_stages.append("fusion_dna_sv_link")
        review = (
            "mkdir -p {outdir}/branches/fusion/consensus && "
            "cp {outdir}/branches/fusion/dna_sv_linked/fusion_consensus.tsv "
            "{outdir}/branches/fusion/consensus/fusion_consensus.tsv && "
            "cp {outdir}/branches/fusion/intermediates/fusion_caller_availability.tsv "
            "{outdir}/branches/fusion/consensus/fusion_caller_availability.tsv"
        )
        stage(lines, "fusion_cross_validation", command=review, outputs={"fusion_consensus": "{outdir}/branches/fusion/consensus/fusion_consensus.tsv", "fusion_caller_availability": "{outdir}/branches/fusion/consensus/fusion_caller_availability.tsv"}, required=True, depends=["fusion_dna_sv_link"])
    generated_star_sj = "{outdir}/rna/star/SJ.out.tab" if args.rna_fastq1 and args.rna_fastq2 else ""
    splice_rna_bam = require(args.splice_rna_bam, "splice RNA BAM") if args.splice_rna_bam else rna_bam
    splice_star_sj = (
        require(args.splice_star_sj, "splice STAR SJ.out.tab") if args.splice_star_sj
        else require(args.star_sj, "STAR SJ.out.tab") if args.star_sj
        else generated_star_sj
    )
    matched_normal_rna_bam = require(args.matched_normal_rna_bam, "matched-normal RNA BAM") if args.matched_normal_rna_bam else ""
    matched_normal_star_sj = require(args.matched_normal_star_sj, "matched-normal STAR SJ.out.tab") if args.matched_normal_star_sj else ""
    normal_junction_sqlite = ""
    if args.normal_junction_sqlite:
        normal_junction_sqlite = require(args.normal_junction_sqlite, "normal junction SQLite index")
    elif args.normal_junctions:
        adjacent_index = Path(str(args.normal_junctions) + ".sqlite")
        if adjacent_index.is_file() and adjacent_index.stat().st_size > 0:
            normal_junction_sqlite = str(adjacent_index.resolve())

    if args.junctions and args.star_sj:
        raise SystemExit("Use only one primary splice evidence input: --junctions or --star-sj")
    primary_splice_source = args.junctions or args.star_sj or generated_star_sj
    if primary_splice_source and (args.snaf or args.splicemutr):
        if args.star_sj:
            primary_arg = f"--star-sj {q(require(args.star_sj, 'STAR SJ.out.tab'))}"
        elif args.junctions:
            primary_arg = f"--junctions {q(require(args.junctions, 'junctions'))}"
        else:
            primary_arg = f"--star-sj {q(generated_star_sj)}"
        command = f"PYTHONPATH={q(root / 'src')} {q(sys.executable)} {q(root / 'scripts/normalize_rna_fusion_splice.py')} --sample-id {q(args.sample_id)} --profile {q(profile)} {primary_arg} --candidate-only"
        if args.gencode_gtf:
            command += f" --annotation-gtf {q(require(args.gencode_gtf, 'matched GENCODE GTF'))}"
        if args.snaf: command += f" --snaf {require(args.snaf, 'SNAF')}"
        if args.splicemutr: command += f" --splicemutr {require(args.splicemutr, 'SpliceMutr')}"
        if args.normal_junctions: command += f" --normal-junctions {require(args.normal_junctions, 'normal junctions')}"
        command += " --outdir {outdir}/branches/splice/intermediates"
        splicemutr_path = Path(args.splicemutr) if args.splicemutr else None
        splice_peptides_for_filter = "{outdir}/branches/splice/intermediates/raw_peptides.tsv"
        if splicemutr_path and splicemutr_path.is_dir():
            corrected_glob = str(splicemutr_path / "formed_transcripts" / "**" / "*_data_splicemutr_cp_corrected.txt")
            command += (
                f" && PYTHONPATH={q(root / 'src')} {q(sys.executable)} "
                f"{q(root / 'scripts/rebuild_splice_origins_from_splicemutr.py')}"
                f" --sample-id {q(args.sample_id)} --genome-build GRCh38"
                " --candidates {outdir}/branches/splice/intermediates/raw_peptides.tsv"
                f" --splicemutr-glob {q(corrected_glob)}"
                " --outdir {outdir}/branches/splice/intermediates/formal_origins"
            )
            splice_peptides_for_filter = "{outdir}/branches/splice/intermediates/formal_origins/raw_peptides.formal_origins.tsv"
        splice_events_for_filter = "{outdir}/branches/splice/intermediates/raw_events.tsv"
        if splice_rna_bam and splice_star_sj:
            qc_outdir = "{outdir}/branches/splice/intermediates/star_bam_qc"
            qc_command = (
                f" && PYTHONPATH={q(root / 'src')} {q(sys.executable)} "
                f"{q(root / 'scripts/build_splice_junction_qc_from_star_bam.py')}"
                f" --events {splice_events_for_filter} --peptides {splice_peptides_for_filter}"
                f" --star-sj {q(splice_star_sj)} --rna-bam {q(splice_rna_bam)}"
                " --splice-consensus {outdir}/branches/splice/intermediates/splice_consensus.tsv"
                f" --samtools {q(args.samtools_executable)} --outdir {qc_outdir}"
            )
            if normal_junction_sqlite:
                qc_command += f" --normal-junction-sqlite {q(normal_junction_sqlite)}"
            if matched_normal_rna_bam:
                qc_command += (
                    f" --matched-normal-rna-bam {q(matched_normal_rna_bam)}"
                    f" --matched-normal-star-sj {q(matched_normal_star_sj)}"
                )
            command += qc_command
            splice_events_for_filter = f"{qc_outdir}/raw_events.enriched.tsv"
            splice_peptides_for_filter = f"{qc_outdir}/raw_peptides.enriched.tsv"
        command += (
            f" && PYTHONPATH={q(root / 'src')} {q(sys.executable)} {q(root / 'scripts/filter_splice_production_candidates.py')}"
            f" --events {splice_events_for_filter}"
            f" --peptides {splice_peptides_for_filter}"
            " --consensus {outdir}/branches/splice/intermediates/splice_consensus.tsv"
            " --outdir {outdir}/branches/splice/production_selected"
            " --min-length 8 --max-length 12 --max-source-binding-rank 2.0"
        )
        stage(lines, "splice_candidates", source="SpliceConsensus", command=command, outputs={"raw_events": "{outdir}/branches/splice/production_selected/raw_events.tsv", "raw_peptides": "{outdir}/branches/splice/production_selected/raw_peptides.tsv", "production_filter_summary": "{outdir}/branches/splice/production_selected/production_filter_summary.json", "formal_origin_summary": "{outdir}/branches/splice/intermediates/formal_origins/rebuild_summary.json", "junction_read_qc": "{outdir}/branches/splice/intermediates/star_bam_qc/splice_junction_qc.enriched.tsv"}, depends=list(dict.fromkeys(hla_dependency + rna_bam_dependency)))
        candidate_stages.append("splice_candidates")
    if not candidate_stages: raise SystemExit("At least one SNV/Fusion/Splice candidate source is required")
    lines += [
        "", "[evidence]", f"purity = {q(purity)}", f"cnv = {q(cnv)}",
        f"hla_loh = {q(hla_loh)}",
        f"evidence_consensus_rules = {q(evidence_consensus_rules)}",
    ]
    if sample_identity:
        lines.append(f"sample_identity = {q(sample_identity)}")
    if rna_vaf:
        lines.append(f"rna_vaf = {q(rna_vaf)}")
    if netmhcstabpan_evidence:
        lines.append(f"netmhcstabpan = {q(netmhcstabpan_evidence)}")
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
