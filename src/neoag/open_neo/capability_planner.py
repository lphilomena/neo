from __future__ import annotations

import json
import os
import re
import sys
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from neoag.controlled_execution.io_utils import load_limited_yaml, write_json, write_tsv
from neoag.production_runner import load_production_manifest

from .rna_fusion_splice_profile import generate_rna_fusion_splice_manifest


TOOL_EXECUTABLES = {
    "bwa": "bwa",
    "samtools": "samtools",
    "gatk": "gatk",
    "bam_matcher": "bam-matcher",
    "optitype": "optitype",
    "hla_la": "HLA-LA.pl",
    "spechla": "spechla",
    "facets": "runFACETS.R",
    "sequenza": "sequenza-utils",
    "purple": "purple",
    "ascat": "ascat.R",
    "lohhla": "LOHHLA",
    "star": "STAR",
    "salmon": "salmon",
    "rsem": "rsem-calculate-expression",
    "easyfuse": "easyfuse",
    "star_fusion": "STAR-Fusion",
    "fusioncatcher": "fusioncatcher-neoag",
    "arriba": "arriba",
    "regtools": "regtools",
    "snaf": "snaf",
    "splicemutr": "splicemutr",
    "netmhcpan": "netMHCpan",
    "mhcflurry": "mhcflurry-predict",
    "netmhcstabpan": "netMHCstabpan",
    "netchop": "netChop",
    "prime": "PRIME",
    "bigmhc_im": "bigmhc_predict",
    "deepimmuno": "deepimmuno-cnn.py",
}

REFERENCE_ALIASES = {
    "reference_fasta": ("reference_fasta", "reference.fasta", "fasta"),
    "gencode_gtf": ("gencode_gtf", "reference.gtf", "gencode.gtf", "gtf"),
    "vep_cache": ("vep_cache", "vep.cache"),
    "vep_plugins": ("vep_plugins", "vep.plugins", "vep_plugins_dir"),
    "facets_snp_vcf": ("facets_snp_vcf", "common_snp", "facets.vcf"),
    "bam_matcher_loci": ("bam_matcher_loci", "bam-matcher.loci", "sample_identity_vcf"),
    "sequenza_fasta": ("sequenza_fasta", "sequenza_reference_fasta", "sequenza.reference_fasta"),
    "sequenza_gc_wiggle": ("sequenza_gc_wiggle", "gc_wiggle"),
    "hla_la_graph": ("hla_la_graph", "prg_mhc", "hla-la.graph"),
    "spechla_db": ("spechla_db", "spechla.db"),
    "lohhla_reference": ("lohhla_reference",),
    "purple_reference": ("purple_reference",),
    "star_index": ("star_index",),
    "ctat_genome_lib": ("ctat_genome_lib", "ctat"),
    "easyfuse_ref": ("easyfuse_ref",),
    "fusioncatcher_ref": ("fusioncatcher_ref",),
    "salmon_index": ("salmon_index",),
    "tx2gene": ("tx2gene", "salmon_tx2gene"),
    "rsem_reference": ("rsem_reference",),
    "normal_expression": ("normal_expression",),
    "normal_hla_ligands": ("normal_hla_ligands", "normal_ligandome"),
    "reference_proteome": ("reference_proteome", "normal_proteome"),
    "normal_junctions": ("normal_junctions",),
}

INPUT_REFERENCE_KEYS = set(REFERENCE_ALIASES) | {"normal_readthrough", "snaf_workflow", "snaf_db", "snaf_python", "altanalyze_image", "splicemutr_workflow", "splicemutr_config", "splicemutr_samples"}
RAW_INPUT_KEYS = {
    "tumor_dna_bam", "normal_dna_bam", "tumor_rna_bam",
    "tumor_dna_fastq", "normal_dna_fastq", "tumor_rna_fastq",
}

DEFAULT_PRODUCTION_PROFILE = "profiles/sarcoma_rna_supported_v2_provisional.toml"
DEFAULT_EVIDENCE_CONSENSUS_RULES = "configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml"
WRAPPER_REQUIRED_PRESENTATION = ["netmhcpan", "mhcflurry", "netmhcstabpan", "netchop"]


@dataclass
class CapabilityDecision:
    domain: str
    tool: str
    status: str
    reason: str
    stage: str = ""
    input_compatible: str = "yes"
    required: str = "false"
    executable: str = ""
    references: str = ""


@dataclass
class AutomaticPlan:
    status: str
    manifest: str
    decisions: list[CapabilityDecision] = field(default_factory=list)
    missing_required: list[str] = field(default_factory=list)
    selected_tools: list[str] = field(default_factory=list)
    routes: list[str] = field(default_factory=list)
    outputs: dict[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "manifest": self.manifest,
            "decisions": [row.__dict__ for row in self.decisions],
            "missing_required": self.missing_required,
            "selected_tools": self.selected_tools,
            "routes": self.routes,
            "outputs": self.outputs,
        }


def is_automatic_production_candidate(inputs: dict[str, Any]) -> bool:
    return not inputs.get("production_manifest") and any(inputs.get(key) for key in RAW_INPUT_KEYS)


def _pair(value: Any) -> tuple[str, str] | None:
    values = value if isinstance(value, list) else ([value] if value else [])
    values = [str(item) for item in values if str(item)]
    return (values[0], values[1]) if len(values) >= 2 else None


def _flatten_paths(data: Any, prefix: str = "") -> dict[str, str]:
    paths: dict[str, str] = {}
    if isinstance(data, dict):
        if isinstance(data.get("path"), str):
            paths[prefix.lower()] = os.path.expandvars(os.path.expanduser(str(data["path"])))
        for key, value in data.items():
            if key in {"path", "sha256"}:
                continue
            child = f"{prefix}.{key}" if prefix else str(key)
            paths.update(_flatten_paths(value, child))
    return paths


def _manifest_tools(path: str | Path | None) -> dict[str, dict[str, Any]]:
    if not path or not Path(path).is_file():
        return {}
    data = load_limited_yaml(path)
    values = data.get("tools") if isinstance(data, dict) else {}
    return {str(key).lower().replace("-", "_"): value for key, value in (values or {}).items() if isinstance(value, dict)}


def _resolve_references(inputs: dict[str, Any], manifest: str | Path | None) -> dict[str, str]:
    result = {key: str(inputs.get(key) or "") for key in INPUT_REFERENCE_KEYS if inputs.get(key)}
    data = load_limited_yaml(manifest) if manifest and Path(manifest).is_file() else {}
    flattened = _flatten_paths(data)
    reference_root = os.environ.get("OPEN_NEO_REFERENCE_ROOT") or os.environ.get("NEOAG_REF_BUNDLE") or ""
    if reference_root:
        flattened = {
            key: value.replace("${OPEN_NEO_REFERENCE_ROOT}", reference_root).replace("$OPEN_NEO_REFERENCE_ROOT", reference_root)
            for key, value in flattened.items()
        }
    for name, aliases in REFERENCE_ALIASES.items():
        if result.get(name):
            continue
        matches = [(key, value) for key, value in flattened.items() if any(alias in key for alias in aliases)]
        if matches:
            result[name] = sorted(matches, key=lambda item: (len(item[0]), item[0]))[0][1]
    return result


def _tool_info(name: str, tools: dict[str, dict[str, Any]]) -> tuple[bool, str, str]:
    cfg = tools.get(name) or tools.get(name.replace("_", "-")) or {}
    executable = str(cfg.get("executable") or TOOL_EXECUTABLES.get(name) or name)
    path = executable if Path(executable).is_file() else shutil.which(executable)
    template = str(cfg.get("command_template") or "")
    return bool(path or template), str(path or executable), template


def _toml(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    if isinstance(value, list):
        return "[" + ", ".join(_toml(item) for item in value) + "]"
    return json.dumps(str(value), ensure_ascii=False)


def _resource_profile(name: str, command: str, thread_ceiling: int) -> tuple[int, float, str]:
    """Return conservative scheduler reservations and align explicit thread flags."""
    key = name.lower()
    profiles = (
        (("easyfuse",), 8, 64.0),
        (("rna_alignment", "star_fusion", "fusioncatcher", "arriba"), 8, 48.0),
        (("rsem",), 4, 32.0),
        (("rna_expression", "salmon"), 4, 16.0),
        (("splicemutr", "snaf"), 6, 48.0),
        (("sequenza", "facets"), 4, 24.0),
        (("purple",), 6, 48.0),
        (("hla_loh",), 6, 40.0),
        (("hla", "optitype", "spechla"), 4, 20.0),
        (("alignment", "mutect", "bwa"), 8, 32.0),
        (("fastq_qc",), 4, 8.0),
        (("fastq_merge",), 2, 8.0),
    )
    cpus, memory_gb = 2, 8.0
    for tokens, requested_cpus, requested_memory in profiles:
        if any(token in key for token in tokens):
            cpus, memory_gb = requested_cpus, requested_memory
            break
    cpus = max(1, min(cpus, max(1, thread_ceiling)))
    command = re.sub(r"(--threads(?:=|\s+))\d+", lambda match: match.group(1) + str(cpus), command)
    command = re.sub(r"(--cores(?:=|\s+))\d+", lambda match: match.group(1) + str(cpus), command)
    return cpus, memory_gb, command


def _render_manifest(path: Path, run: dict[str, Any], stages: dict[str, dict[str, Any]], evidence: dict[str, Any]) -> None:
    lines = ["# Generated by Open-Neo capability-aware planner.", "# Review the DAG before approved execution.", "", "[run]"]
    lines.extend(f"{key} = {_toml(value)}" for key, value in run.items() if key != "required_tool_groups")
    for domain, rule in (run.get("required_tool_groups") or {}).items():
        lines += ["", f"[run.required_tool_groups.{domain}]"]
        lines.extend(f"{key} = {_toml(value)}" for key, value in rule.items())
    for name, spec in stages.items():
        lines += ["", f"[stages.{name}]"]
        for key in ("required", "source", "depends_on", "cpus", "memory_gb", "command"):
            if key in spec and spec[key] is not None and spec[key] != "":
                lines.append(f"{key} = {_toml(spec[key])}")
        lines.append(f"[stages.{name}.outputs]")
        lines.extend(f"{key} = {_toml(value)}" for key, value in (spec.get("outputs") or {}).items())
    lines += ["", "[evidence]"]
    lines.extend(f"{key} = {_toml(value)}" for key, value in evidence.items() if value)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _template_command(template: str, values: dict[str, Any]) -> str:
    if not template:
        return ""
    try:
        return template.format_map({key: str(value) for key, value in values.items()})
    except KeyError:
        return ""


def _write_hla_file(path: Path, alleles: list[str]) -> str:
    normalized = []
    for raw in alleles:
        value = str(raw).strip()
        if value and not value.upper().startswith("HLA-"):
            value = "HLA-" + value
        if value and value not in normalized:
            normalized.append(value)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(normalized) + "\n", encoding="utf-8")
    return str(path)


def build_automatic_production_plan(
    inputs: dict[str, Any],
    manifest_path: str | Path,
    *,
    project_root: str | Path,
    outdir: str | Path,
    tools_manifest: str | Path | None = None,
    reference_manifest: str | Path | None = None,
    policy: str = "all-available",
) -> AutomaticPlan:
    root = Path(project_root).resolve()
    target = Path(outdir).resolve()
    manifest = Path(manifest_path).resolve()
    tools = _manifest_tools(tools_manifest)
    refs = _resolve_references(inputs, reference_manifest)
    decisions: list[CapabilityDecision] = []
    stages: dict[str, dict[str, Any]] = {}
    evidence: dict[str, Any] = {}
    selected: list[str] = []
    missing_required: list[str] = []
    routes: list[str] = []
    sample_id = str(inputs.get("sample_id") or inputs.get("case_id") or "SAMPLE001")
    tumor_id = str(inputs.get("tumor_sample_id") or inputs.get("tumor_id") or f"{sample_id}_T")
    normal_id = str(inputs.get("normal_sample_id") or inputs.get("normal_id") or f"{sample_id}_N")
    threads = int(inputs.get("threads") or inputs.get("rna_threads") or 16)
    samtools_executable = str(inputs.get("samtools_executable") or "")

    minimal_tools = {"bwa", "samtools", "gatk", "bam_matcher", "optitype", "facets", "netmhcpan", "mhcflurry", "star", "salmon", "rsem", "easyfuse", "regtools"}
    balanced_exclusions = {"ascat", "deepimmuno", "netmhcstabpan"}

    def permitted(tool: str) -> bool:
        if policy == "minimal":
            return tool in minimal_tools
        if policy == "balanced":
            return tool not in balanced_exclusions
        return True

    def decide(domain: str, tool: str, status: str, reason: str, *, stage: str = "", required: bool = False, executable: str = "", references: list[str] | None = None, compatible: bool = True) -> None:
        decisions.append(CapabilityDecision(
            domain, tool, status, reason, stage,
            "yes" if compatible else "no", "true" if required else "false",
            executable, ",".join(references or []),
        ))
        if status == "SELECTED" and tool not in selected:
            selected.append(tool)

    def add_stage(name: str, *, command: str, outputs: dict[str, str], required: bool = False, source: str = "", depends: list[str] | None = None) -> None:
        python_bin_dir = Path(sys.executable).resolve().parent
        cpus, memory_gb, command = _resource_profile(name, command, threads)
        runtime_command = f'PATH="{python_bin_dir}:$PATH" {command}'
        stages[name] = {
            "required": required, "command": runtime_command, "outputs": outputs,
            "cpus": cpus, "memory_gb": memory_gb,
        }
        if source:
            stages[name]["source"] = source
        if depends:
            stages[name]["depends_on"] = depends

    tumor_bam = str(inputs.get("tumor_dna_bam") or "")
    normal_bam = str(inputs.get("normal_dna_bam") or "")
    tumor_dna_fastq = _pair(inputs.get("tumor_dna_fastq"))
    normal_dna_fastq = _pair(inputs.get("normal_dna_fastq"))
    tumor_rna_fastq = _pair(inputs.get("tumor_rna_fastq"))
    if not tumor_rna_fastq:
        r1_values = inputs.get("rna_fastq1")
        r2_values = inputs.get("rna_fastq2")
        r1_values = r1_values if isinstance(r1_values, list) else ([r1_values] if r1_values else [])
        r2_values = r2_values if isinstance(r2_values, list) else ([r2_values] if r2_values else [])
        if r1_values and r2_values:
            tumor_rna_fastq = (str(r1_values[0]), str(r2_values[0]))
    tumor_rna_bam = str(inputs.get("tumor_rna_bam") or "")
    somatic_vcf = str(inputs.get("somatic_vcf") or "")
    hla_file = str(inputs.get("hla_file") or "")
    hla_alleles = [str(value) for value in inputs.get("hla_alleles") or []]
    invalid_hla_file = hla_file if hla_file and not Path(hla_file).is_file() else ""
    if invalid_hla_file:
        hla_file = ""

    if not tumor_bam and tumor_dna_fastq:
        available, executable, _ = _tool_info("bwa", tools)
        samtools, _, _ = _tool_info("samtools", tools)
        if available and samtools and refs.get("reference_fasta"):
            tumor_bam = "{outdir}/dna/tumor/tumor.sorted.bam"
            add_stage(
                "tumor_dna_alignment", required=True,
                command=f"bash {root / 'scripts/run_dna_fastq_to_bam.sh'} --fastq1 {tumor_dna_fastq[0]} --fastq2 {tumor_dna_fastq[1]} --reference {refs['reference_fasta']} --sample-id tumor --threads {threads} --outdir {{outdir}}/dna/tumor",
                outputs={"tumor_bam": tumor_bam, "tumor_bai": tumor_bam + ".bai"},
            )
            decide("dna_alignment", "bwa+samtools", "SELECTED", "tumor DNA FASTQ pair and GRCh38 reference are available", stage="tumor_dna_alignment", required=True, executable=executable, references=["reference_fasta"])
        else:
            missing_required.append("tumor_dna_alignment")
            decide("dna_alignment", "bwa+samtools", "BLOCKED", "tumor DNA FASTQ needs BWA, samtools and reference_fasta", required=True)
    if not normal_bam and normal_dna_fastq:
        available, executable, _ = _tool_info("bwa", tools)
        samtools, _, _ = _tool_info("samtools", tools)
        if available and samtools and refs.get("reference_fasta"):
            normal_bam = "{outdir}/dna/normal/normal.sorted.bam"
            add_stage(
                "normal_dna_alignment", required=True,
                command=f"bash {root / 'scripts/run_dna_fastq_to_bam.sh'} --fastq1 {normal_dna_fastq[0]} --fastq2 {normal_dna_fastq[1]} --reference {refs['reference_fasta']} --sample-id normal --threads {threads} --outdir {{outdir}}/dna/normal",
                outputs={"normal_bam": normal_bam, "normal_bai": normal_bam + ".bai"},
            )
            decide("dna_alignment", "bwa+samtools", "SELECTED", "normal DNA FASTQ pair and GRCh38 reference are available", stage="normal_dna_alignment", required=True, executable=executable, references=["reference_fasta"])
        else:
            missing_required.append("normal_dna_alignment")
            decide("dna_alignment", "bwa+samtools", "BLOCKED", "normal DNA FASTQ needs BWA, samtools and reference_fasta", required=True)

    alignment_deps = [name for name in ("tumor_dna_alignment", "normal_dna_alignment") if name in stages]
    for label, bam in (("tumor", tumor_bam), ("normal", normal_bam)):
        if not bam or "{outdir}" in bam or Path(bam + ".bai").exists() or Path(bam).with_suffix(".bai").exists():
            continue
        samtools_available, samtools_exe, _ = _tool_info("samtools", tools)
        if samtools_available:
            stage = f"{label}_bam_index"
            add_stage(
                stage,
                required=True,
                command=f"{samtools_exe or 'samtools'} index -@ {threads} {bam}",
                outputs={f"{label}_bai": bam + ".bai"},
            )
            alignment_deps.append(stage)
            decide("dna_alignment", "samtools", "SELECTED", f"{label} BAM index is missing and will be generated", stage=stage, required=True, executable=samtools_exe)
        else:
            missing_required.append(f"{label}_bam_index")
            decide("dna_alignment", "samtools", "BLOCKED", f"{label} BAM index is missing and samtools is unavailable", required=True)
    bam_pair = bool(tumor_bam and normal_bam)
    paired_analysis_deps = list(alignment_deps)
    if bam_pair:
        matcher_available, matcher_exe, matcher_template = _tool_info("bam_matcher", tools)
        matcher_loci = refs.get("bam_matcher_loci")
        values = {
            "tumor_bam": tumor_bam, "normal_bam": normal_bam,
            "reference_fasta": refs.get("reference_fasta", ""),
            "bam_matcher_loci": matcher_loci or "",
            "outdir": "{outdir}/sample_identity/bam_matcher",
            "sample_id": sample_id,
        }
        matcher_command = _template_command(matcher_template, values)
        if matcher_command and matcher_exe:
            matcher_command = f"PATH={Path(matcher_exe).parent}:$PATH {matcher_command}"
        if matcher_available and refs.get("reference_fasta") and matcher_loci:
            matcher_command = matcher_command or (
                f"bash {root / 'scripts/run_bam_matcher_pair.sh'} --bam1 {normal_bam} --bam2 {tumor_bam} "
                f"--reference {refs['reference_fasta']} --loci {matcher_loci} --outdir {{outdir}}/sample_identity/bam_matcher"
            )
            add_stage(
                "sample_identity_bam_matcher", required=True, command=matcher_command,
                outputs={
                    "sample_identity": "{outdir}/sample_identity/bam_matcher/sample_identity.tsv",
                    "raw_report": "{outdir}/sample_identity/bam_matcher/bam_matcher.short.tsv",
                },
                depends=alignment_deps,
            )
            paired_analysis_deps.append("sample_identity_bam_matcher")
            decide(
                "sample_identity", "bam_matcher", "SELECTED",
                "tumor-normal BAM pair, GRCh38 FASTA and identity SNP panel are available",
                stage="sample_identity_bam_matcher", required=True, executable=matcher_exe,
                references=["reference_fasta", "bam_matcher_loci"],
            )
        elif matcher_available and not matcher_loci:
            decide(
                "sample_identity", "bam_matcher", "UNASSESSED",
                "BAM-matcher is installed but bam_matcher_loci is not declared; bundled hg19 loci are forbidden for GRCh38",
                executable=matcher_exe, references=["bam_matcher_loci"],
            )
        else:
            decide(
                "sample_identity", "bam_matcher", "UNAVAILABLE",
                "BAM-matcher is optional; install it and declare a build-matched identity SNP panel to enable genotype identity QC",
                executable=matcher_exe, references=["reference_fasta", "bam_matcher_loci"],
            )
    typing_input_bam = normal_bam or tumor_bam

    if invalid_hla_file:
        decide("hla_typing", "provided_hla", "INVALID", f"provided HLA file is missing: {invalid_hla_file}", required=True)
    if hla_alleles and not hla_file:
        hla_file = _write_hla_file(manifest.with_name("provided_hla.txt"), hla_alleles)
    if hla_file:
        decide("hla_typing", "provided_hla", "REUSED", "HLA alleles/file supplied by user", required=True)
    else:
        hla_results: list[str] = []
        if normal_dna_fastq or tumor_dna_fastq:
            pair = normal_dna_fastq or tumor_dna_fastq
            available, executable, template = _tool_info("optitype", tools)
            command = _template_command(template, {"fastq1": pair[0], "fastq2": pair[1], "sample_id": sample_id, "outdir": "{outdir}/hla/optitype", "threads": threads})
            if available and not command:
                command = f"optitype run -i {pair[0]} -i {pair[1]} --dna -o {{outdir}}/hla/optitype --solver cbc --threads {threads}"
            if available and command:
                add_stage("hla_optitype", command=command, outputs={
                    "result_dir": "{outdir}/hla/optitype",
                    "complete_marker": "{outdir}/hla/optitype/.complete",
                }, depends=alignment_deps)
                hla_results.append("{outdir}/hla/optitype")
                decide("hla_typing", "optitype", "SELECTED", "paired DNA FASTQ is compatible", stage="hla_optitype", executable=executable)
            else:
                decide("hla_typing", "optitype", "UNAVAILABLE", "paired DNA FASTQ is present but OptiType is unavailable")
        if typing_input_bam:
            optitype_available, optitype_executable, optitype_template = _tool_info("optitype", tools)
            if "hla_optitype" not in stages and permitted("optitype"):
                optitype_command = _template_command(optitype_template, {
                    "bam": typing_input_bam, "sample_id": normal_id if normal_bam else tumor_id,
                    "outdir": "{outdir}/hla/optitype", "threads": threads,
                })
                if optitype_available and not optitype_command:
                    samtools_arg = f" --samtools-bin {samtools_executable}" if samtools_executable else ""
                    optitype_command = (
                        f"bash {root / 'scripts/run_optitype_sample.sh'} --bam {typing_input_bam} "
                        f"--sample-id {normal_id if normal_bam else tumor_id} --threads {threads} "
                        "--outdir {outdir}/hla/optitype "
                        f"--optitype-bin {optitype_executable}{samtools_arg}"
                    )
                if optitype_available and optitype_command:
                    add_stage("hla_optitype", command=optitype_command, outputs={
                        "result_dir": "{outdir}/hla/optitype",
                        "complete_marker": "{outdir}/hla/optitype/.complete",
                    }, depends=alignment_deps)
                    hla_results.append("{outdir}/hla/optitype")
                    decide("hla_typing", "optitype", "SELECTED", "repository-owned BAM-to-HLA-read runner is available", stage="hla_optitype", executable=optitype_executable)
                else:
                    decide("hla_typing", "optitype", "UNAVAILABLE", "OptiType is unavailable for BAM typing", executable=optitype_executable)
            for tool, graph_ref in (("hla_la", "hla_la_graph"), ("spechla", "spechla_db")):
                available, executable, template = _tool_info(tool, tools)
                if tool == "spechla" and tools.get("spechla") and (root / "scripts/run_spechla_sample.sh").is_file():
                    available = True
                if not permitted(tool):
                    decide("hla_typing", tool, "POLICY_SKIPPED", f"excluded by {policy} policy", executable=executable)
                    continue
                values = {"bam": typing_input_bam, "sample_id": sample_id, "outdir": f"{{outdir}}/hla/{tool}", "threads": threads, **refs}
                command = _template_command(template, values)
                if available and not command and refs.get(graph_ref):
                    runner = "run_hla_la_sample.sh" if tool == "hla_la" else "run_spechla_sample.sh"
                    command = (
                        f"bash {root / 'scripts' / runner} --bam {typing_input_bam} "
                        f"--sample-id {normal_id if normal_bam else tumor_id} --threads {threads} "
                        f"--outdir {{outdir}}/hla/{tool}"
                    )
                if available and command and (not graph_ref or refs.get(graph_ref)):
                    stage = f"hla_{tool}"
                    outputs = {"result_dir": f"{{outdir}}/hla/{tool}"}
                    if tool == "spechla":
                        outputs["complete_marker"] = "{outdir}/hla/spechla/.complete"
                    add_stage(stage, command=command, outputs=outputs, depends=alignment_deps)
                    hla_results.append(f"{{outdir}}/hla/{tool}")
                    decide("hla_typing", tool, "SELECTED", "BAM input and a validated template or repository-owned runner are available", stage=stage, executable=executable, references=[graph_ref])
                elif available and not command:
                    decide("hla_typing", tool, "TEMPLATE_REQUIRED", "tool exists but no validated sample-level command_template is declared", executable=executable, references=[graph_ref])
                else:
                    decide("hla_typing", tool, "UNAVAILABLE", "tool or required reference is unavailable", executable=executable, references=[graph_ref], compatible=bool(typing_input_bam))
        hla_file = "{outdir}/hla/consensus/hla_consensus.txt"
        compare_args = " ".join(f"--result-dir {value}" for value in hla_results)
        if hla_results:
            add_stage(
                "hla_consensus", required=True, depends=[name for name in stages if name.startswith("hla_")],
                command=f"PYTHONPATH={root / 'src'} python -m neoag.agent_skills.hla_typing_compare {compare_args} --sample-id {sample_id} --outdir {{outdir}}/hla/consensus && PYTHONPATH={root / 'src'} python {root / 'scripts/hla_consensus_to_file.py'} --consensus {{outdir}}/hla/consensus/hla_typing_consensus.tsv --output {hla_file}",
                outputs={"hla_file": hla_file, "hla_consensus": "{outdir}/hla/consensus/hla_typing_consensus.tsv"},
            )
        else:
            missing_required.append("hla_typing")
            decide("hla_typing", "consensus", "BLOCKED", "no supplied HLA and no executable HLA caller with a validated runner", required=True)

    if not somatic_vcf and bam_pair:
        available, executable, template = _tool_info("gatk", tools)
        values = {"tumor_bam": tumor_bam, "normal_bam": normal_bam, "sample_id": sample_id, "outdir": "{outdir}/dna/mutect2", "threads": threads, **refs}
        command = _template_command(template, values)
        if available and refs.get("reference_fasta"):
            command = command or f"bash {root / 'scripts/run_mutect2_tumor_normal.sh'} --tumor-bam {tumor_bam} --normal-bam {normal_bam} --reference {refs['reference_fasta']} --sample-id {sample_id} --outdir {{outdir}}/dna/mutect2"
            somatic_vcf = "{outdir}/dna/mutect2/somatic.pass.vcf.gz"
            add_stage("somatic_variant_calling", required=True, command=command, outputs={"somatic_vcf": somatic_vcf, "somatic_vcf_index": somatic_vcf + ".tbi"}, depends=paired_analysis_deps)
            decide("somatic_variant", "gatk_mutect2", "SELECTED", "tumor-normal BAM pair and reference are available", stage="somatic_variant_calling", required=True, executable=executable, references=["reference_fasta"])
        else:
            missing_required.append("somatic_vcf")
            decide("somatic_variant", "gatk_mutect2", "BLOCKED", "no somatic VCF and Mutect2/reference are unavailable", required=True)
    elif somatic_vcf:
        decide("somatic_variant", "provided_vcf", "REUSED", "somatic VCF supplied")

    hla_dependency = ["hla_consensus"] if "hla_consensus" in stages else []
    if somatic_vcf and hla_file:
        deps = list(dict.fromkeys(
            (["somatic_variant_calling"] if "somatic_variant_calling" in stages else [])
            + paired_analysis_deps + hla_dependency
        ))
        command = (
            f"PYTHONPATH={root / 'src'} python {root / 'scripts/run_candidate_upstream.py'} --mode snv --input {somatic_vcf} "
            f"--hla-file {hla_file} --sample-id {sample_id} --outdir {{outdir}}/branches/snv/upstream"
        )
        for flag, value in (
            ("--reference-fasta", refs.get("reference_fasta")),
            ("--vep-cache", refs.get("vep_cache")),
            ("--vep-plugins", refs.get("vep_plugins") or os.environ.get("NEOAG_VEP_PLUGINS")),
            ("--normal-proteome", refs.get("reference_proteome")),
        ):
            if value:
                command += f" {flag} {value}"
        add_stage(
            "snv_indel_candidates", required=True, source="SNV_INDEL", command=command,
            outputs={"raw_events": "{outdir}/branches/snv/upstream/parsed/raw_events.tsv", "raw_peptides": "{outdir}/branches/snv/upstream/parsed/raw_peptides.tsv"}, depends=deps,
        )
        routes.append("snv_indel")

    if bam_pair:
        purity_dirs: list[str] = []
        facets_available, facets_exe, _ = _tool_info("facets", tools)
        facets_ref = refs.get("facets_snp_vcf")
        if facets_available and facets_ref:
            command = (
                f"FACETS_MODE=omni2p5 FACETS_SNP_VCF={facets_ref} FACETS_CVAL_PRE=50 "
                f"FACETS_CVAL_PROC=300 FACETS_MIN_NHET=10 FACETS_TARGET_ROWS=1000000 "
                f"PATIENT_ID={sample_id} TUMOR_BAM={tumor_bam} NORMAL_BAM={normal_bam} "
                f"OUTDIR={{outdir}}/purity/facets bash {root / 'scripts/run_facets_sample.sh'}"
            )
            add_stage(
                "purity_facets",
                command=command,
                outputs={
                    "result_dir": "{outdir}/purity/facets",
                    "summary": "{outdir}/purity/facets/facets_omni2p5_summary.tsv",
                },
                depends=paired_analysis_deps,
            )
            purity_dirs.append("{outdir}/purity/facets")
            decide("purity_cnv", "facets", "SELECTED", "BAM pair and omni2p5 reference are available; robust parameters 50/300/10/1000000 are fixed", stage="purity_facets", executable=facets_exe, references=["facets_snp_vcf"])
        else:
            decide("purity_cnv", "facets", "UNAVAILABLE", "FACETS executable or SNP reference is missing", references=["facets_snp_vcf"])

        seq_available, seq_exe, _ = _tool_info("sequenza", tools)
        if not permitted("sequenza"):
            decide("purity_cnv", "sequenza", "POLICY_SKIPPED", f"excluded by {policy} policy", executable=seq_exe)
        elif seq_available and (refs.get("sequenza_fasta") or refs.get("reference_fasta")) and refs.get("sequenza_gc_wiggle"):
            seq_ref = refs.get("sequenza_fasta") or refs["reference_fasta"]
            gc = refs["sequenza_gc_wiggle"]
            command = f"SAMPLE_ID={sample_id} TUMOR_BAM={tumor_bam} NORMAL_BAM={normal_bam} REF_FASTA={seq_ref} GC_WIGGLE={gc} OUTDIR={{outdir}}/purity/sequenza BIN_WINDOW=${{SEQUENZA_BIN_WINDOW:-500}} bash {root / 'scripts/run_sequenza_sample_by_chrom.sh'}"
            summary = f"{{outdir}}/purity/sequenza/sequenza_fit/{sample_id}.sequenza_summary.tsv"
            add_stage("purity_sequenza", command=command, outputs={"purity": summary}, depends=paired_analysis_deps)
            purity_dirs.append("{outdir}/purity/sequenza")
            decide("purity_cnv", "sequenza", "SELECTED", "BAM pair and chr-prefixed Sequenza reference FASTA are available; resume-safe binning defaults to SEQUENZA_BIN_WINDOW or 500", stage="purity_sequenza", executable=seq_exe, references=["sequenza_fasta", "sequenza_gc_wiggle"])
        else:
            decide("purity_cnv", "sequenza", "UNAVAILABLE", "Sequenza or reference FASTA is missing", references=["reference_fasta"])

        for tool in ("purple", "ascat"):
            available, executable, template = _tool_info(tool, tools)
            if not permitted(tool):
                decide("purity_cnv", tool, "POLICY_SKIPPED", f"excluded by {policy} policy", executable=executable)
                continue
            command = _template_command(template, {
                "tumor_bam": tumor_bam, "normal_bam": normal_bam,
                "sample_id": sample_id, "tumor_id": tumor_id, "normal_id": normal_id,
                "outdir": f"{{outdir}}/purity/{tool}", **refs,
            })
            if tool == "purple" and available and refs.get("purple_reference") and not command:
                required_hmf_tools = ("amber", "cobalt", "purple")

                def valid_hmf_env(candidate: Path) -> bool:
                    return all((candidate / "bin" / name).is_file() for name in required_hmf_tools)

                executable_path = Path(executable).expanduser()
                inferred_root = executable_path.parent
                inferred_root = inferred_root.parent if inferred_root.name == "bin" else inferred_root
                candidates = []
                if os.environ.get("HMF_ENV"):
                    candidates.append(Path(os.environ["HMF_ENV"]).expanduser())
                if os.environ.get("HMFTOOLS_HOME"):
                    candidates.append(Path(os.environ["HMFTOOLS_HOME"]).expanduser() / ".conda")
                candidates.extend([
                    inferred_root,
                    inferred_root / ".conda",
                    inferred_root / "tools" / "HMFTOOLS" / ".conda",
                    root / "tools" / "HMFTOOLS" / ".conda",
                ])
                hmf_env = next((candidate for candidate in candidates if valid_hmf_env(candidate)), candidates[0])
                command = (
                    f"HMF_ENV={hmf_env} "
                    f"HMFTOOLS_REFERENCE_ROOT={refs['purple_reference']} "
                    f"HMFTOOLS_REFERENCE_FASTA={refs.get('reference_fasta', '')} "
                    f"bash {root / 'scripts/run_purple_sample.sh'} --sample-id {sample_id} "
                    f"--tumor-id {tumor_id} --normal-id {normal_id} --tumor-bam {tumor_bam} "
                    f"--normal-bam {normal_bam} --threads {threads} --outdir {{outdir}}/purity/purple"
                )
            if available and command:
                stage = f"purity_{tool}"
                outputs = {"result_dir": f"{{outdir}}/purity/{tool}"}
                if tool == "purple":
                    outputs["complete_marker"] = "{outdir}/purity/purple/.complete"
                add_stage(stage, command=command, outputs=outputs, depends=paired_analysis_deps)
                purity_dirs.append(f"{{outdir}}/purity/{tool}")
                decide("purity_cnv", tool, "SELECTED", "validated command_template or repository-owned runner is available", stage=stage, executable=executable, references=["purple_reference"] if tool == "purple" else [])
            elif available:
                decide("purity_cnv", tool, "TEMPLATE_REQUIRED", "installed tool has no validated sample-level command_template", executable=executable)
            else:
                decide("purity_cnv", tool, "UNAVAILABLE", "tool unavailable")

        if purity_dirs:
            args_dirs = " ".join(f"--result-dir {value}" for value in purity_dirs)
            dependencies = [name for name in stages if name.startswith("purity_")]
            add_stage(
                "purity_consensus",
                command=f"PYTHONPATH={root / 'src'} python -m neoag.agent_skills.purity_cnv_review {args_dirs} --sample-id {sample_id} --outdir {{outdir}}/purity/consensus",
                outputs={
                    "result_dir": "{outdir}/purity/consensus",
                    "purity": "{outdir}/purity/consensus/recommended_purity.tsv",
                    "cnv": "{outdir}/purity/consensus/recommended_cnv_segments.tsv",
                    "purity_consensus": "{outdir}/purity/consensus/purity_cnv_consensus.tsv",
                    "purity_tool_results": "{outdir}/purity/consensus/purity_cnv_tool_summary.tsv",
                    "hla_6p21_consensus": "{outdir}/purity/consensus/hla_6p21_cnv_consensus.tsv",
                },
                depends=dependencies,
            )
            evidence["purity"] = "{outdir}/purity/consensus/recommended_purity.tsv"
            evidence["cnv"] = "{outdir}/purity/consensus/recommended_cnv_segments.tsv"

        loh_available, loh_exe, _ = _tool_info("lohhla", tools)
        spechla_available, spechla_exe, _ = _tool_info("spechla", tools)
        if tools.get("spechla") and (root / "scripts/run_spechla_sample.sh").is_file():
            spechla_available = True
        selected_loh_tools: list[str] = []
        if permitted("lohhla") and loh_available:
            selected_loh_tools.append("lohhla")
        elif not permitted("lohhla"):
            decide("hla_loh", "lohhla", "POLICY_SKIPPED", f"excluded by {policy} policy", executable=loh_exe)
        else:
            decide("hla_loh", "lohhla", "UNAVAILABLE", "LOHHLA is unavailable", executable=loh_exe)
        if permitted("spechla") and spechla_available and refs.get("spechla_db"):
            selected_loh_tools.append("spechla")
        elif not permitted("spechla"):
            decide("hla_loh", "spechla", "POLICY_SKIPPED", f"excluded by {policy} policy", executable=spechla_exe)
        else:
            decide("hla_loh", "spechla", "UNAVAILABLE", "SpecHLA or spechla_db is unavailable", executable=spechla_exe, references=["spechla_db"])

        if selected_loh_tools and hla_file and purity_dirs:
            stage = "hla_loh_multi_tool"
            deps = list(dict.fromkeys(paired_analysis_deps + hla_dependency + ["purity_consensus"]))
            command = (
                f"PYTHONPATH={root / 'src'} bash {root / 'scripts/run_hla_loh_multi_tool.sh'} "
                f"--sample-id {sample_id} --tumor-id {tumor_id} --normal-id {normal_id} "
                f"--tumor-bam {tumor_bam} --normal-bam {normal_bam} --hla-file {hla_file} "
                f"--purity-tsv {{outdir}}/purity/consensus/recommended_purity.tsv "
                f"--tools {','.join(selected_loh_tools)} --threads {threads} --outdir {{outdir}}/hla_loh"
            )
            outputs = {
                "result_dir": "{outdir}/hla_loh",
                "hla_loh": "{outdir}/hla_loh/recommended_hla_loh.tsv",
                "hla_loh_consensus": "{outdir}/hla_loh/hla_loh_consensus.tsv",
            }
            if "lohhla" in selected_loh_tools:
                outputs["lohhla_hla_loh"] = "{outdir}/hla_loh/lohhla/hla_loh.tsv"
            if "spechla" in selected_loh_tools:
                outputs["spechla_hla_loh"] = "{outdir}/hla_loh/spechla/hla_loh.tsv"
            add_stage(stage, command=command, outputs=outputs, depends=deps)
            evidence["hla_loh"] = "{outdir}/hla_loh/recommended_hla_loh.tsv"
            for tool in selected_loh_tools:
                decide(
                    "hla_loh", tool, "SELECTED",
                    "scheduled with recommended cross-tool purity/ploidy and allele-level consensus",
                    stage=stage, executable=loh_exe if tool == "lohhla" else spechla_exe,
                    references=["lohhla_reference"] if tool == "lohhla" else ["spechla_db"],
                )
        elif selected_loh_tools:
            decide("hla_loh", "consensus", "BLOCKED", "HLA file or recommended purity/CNV stage is unavailable")

    if tumor_rna_fastq:
        routes.extend(["rna_expression", "fusion", "splice"])
        rna_inputs = {**inputs, **refs, "hla_file": hla_file, "sample_id": sample_id}
        submanifest = manifest.with_name("rna_fusion_splice.generated.toml")
        generated = generate_rna_fusion_splice_manifest(rna_inputs, submanifest, project_root=root, outdir=target)
        rna_cfg = load_production_manifest(submanifest)
        for name, spec in (rna_cfg.get("stages") or {}).items():
            if name not in stages:
                stages[name] = spec
        evidence.update(rna_cfg.get("evidence") or {})
        for name, tool in (
            ("rna_expression", "salmon"), ("rsem_expression_crosscheck", "rsem"),
            ("easyfuse_discovery", "easyfuse"),
            ("star_fusion_discovery", "star_fusion"), ("fusioncatcher_discovery", "fusioncatcher"),
            ("arriba_discovery", "arriba"),
            ("junction_extraction", "regtools"), ("snaf_discovery", "snaf"),
            ("splicemutr_discovery", "splicemutr"),
        ):
            spec = stages.get(name) or {}
            command_text = str(spec.get("command") or "")
            decision_tool = tool
            if name == "rna_expression" and "run_rsem_fastq_to_tpm.sh" in command_text:
                decision_tool = "rsem"
            if not permitted(decision_tool):
                spec["command"] = ""
                decide("rna" if name in {"rna_expression", "rsem_expression_crosscheck"} else "fusion_splice", decision_tool, "POLICY_SKIPPED", f"excluded by {policy} policy")
                continue
            status = "SELECTED" if spec.get("command") else "UNASSESSED"
            decide("rna" if name in {"rna_expression", "rsem_expression_crosscheck"} else "fusion_splice", decision_tool, status, "RNA FASTQ profile selected" if status == "SELECTED" else "runner/reference/workflow not configured", stage=name if status == "SELECTED" else "")
        if somatic_vcf and "rna_alignment" in stages:
            add_stage(
                "rna_alt_vaf",
                command=(
                    f"PYTHONPATH={root / 'src'} python {root / 'scripts/rna_allele_counts_pysam.py'} "
                    f"--somatic-vcf {somatic_vcf} --rna-bam {{outdir}}/rna/star/Aligned.sortedByCoord.out.bam "
                    "--output-tsv {outdir}/rna/rna_alt_vaf.tsv"
                ),
                outputs={"rna_vaf": "{outdir}/rna/rna_alt_vaf.tsv"},
                depends=["rna_alignment"],
            )
            evidence["rna_vaf"] = "{outdir}/rna/rna_alt_vaf.tsv"
            decide("rna", "rna_alt_vaf", "SELECTED", "STAR RNA BAM and supplied somatic VCF are available", stage="rna_alt_vaf")
    elif tumor_rna_bam:
        routes.append("rna_evidence")
        if somatic_vcf:
            add_stage(
                "rna_alt_vaf", command=f"PYTHONPATH={root / 'src'} python {root / 'scripts/rna_allele_counts_pysam.py'} --somatic-vcf {somatic_vcf} --rna-bam {tumor_rna_bam} --output-tsv {{outdir}}/rna/rna_alt_vaf.tsv",
                outputs={"rna_vaf": "{outdir}/rna/rna_alt_vaf.tsv"}, depends=["somatic_variant_calling"] if "somatic_variant_calling" in stages else None,
            )
            evidence["rna_vaf"] = "{outdir}/rna/rna_alt_vaf.tsv"
            decide("rna", "rna_alt_vaf", "SELECTED", "RNA BAM and somatic VCF are available", stage="rna_alt_vaf")
        reg_available, reg_exe, _ = _tool_info("regtools", tools)
        if reg_available:
            add_stage("junction_extraction", command=f"bash {root / 'scripts/run_regtools_junctions.sh'} --bam {tumor_rna_bam} --sample-id {sample_id} --out {{outdir}}/rna/regtools_junctions.tsv", outputs={"junctions": "{outdir}/rna/regtools_junctions.tsv"})
            evidence["rna_junction_tsv"] = "{outdir}/rna/regtools_junctions.tsv"
            decide("rna", "regtools", "SELECTED", "RNA BAM is available", stage="junction_extraction", executable=reg_exe)

    for direct_key, mode, source in (
        ("fusion_tsv", "fusion", "Fusion"),
        ("splice_junction_tsv", "splice", "Splice"),
    ):
        value = str(inputs.get(direct_key) or "")
        stage_name = f"provided_{mode}_candidates"
        if value and stage_name not in stages and hla_file:
            command = f"PYTHONPATH={root / 'src'} python {root / 'scripts/run_candidate_upstream.py'} --mode {mode} --input {value} --hla-file {hla_file} --sample-id {sample_id} --outdir {{outdir}}/branches/{mode}/provided"
            add_stage(stage_name, command=command, source=source, outputs={"raw_events": f"{{outdir}}/branches/{mode}/provided/parsed/raw_events.tsv", "raw_peptides": f"{{outdir}}/branches/{mode}/provided/parsed/raw_peptides.tsv"}, depends=hla_dependency)
            routes.append(mode)

    presentation_predictors = []
    required_predictors = []
    for tool in (*WRAPPER_REQUIRED_PRESENTATION, "prime", "bigmhc_im", "deepimmuno"):
        available, executable, _ = _tool_info(tool, tools)
        if not permitted(tool):
            decide("presentation", tool, "POLICY_SKIPPED", f"excluded by {policy} policy", executable=executable)
            continue
        if available:
            presentation_predictors.append(tool)
            if tool in WRAPPER_REQUIRED_PRESENTATION:
                required_predictors.append(tool)
            decide("presentation", tool, "SELECTED", "predictor available and compatible with peptide-HLA candidates", executable=executable)
        else:
            decide("presentation", tool, "UNAVAILABLE", "predictor not available")
    for tool in WRAPPER_REQUIRED_PRESENTATION:
        if tool not in presentation_predictors:
            missing_required.append(tool)

    for stage_name, stage_spec in stages.items():
        if "cpus" not in stage_spec or "memory_gb" not in stage_spec:
            cpus, memory_gb, command = _resource_profile(
                stage_name, str(stage_spec.get("command") or ""), threads,
            )
            stage_spec["cpus"] = cpus
            stage_spec["memory_gb"] = memory_gb
            stage_spec["command"] = command

    expected_sources = [str(spec.get("source")) for spec in stages.values() if spec.get("source")]
    profile_value = str(inputs.get("profile") or (root / DEFAULT_PRODUCTION_PROFILE))
    evidence_consensus_rules = str(inputs.get("evidence_consensus_rules") or (root / DEFAULT_EVIDENCE_CONSENSUS_RULES))
    run = {
        "sample_id": sample_id,
        "profile": profile_value,
        "outdir": str(target),
        "hla_file": hla_file,
        "hla_alleles": hla_alleles,
        "tools_stub": False,
        "immunogenicity_stub": False,
        "expected_peptide_sources": list(dict.fromkeys(expected_sources)),
        "presentation_predictors": presentation_predictors or WRAPPER_REQUIRED_PRESENTATION,
        "required_presentation_predictors": required_predictors or WRAPPER_REQUIRED_PRESENTATION,
        "reports": "patient,technical",
        "automatic_tool_policy": policy,
        "total_cpus": int(inputs.get("total_cpus") or threads),
        "total_memory_gb": float(inputs.get("total_memory_gb") or 0),
        "max_parallel_stages": int(inputs.get("max_parallel_stages") or 3),
        "default_stage_cpus": 1,
        "default_stage_memory_gb": 4.0,
    }
    if bam_pair:
        run["required_tool_groups"] = {
            "purity_cnv": {"tools": ["facets", "sequenza", "purple"], "min_successful": 2, "require_all_declared": True},
            "hla_loh": {"tools": ["lohhla", "spechla"], "min_successful": 2, "require_all_declared": True},
        }
    if tumor_rna_fastq:
        run["rna_profile"] = "rna_fusion_splice_v1"
    evidence["evidence_consensus_rules"] = evidence_consensus_rules
    for key in ("normal_expression", "normal_hla_ligands", "reference_proteome", "normal_junctions"):
        if refs.get(key):
            evidence[key] = refs[key]
    _render_manifest(manifest, run, stages, evidence)

    status = "BLOCKED" if missing_required else ("PARTIAL" if any(row.status in {"UNASSESSED", "TEMPLATE_REQUIRED"} for row in decisions) else "READY")
    outputs = {
        "automatic_production_manifest": str(manifest),
        "capability_decisions": str(manifest.parent / "capability_decisions.tsv"),
        "capability_plan": str(manifest.parent / "capability_plan.json"),
    }
    write_tsv(outputs["capability_decisions"], [row.__dict__ for row in decisions])
    plan = AutomaticPlan(status, str(manifest), decisions, sorted(set(missing_required)), selected, list(dict.fromkeys(routes)), outputs)
    write_json(outputs["capability_plan"], plan.to_dict())
    return plan
