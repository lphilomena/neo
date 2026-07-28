from __future__ import annotations

import json
import gzip
import os
import shlex
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from neoag.controlled_execution.io_utils import load_limited_yaml, write_json, write_tsv


TOOL_ENV_PATHS = {
    "hla_la": ("HLALA_BIN", "HLA_LA_BIN"),
    "spechla": ("SPECHLA_BIN",),
    "purple": ("PURPLE_BIN",),
    "ascat": ("ASCAT_BIN",),
    "snaf": ("SNAF_BIN",),
    "splicemutr": ("SPLICEMUTR_BIN",),
    "bam_matcher": ("BAM_MATCHER_BIN",),
    "netmhcpan": ("NEOAG_NETMHCPAN_BIN",),
    "netmhcstabpan": ("NETMHCSTABPAN_BIN",),
}

REFERENCE_ENV_PATHS = {
    "reference_fasta": ("NEOAG_REFERENCE_FASTA",),
    "gencode_gtf": ("NEOAG_GENCODE_GTF",),
    "vep_cache": ("NEOAG_VEP_CACHE",),
    "hla_la_graph": ("HLALA_GRAPH", "HLA_LA_GRAPH"),
    "spechla_db": ("SPECHLA_DB",),
    "lohhla_reference": ("LOHHLA_HOME",),
    "facets_snp_vcf": ("FACETS_SNP_VCF", "COMMON_SNP_VCF"),
    "sequenza_gc_wiggle": ("GC_WIGGLE",),
    "purple_reference": ("PURPLE_REFERENCE",),
    "bam_matcher_loci": ("BAM_MATCHER_LOCI",),
    "snaf_workflow": ("SNAF_WORKFLOW",),
    "snaf_db": ("NEOAG_SNAF_DB", "SNAF_DB"),
    "snaf_python": ("SNAF_PYTHON",),
    "splicemutr_workflow": ("SPLICEMUTR_WORKFLOW",),
    "ascat_loci": ("ASCAT_LOCI",),
    "ascat_alleles": ("ASCAT_ALLELES",),
}

REFERENCE_CANDIDATES = {
    "reference_fasta": ("data/ref/hg38/Homo_sapiens_assembly38.fasta", "GRCh38/fasta/GRCh38.fa"),
    "gencode_gtf": ("data/ref/hg38/gencode.gtf", "GRCh38/gencode/gencode.v44.annotation.gtf.gz"),
    "vep_cache": ("data/vep/homo_sapiens/105_GRCh38", "GRCh38/vep/vep_115_GRCh38"),
    "hla_la_graph": ("data/hla/PRG_MHC_GRCh38_withIMGT", "hla/hla-la/graph/PRG_MHC_GRCh38_withIMGT"),
    "spechla_db": ("data/hla/spechla_db", "hla/spechla/db"),
    "facets_snp_vcf": ("data/facets/reference/common_snp.hg38.vcf.gz", "cnv/facets/common_snp.vcf.gz"),
    "sequenza_gc_wiggle": ("data/sequenza/reference/Homo_sapiens.GRCh38.dna.primary_assembly.chr.gc50.wig.gz",),
    "purple_reference": ("data/hmf/purple_reference", "cnv/purple"),
    "bam_matcher_loci": ("data/sample_identity/bam_matcher.common_snps.hg38.vcf", "GRCh38/sample_identity/bam_matcher.common_snps.vcf"),
    "snaf_db": ("data/snaf/reference/data", "splice/snaf/reference/data"),
    "ascat_loci": ("data/ascat/G1000_loci_hg38.txt", "cnv/ascat/loci/G1000_loci_hg38.txt"),
    "ascat_alleles": ("data/ascat/G1000_alleles_hg38.txt", "cnv/ascat/loci/G1000_alleles_hg38.txt"),
}

SPECIAL_REQUIREMENTS = {
    "hla_la": ("hla_la_graph",),
    "spechla": ("spechla_db",),
    "purple": ("purple_reference",),
    "ascat": ("ascat_loci", "ascat_alleles"),
    "snaf": ("snaf_db",),
    "splicemutr": ("splicemutr_workflow",),
    "bam_matcher": ("bam_matcher_loci", "reference_fasta"),
}

BUILTIN_SAMPLE_RUNNERS = {"bam_matcher", "facets", "sequenza", "lohhla", "gatk", "optitype"}


@dataclass
class AutoConfigResult:
    status: str
    tools_manifest: str
    reference_manifest: str
    outputs: dict[str, str] = field(default_factory=dict)
    rows: list[dict[str, str]] = field(default_factory=list)


def _dump_yaml(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import yaml  # type: ignore

        text = yaml.safe_dump(data, sort_keys=False, allow_unicode=True)
    except Exception:
        text = json.dumps(data, ensure_ascii=False, indent=2)
    path.write_text(text, encoding="utf-8")


def _expand(value: Any) -> str:
    return os.path.expandvars(os.path.expanduser(str(value or "")))


def _existing(value: Any) -> str:
    path = Path(_expand(value)) if value else None
    return str(path.resolve()) if path and path.exists() else ""


def _resolve_executable(name: str, spec: dict[str, Any], roots: list[Path]) -> str:
    declared = str(spec.get("executable") or spec.get("path") or "")
    if declared:
        if Path(_expand(declared)).is_file():
            return str(Path(_expand(declared)).resolve())
        found = shutil.which(declared)
        if found:
            return found
    for env_name in TOOL_ENV_PATHS.get(name, ()):
        if _existing(os.environ.get(env_name)):
            return _existing(os.environ[env_name])
    candidates = [declared, name, name.replace("_", "-"), name.replace("_", "")]
    for root in roots:
        for candidate in candidates:
            if not candidate:
                continue
            for path in (root / "bin" / candidate, root / "tools" / name / candidate):
                if path.is_file() and os.access(path, os.X_OK):
                    return str(path.resolve())
        envs = root / "envs"
        if envs.is_dir():
            for candidate in candidates:
                matches = sorted(envs.glob(f"*/bin/{candidate}")) if candidate else []
                if matches:
                    return str(matches[0].resolve())
    return ""


def _resolve_reference(name: str, spec: dict[str, Any], roots: list[Path]) -> str:
    declared = _existing(spec.get("path"))
    if declared:
        return declared
    for env_name in REFERENCE_ENV_PATHS.get(name, ()):
        value = _existing(os.environ.get(env_name))
        if value:
            return value
    for root in roots:
        for relative in REFERENCE_CANDIDATES.get(name, ()):
            candidate = root / relative
            if candidate.exists():
                return str(candidate.resolve())
    return ""


def _reference_build_status(name: str, path: str, genome_build: str) -> tuple[str, str]:
    """Return a conservative build compatibility result for build-sensitive assets."""
    if not path or not genome_build:
        return "UNKNOWN", "genome build could not be verified"
    expected = genome_build.lower().replace("grch", "hg")
    haystack = path.lower()
    if name in {"reference_fasta", "bam_matcher_loci", "facets_snp_vcf", "ascat_loci", "ascat_alleles"}:
        if "hg19" in haystack or "grch37" in haystack:
            return ("MATCH", "path identifies GRCh37/hg19") if expected == "hg19" else ("MISMATCH", "GRCh37/hg19 asset cannot be used with GRCh38")
        if "hg38" in haystack or "grch38" in haystack or "assembly38" in haystack:
            return ("MATCH", "path identifies GRCh38/hg38") if expected == "hg38" else ("MISMATCH", "GRCh38/hg38 asset does not match requested build")
    if name == "bam_matcher_loci" and Path(path).is_file():
        opener = gzip.open if path.endswith(".gz") else open
        try:
            with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
                header = "".join(handle.readline() for _ in range(30)).lower()
            if "grch37" in header or "hg19" in header:
                return ("MATCH", "header identifies GRCh37/hg19") if expected == "hg19" else ("MISMATCH", "loci header identifies GRCh37/hg19")
            if "grch38" in header or "hg38" in header:
                return ("MATCH", "header identifies GRCh38/hg38") if expected == "hg38" else ("MISMATCH", "loci header identifies GRCh38/hg38")
        except OSError:
            return "UNKNOWN", "reference exists but its build metadata could not be read"
    return "UNKNOWN", "reference exists; build is not encoded in its path or header"


def _smoke(name: str, executable: str, spec: dict[str, Any], *, enabled: bool, timeout: int = 20) -> dict[str, str]:
    if not executable:
        return {"tool": name, "status": "NOT_RUN", "command": "", "message": "executable unavailable"}
    command_text = str(spec.get("smoke_command") or "")
    if command_text:
        command = shlex.split(command_text)
        command[0] = executable
    else:
        command = [executable, "--help"]
    if not enabled:
        return {"tool": name, "status": "PLANNED", "command": shlex.join(command), "message": "smoke disabled in plan mode"}
    try:
        proc = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=timeout, check=False)
        status = "PASS" if proc.returncode == 0 else "WARN"
        message = "" if proc.returncode == 0 else f"exit_code={proc.returncode}"
    except (OSError, subprocess.TimeoutExpired) as exc:
        status, message = "WARN", str(exc)
    return {"tool": name, "status": status, "command": shlex.join(command), "message": message}


def _template_for(name: str, executable: str, refs: dict[str, str], project_root: Path) -> str:
    if name == "hla_la" and executable and refs.get("hla_la_graph"):
        return (
            f"{shlex.quote(executable)} --BAM {{bam}} --graph {{hla_la_graph}} --sampleID {{sample_id}} "
            "--maxThreads {threads} --workingDir {outdir}"
        )
    if name == "bam_matcher" and (project_root / "scripts/run_bam_matcher_pair.sh").is_file():
        return (
            f"bash {shlex.quote(str(project_root / 'scripts/run_bam_matcher_pair.sh'))} "
            "--bam1 {normal_bam} --bam2 {tumor_bam} --reference {reference_fasta} "
            "--loci {bam_matcher_loci} --outdir {outdir}"
        )
    return ""


def configure_machine(
    *,
    project_root: str | Path,
    tools_manifest: str | Path,
    reference_manifest: str | Path,
    outdir: str | Path,
    tools_root: str | Path,
    reference_root: str | Path,
    licensed_root: str | Path,
    run_smoke: bool = False,
    publish_local: bool = False,
) -> AutoConfigResult:
    project = Path(project_root).resolve()
    output = Path(outdir).resolve()
    tool_data = load_limited_yaml(tools_manifest)
    ref_data = load_limited_yaml(reference_manifest)
    tools = tool_data.setdefault("tools", {})
    references = ref_data.setdefault("references", {})
    genome_build = str(ref_data.get("genome_build") or "")
    search_tool_roots = [project, Path(tools_root).resolve(), Path(licensed_root).resolve()]
    search_ref_roots = [Path(reference_root).resolve(), project]

    resolved_refs: dict[str, str] = {}
    rows: list[dict[str, str]] = []
    fixes: list[str] = []
    for name, spec_value in references.items():
        spec = spec_value if isinstance(spec_value, dict) else {"path": spec_value}
        path = _resolve_reference(str(name), spec, search_ref_roots)
        if path:
            spec["path"] = path
            build_status, build_reason = _reference_build_status(str(name), path, genome_build)
            if build_status == "MISMATCH":
                status, reason = "BUILD_MISMATCH", build_reason
                fixes.append(f"Reference `{name}`: {build_reason}; provide a {genome_build or 'matching-build'} asset.")
            else:
                status, reason = "CONFIGURED", f"reference path exists; {build_reason}"
        else:
            build_status = "UNKNOWN"
            status, reason = "MISSING", "reference path was not found in manifest, environment or standard layout"
            fixes.append(f"Reference `{name}`: provide a verified path under the reference root or set its documented environment variable.")
        references[name] = spec
        resolved_refs[str(name)] = path if status == "CONFIGURED" else ""
        rows.append({"component_type": "reference", "component": str(name), "status": status, "resolved_path": path, "requirements": "", "template_status": "N/A", "build_status": build_status, "reason": reason})

    smoke_rows: list[dict[str, str]] = []
    templates: dict[str, dict[str, str]] = {}
    for name, spec_value in tools.items():
        name = str(name)
        spec = spec_value if isinstance(spec_value, dict) else {"executable": spec_value}
        executable = _resolve_executable(name, spec, search_tool_roots)
        if executable:
            spec["executable"] = executable
        required_refs = SPECIAL_REQUIREMENTS.get(name, ())
        missing_refs = [ref for ref in required_refs if not resolved_refs.get(ref)]
        existing_template = str(spec.get("command_template") or "")
        generated_template = existing_template or _template_for(name, executable, resolved_refs, project)
        needs_template = name in SPECIAL_REQUIREMENTS and name not in BUILTIN_SAMPLE_RUNNERS
        if generated_template:
            spec["command_template"] = generated_template
            templates[name] = {"command_template": generated_template, "source": "preserved" if existing_template else "generated"}
        if not executable:
            status, reason = "UNAVAILABLE", "executable not found"
        elif missing_refs:
            status, reason = "PARTIAL", "missing references: " + ",".join(missing_refs)
        elif needs_template and not generated_template:
            status, reason = "PARTIAL", "sample-level command template requires manual workflow confirmation"
        else:
            status, reason = "CONFIGURED", "executable, references and runner/template are available"
        if status != "CONFIGURED":
            fixes.append(f"Tool `{name}`: {reason}.")
        tools[name] = spec
        smoke = _smoke(name, executable, spec, enabled=run_smoke)
        smoke_rows.append(smoke)
        rows.append({
            "component_type": "tool", "component": name, "status": status,
            "resolved_path": executable, "requirements": ",".join(required_refs),
            "template_status": "CONFIGURED" if generated_template else ("BUILTIN" if name in BUILTIN_SAMPLE_RUNNERS else "NOT_CONFIGURED"),
            "build_status": "N/A",
            "reason": reason,
        })

    manifests = output / "manifests"
    local_tools = manifests / "tools_manifest.configured.yaml"
    local_refs = manifests / "reference_manifest.configured.yaml"
    command_templates = manifests / "command_templates.yaml"
    _dump_yaml(local_tools, tool_data)
    _dump_yaml(local_refs, ref_data)
    _dump_yaml(command_templates, {"schema_version": 1, "command_templates": templates})
    status_path = output / "configuration_status.tsv"
    smoke_path = output / "smoke_tests.tsv"
    fixes_path = output / "recommended_fixes.md"
    write_tsv(status_path, rows)
    write_tsv(smoke_path, smoke_rows)
    fixes_path.write_text("# Recommended configuration fixes\n\n" + ("\n".join(f"- {item}" for item in fixes) if fixes else "No configuration fixes are currently required.") + "\n", encoding="utf-8")
    outputs = {
        "configured_tools_manifest": str(local_tools),
        "configured_reference_manifest": str(local_refs),
        "command_templates": str(command_templates),
        "configuration_status": str(status_path),
        "configuration_smoke_tests": str(smoke_path),
        "configuration_recommended_fixes": str(fixes_path),
    }
    if publish_local:
        local_dir = project / "configs" / "local"
        local_dir.mkdir(parents=True, exist_ok=True)
        for source, name in ((local_tools, "tools_manifest.generated.yaml"), (local_refs, "reference_manifest.generated.yaml"), (command_templates, "command_templates.generated.yaml")):
            target = local_dir / name
            shutil.copy2(source, target)
            outputs[f"published_{name.replace('.', '_')}"] = str(target)
    configured = sum(row["status"] == "CONFIGURED" for row in rows if row["component_type"] == "tool")
    available = sum(bool(row["resolved_path"]) for row in rows if row["component_type"] == "tool")
    status = "READY" if available and configured == available else "PARTIAL"
    write_json(output / "configuration_summary.json", {"status": status, "configured_tools": configured, "available_tools": available, "outputs": outputs})
    outputs["configuration_summary"] = str(output / "configuration_summary.json")
    return AutoConfigResult(status, str(local_tools), str(local_refs), outputs, rows)
