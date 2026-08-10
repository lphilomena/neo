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
    "nextflow": ("NEXTFLOW_BIN",),
    "star": ("STAR_BIN",),
    "regtools": ("REGTOOLS_BIN",),
    "salmon": ("SALMON_BIN",),
    "rsem": ("RSEM_CALCULATE_EXPRESSION_BIN",),
    "vep": ("VEP_BIN",),
    "gatk": ("GATK_BIN",),
    "optitype": ("OPTITYPE_BIN",),
    "hla_la": ("HLALA_BIN", "HLA_LA_BIN"),
    "spechla": ("SPECHLA_BIN",),
    "lohhla": ("LOHHLA_BIN",),
    "facets": ("FACETS_BIN",),
    "sequenza": ("SEQUENZA_BIN",),
    "purple": ("PURPLE_BIN",),
    "ascat": ("ASCAT_BIN",),
    "snaf": ("SNAF_BIN",),
    "splicemutr": ("SPLICEMUTR_BIN",),
    "bam_matcher": ("BAM_MATCHER_BIN",),
    "netmhcpan": ("NEOAG_NETMHCPAN_BIN",),
    "netmhcstabpan": ("NETMHCSTABPAN_BIN",),
    "prime": ("PRIME_BIN",),
    "mixmhcpred": ("MIXMHCPRED_BIN",),
    "bigmhc_im": ("BIGMHC_BIN",),
}

# Manifest names, upstream executable names and project wrappers intentionally
# differ. Keep the mapping explicit so a copied release can be configured
# without relying on the old machine's PATH.
TOOL_EXECUTABLE_CANDIDATES = {
    "nextflow": ("nextflow", "neoag-nextflow"),
    "star": ("STAR",),
    "regtools": ("regtools", "regtools-neoag"),
    "salmon": ("salmon",),
    "rsem": ("rsem-calculate-expression",),
    "vep": ("vep", "vep-neoag"),
    "gatk": ("gatk",),
    "java": ("java",),
    "optitype": ("optitype", "OptiTypePipeline.py"),
    "pvacseq": ("pvacseq", "pvacseq-neoag7"),
    "pvacfuse": ("pvacfuse", "pvacfuse-neoag7"),
    "pvacsplice": ("pvacsplice", "pvacsplice-neoag7"),
    "pvacbind": ("pvacbind", "pvacbind-neoag7"),
    "mhcflurry": ("mhcflurry-predict",),
    "hla_la": ("HLA-LA.pl", "HLA-LA"),
    "spechla": ("SpecHLA.py", "SpecHLA"),
    "lohhla": ("LOHHLA", "LOHHLA.pl"),
    "facets": ("runFACETS.R", "snp-pileup"),
    "sequenza": ("sequenza-utils",),
    "ascat": ("ascat.R", "ascat-v3"),
    "purple": ("purple",),
    "amber": ("amber",),
    "cobalt": ("cobalt",),
    "easyfuse": ("easyfuse", "easyfuse-neoag"),
    "star_fusion": ("STAR-Fusion", "star-fusion-neoag"),
    "fusioncatcher": ("fusioncatcher", "fusioncatcher-neoag"),
    "spladder": ("spladder", "spladder-neoag"),
    "irfinder_s": ("IRFinder", "irfinder-s-neoag"),
    "immunopepper": ("immunopepper", "immunopepper-neoag"),
    "splicemutr": ("SpliceMutr", "splicemutr-neoag"),
    "snaf": ("snaf", "snaf-neoag"),
    "bigmhc_im": ("bigmhc_predict",),
    "pyclone": ("pyclone", "pyclone-vi"),
    "polysolver": ("run-polysolver", "shell_call_hla_type"),
    "bam_matcher": ("bam-matcher", "bam_matcher.py"),
    "netmhcpan": ("netMHCpan",),
    "netmhcstabpan": ("netMHCstabpan",),
    "prime": ("PRIME",),
    "mixmhcpred": ("MixMHCpred",),
}

TOOL_RELATIVE_CANDIDATES = {
    "netmhcpan": ("data/predictors/netMHCpan/netMHCpan",),
    "prime": ("data/predictors/prime/PRIME",),
    "mixmhcpred": ("data/predictors/mixMHCpred_install/MixMHCpred",),
    "lohhla": ("data/tools/lohhla/LOHHLA",),
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
    "normal_proteome": ("NEOAG_NORMAL_PROTEOME",),
    "normal_ligandome": ("NEOAG_NORMAL_LIGANDOME",),
    "normal_junctions": ("NEOAG_NORMAL_JUNCTIONS",),
    "normal_junctions_liver": ("NEOAG_NORMAL_JUNCTIONS_LIVER",),
    "ctat_genome_lib": ("CTAT_GENOME_LIB",),
    "easyfuse_ref": ("EASYFUSE_REF",),
    "salmon_index": ("SALMON_INDEX",),
    "salmon_tx2gene": ("SALMON_TX2GENE", "TX2GENE_TSV"),
    "rsem_reference": ("RSEM_REFERENCE",),
    "bigmhc_models": ("BIGMHC_MODELS",),
}

REFERENCE_CANDIDATES = {
    "reference_fasta": ("data/ref/hg38/Homo_sapiens_assembly38.fasta", "GRCh38/fasta/GRCh38.fa"),
    "gencode_gtf": ("data/ref/hg38/gencode.gtf", "GRCh38/gencode/gencode.v44.annotation.gtf.gz"),
    "vep_cache": ("data/vep/homo_sapiens/105_GRCh38", "GRCh38/vep/vep_115_GRCh38"),
    "hla_la_graph": ("data/hla/PRG_MHC_GRCh38_withIMGT", "hla/hla-la/graph/PRG_MHC_GRCh38_withIMGT"),
    "spechla_db": ("data/hla/spechla/db", "data/hla/spechla_db", "hla/spechla/db"),
    "facets_snp_vcf": ("data/facets/reference/common_snp.hg38.vcf.gz", "cnv/facets/common_snp.vcf.gz"),
    "facets_common_snp_vcf": ("data/facets/reference/common_snp.hg38.vcf.gz", "cnv/facets/common_snp.vcf.gz"),
    "sequenza_gc_wiggle": ("data/sequenza/reference/Homo_sapiens.GRCh38.dna.primary_assembly.chr.gc50.wig.gz",),
    "purple_reference": ("data/hmf/purple_reference", "cnv/purple"),
    "bam_matcher_loci": ("data/sample_identity/bam_matcher.common_snps.hg38.vcf", "GRCh38/sample_identity/bam_matcher.common_snps.vcf"),
    "snaf_db": ("data/snaf/reference/data", "splice/snaf/reference/data"),
    "ascat_loci": ("data/ascat/G1000_loci_hg38.txt", "cnv/ascat/loci/G1000_loci_hg38.txt"),
    "ascat_alleles": ("data/ascat/G1000_alleles_hg38.txt", "cnv/ascat/loci/G1000_alleles_hg38.txt"),
    "normal_proteome": ("data/normal/proteome/Homo_sapiens.GRCh38.pep.all.fa",),
    "normal_ligandome": ("data/normal/ligandome/normal_ms_ligands.tsv",),
    "normal_junctions": ("data/normal/junctions/normal_junctions.GRCh38.tsv.gz",),
    "normal_junctions_liver": ("data/normal/junctions/gtex_v8_liver.GRCh38.tsv.gz",),
    "ctat_genome_lib": (
        "data/ref/ctat/current/ctat_genome_lib_build_dir",
        "data/ref/ctat/current",
        "data/easyfuse/easyfuse_ref_v4/starfusion_index",
    ),
    "easyfuse_ref": ("data/easyfuse/easyfuse_ref_v4", "data/ref/ctat/current"),
    "star_index": (
        "data/ref/ctat/current/ctat_genome_lib_build_dir/ref_genome.fa.star.idx",
        "data/easyfuse/easyfuse_ref_v4/star_index",
        "data/easyfuse/easyfuse_ref_v4/starfusion_index/ref_genome.fa.star.idx",
    ),
    "salmon_index": ("data/rna/salmon_index", "data/rna/gencode_v49/salmon_index"),
    "salmon_tx2gene": ("data/rna/tx2gene.tsv", "data/rna/gencode_v49/tx2gene.tsv"),
    "tx2gene": ("data/rna/tx2gene.tsv", "data/rna/gencode_v49/tx2gene.tsv"),
    "rsem_reference": ("data/rna/rsem_reference",),
    "bigmhc_models": ("data/predictors/bigmhc/models",),
    "lohhla_reference": ("data/lohhla/polysolver", "data/tools/lohhla"),
    "iedb_human_ms_ligands_detail": ("data/normal/iedb_mhc_ligand/2026-07-14/build/iedb_human_ms_ligands_detail.tsv",),
    "iedb_human_ms_ligands": ("data/normal/iedb_mhc_ligand/2026-07-14/build/iedb_human_ms_ligands.tsv",),
    "iedb_strict_normal_ligands": ("data/normal/iedb_mhc_ligand/2026-07-14/build/iedb_human_normal_direct_ex_vivo_ligands.tsv",),
    "cancer_gene_list": ("data/annotation/cancer_gene/oncokb_v7.3_2026-06-25/cancerGeneList.tsv",),
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


def _dump_yaml_scalar(value: Any) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return str(value)
    text = str(value)
    if text == "" or any(ch in text for ch in ":#{}[]&*!|>'\"%@`,\n") or text.strip() != text:
        return json.dumps(text, ensure_ascii=False)
    return text


def _dump_yaml_lines(data: Any, indent: int = 0) -> list[str]:
    pad = "  " * indent
    if isinstance(data, dict):
        if not data:
            return [f"{pad}{{}}"]
        lines: list[str] = []
        for key, value in data.items():
            key_text = str(key)
            if isinstance(value, dict):
                lines.append(f"{pad}{key_text}:")
                lines.extend(_dump_yaml_lines(value, indent + 1) if value else [f"{pad}  {{}}"])
            elif isinstance(value, list):
                lines.append(f"{pad}{key_text}:")
                if not value:
                    lines.append(f"{pad}  []")
                else:
                    for item in value:
                        if isinstance(item, (dict, list)):
                            lines.append(f"{pad}-")
                            lines.extend(_dump_yaml_lines(item, indent + 1))
                        else:
                            lines.append(f"{pad}- {_dump_yaml_scalar(item)}")
            else:
                lines.append(f"{pad}{key_text}: {_dump_yaml_scalar(value)}")
        return lines
    if isinstance(data, list):
        if not data:
            return [f"{pad}[]"]
        lines = []
        for item in data:
            if isinstance(item, (dict, list)):
                lines.append(f"{pad}-")
                lines.extend(_dump_yaml_lines(item, indent + 1))
            else:
                lines.append(f"{pad}- {_dump_yaml_scalar(item)}")
        return lines
    return [f"{pad}{_dump_yaml_scalar(data)}"]


def _dump_yaml(path: Path, data: dict[str, Any]) -> None:
    """Write a YAML-looking manifest.

    Prefer PyYAML. Without it, emit indentation YAML (not JSON) so a later
    ``load_limited_yaml`` pass without PyYAML still parses nested ``references``.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        import yaml  # type: ignore

        text = yaml.safe_dump(data, sort_keys=False, allow_unicode=True)
    except Exception:
        text = "\n".join(_dump_yaml_lines(data)) + "\n"
    path.write_text(text, encoding="utf-8")


def _expand(value: Any) -> str:
    return os.path.expandvars(os.path.expanduser(str(value or "")))


def _existing(value: Any) -> str:
    path = Path(_expand(value)) if value else None
    if not path:
        return ""
    try:
        return str(path.resolve()) if path.exists() else ""
    except OSError:
        # A manifest may contain an otherwise valid path owned by another
        # account or copied from another machine. Treat it as unavailable and
        # continue portable-root discovery instead of aborting installation.
        return ""


def _resolve_executable(name: str, spec: dict[str, Any], roots: list[Path]) -> str:
    declared = str(spec.get("executable") or spec.get("path") or "")
    if declared:
        declared_path = _existing(declared)
        if declared_path and Path(declared_path).is_file():
            return declared_path
        found = shutil.which(declared)
        if found:
            return found
    for env_name in TOOL_ENV_PATHS.get(name, ()):
        if _existing(os.environ.get(env_name)):
            return _existing(os.environ[env_name])
    candidates = list(dict.fromkeys([
        declared,
        *TOOL_EXECUTABLE_CANDIDATES.get(name, ()),
        name,
        name.replace("_", "-"),
        name.replace("_", ""),
    ]))
    for root in roots:
        for relative in TOOL_RELATIVE_CANDIDATES.get(name, ()):
            resolved = _existing(root / relative)
            if resolved and Path(resolved).is_file() and os.access(resolved, os.X_OK):
                return resolved
        for candidate in candidates:
            if not candidate:
                continue
            for path in (root / "bin" / candidate, root / "tools" / name / candidate):
                resolved = _existing(path)
                if resolved and Path(resolved).is_file() and os.access(resolved, os.X_OK):
                    return resolved
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
            resolved = _existing(candidate)
            if resolved:
                return resolved
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
    reference_root_path = Path(reference_root).resolve()
    search_tool_roots = [project, Path(tools_root).resolve(), Path(licensed_root).resolve(), reference_root_path]
    # Doctor/tier checks resolve tools via configured absolute paths. Include the
    # active conda base so STAR/gatk/etc. installed into envs/neoag-* are found
    # even when not on PATH (Skill1 standard puts STAR in neoag-splicemutr/fusion).
    conda_base = Path(os.environ.get("NEOAG_CONDA_BASE") or "").expanduser()
    if conda_base.is_dir():
        search_tool_roots.append(conda_base)
    search_ref_roots = [reference_root_path, project]

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
