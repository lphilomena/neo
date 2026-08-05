from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path
from typing import Any

from neoag.controlled_execution.doctor import CheckRow, DoctorResult, run_doctor
from neoag.controlled_execution.io_utils import (
    load_limited_yaml,
    markdown_table,
    now_iso,
    sha256_file,
    write_json,
    write_tsv,
)

from .contracts import MacroResult, MacroStep
from .auto_config import configure_machine
from .errors import FailureCode
from .state import RunLayout, audit, new_run_id, safe_identifier, update_case_state


DEFAULT_ASSET_SOURCE_HOST = "na@10.200.50.134"
DEFAULT_ASSET_SOURCE_ROOT = "/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git"


def _apply_default_asset_source(args: dict[str, Any]) -> dict[str, Any]:
    configured = dict(args)
    if not configured.get("asset_source_host"):
        configured["asset_source_host"] = os.environ.get(
            "OPEN_NEO_ASSET_SOURCE_HOST", DEFAULT_ASSET_SOURCE_HOST,
        )
    if not configured.get("asset_source_root"):
        configured["asset_source_root"] = os.environ.get(
            "OPEN_NEO_ASSET_SOURCE_ROOT", DEFAULT_ASSET_SOURCE_ROOT,
        )
    return configured


TIER_REQUIRED_TOOLS = {
    "review": ["python", "neoag", "neoag-skill"],
    "core": ["python", "neoag", "neoag-skill"],
    "prediction": ["python", "neoag", "neoag-skill", "vep", "netmhcpan", "mhcflurry", "pvacseq"],
    "full": [
        "python", "neoag", "neoag-skill", "vep", "netmhcpan", "mhcflurry",
        "pvacseq", "java", "nextflow", "bwa", "star", "gatk",
    ],
}

TIER_TOOL_GROUPS = {
    "prediction": {
        "immunogenicity": ["prime", "bigmhc_im", "deepimmuno"],
    },
    "full": {
        "immunogenicity": ["prime", "bigmhc_im", "deepimmuno"],
        "hla_typing": ["optitype", "spechla", "hla_la"],
        "hla_loh": ["lohhla", "spechla"],
        "purity_cnv": ["facets", "purple", "ascat", "sequenza"],
        "fusion": ["easyfuse", "arriba", "star_fusion", "fusioncatcher"],
        "splice": ["snaf", "splicemutr"],
        "ccf": ["pyclone_vi", "facets", "purple"],
        "sample_identity": ["bam_matcher"],
    },
}

TIER_REQUIRED_REFERENCES = {
    "review": [],
    "core": [],
    "prediction": [
        "reference_fasta", "reference_fasta_fai", "reference_fasta_dict",
        "gencode_gtf", "vep_cache", "normal_proteome",
    ],
    "full": [
        "reference_fasta", "reference_fasta_fai", "reference_fasta_dict",
        "gencode_gtf", "vep_cache", "normal_proteome", "normal_ligandome",
        "normal_junctions", "star_index", "ctat_genome_lib", "salmon_index",
        "tx2gene",
    ],
}

TIER_REFERENCE_GROUPS = {
    "full": {
        "hla_typing_reference": ["hla_reference", "spechla_db", "hla_la_graph"],
        "hla_loh_reference": ["lohhla_reference", "spechla_db"],
        "purity_reference": ["facets_snp_vcf", "sequenza_gc_wiggle", "purple_reference"],
        "sample_identity_reference": ["bam_matcher_loci"],
    },
}

REFERENCE_ALIASES = {
    "reference_fasta": ["reference_fasta", "reference.fasta", "fasta"],
    "gencode_gtf": ["gencode_gtf", "gencode.gtf", "reference.gtf", "gtf"],
    "vep_cache": ["vep_cache", "vep.cache"],
    "normal_proteome": ["normal_proteome"],
    "normal_ligandome": ["normal_ligandome", "normal_hla_ligands"],
    "normal_junctions": ["normal_junctions"],
    "star_index": ["star_index"],
    "ctat_genome_lib": ["ctat_genome_lib", "ctat"],
    "salmon_index": ["salmon_index"],
    "tx2gene": ["tx2gene"],
    "hla_reference": ["hla_reference"],
    "spechla_db": ["spechla_db", "spechla.db"],
    "hla_la_graph": ["hla_la_graph", "hla-la.graph", "prg_mhc"],
    "lohhla_reference": ["lohhla_reference"],
    "facets_snp_vcf": ["facets_snp_vcf", "facets.vcf", "common_snp"],
    "sequenza_gc_wiggle": ["sequenza_gc_wiggle", "gc_wiggle"],
    "purple_reference": ["purple_reference"],
    "bam_matcher_loci": ["bam_matcher_loci", "sample_identity_vcf"],
    "ascat_loci": ["ascat_loci"],
    "ascat_alleles": ["ascat_alleles"],
    "snaf_workflow": ["snaf_workflow"],
    "splicemutr_workflow": ["splicemutr_workflow"],
}

OK_STATUSES = {"OK", "INFO"}
EXECUTION_MODES = {"repair", "install", "resume"}


def _rewrite_asset_manifest(
    source: Path,
    output: Path,
    *,
    tools_root: Path,
    reference_root: Path,
    licensed_root: Path,
    source_root: Path | None,
) -> Path:
    comments: list[str] = []
    data_lines: list[str] = []
    for line in source.read_text(encoding="utf-8").splitlines():
        if not data_lines and (not line.strip() or line.lstrip().startswith("#")):
            comments.append(line)
        else:
            data_lines.append(line)
    if not data_lines:
        raise ValueError(f"asset manifest has no TSV header: {source}")
    reader = csv.DictReader(data_lines, delimiter="\t")
    required_columns = {"asset_name", "source_path", "target_path"}
    if not required_columns.issubset(set(reader.fieldnames or [])):
        raise ValueError(f"asset manifest is missing required columns: {source}")
    rows = []
    for row in reader:
        target = str(row.get("target_path") or "")
        target = target.replace("/srv/neoag-tools", str(tools_root), 1)
        target = target.replace("/srv/neoag-licensed", str(licensed_root), 1)
        target = target.replace("/srv/neoag-assets/install", str(reference_root), 1)
        source_path = str(row.get("source_path") or "")
        if source_root is not None:
            source_path = source_path.replace("/srv/neoag-assets/source", str(source_root), 1)
        row["source_path"] = source_path
        row["target_path"] = target
        rows.append(row)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        for line in comments:
            handle.write(line + "\n")
        writer = csv.DictWriter(handle, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return output


def _required_asset_sources_missing(manifest: str | Path, source_host: str) -> list[str]:
    if source_host:
        return []
    lines = [
        line for line in Path(manifest).read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    reader = csv.DictReader(lines, delimiter="\t")
    missing = []
    for row in reader:
        required = str(row.get("required") or "1").strip().lower() not in {"0", "no", "false", "optional"}
        source = str(row.get("source_path") or "")
        if required and ":" not in source and not Path(source).exists():
            missing.append(f"{row.get('asset_name')}={source}")
    return missing


def _environment_inventory(project_root: Path) -> list[dict[str, Any]]:
    commands = {
        "python": sys.executable,
        "java": "java",
        "docker": "docker",
        "apptainer": "apptainer",
        "nextflow": "nextflow",
        "git": "git",
    }
    rows = []
    for name, command in commands.items():
        path = command if Path(command).is_absolute() and Path(command).exists() else shutil.which(command)
        rows.append({"component": name, "status": "FOUND" if path else "MISSING", "path": str(path or "")})
    disk_root = project_root if project_root.exists() else project_root.parent
    disk = shutil.disk_usage(disk_root)
    rows.append({"component": "disk_free_bytes", "status": "INFO", "path": str(disk.free)})
    rows.append({"component": "python_version", "status": "INFO", "path": platform.python_version()})
    rows.append({"component": "platform", "status": "INFO", "path": platform.platform()})
    return rows


def _safe_release_members(archive: tarfile.TarFile, destination: Path) -> list[tarfile.TarInfo]:
    safe: list[tarfile.TarInfo] = []
    root = destination.resolve()
    for member in archive.getmembers():
        name = member.name.replace("\\", "/")
        target = (destination / name).resolve()
        try:
            target.relative_to(root)
        except ValueError as exc:
            raise ValueError(f"unsafe archive path: {member.name}") from exc
        if name.startswith("/") or member.issym() or member.islnk() or not (member.isfile() or member.isdir()):
            raise ValueError(f"unsupported archive member: {member.name}")
        safe.append(member)
    return safe


def _find_staged_project_root(destination: Path) -> Path:
    candidates = [path.parent for path in destination.rglob("pyproject.toml")]
    candidates += [path.parent for path in destination.rglob("setup.py") if path.parent not in candidates]
    if not candidates:
        raise ValueError("release archive has no pyproject.toml or setup.py")
    depths = {path: len(path.relative_to(destination).parts) for path in candidates}
    shallowest = min(depths.values())
    best = sorted(path for path, depth in depths.items() if depth == shallowest)
    if len(best) != 1:
        raise ValueError("release archive has multiple candidate project roots")
    return best[0]


def _stage_release(archive_path: Path, sha256: str, layout: RunLayout) -> tuple[Path, dict[str, str], bool]:
    stage_root = layout.root / "release_staging" / sha256[:16]
    manifest_path = layout.root / "release_staging.json"
    previous = {}
    if manifest_path.is_file():
        try:
            previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        except Exception:
            previous = {}
    previous_root = Path(str(previous.get("project_root") or ""))
    if previous.get("archive_sha256") == sha256 and previous_root.is_dir():
        return previous_root, {"release_staging": str(manifest_path), "staged_project_root": str(previous_root)}, True
    stage_root.mkdir(parents=True, exist_ok=True)
    try:
        with tarfile.open(archive_path, "r:*") as archive:
            members = _safe_release_members(archive, stage_root)
            for member in members:
                target = stage_root / member.name
                if member.isdir():
                    target.mkdir(parents=True, exist_ok=True)
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                source = archive.extractfile(member)
                if source is None:
                    raise ValueError(f"could not read archive member: {member.name}")
                with source, target.open("wb") as handle:
                    shutil.copyfileobj(source, handle)
                target.chmod(member.mode & 0o777)
    except (tarfile.TarError, OSError, ValueError) as exc:
        raise ValueError(f"release staging failed: {exc}") from exc
    project_root = _find_staged_project_root(stage_root)
    payload = {
        "schema_version": "open-neo-release-staging-v1",
        "archive": str(archive_path.resolve()),
        "archive_sha256": sha256,
        "staged_at": now_iso(),
        "staging_root": str(stage_root.resolve()),
        "project_root": str(project_root.resolve()),
    }
    write_json(manifest_path, payload)
    return project_root, {"release_staging": str(manifest_path), "staged_project_root": str(project_root)}, False


def _write_local_manifest_templates(layout: RunLayout, project_root: Path, args: dict[str, Any]) -> dict[str, str]:
    deploy_root = Path(str(args.get("deploy_root") or "/opt/neoag"))
    tools_root = Path(str(args.get("tools_root") or deploy_root / "env_tool"))
    reference_root = Path(str(args.get("reference_root") or deploy_root / "refs"))
    licensed_root = Path(str(args.get("licensed_root") or deploy_root / "licensed_tools"))
    tools = layout.manifests / "tools_manifest.local.yaml"
    refs = layout.manifests / "reference_manifest.local.yaml"
    tools.write_text(
        """manifest_version: 2
tools:
  bwa:
    executable: bwa
    mode: conda_or_container
  star:
    executable: STAR
    mode: conda_or_container
  gatk:
    executable: gatk
    mode: conda_or_container
  bam_matcher:
    executable: bam-matcher
    mode: conda
    recommended: true
  vep:
    executable: vep
    mode: conda_or_container
  netmhcpan:
    executable: netMHCpan
    mode: local_license
    license_required: true
    distributable: false
  netmhcstabpan:
    executable: NetMHCstabpan
    mode: local_license
    license_required: true
    distributable: false
  mhcflurry:
    executable: mhcflurry-predict
    mode: conda
  pvacseq:
    executable: pvacseq
    mode: conda
  prime:
    executable: PRIME
    mode: local
  bigmhc_im:
    executable: bigmhc_predict
    mode: conda
  deepimmuno:
    executable: deepimmuno-cnn.py
    mode: conda
  optitype:
    executable: OptiTypePipeline.py
    mode: conda
  hla_la:
    executable: HLA-LA.pl
    mode: local
  spechla:
    executable: SpecHLA
    mode: local
  lohhla:
    executable: LOHHLA
    mode: local_or_container
  facets:
    executable: runFACETS.R
    mode: conda_or_container
  purple:
    executable: purple
    mode: local_or_container
  ascat:
    executable: ascat.R
    mode: conda
  sequenza:
    executable: sequenza-utils
    mode: conda
  arriba:
    executable: arriba
    mode: conda_or_container
  star_fusion:
    executable: STAR-Fusion
    mode: conda_or_container
  fusioncatcher:
    executable: fusioncatcher
    mode: conda_or_container
  easyfuse:
    executable: easyfuse
    mode: conda_or_container
  snaf:
    executable: snaf
    mode: conda
  splicemutr:
    executable: splicemutr
    mode: conda
  pyclone_vi:
    executable: pyclone-vi
    mode: conda
""",
        encoding="utf-8",
    )
    refs_text = """manifest_version: 2
genome_build: GRCh38
references:
  reference_fasta:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/ref/hg38/Homo_sapiens_assembly38.fasta'
  gencode_gtf:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/ref/hg38/gencode.gtf'
  vep_cache:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/vep/homo_sapiens/105_GRCh38'
  normal_proteome:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/normal/proteome/Homo_sapiens.GRCh38.pep.all.fa'
  normal_ligandome:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/normal/ligandome/normal_ms_ligands.tsv'
  normal_junctions:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/normal/junctions/normal_junctions.tsv'
  star_index:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/rna/star_index'
  ctat_genome_lib:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/ctat/current'
  salmon_index:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/rna/salmon_index'
  tx2gene:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/rna/tx2gene.tsv'
  hla_reference:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/hla'
  spechla_db:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/hla/spechla_db'
  hla_la_graph:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/hla/PRG_MHC_GRCh38_withIMGT'
  lohhla_reference:
    path: '${OPEN_NEO_TOOLS_ROOT}/tools/lohhla'
  facets_snp_vcf:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/facets/reference/common_snp.hg38.vcf.gz'
  sequenza_gc_wiggle:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/sequenza/reference/Homo_sapiens.GRCh38.dna.primary_assembly.chr.gc50.wig.gz'
  purple_reference:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/hmf/purple_reference'
  bam_matcher_loci:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/sample_identity/bam_matcher.common_snps.hg38.vcf'
  ascat_loci:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/ascat/G1000_loci_hg38.txt'
  ascat_alleles:
    path: '${OPEN_NEO_REFERENCE_ROOT}/data/ascat/G1000_alleles_hg38.txt'
  snaf_workflow:
    path: '${OPEN_NEO_REFERENCE_ROOT}/workflows/snaf.workflow.yaml'
  splicemutr_workflow:
    path: '${OPEN_NEO_REFERENCE_ROOT}/workflows/splicemutr.workflow.yaml'
"""
    refs_text = refs_text.replace("${OPEN_NEO_REFERENCE_ROOT}", str(reference_root))
    refs_text = refs_text.replace("${OPEN_NEO_TOOLS_ROOT}", str(tools_root))
    refs.write_text(refs_text, encoding="utf-8")
    paths = layout.manifests / "paths.env"
    paths.write_text(
        "\n".join([
            f"export OPEN_NEO_PROJECT_ROOT={shlex.quote(str(project_root))}",
            f"export OPEN_NEO_TOOLS_ROOT={shlex.quote(str(tools_root))}",
            f"export OPEN_NEO_REFERENCE_ROOT={shlex.quote(str(reference_root))}",
            f"export OPEN_NEO_LICENSED_ROOT={shlex.quote(str(licensed_root))}",
            *([f"export NEOAG_CONDA_BASE={shlex.quote(str(args['conda_base']))}"] if args.get("conda_base") else []),
            "",
        ]),
        encoding="utf-8",
    )
    outputs = {
        "tools_manifest_template": str(tools),
        "reference_manifest_template": str(refs),
        "paths_env": str(paths),
    }
    source_manifest = Path(str(args.get("asset_manifest") or project_root / "configs/assets/production_assets.tsv"))
    if source_manifest.is_file():
        asset_manifest = _rewrite_asset_manifest(
            source_manifest,
            layout.manifests / "production_assets.local.tsv",
            tools_root=tools_root,
            reference_root=reference_root,
            licensed_root=licensed_root,
            source_root=Path(str(args["asset_source_root"])).resolve() if args.get("asset_source_root") else None,
        )
        outputs["asset_manifest_local"] = str(asset_manifest)
    return outputs


def _tool_statuses(doctor_rows: list[CheckRow]) -> dict[str, str]:
    statuses = {str(row.name).lower(): str(row.status) for row in doctor_rows if row.category in {"tool", "core"}}
    if statuses.get("python_import_neoag") == "OK":
        statuses["neoag"] = "OK"
        try:
            import neoag.skill_taxonomy.cli  # noqa: F401
            skill_module_ok = True
        except Exception:
            skill_module_ok = False
        if skill_module_ok or shutil.which("neoag-skill") or shutil.which("neoag"):
            statuses["neoag-skill"] = "OK"
    return statuses


def _flatten_manifest_paths(data: Any, prefix: str = "") -> dict[str, str]:
    paths: dict[str, str] = {}
    if isinstance(data, dict):
        if isinstance(data.get("path"), str):
            paths[prefix.lower()] = str(data["path"])
        for key, value in data.items():
            if key in {"path", "sha256"}:
                continue
            child = f"{prefix}.{key}" if prefix else str(key)
            paths.update(_flatten_manifest_paths(value, child))
    elif isinstance(data, str) and ("/" in data or data.startswith("$") or data.startswith(".")):
        paths[prefix.lower()] = data
    return paths


def _declared_reference(name: str, paths: dict[str, str]) -> str:
    aliases = REFERENCE_ALIASES.get(name, [name])
    candidates = [(key, value) for key, value in paths.items() if any(alias in key for alias in aliases)]
    if not candidates:
        return ""
    return sorted(candidates, key=lambda item: (len(item[0]), item[0]))[0][1]


def _safe_path_exists(path: Path) -> bool:
    try:
        return path.exists()
    except OSError:
        return False


def _safe_path_is_file(path: Path) -> bool:
    try:
        return path.is_file()
    except OSError:
        return False


def _reference_status(name: str, paths: dict[str, str]) -> tuple[str, str]:
    if name in {"reference_fasta_fai", "reference_fasta_dict"}:
        fasta = _declared_reference("reference_fasta", paths)
        if not fasta:
            return "MISSING", "reference_fasta not declared"
        expanded = Path(os.path.expandvars(os.path.expanduser(fasta)))
        if name == "reference_fasta_fai":
            sidecars = [Path(str(expanded) + ".fai")]
        else:
            sidecars = [expanded.with_suffix(".dict"), Path(str(expanded).rsplit(".fasta", 1)[0].rsplit(".fa", 1)[0] + ".dict")]
        existing = next((path for path in sidecars if _safe_path_is_file(path)), None)
        return ("OK", str(existing)) if existing else ("MISSING", str(sidecars[0]))
    value = _declared_reference(name, paths)
    if not value:
        return "MISSING", "not declared"
    expanded = Path(os.path.expandvars(os.path.expanduser(value)))
    return ("OK", str(expanded)) if _safe_path_exists(expanded) else ("MISSING", str(expanded))


def _assess_tier(tier: str, doctor_rows: list[CheckRow], reference_manifest: str | Path | None) -> tuple[str, list[dict[str, str]]]:
    tool_status = _tool_statuses(doctor_rows)
    rows: list[dict[str, str]] = []
    for name in TIER_REQUIRED_TOOLS[tier]:
        status = tool_status.get(name, "MISSING")
        rows.append({"kind": "tool", "requirement": name, "required": "true", "status": "OK" if status in OK_STATUSES else "MISSING", "evidence": status})
    for group, members in TIER_TOOL_GROUPS.get(tier, {}).items():
        available = [name for name in members if tool_status.get(name) in OK_STATUSES]
        rows.append({"kind": "tool_group", "requirement": group, "required": "true", "status": "OK" if available else "MISSING", "evidence": ",".join(available) or "any of: " + ",".join(members)})

    ref_data: dict[str, Any] = {}
    if reference_manifest:
        try:
            ref_data = load_limited_yaml(reference_manifest)
        except Exception:
            ref_data = {}
    paths = _flatten_manifest_paths(ref_data)
    for name in TIER_REQUIRED_REFERENCES[tier]:
        status, evidence = _reference_status(name, paths)
        rows.append({"kind": "reference", "requirement": name, "required": "true", "status": status, "evidence": evidence})
    for group, members in TIER_REFERENCE_GROUPS.get(tier, {}).items():
        available = []
        evidence = []
        for name in members:
            status, detail = _reference_status(name, paths)
            evidence.append(f"{name}={status}:{detail}")
            if status == "OK":
                available.append(name)
        rows.append({"kind": "reference_group", "requirement": group, "required": "true", "status": "OK" if available else "MISSING", "evidence": ";".join(evidence)})

    missing = [row for row in rows if row["required"] == "true" and row["status"] != "OK"]
    core_missing = any(row["requirement"] in {"python", "neoag", "neoag-skill"} for row in missing)
    return ("BLOCKED" if core_missing else "PARTIAL" if missing else "READY"), rows


def _deployment_command(args: dict[str, Any], project_root: Path, layout: RunLayout, *, execute: bool) -> list[str]:
    deploy_root = Path(str(args.get("deploy_root") or "/opt/neoag"))
    script = project_root / ".agents/skills/neoag-remote-deploy/scripts/16_install_new_machine.sh"
    tier = str(args.get("deployment_tier") or "core").lower()
    default_profile = "standard" if tier in {"prediction", "full"} else "minimal"
    installer_profile = str(args.get("installer_profile") or default_profile)
    command = [
        "bash", str(script),
        "--project-root", str(project_root),
        "--tools-root", str(args.get("tools_root") or deploy_root / "env_tool"),
        "--reference-root", str(args.get("reference_root") or deploy_root / "refs"),
        "--licensed-root", str(args.get("licensed_root") or deploy_root / "licensed_tools"),
        "--outdir", str(layout.root / "deployment"),
        "--" + installer_profile,
    ]
    conda_base = args.get("conda_base") or os.environ.get("NEOAG_CONDA_BASE", "")
    if conda_base:
        command += ["--conda-base", str(conda_base)]
    asset_manifest = args.get("asset_manifest")
    deployment_reference_manifest = args.get("deployment_reference_manifest") or args.get("reference_manifest")
    if asset_manifest:
        command += ["--asset-manifest", str(asset_manifest)]
    if deployment_reference_manifest:
        command += ["--reference-manifest", str(deployment_reference_manifest)]
    if args.get("asset_source_host"):
        command += ["--asset-source-host", str(args["asset_source_host"])]
    if bool(args.get("allow_download", False)):
        command.append("--allow-download")
    if bool(args.get("install_claude_code", False)):
        command += ["--claude-code", "--claude-code-channel", str(args.get("claude_code_channel") or "stable")]
    if bool(args.get("no_sync_assets", False)):
        command.append("--no-sync-assets")
    if execute:
        command.append("--execute")
    return command


def _collect_claude_code_status(
    layout: RunLayout, *, requested: bool, executed: bool,
) -> tuple[dict[str, str], dict[str, str]]:
    report = layout.root / "deployment/readme_tools/claude_code/claude_code_install_report.md"
    row = {
        "requested": "true" if requested else "false",
        "status": "NOT_REQUESTED",
        "version": "",
        "binary": "",
        "report": "",
    }
    if requested and not executed:
        row["status"] = "PLANNED"
    elif requested and not report.is_file():
        row["status"] = "MISSING_REPORT"
    elif requested:
        text = report.read_text(encoding="utf-8", errors="replace")
        values: dict[str, str] = {}
        for line in text.splitlines():
            for label, key in (("Binary", "binary"), ("Version", "version")):
                prefix = f"{label}: `"
                if line.startswith(prefix) and line.endswith("`"):
                    values[key] = line[len(prefix):-1]
        row.update(values)
        row["report"] = str(report)
        if row["binary"] and row["version"] and row["version"] != "UNASSESSED":
            row["status"] = "READY"
        else:
            row["status"] = "UNVERIFIED"

    status_tsv = layout.root / "claude_code_status.tsv"
    status_json = layout.root / "claude_code_status.json"
    write_tsv(status_tsv, [row])
    write_json(status_json, row)
    outputs = {
        "claude_code_status": str(status_tsv),
        "claude_code_status_json": str(status_json),
    }
    if row["report"]:
        outputs["claude_code_install_report"] = row["report"]
    return row, outputs


def _command_hash(command: list[str]) -> str:
    return hashlib.sha256(json.dumps(command, ensure_ascii=True).encode("utf-8")).hexdigest()


def _run_deployment(command: list[str], project_root: Path, log_path: Path, checkpoint_path: Path, *, timeout: int, resume: bool) -> tuple[str, str, str]:
    command_hash = _command_hash(command)
    if resume and checkpoint_path.is_file():
        try:
            checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        except Exception:
            checkpoint = {}
        if checkpoint.get("status") == "PASS" and checkpoint.get("command_hash") == command_hash:
            return "REUSED", str(log_path), ""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    started = now_iso()
    write_json(checkpoint_path, {"status": "RUNNING", "command_hash": command_hash, "started_at": started, "command": command})
    try:
        proc = subprocess.run(
            command, cwd=project_root, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, check=False, timeout=timeout,
        )
        output = proc.stdout or ""
        status = "PASS" if proc.returncode == 0 else "FAILED"
        failure = "" if proc.returncode == 0 else f"installer exit code {proc.returncode}"
    except subprocess.TimeoutExpired as exc:
        output = (exc.stdout or "") if isinstance(exc.stdout, str) else (exc.stdout or b"").decode("utf-8", errors="replace")
        status = "FAILED"
        failure = f"installation timed out after {timeout} seconds"
    log_path.write_text(output, encoding="utf-8")
    write_json(checkpoint_path, {
        "status": status, "command_hash": command_hash, "started_at": started,
        "finished_at": now_iso(), "failure_reason": failure, "log": str(log_path),
    })
    return status, str(log_path), failure


def _doctor_delta(before: DoctorResult | None, after: DoctorResult, path: Path) -> str:
    old = {(row.category, row.name): row.status for row in before.rows} if before else {}
    rows = []
    for row in after.rows:
        prior = old.get((row.category, row.name), "NOT_CHECKED")
        rows.append({
            "category": row.category, "name": row.name, "before": prior,
            "after": row.status, "changed": "true" if prior != row.status else "false",
        })
    write_tsv(path, rows)
    return str(path)


def _run_doctor(args: dict[str, Any], project_root: Path, outdir: Path, tier: str, *, allow_execute: bool, smoke: bool) -> DoctorResult:
    return run_doctor(
        project_root=project_root,
        outdir=outdir,
        tools_manifest=args.get("tools_manifest"),
        reference_manifest=args.get("reference_manifest"),
        sample_manifest=args.get("sample_manifest"),
        profile=str(args.get("profile") or "local"),
        run_demo=bool(args.get("run_demo", smoke and tier in {"core", "prediction", "full"})),
        run_pytest=bool(args.get("run_pytest", False)),
        run_nextflow=bool(args.get("run_nextflow", False)),
        mini_smoke=bool(args.get("mini_smoke", smoke and tier in {"prediction", "full"})),
        release_audit=bool(args.get("release_audit", True)),
        allow_execute=allow_execute,
    )


def run_install_check(args: dict[str, Any]) -> dict[str, Any]:
    args = _apply_default_asset_source(args)
    case_id = safe_identifier(str(args.get("case_id") or "INSTALL"))
    mode = str(args.get("mode") or "verify").lower()
    tier = str(args.get("deployment_tier") or "core").lower()
    if tier not in TIER_REQUIRED_TOOLS:
        tier = "core"
    layout = RunLayout.create(args.get("outdir") or f"work/open-neo-install-check/{case_id}")
    result = MacroResult(
        "open-neo-install-check", case_id, new_run_id(case_id, "install"), mode,
        approval_required=mode in EXECUTION_MODES, approved=bool(args.get("approved", False)),
    )

    archive = Path(str(args.get("release_tarball") or "")) if args.get("release_tarball") else None
    project_root = Path(args.get("project_root") or ".").resolve()
    archive_hash = ""
    if archive:
        if not args.get("sha256"):
            result.blocking_issues.append(FailureCode.CHECKSUM_REQUIRED.value)
            result.steps.append(MacroStep("01", "release-checksum", "BLOCKED", "A release archive requires --sha256", failure_code=FailureCode.CHECKSUM_REQUIRED.value))
            result.finish("BLOCKED").write(layout.skill_result)
            return result.to_dict()
        if not archive.is_file():
            result.blocking_issues.append(FailureCode.CHECKSUM_FAILED.value)
            result.steps.append(MacroStep("01", "release-checksum", "FAILED", f"file not found: {archive}", failure_code=FailureCode.CHECKSUM_FAILED.value))
            result.finish("BLOCKED").write(layout.skill_result)
            return result.to_dict()
        archive_hash = sha256_file(archive)
        if archive_hash.lower() != str(args["sha256"]).lower():
            result.blocking_issues.append(FailureCode.CHECKSUM_FAILED.value)
            result.steps.append(MacroStep("01", "release-checksum", "FAILED", f"observed={archive_hash}", failure_code=FailureCode.CHECKSUM_FAILED.value))
            result.finish("BLOCKED").write(layout.skill_result)
            return result.to_dict()
        try:
            project_root, staging_outputs, reused = _stage_release(archive, archive_hash, layout)
        except ValueError as exc:
            result.blocking_issues.append(FailureCode.RELEASE_STAGING_FAILED.value)
            result.steps.append(MacroStep("01", "release-staging", "FAILED", str(exc), failure_code=FailureCode.RELEASE_STAGING_FAILED.value))
            result.finish("BLOCKED").write(layout.skill_result)
            return result.to_dict()
        result.outputs.update(staging_outputs)
        result.steps.append(MacroStep("01", "release-staging", "REUSED" if reused else "PASS", f"sha256={archive_hash}", outputs=staging_outputs))

    if not ((project_root / "pyproject.toml").is_file() or (project_root / "setup.py").is_file()):
        result.blocking_issues.append(FailureCode.PROJECT_ROOT_INVALID.value)
        result.steps.append(MacroStep("01", "project-root", "BLOCKED", f"not a project root: {project_root}", failure_code=FailureCode.PROJECT_ROOT_INVALID.value))
        result.finish("BLOCKED").write(layout.skill_result)
        return result.to_dict()

    audit(layout, "install_check.start", "START", mode=mode, tier=tier, project_root=str(project_root), release_sha256=archive_hash)
    inventory = _environment_inventory(project_root)
    write_tsv(layout.root / "environment_inventory.tsv", inventory)
    result.steps.append(MacroStep("02", "environment-preflight", "PASS", outputs={"environment_inventory": str(layout.root / "environment_inventory.tsv")}))
    result.outputs["environment_inventory"] = str(layout.root / "environment_inventory.tsv")

    templates = _write_local_manifest_templates(layout, project_root, args)
    result.outputs.update(templates)
    result.steps.append(MacroStep("03", "local-manifest-templates", "PASS", outputs=templates))
    if not args.get("tools_manifest"):
        args["tools_manifest"] = templates["tools_manifest_template"]
    if not args.get("reference_manifest"):
        args["reference_manifest"] = templates["reference_manifest_template"]
    if not args.get("deployment_reference_manifest"):
        args["deployment_reference_manifest"] = args["reference_manifest"]
    if templates.get("asset_manifest_local"):
        args["asset_manifest"] = templates["asset_manifest_local"]

    deploy_root = Path(str(args.get("deploy_root") or "/opt/neoag"))
    auto_before = configure_machine(
        project_root=project_root,
        tools_manifest=args["tools_manifest"],
        reference_manifest=args["reference_manifest"],
        outdir=layout.root / "auto_configuration_before",
        tools_root=args.get("tools_root") or deploy_root / "env_tool",
        reference_root=args.get("reference_root") or deploy_root / "refs",
        licensed_root=args.get("licensed_root") or deploy_root / "licensed_tools",
        run_smoke=mode == "verify",
        publish_local=False,
    )
    args["tools_manifest"] = auto_before.tools_manifest
    args["reference_manifest"] = auto_before.reference_manifest
    result.outputs.update({f"auto_before_{key}": value for key, value in auto_before.outputs.items()})
    result.steps.append(MacroStep("03b", "automatic-machine-configuration", auto_before.status, outputs=auto_before.outputs))

    if mode in EXECUTION_MODES and not result.approved:
        result.blocking_issues.append(FailureCode.APPROVAL_REQUIRED.value)
        result.steps.append(MacroStep("04", "approval-gate", "APPROVAL_REQUIRED", "Installation, repair or resume requires explicit approval", failure_code=FailureCode.APPROVAL_REQUIRED.value))
        result.finish("APPROVAL_REQUIRED").write(layout.skill_result)
        return result.to_dict()

    if mode in EXECUTION_MODES and bool(args.get("install_claude_code", False)) and not bool(args.get("allow_download", False)):
        result.blocking_issues.append(FailureCode.APPROVAL_REQUIRED.value)
        result.steps.append(MacroStep(
            "04", "claude-code-download-approval", "APPROVAL_REQUIRED",
            "Claude Code installation requires explicit --allow-download approval",
            failure_code=FailureCode.APPROVAL_REQUIRED.value,
        ))
        result.finish("APPROVAL_REQUIRED").write(layout.skill_result)
        return result.to_dict()

    if mode in EXECUTION_MODES and not bool(args.get("no_sync_assets", False)):
        if not args.get("asset_manifest"):
            result.blocking_issues.append(FailureCode.ASSET_SOURCE_UNCONFIGURED.value)
            result.steps.append(MacroStep("04", "asset-source", "BLOCKED", "No asset manifest is available", failure_code=FailureCode.ASSET_SOURCE_UNCONFIGURED.value))
            result.finish("BLOCKED").write(layout.skill_result)
            return result.to_dict()
        missing_sources = _required_asset_sources_missing(args["asset_manifest"], str(args.get("asset_source_host") or ""))
        if missing_sources:
            detail = "Required asset sources are unavailable: " + "; ".join(missing_sources[:8])
            result.blocking_issues.append(FailureCode.ASSET_SOURCE_UNCONFIGURED.value)
            result.steps.append(MacroStep("04", "asset-source", "BLOCKED", detail, failure_code=FailureCode.ASSET_SOURCE_UNCONFIGURED.value))
            result.finish("BLOCKED").write(layout.skill_result)
            return result.to_dict()

    before: DoctorResult | None = None
    deployment_failure = ""
    if mode in EXECUTION_MODES:
        before = _run_doctor(args, project_root, layout.root / "doctor_before", tier, allow_execute=False, smoke=False)
        result.outputs.update({f"doctor_before_{key}": value for key, value in before.outputs.items()})
        result.steps.append(MacroStep("04", "doctor-before-install", before.status, outputs=before.outputs))

    deploy_command = _deployment_command(args, project_root, layout, execute=mode in EXECUTION_MODES)
    write_json(layout.root / "deployment_command.json", {"command": deploy_command, "execute": mode in EXECUTION_MODES})
    result.outputs["deployment_command"] = str(layout.root / "deployment_command.json")
    if mode in EXECUTION_MODES:
        if not Path(deploy_command[1]).is_file():
            deployment_failure = f"installer not found: {deploy_command[1]}"
            result.steps.append(MacroStep("05", "portable-deployment", "BLOCKED", deployment_failure, failure_code=FailureCode.CORE_INSTALL_FAILED.value))
        else:
            status, deploy_log, deployment_failure = _run_deployment(
                deploy_command, project_root, layout.logs / "portable_deployment.log",
                layout.root / "deployment_checkpoint.json",
                timeout=int(args.get("install_timeout") or 7200), resume=mode == "resume",
            )
            failure_code = ""
            if deployment_failure:
                failure_code = FailureCode.INSTALL_TIMEOUT.value if "timed out" in deployment_failure else FailureCode.CORE_INSTALL_FAILED.value
            result.steps.append(MacroStep("05", "portable-deployment", status, deployment_failure, outputs={"log": deploy_log, "checkpoint": str(layout.root / "deployment_checkpoint.json")}, failure_code=failure_code))
            result.outputs.update({"deployment_log": deploy_log, "deployment_checkpoint": str(layout.root / "deployment_checkpoint.json")})
    else:
        result.steps.append(MacroStep("05", "portable-deployment", "PLANNED", "No installation command executed", outputs={"command": str(layout.root / "deployment_command.json")}))

    claude_requested = bool(args.get("install_claude_code", False))
    claude_status, claude_outputs = _collect_claude_code_status(
        layout, requested=claude_requested, executed=mode in EXECUTION_MODES,
    )
    result.outputs.update(claude_outputs)
    if claude_requested:
        claude_step_status = {
            "READY": "PASS", "PLANNED": "PLANNED", "UNVERIFIED": "FAILED", "MISSING_REPORT": "FAILED",
        }[claude_status["status"]]
        result.steps.append(MacroStep(
            "05c", "claude-code-readiness", claude_step_status,
            f"status={claude_status['status']}; version={claude_status['version'] or 'unknown'}",
            outputs=claude_outputs,
            failure_code=FailureCode.CORE_INSTALL_FAILED.value if claude_step_status == "FAILED" else "",
        ))
        if mode in EXECUTION_MODES and claude_status["status"] != "READY" and not deployment_failure:
            deployment_failure = f"Claude Code verification failed: {claude_status['status']}"

    if mode in EXECUTION_MODES and not deployment_failure:
        auto_final = configure_machine(
            project_root=project_root,
            tools_manifest=args["tools_manifest"],
            reference_manifest=args["reference_manifest"],
            outdir=layout.root / "auto_configuration",
            tools_root=args.get("tools_root") or deploy_root / "env_tool",
            reference_root=args.get("reference_root") or deploy_root / "refs",
            licensed_root=args.get("licensed_root") or deploy_root / "licensed_tools",
            run_smoke=True,
            publish_local=True,
        )
        args["tools_manifest"] = auto_final.tools_manifest
        args["reference_manifest"] = auto_final.reference_manifest
        result.outputs.update(auto_final.outputs)
        result.steps.append(MacroStep("05b", "automatic-machine-configuration-after-install", auto_final.status, outputs=auto_final.outputs))
    else:
        auto_final = auto_before
        result.outputs.update(auto_final.outputs)

    doctor = _run_doctor(args, project_root, layout.root / "doctor", tier, allow_execute=mode != "plan", smoke=True)
    result.steps.append(MacroStep("06", "doctor-after-install", doctor.status, outputs=doctor.outputs))
    result.outputs.update({f"doctor_{key}": value for key, value in doctor.outputs.items()})
    delta_path = _doctor_delta(before, doctor, layout.root / "deployment_delta.tsv")
    result.outputs["deployment_delta"] = delta_path

    tier_status, requirement_rows = _assess_tier(tier, doctor.rows, args.get("reference_manifest"))
    write_tsv(layout.root / "tier_requirements.tsv", requirement_rows)
    result.outputs["tier_requirements"] = str(layout.root / "tier_requirements.tsv")
    missing = [row for row in requirement_rows if row["required"] == "true" and row["status"] != "OK"]

    final_status = tier_status
    if deployment_failure:
        final_status = "FAILED"
        result.blocking_issues.append(result.steps[-2].failure_code or FailureCode.CORE_INSTALL_FAILED.value)
    elif doctor.status == "UNSAFE":
        final_status = "UNSAFE"
    elif doctor.status == "BLOCKED" and tier_status != "READY":
        final_status = "BLOCKED"
    elif doctor.status in {"PARTIAL", "BLOCKED"} and tier_status == "READY":
        result.warnings.append(f"doctor_status={doctor.status}; optional tools or references outside tier are incomplete")
    result.warnings.extend(f"tier_missing:{row['kind']}:{row['requirement']}" for row in missing)

    report_rows = [{
        "deployment_tier": tier, "tier_status": final_status,
        "doctor_status": doctor.status, "missing_required_count": len(missing),
        "project_root": str(project_root), "release_sha256": archive_hash,
        "claude_code_requested": claude_status["requested"],
        "claude_code_status": claude_status["status"],
        "claude_code_version": claude_status["version"],
        "claude_code_binary": claude_status["binary"],
        "claude_code_report": claude_status["report"],
    }]
    write_tsv(layout.root / "deployment_status.tsv", report_rows)
    md = [
        "# Open-Neo installation and environment check", "",
        f"- Deployment tier: **{tier}**",
        f"- Tier status: **{final_status}**",
        f"- Doctor status: **{doctor.status}**",
        f"- Project root: `{project_root}`",
        f"- Release SHA256: `{archive_hash or 'source-checkout'}`", "",
        "## Tier requirements", "", markdown_table(requirement_rows, max_rows=120), "",
        "## Environment inventory", "", markdown_table(inventory, max_rows=30), "",
        "## Installation delta", "", f"See `{delta_path}`.", "",
        "## Claude Code", "",
        markdown_table([claude_status]), "",
        "## Boundary", "",
        "READY means all required tools, capability groups and reference assets for the requested tier were observed. Licensed tools remain machine-local and are never redistributed.",
    ]
    (layout.root / "deployment_report.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    result.outputs.update({"deployment_report": str(layout.root / "deployment_report.md"), "deployment_status": str(layout.root / "deployment_status.tsv")})
    result.provenance = {
        "python": platform.python_version(), "project_root": str(project_root),
        "deployment_tier": tier, "release_sha256": archive_hash,
    }
    update_case_state(layout, case_id=case_id, current_intent="install_check", deployment_tier=tier, status=final_status, outputs=result.outputs)
    audit(layout, "install_check.finish", final_status, missing=missing, deployment_failure=deployment_failure)
    result.finish(final_status).write(layout.skill_result)
    return result.to_dict()
