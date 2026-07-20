from __future__ import annotations

import argparse
import os
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

from .audit import AuditLogger, run_command
from .io_utils import ensure_dir, markdown_table, now_iso, write_json, write_tsv
from .manifests import load_manifest, manifest_paths
from .release_audit import scan_release_boundary, write_release_audit

DEFAULT_TOOL_EXECUTABLES = {
    "python": sys.executable,
    "neoag": "neoag",
    "neoag-agent": "neoag-agent",
    "neoag-llm-agent": "neoag-llm-agent",
    "nextflow": "nextflow",
    "java": "java",
    "bwa": "bwa",
    "star": "STAR",
    "gatk": "gatk",
    "vep": "vep",
    "netmhcpan": "netMHCpan",
    "netmhcstabpan": "NetMHCstabpan",
    "mhcflurry": "mhcflurry-predict",
    "prime": "PRIME",
    "bigmhc_im": "bigmhc_predict",
    "deepimmuno": "deepimmuno-cnn.py",
    "optitype": "OptiTypePipeline.py",
    "hla_la": "HLA-LA.pl",
    "spechla": "SpecHLA",
    "lohhla": "LOHHLA",
    "facets": "runFACETS.R",
    "purple": "purple",
    "ascat": "ascat.R",
    "sequenza": "sequenza-utils",
    "arriba": "arriba",
    "star_fusion": "STAR-Fusion",
    "fusioncatcher": "fusioncatcher",
    "easyfuse": "easyfuse",
    "pyclone": "pyclone",
    "pyclone_vi": "pyclone-vi",
}

CRITICAL_REFERENCES = ["reference_fasta", "gencode_gtf", "vep_cache"]
LICENSE_RESTRICTED = {"netmhcpan", "netmhcstabpan"}


@dataclass
class CheckRow:
    category: str
    name: str
    status: str
    message: str = ""
    path: str = ""
    blocking: str = "false"


@dataclass
class DoctorResult:
    status: str
    checked_at: str
    profile: str
    rows: list[CheckRow] = field(default_factory=list)
    outputs: dict[str, str] = field(default_factory=dict)
    summary: dict[str, int] = field(default_factory=dict)


def _row(category: str, name: str, status: str, message: str = "", path: str = "", blocking: bool = False) -> CheckRow:
    return CheckRow(category, name, status, message, path, "true" if blocking else "false")


def _flatten_tools(tools_manifest: dict[str, Any]) -> dict[str, str]:
    out = dict(DEFAULT_TOOL_EXECUTABLES)
    tools = tools_manifest.get("tools") if isinstance(tools_manifest.get("tools"), dict) else tools_manifest
    if isinstance(tools, dict):
        for name, spec in tools.items():
            if isinstance(spec, str):
                out[str(name)] = spec
            elif isinstance(spec, dict):
                exe = spec.get("executable") or spec.get("path") or spec.get("cmd")
                if exe:
                    out[str(name)] = str(exe)
    return out


def _check_executable(name: str, exe: str) -> CheckRow:
    if os.path.sep in exe or exe.startswith("."):
        p = Path(exe)
        if p.exists():
            return _row("tool", name, "OK", "executable path exists", str(p))
        return _row("tool", name, "MISSING", "executable path not found", str(p), blocking=name in {"python"})
    found = shutil.which(exe)
    if found:
        return _row("tool", name, "OK", "found on PATH", found)
    severity = "LICENSE_REQUIRED" if name in LICENSE_RESTRICTED else "MISSING"
    msg = "licensed tool; install locally and record path in tools_manifest" if name in LICENSE_RESTRICTED else "not found on PATH"
    return _row("tool", name, severity, msg, exe, blocking=name in {"python"})


def _check_manifest_paths(data: dict[str, Any], category: str) -> list[CheckRow]:
    rows: list[CheckRow] = []
    for key, val in manifest_paths(data):
        if val.startswith("docker://") or val.startswith("oras://") or ":" in val and not val.startswith("/") and not val.startswith("."):
            rows.append(_row(category, key, "INFO", "non-local reference/image; not checked as a filesystem path", val))
            continue
        p = Path(os.path.expandvars(os.path.expanduser(val)))
        if p.exists():
            rows.append(_row(category, key, "OK", "path exists", str(p)))
        else:
            blocking = any(k in key.lower() for k in CRITICAL_REFERENCES)
            rows.append(_row(category, key, "MISSING", "path does not exist", str(p), blocking=blocking))
    return rows


def _walk_manifest_items(data: Any, prefix: str = "") -> list[tuple[str, Any]]:
    rows: list[tuple[str, Any]] = []
    if isinstance(data, dict):
        for k, v in data.items():
            key = f"{prefix}.{k}" if prefix else str(k)
            rows.extend(_walk_manifest_items(v, key))
    elif isinstance(data, list):
        for i, v in enumerate(data):
            rows.extend(_walk_manifest_items(v, f"{prefix}[{i}]"))
    else:
        rows.append((prefix, data))
    return rows


def _manifest_value(data: dict[str, Any], *needles: str) -> str | None:
    lowered = [n.lower() for n in needles]
    for key, val in _walk_manifest_items(data):
        if not isinstance(val, str):
            continue
        key_l = key.lower()
        if all(n in key_l for n in lowered):
            return val
    return None


def _path_exists(value: str | None) -> tuple[Path | None, bool]:
    if not value or value.startswith(("docker://", "oras://")):
        return None, False
    p = Path(os.path.expandvars(os.path.expanduser(value)))
    return p, p.exists()


def _check_expected_sidecar(path_value: str | None, suffix: str, name: str, *, blocking: bool = False) -> CheckRow:
    p, exists = _path_exists(path_value)
    if not p:
        return _row("reference_specific", name, "INFO", "reference path not declared", "")
    sidecar = Path(str(p) + suffix)
    if exists and sidecar.exists():
        return _row("reference_specific", name, "OK", "sidecar exists", str(sidecar))
    if exists:
        return _row("reference_specific", name, "MISSING", "main file exists but sidecar is missing", str(sidecar), blocking=blocking)
    return _row("reference_specific", name, "MISSING", "main reference file is missing", str(p), blocking=blocking)


def _check_reference_specifics(ref_data: dict[str, Any]) -> list[CheckRow]:
    rows: list[CheckRow] = []
    fasta = _manifest_value(ref_data, "reference", "fasta") or _manifest_value(ref_data, "fasta")
    vep_cache = _manifest_value(ref_data, "vep", "cache")
    gtf = _manifest_value(ref_data, "gencode", "gtf") or _manifest_value(ref_data, "gtf")
    facets_vcf = _manifest_value(ref_data, "facets", "vcf") or _manifest_value(ref_data, "snp", "vcf")

    rows.append(_check_expected_sidecar(fasta, ".fai", "reference_fasta.fai", blocking=True))
    if fasta:
        p = Path(os.path.expandvars(os.path.expanduser(fasta)))
        dict_candidates = [p.with_suffix(".dict"), Path(str(p).rsplit(".fasta", 1)[0].rsplit(".fa", 1)[0] + ".dict")]
        if p.exists() and any(x.exists() for x in dict_candidates):
            rows.append(_row("reference_specific", "reference_fasta.dict", "OK", "sequence dictionary exists", str(next(x for x in dict_candidates if x.exists()))))
        elif p.exists():
            rows.append(_row("reference_specific", "reference_fasta.dict", "MISSING", "sequence dictionary is missing", str(dict_candidates[0]), blocking=True))

    p, exists = _path_exists(gtf)
    if p:
        rows.append(_row("reference_specific", "gencode_gtf", "OK" if exists else "MISSING", "GTF exists" if exists else "GTF missing", str(p), blocking=not exists))

    p, exists = _path_exists(vep_cache)
    if p:
        if exists:
            has_species = p.name.lower() == "homo_sapiens" or any(x.name.lower() == "homo_sapiens" for x in p.parents) or (p.is_dir() and any(x.is_dir() and "grch38" in x.name.lower() for x in p.iterdir()))
            rows.append(_row("reference_specific", "vep_cache_layout", "OK" if has_species else "WARN", "VEP cache path exists" if has_species else "VEP cache exists but species/build layout was not recognized", str(p)))
        else:
            rows.append(_row("reference_specific", "vep_cache_layout", "MISSING", "VEP cache path missing", str(p), blocking=True))
    else:
        rows.append(_row("reference_specific", "vep_cache_layout", "INFO", "vep_cache not declared", ""))

    rows.append(_check_expected_sidecar(facets_vcf, ".tbi", "facets_snp_vcf.tbi"))

    for label, needles in {
        "lohhla_reference": ("lohhla",),
        "ascat_loci": ("ascat", "loci"),
        "ascat_alleles": ("ascat", "alleles"),
        "arriba_reference": ("arriba",),
        "hla_reference": ("hla", "reference"),
    }.items():
        value = _manifest_value(ref_data, *needles)
        p, exists = _path_exists(value)
        if p:
            rows.append(_row("reference_specific", label, "OK" if exists else "MISSING", "declared path exists" if exists else "declared path missing", str(p)))
        else:
            rows.append(_row("reference_specific", label, "INFO", "not declared in reference manifest", ""))
    return rows


def _check_tool_specifics(tools_data: dict[str, Any]) -> list[CheckRow]:
    rows: list[CheckRow] = []
    for name in ["neoag-doctor", "neoag-release-audit", "neoag-pipeline-full", "neoag-gateway"]:
        found = shutil.which(name)
        rows.append(_row("entrypoint", name, "OK" if found else "INFO", "console script found" if found else "console script not on PATH; run pip install -e . to expose it", found or name))

    for env_name, label in [("NETMHCPAN_HOME", "netmhcpan_home"), ("NETMHCpan", "netmhcpan_env_legacy"), ("MHCFLURRY_DATA_DIR", "mhcflurry_data_dir"), ("SPECHLA_HOME", "spechla_home"), ("HLA_LA_HOME", "hla_la_home"), ("PURPLE_HOME", "purple_home")]:
        value = os.environ.get(env_name)
        if value:
            p = Path(os.path.expandvars(os.path.expanduser(value)))
            rows.append(_row("tool_specific", label, "OK" if p.exists() else "MISSING", f"{env_name} set" if p.exists() else f"{env_name} set but path missing", str(p)))
        else:
            rows.append(_row("tool_specific", label, "INFO", f"{env_name} not set; manifest executable/path may still be sufficient", ""))

    tools = tools_data.get("tools") if isinstance(tools_data.get("tools"), dict) else tools_data
    if isinstance(tools, dict):
        for name in ["vep", "netmhcpan", "mhcflurry", "lohhla", "facets", "ascat", "arriba", "easyfuse", "spechla", "hla_la"]:
            spec = tools.get(name)
            if isinstance(spec, dict):
                mode = spec.get("mode", "")
                if spec.get("license_required"):
                    rows.append(_row("tool_specific", f"{name}.license", "LICENSE_REQUIRED", "licensed/restricted tool; verify local license and do not redistribute binaries", str(mode)))
                elif mode:
                    rows.append(_row("tool_specific", f"{name}.mode", "INFO", f"declared install mode: {mode}", ""))
            elif spec is None:
                rows.append(_row("tool_specific", f"{name}.manifest", "INFO", "tool not declared in tools manifest; default PATH check was used", ""))
    return rows


def _check_workflow_readiness(project_root: Path, profile: str) -> list[CheckRow]:
    rows: list[CheckRow] = []
    for dirname in ["bin", "work", "results"]:
        p = project_root / dirname
        if not p.exists():
            if dirname in {"work", "results"}:
                try:
                    p.mkdir(parents=True, exist_ok=True)
                except Exception as exc:
                    rows.append(_row("workflow", f"{dirname}_write", "MISSING", f"directory missing and could not be created: {exc}", str(p), blocking=True))
                    continue
            else:
                rows.append(_row("workflow", dirname, "MISSING", "directory not found", str(p), blocking=dirname == "bin"))
                continue
        if dirname == "bin":
            scripts = [x for x in p.iterdir() if x.is_file()]
            non_exec = [x.name for x in scripts if x.suffix in {"", ".sh", ".py"} and not os.access(x, os.X_OK)]
            status = "OK" if not non_exec else "WARN"
            msg = "bin scripts executable" if not non_exec else "non-executable bin files: " + ", ".join(non_exec[:10])
            rows.append(_row("workflow", "bin_permissions", status, msg, str(p)))
        else:
            test = p / ".neoag_doctor_write_test"
            try:
                test.write_text("ok\n", encoding="utf-8")
                test.unlink()
                rows.append(_row("workflow", f"{dirname}_write", "OK", "directory writable", str(p)))
            except Exception as exc:
                rows.append(_row("workflow", f"{dirname}_write", "FAIL", str(exc), str(p), blocking=True))
    for prof in ["local", "slurm", profile]:
        candidates = [
            project_root / "conf" / f"{prof}.config",
            project_root / "conf" / f"{prof}.conf",
            project_root / "profiles" / f"{prof}.config",
            project_root / "nextflow.config",
        ]
        if any(x.exists() for x in candidates):
            rows.append(_row("workflow", f"profile_{prof}", "OK", "profile or nextflow.config found", str(next(x for x in candidates if x.exists()))))
        else:
            rows.append(_row("workflow", f"profile_{prof}", "INFO", "profile not found; required only for this execution backend", ",".join(str(x) for x in candidates)))
    return rows


def _run_version_checks(project_root: Path, logger: AuditLogger, *, allow_execute: bool) -> list[CheckRow]:
    rows: list[CheckRow] = []
    commands = {
        "java_version": ["java", "-version"],
        "nextflow_version": ["nextflow", "-version"],
        "gatk_mutect2_help": ["gatk", "Mutect2", "--help"],
    }
    for name, cmd in commands.items():
        if not shutil.which(cmd[0]):
            rows.append(_row("workflow" if name != "gatk_mutect2_help" else "tool_smoke", name, "MISSING", f"{cmd[0]} not found on PATH", cmd[0], blocking=name == "java_version"))
            continue
        res = run_command(cmd, cwd=project_root, timeout=60, logger=logger, risk_level="LOW", allow_execute=allow_execute)
        if res.get("dry_run"):
            rows.append(_row("workflow" if name != "gatk_mutect2_help" else "tool_smoke", name, "DRY_RUN", " ".join(cmd), " ".join(cmd)))
        else:
            output = (res.get("stdout") or res.get("stderr") or "")[-1000:]
            rows.append(_row("workflow" if name != "gatk_mutect2_help" else "tool_smoke", name, "OK" if res.get("ok") else "FAIL", output, " ".join(cmd), blocking=False))
    return rows


def _run_tool_mini_smoke(project_root: Path, outdir: Path, logger: AuditLogger, ref_data: dict[str, Any], *, allow_execute: bool) -> list[CheckRow]:
    rows: list[CheckRow] = []
    od = ensure_dir(outdir)
    # NetMHCpan smoke: one peptide and one common allele. Licensed installs may still fail cleanly.
    if shutil.which("netMHCpan"):
        pep = od / "netmhcpan.pep.txt"
        pep.write_text("SYFPEITHI\n", encoding="utf-8")
        res = run_command(["netMHCpan", "-p", str(pep), "-a", "HLA-A02:01"], cwd=project_root, timeout=120, logger=logger, risk_level="LOW", allow_execute=allow_execute)
        rows.append(_row("smoke", "netmhcpan_1peptide_1hla", "OK" if res.get("ok") else ("DRY_RUN" if res.get("dry_run") else "FAIL"), (res.get("stdout") or res.get("stderr") or "")[-1000:], "netMHCpan -p peptide -a HLA-A02:01"))
    else:
        rows.append(_row("smoke", "netmhcpan_1peptide_1hla", "MISSING", "netMHCpan not found", "netMHCpan"))
    # MHCflurry model load via predictor help/import-light command.
    if shutil.which("mhcflurry-predict"):
        res = run_command(["mhcflurry-predict", "--help"], cwd=project_root, timeout=60, logger=logger, risk_level="LOW", allow_execute=allow_execute)
        rows.append(_row("smoke", "mhcflurry_model_cli", "OK" if res.get("ok") else ("DRY_RUN" if res.get("dry_run") else "FAIL"), (res.get("stdout") or res.get("stderr") or "")[-1000:], "mhcflurry-predict --help"))
    else:
        rows.append(_row("smoke", "mhcflurry_model_cli", "MISSING", "mhcflurry-predict not found", "mhcflurry-predict"))
    # VEP offline cache tiny VCF smoke.
    vep_cache = _manifest_value(ref_data, "vep", "cache")
    fasta = _manifest_value(ref_data, "reference", "fasta") or _manifest_value(ref_data, "fasta")
    if shutil.which("vep") and vep_cache:
        tiny = od / "vep_tiny.vcf"
        tiny.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t10000\t.\tA\tG\t.\t.\t.\n", encoding="utf-8")
        cmd = ["vep", "--offline", "--cache", "--dir_cache", str(Path(vep_cache).parent.parent if Path(vep_cache).name.lower().startswith("115") else Path(vep_cache)), "--input_file", str(tiny), "--output_file", str(od / "vep_tiny.out"), "--force_overwrite"]
        if fasta:
            cmd += ["--fasta", fasta]
        res = run_command(cmd, cwd=project_root, timeout=180, logger=logger, risk_level="LOW", allow_execute=allow_execute)
        rows.append(_row("smoke", "vep_tiny_offline_cache", "OK" if res.get("ok") else ("DRY_RUN" if res.get("dry_run") else "FAIL"), (res.get("stdout") or res.get("stderr") or "")[-1000:], "vep tiny offline cache"))
    else:
        rows.append(_row("smoke", "vep_tiny_offline_cache", "MISSING", "vep or vep_cache not available", ""))
    # FACETS/LOHHLA config-level smoke.
    rows.append(_row("smoke", "facets_fixture_ready", "OK" if shutil.which("snp-pileup") or shutil.which("runFACETS.R") else "MISSING", "snp-pileup/runFACETS.R availability checked", "snp-pileup/runFACETS.R"))
    rows.append(_row("smoke", "lohhla_reference_ready", "OK" if _manifest_value(ref_data, "lohhla") else "MISSING", "LOHHLA reference declared" if _manifest_value(ref_data, "lohhla") else "LOHHLA reference not declared", str(_manifest_value(ref_data, "lohhla") or "")))
    for label, needles in {"arriba_reference_ready": ("arriba",), "star_fusion_reference_ready": ("star", "fusion")}.items():
        value = _manifest_value(ref_data, *needles)
        p, exists = _path_exists(value)
        rows.append(_row("smoke", label, "OK" if exists else ("MISSING" if value else "INFO"), "reference path exists" if exists else "reference not declared or missing", str(p or value or "")))
    return rows


def _write_recommended_fixes(rows: list[dict[str, Any]], outdir: Path) -> Path:
    lines = ["# Recommended fixes", ""]
    problematic = [r for r in rows if r.get("status") in {"MISSING", "FAIL", "WARN", "LICENSE_REQUIRED", "UNSAFE"}]
    if not problematic:
        lines.append("No blocking or warning-level issues detected.")
    for r in problematic:
        name = r.get("name", "")
        status = r.get("status", "")
        path = r.get("path", "")
        msg = r.get("message", "")
        fix = "Review the manifest entry and installation path."
        if "license" in name or status == "LICENSE_REQUIRED":
            fix = "Install/configure this licensed tool locally; do not redistribute binaries or licensed data."
        elif "reference" in r.get("category", "") or "fasta" in name or "cache" in name:
            fix = "Download or symlink the required reference file/cache and update reference_manifest."
        elif r.get("category") == "tool":
            fix = "Install the tool or set its executable path in tools_manifest."
        elif r.get("category") == "workflow":
            fix = "Check permissions, profile configuration, Java/Nextflow install, or writable work/results directories."
        elif r.get("category") == "release":
            fix = "Remove private paths, patient hints, generated caches, or runtime artifacts from the release boundary."
        lines.extend([f"## {name} [{status}]", "", f"- Message: {msg}", f"- Path: {path}", f"- Suggested fix: {fix}", ""])
    p = outdir / "recommended_fixes.md"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def _run_mini_smoke(project_root: Path, outdir: Path, logger: AuditLogger) -> list[CheckRow]:
    rows: list[CheckRow] = []
    modules = [
        "neoag.cli",
        "neoag.controlled_execution.doctor",
        "neoag.controlled_execution.release_audit",
        "neoag.controlled_execution.pipeline_runner",
    ]
    for module in modules:
        try:
            __import__(module)
            rows.append(_row("mini_smoke", module, "OK", "module importable"))
        except Exception as exc:
            rows.append(_row("mini_smoke", module, "FAIL", str(exc), blocking=True))
    for label, cmd in {
        "neoag_cli_help": [sys.executable, "-m", "neoag.cli", "--help"],
        "controlled_execution_doctor_help": [sys.executable, "-m", "neoag.controlled_execution.doctor", "--help"],
        "controlled_execution_release_audit_help": [sys.executable, "-m", "neoag.controlled_execution.release_audit", "--help"],
    }.items():
        res = run_command(cmd, cwd=project_root, timeout=60, logger=logger, risk_level="LOW", allow_execute=True)
        rows.append(_row("mini_smoke", label, "OK" if res.get("ok") else "FAIL", (res.get("stdout") or res.get("stderr") or "")[-1000:], " ".join(cmd), blocking=not res.get("ok")))
    return rows


def _run_core_smoke(project_root: Path, outdir: Path, logger: AuditLogger, *, run_demo: bool, run_pytest: bool, run_nextflow: bool, allow_execute: bool) -> list[CheckRow]:
    rows: list[CheckRow] = []
    if run_pytest:
        res = run_command([sys.executable, "-m", "pytest", "-q"], cwd=project_root, timeout=300, logger=logger, risk_level="LOW", allow_execute=allow_execute)
        rows.append(_row("smoke", "pytest", "OK" if res.get("ok") else ("DRY_RUN" if res.get("dry_run") else "FAIL"), (res.get("stdout") or res.get("stderr") or "")[-1000:]))
    if run_demo:
        demo_out = outdir / "run_demo"
        res = run_command([sys.executable, "-m", "neoag.cli", "run-demo", "--sample-id", "DOCTOR", "--outdir", str(demo_out)], cwd=project_root, timeout=300, logger=logger, risk_level="LOW", allow_execute=allow_execute)
        rows.append(_row("smoke", "run-demo", "OK" if res.get("ok") else ("DRY_RUN" if res.get("dry_run") else "FAIL"), (res.get("stdout") or res.get("stderr") or "")[-1000:], str(demo_out)))
    if run_nextflow:
        wrapper = project_root / "bin" / "neoag-nextflow"
        if not wrapper.exists():
            rows.append(_row("smoke", "nextflow_fixture", "MISSING", "bin/neoag-nextflow not found", str(wrapper)))
        else:
            res = run_command([str(wrapper), "run", "workflows/main.nf", "--pvac_files", "data/fixtures/pvacseq_aggregated.tsv", "--outdir", str(outdir / "nextflow_demo")], cwd=project_root, timeout=600, logger=logger, risk_level="MEDIUM", allow_execute=allow_execute)
            rows.append(_row("smoke", "nextflow_fixture", "OK" if res.get("ok") else ("DRY_RUN" if res.get("dry_run") else "FAIL"), (res.get("stdout") or res.get("stderr") or "")[-1000:]))
    return rows


def run_doctor(
    *,
    project_root: str | Path = ".",
    outdir: str | Path,
    tools_manifest: str | Path | None = None,
    reference_manifest: str | Path | None = None,
    sample_manifest: str | Path | None = None,
    profile: str = "local",
    run_demo: bool = False,
    run_pytest: bool = False,
    run_nextflow: bool = False,
    mini_smoke: bool = False,
    release_audit: bool = True,
    allow_execute: bool = True,
) -> DoctorResult:
    root = Path(project_root).resolve()
    od = ensure_dir(outdir)
    logger = AuditLogger(od / "audit_log.jsonl")
    rows: list[CheckRow] = []
    logger.log("doctor.start", "START", f"project_root={root}")

    # Core importability and package entrypoints.
    try:
        import neoag  # noqa: F401
        rows.append(_row("core", "python_import_neoag", "OK", "neoag importable"))
    except Exception as exc:
        rows.append(_row("core", "python_import_neoag", "FAIL", str(exc), blocking=True))
    for name in ["python", "neoag", "neoag-agent", "neoag-llm-agent"]:
        rows.append(_check_executable(name, DEFAULT_TOOL_EXECUTABLES[name]))
    rows.extend(_check_workflow_readiness(root, profile))
    rows.extend(_run_version_checks(root, logger, allow_execute=allow_execute))

    # Manifest loading.
    for label, path in [("tools_manifest", tools_manifest), ("reference_manifest", reference_manifest), ("sample_manifest", sample_manifest)]:
        res = load_manifest(path)
        if path:
            rows.append(_row("manifest", label, "OK" if res.exists and not res.error else "FAIL", res.error or "loaded", res.path, blocking=(label == "reference_manifest" and bool(res.error))))
    tools_data = load_manifest(tools_manifest).data if tools_manifest else {}
    ref_data = load_manifest(reference_manifest).data if reference_manifest else {}
    sample_data = load_manifest(sample_manifest).data if sample_manifest else {}

    # Tool checks.
    for name, exe in sorted(_flatten_tools(tools_data).items()):
        if name in {"python", "neoag", "neoag-agent", "neoag-llm-agent"}:
            continue
        rows.append(_check_executable(name, exe))
    rows.extend(_check_tool_specifics(tools_data))

    # Reference/sample path checks.
    rows.extend(_check_manifest_paths(ref_data, "reference"))
    rows.extend(_check_manifest_paths(sample_data, "sample_input"))
    rows.extend(_check_reference_specifics(ref_data))

    # Minimal reference sanity.
    if reference_manifest and not ref_data:
        rows.append(_row("reference", "reference_manifest_content", "FAIL", "manifest could not be parsed or is empty", str(reference_manifest), blocking=True))

    # Release boundary audit.
    if release_audit:
        audit_result = scan_release_boundary(root)
        write_release_audit(audit_result, od / "release_audit")
        rows.append(_row("release", "boundary", audit_result["status"], f"private_path_hits={audit_result['summary']['private_path_hits']}; cache_hits={audit_result['summary']['cache_hits']}", str(root), blocking=audit_result["status"] == "UNSAFE"))

    # Optional smoke tests.
    if mini_smoke:
        rows.extend(_run_mini_smoke(root, od / "mini_smoke", logger))
        rows.extend(_run_tool_mini_smoke(root, od / "tool_smoke", logger, ref_data, allow_execute=allow_execute))
    rows.extend(_run_core_smoke(root, od / "smoke", logger, run_demo=run_demo, run_pytest=run_pytest, run_nextflow=run_nextflow, allow_execute=allow_execute))

    counts: dict[str, int] = {}
    for r in rows:
        counts[r.status] = counts.get(r.status, 0) + 1
    if any(r.status in {"FAIL"} and r.blocking == "true" for r in rows):
        status = "BLOCKED"
    elif any(r.status in {"UNSAFE"} for r in rows):
        status = "UNSAFE"
    elif any(r.status in {"MISSING", "LICENSE_REQUIRED", "FAIL", "WARN"} for r in rows):
        status = "PARTIAL"
    else:
        status = "READY"
    result = DoctorResult(status=status, checked_at=now_iso(), profile=profile, rows=rows, summary=counts)
    outputs = write_doctor_outputs(result, od)
    result.outputs = outputs
    logger.log("doctor.finish", status, metadata={"summary": counts, "outputs": outputs})
    return result


def write_doctor_outputs(result: DoctorResult, outdir: str | Path) -> dict[str, str]:
    od = ensure_dir(outdir)
    rows = [r.__dict__ for r in result.rows]
    fields = ["category", "name", "status", "message", "path", "blocking"]
    write_tsv(od / "doctor_checks.tsv", rows, fields)
    tool_rows = [r for r in rows if r.get("category") in {"tool", "tool_specific", "entrypoint"}]
    reference_rows = [r for r in rows if r.get("category") in {"reference", "reference_specific", "sample_input"}]
    smoke_rows = [r for r in rows if r.get("category") in {"smoke", "mini_smoke", "tool_smoke"}]
    write_tsv(od / "tool_status.tsv", tool_rows, fields)
    write_tsv(od / "reference_status.tsv", reference_rows, fields)
    write_tsv(od / "smoke_tests.tsv", smoke_rows, fields)
    blocking = [r for r in rows if r.get("blocking") == "true" and r.get("status") not in {"OK", "INFO"}]
    write_tsv(od / "blocking_issues.tsv", blocking, fields)
    fixes = _write_recommended_fixes(rows, od)
    write_json(od / "doctor_status.json", {"status": result.status, "checked_at": result.checked_at, "profile": result.profile, "summary": result.summary, "checks": rows})
    md = [
        "# NeoAg Doctor",
        "",
        f"Status: **{result.status}**",
        f"Checked at: {result.checked_at}",
        "",
        "## Summary",
        markdown_table([result.summary]),
        "",
        "## Blocking issues",
        markdown_table(blocking),
        "",
        "## Tool status",
        markdown_table(tool_rows, fields, max_rows=80),
        "",
        "## Reference status",
        markdown_table(reference_rows, fields, max_rows=80),
        "",
        "## Smoke tests",
        markdown_table(smoke_rows, fields, max_rows=80),
        "",
        "## All checks",
        markdown_table(rows, fields, max_rows=240),
    ]
    (od / "doctor_summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    return {
        "doctor_status": str(od / "doctor_status.json"),
        "doctor_checks": str(od / "doctor_checks.tsv"),
        "tool_status": str(od / "tool_status.tsv"),
        "reference_status": str(od / "reference_status.tsv"),
        "smoke_tests": str(od / "smoke_tests.tsv"),
        "doctor_summary": str(od / "doctor_summary.md"),
        "blocking_issues": str(od / "blocking_issues.tsv"),
        "recommended_fixes": str(fixes),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="NeoAg Doctor: controlled-execution read-only health and release-boundary check")
    ap.add_argument("--project-root", default=".")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tools-manifest")
    ap.add_argument("--reference-manifest")
    ap.add_argument("--sample-manifest")
    ap.add_argument("--profile", default="local")
    ap.add_argument("--run-demo", action="store_true")
    ap.add_argument("--run-pytest", action="store_true")
    ap.add_argument("--run-nextflow", action="store_true")
    ap.add_argument("--mini-smoke", action="store_true")
    ap.add_argument("--skip-release-audit", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)
    res = run_doctor(
        project_root=args.project_root,
        outdir=args.outdir,
        tools_manifest=args.tools_manifest,
        reference_manifest=args.reference_manifest,
        sample_manifest=args.sample_manifest,
        profile=args.profile,
        run_demo=args.run_demo,
        run_pytest=args.run_pytest,
        run_nextflow=args.run_nextflow,
        mini_smoke=args.mini_smoke,
        release_audit=not args.skip_release_audit,
        allow_execute=not args.dry_run,
    )
    print(f"NeoAg Doctor status: {res.status}")
    for k, v in res.outputs.items():
        print(f"  {k}: {v}")
    return 0 if res.status in {"READY", "PARTIAL"} else 2


if __name__ == "__main__":
    raise SystemExit(main())
