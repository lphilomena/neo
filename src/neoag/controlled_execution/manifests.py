from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .io_utils import load_limited_yaml


@dataclass
class ManifestLoadResult:
    path: str
    exists: bool
    data: dict[str, Any] = field(default_factory=dict)
    error: str = ""


def load_manifest(path: str | Path | None) -> ManifestLoadResult:
    if not path:
        return ManifestLoadResult(path="", exists=False, data={}, error="not_provided")
    p = Path(path)
    if not p.is_file():
        return ManifestLoadResult(path=str(path), exists=False, data={}, error="file_not_found")
    try:
        return ManifestLoadResult(path=str(p), exists=True, data=load_limited_yaml(p), error="")
    except Exception as exc:
        return ManifestLoadResult(path=str(p), exists=True, data={}, error=str(exc))


def manifest_paths(data: dict[str, Any]) -> list[tuple[str, str]]:
    """Flatten path-like keys from nested manifest dictionaries."""
    rows: list[tuple[str, str]] = []
    def walk(prefix: str, obj: Any):
        if isinstance(obj, dict):
            for k, v in obj.items():
                p = f"{prefix}.{k}" if prefix else str(k)
                walk(p, v)
        elif isinstance(obj, list):
            for i, v in enumerate(obj):
                walk(f"{prefix}[{i}]", v)
        else:
            if obj is not None and isinstance(obj, str):
                key_l = prefix.lower()
                if any(token in key_l for token in ["path", "file", "fasta", "gtf", "vcf", "bam", "tsv", "cache", "dir", "image", "bed", "proteome", "ligandome"]):
                    rows.append((prefix, obj))
    walk("", data)
    return rows


def infer_sample_id(sample_manifest: dict[str, Any]) -> str:
    for key in ["sample_id", "case_id", "patient_id"]:
        val = sample_manifest.get(key)
        if val:
            return str(val)
    return "UNKNOWN_SAMPLE"


def reference_entry(data: dict[str, Any], key: str) -> dict[str, Any]:
    refs = data.get("references")
    if isinstance(refs, dict) and isinstance(refs.get(key), dict):
        return dict(refs[key])
    if data.get(key):
        return {"path": data.get(key), "required": key in {"reference_fasta", "gencode_gtf", "vep_cache"}}
    return {}


def reference_path(data: dict[str, Any], key: str) -> str:
    entry = reference_entry(data, key)
    return str(entry.get("path") or "")


def reference_required(data: dict[str, Any], key: str) -> bool:
    entry = reference_entry(data, key)
    return bool(entry.get("required"))


VALID_TOOL_MODES = {"container", "apptainer", "local_license", "conda", "conda_or_container", "local_or_container", "path", "system"}


def normalize_tools_manifest(data: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Return tools keyed by tool name, accepting either top-level tools or a tools: block."""
    tools = data.get("tools") if isinstance(data.get("tools"), dict) else data
    out: dict[str, dict[str, Any]] = {}
    if not isinstance(tools, dict):
        return out
    for name, spec in tools.items():
        if not isinstance(spec, dict):
            out[str(name)] = {"executable": str(spec)}
            continue
        item = dict(spec)
        if "path" in item and "executable" not in item:
            item["executable"] = item["path"]
        out[str(name)] = item
    return out


def validate_sample_manifest(data: dict[str, Any]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    def add(field: str, status: str, message: str):
        rows.append({"manifest": "sample", "field": field, "status": status, "message": message})
    if not data.get("sample_id"):
        add("sample_id", "FAIL", "sample_id is required")
    else:
        add("sample_id", "OK", "sample_id present")
    if not isinstance(data.get("inputs", {}), dict):
        add("inputs", "FAIL", "inputs must be a mapping")
    else:
        add("inputs", "OK", "inputs mapping present")
    for section in ["tumor", "normal"]:
        if section in data and not isinstance(data.get(section), dict):
            add(section, "FAIL", f"{section} must be a mapping")
        elif section in data:
            add(section, "OK", f"{section} mapping present")
        else:
            add(section, "INFO", f"{section} section not declared")
    sv = data.get("inputs", {}).get("sv_vcf") if isinstance(data.get("inputs"), dict) else None
    if sv is not None and not isinstance(sv, list):
        add("inputs.sv_vcf", "WARN", "sv_vcf should be a list when multiple callers are supported")
    return rows


def validate_reference_manifest(data: dict[str, Any]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if data.get("manifest_version"):
        rows.append({"manifest": "reference", "field": "manifest_version", "status": "OK", "message": str(data.get("manifest_version"))})
    genome_build = data.get("genome_build")
    rows.append({"manifest": "reference", "field": "genome_build", "status": "OK" if genome_build else "FAIL", "message": "present" if genome_build else "required field missing"})
    required = ["reference_fasta", "gencode_gtf", "vep_cache"]
    for key in required:
        value = reference_path(data, key)
        rows.append({"manifest": "reference", "field": key, "status": "OK" if value else "FAIL", "message": "present" if value else "required field missing"})
    optional = ["normal_proteome", "normal_ligandome", "normal_junctions", "hla_reference", "facets_snp_vcf", "hla_la_graph", "lohhla_reference", "normal_readthrough"]
    for key in optional:
        value = reference_path(data, key)
        rows.append({"manifest": "reference", "field": key, "status": "OK" if value else "INFO", "message": "present" if value else "optional field not declared"})
    refs = data.get("references")
    if isinstance(refs, dict):
        for key, spec in refs.items():
            if isinstance(spec, dict):
                if spec.get("required") and not spec.get("path"):
                    rows.append({"manifest": "reference", "field": f"references.{key}.path", "status": "FAIL", "message": "required reference missing path"})
                for meta in ["source", "version"]:
                    if not spec.get(meta):
                        rows.append({"manifest": "reference", "field": f"references.{key}.{meta}", "status": "WARN", "message": f"{meta} not declared"})
    return rows


def validate_tools_manifest(data: dict[str, Any]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    tools = normalize_tools_manifest(data)
    if not tools:
        return [{"manifest": "tools", "field": "tools", "status": "FAIL", "message": "no tools declared"}]
    for name, spec in tools.items():
        mode = str(spec.get("mode", "system"))
        if mode not in VALID_TOOL_MODES:
            rows.append({"manifest": "tools", "field": f"{name}.mode", "status": "WARN", "message": f"unknown mode: {mode}"})
        else:
            rows.append({"manifest": "tools", "field": f"{name}.mode", "status": "OK", "message": mode})
        if not (spec.get("executable") or spec.get("image")):
            rows.append({"manifest": "tools", "field": name, "status": "WARN", "message": "neither executable nor image declared"})
        if name.lower() in {"netmhcpan", "netmhcstabpan"} and not spec.get("license_required"):
            rows.append({"manifest": "tools", "field": f"{name}.license_required", "status": "WARN", "message": "licensed tool should declare license_required: true"})
    return rows


def validate_manifests(sample: dict[str, Any], reference: dict[str, Any], tools: dict[str, Any]) -> list[dict[str, str]]:
    return validate_sample_manifest(sample) + validate_reference_manifest(reference) + validate_tools_manifest(tools)


def tool_container_digests(tools_data: dict[str, Any]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for name, spec in normalize_tools_manifest(tools_data).items():
        image = str(spec.get("image") or "")
        digest = ""
        if "@sha256:" in image:
            digest = image.split("@sha256:", 1)[1]
        rows.append({"tool": name, "mode": str(spec.get("mode", "")), "image": image, "container_digest": digest, "executable": str(spec.get("executable") or spec.get("path") or "")})
    return rows
