"""Rebuild patient/technical reports from an existing production final/ directory."""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any, Mapping

from .agent_skills.purity_cnv_review import collect_tool_results, consensus as purity_consensus
from .config import load_profile
from .reports_dual import load_report_bundle, make_patient_report, make_technical_report
from .utils import read_tsv, write_json


def _copy_file(source: str | Path | None, dest: Path) -> Path | None:
    if not source:
        return None
    src = Path(source)
    if not src.is_file():
        return None
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.resolve() == src.resolve():
        return dest
    shutil.copy2(src, dest)
    return dest


def _toml_load(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        import tomllib
    except ImportError:  # pragma: no cover
        return {}
    return tomllib.loads(path.read_text(encoding="utf-8"))


def _stage_output(manifest: Mapping[str, Any], stage: str, key: str) -> str:
    return str((((manifest.get("stages") or {}).get(stage) or {}).get("outputs") or {}).get(key) or "")


def materialize_hla_loh_layout(final_dir: Path, *, hla_loh: str | Path | None = None, manifest: Mapping[str, Any] | None = None) -> dict[str, str]:
    """Copy per-tool HLA LOH tables into the layout the patient report scans."""
    production = final_dir.parent
    manifest = dict(manifest or {})
    lohhla = (
        _stage_output(manifest, "hla_loh_lohhla", "lohhla_hla_loh")
        or str(production / "evidence/hla_loh.tsv")
        or ""
    )
    spechla = (
        _stage_output(manifest, "hla_loh_spechla", "spechla_hla_loh")
        or str(production / "evidence/hla_loh.spechla.tsv")
        or ""
    )
    consensus = (
        str(hla_loh or "")
        or _stage_output(manifest, "hla_loh_consensus", "hla_loh_consensus")
        or str(production / "evidence/hla_loh/hla_loh_consensus.tsv")
    )
    copied = {
        "lohhla": str(_copy_file(lohhla, final_dir / "hla_loh/lohhla/hla_loh.tsv") or ""),
        "spechla": str(_copy_file(spechla, final_dir / "hla_loh/spechla/hla_loh.tsv") or ""),
        "consensus": str(_copy_file(consensus, final_dir / "hla_loh_consensus/hla_loh_consensus.tsv") or ""),
    }
    return {key: value for key, value in copied.items() if value}


def _purity_records(final_dir: Path, manifest: Mapping[str, Any], provenance: Mapping[str, Any]) -> tuple[list[dict[str, str]], dict[str, Any]]:
    production = final_dir.parent
    search: list[Path] = []
    for stage, key in (("purity_facets", "facets_result"), ("purity_ascat", "ascat_result"), ("purity_sequenza", "sequenza_result"), ("purity_purple", "purple_result")):
        value = _stage_output(manifest, stage, key)
        if value:
            search.append(Path(value))
    for name in ("facets", "ascat", "sequenza", "purple"):
        record = ((provenance.get("tools") or {}).get(name) or {}) if isinstance(provenance.get("tools"), Mapping) else {}
        if record.get("file"):
            search.append(Path(str(record["file"])))
    search = [path for path in search if path.exists()]
    rows = collect_tool_results(search, sample_id=None) if search else []
    tools = []
    for row in rows:
        if str(row.get("status") or "").upper() == "MISSING":
            continue
        tools.append({
            "tool": str(row.get("tool") or ""),
            "purity": "" if row.get("purity") in {None, ""} else f"{row['purity']}",
            "ploidy": "" if row.get("ploidy") in {None, ""} else f"{row['ploidy']}",
            "status": str(row.get("status") or "ASSESSED"),
            "note": str(row.get("notes") or row.get("parse_method") or ""),
        })
    if not tools:
        return [], {}
    assessed_rows = [row for row in rows if str(row.get("status") or "").upper() != "MISSING"]
    cons = purity_consensus(assessed_rows)
    selected = next((row for row in tools if row.get("purity")), tools[0])
    ploidy = next((row.get("ploidy") or "" for row in tools if row.get("ploidy")), selected.get("ploidy") or "")
    consensus = {
        "recommended_purity": str(cons.get("recommended_purity") or selected.get("purity") or ""),
        "recommended_ploidy": ploidy,
        "selected_tool": "多工具并列" if len(tools) > 1 else selected.get("tool") or "",
        "status": str(cons.get("status") or ("MULTI_TOOL_REVIEW" if len(tools) > 1 else "SINGLE_TOOL_NO_CROSSCHECK")),
        "basis": str(cons.get("interpretation") or ("已并列保留 " + "、".join(row["tool"] for row in tools) + " 结果。")),
    }
    return tools, consensus


def _input_files(final_dir: Path, manifest: Mapping[str, Any], generated: Mapping[str, Any]) -> dict[str, str]:
    production = final_dir.parent
    inputs = dict((generated.get("inputs") or {}) if isinstance(generated.get("inputs"), Mapping) else {})
    files: dict[str, str] = {}
    somatic = _stage_output(manifest, "snv_indel_candidates", "raw_events")
    command = str(((manifest.get("stages") or {}).get("snv_indel_candidates") or {}).get("command") or "")
    marker = "--input "
    if marker in command:
        token = command.split(marker, 1)[1].split(" --", 1)[0].strip().strip('"')
        if token:
            files["somatic_vcf"] = token
    for key, dest in (
        ("expression", "gene_expression"),
        ("transcript_expression", "transcript_expression"),
        ("purity", "purity"),
        ("cnv", "cnv_segments"),
        ("hla_loh", "hla_loh_consensus"),
        ("raw_events", "combined_raw_events"),
        ("raw_peptides", "combined_raw_peptides"),
    ):
        value = str(inputs.get(key) or "")
        if value:
            files[dest] = value
    hla_file = str((manifest.get("run") or {}).get("hla_file") or production / "evidence/hla_typing/recommended_hla.txt")
    if Path(hla_file.replace("{outdir}", str(production))).is_file() or "{outdir}" in hla_file:
        files["hla_consensus"] = hla_file.replace("{outdir}", str(production))
    return {key: value for key, value in files.items() if value}


def enrich_report_provenance(
    final_dir: Path,
    provenance: dict[str, Any],
    *,
    manifest: Mapping[str, Any] | None = None,
    generated: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    production = final_dir.parent
    manifest = dict(manifest or _toml_load(production / "manifest/production.results.toml"))
    generated = dict(generated or _toml_load(production / "run.production.generated.toml"))
    prov = dict(provenance)
    tools, consensus = _purity_records(final_dir, manifest, prov)
    if tools:
        prov["purity_cnv_tools"] = tools
    if consensus:
        prov["purity_cnv_consensus"] = consensus
    inputs = _input_files(final_dir, manifest, generated)
    if inputs:
        prov["input_files"] = {**dict(prov.get("input_files") or {}), **inputs}
        has_tumor = any("tumor" in str(value).lower() or key in {"somatic_vcf", "purity", "cnv"} for key, value in inputs.items())
        has_normal = any("normal" in str(value).lower() or "blood" in str(value).lower() or key in {"hla_consensus", "hla_loh_consensus"} for key, value in inputs.items())
        if has_tumor and has_normal and not prov.get("pairing_status"):
            prov["pairing_status"] = "已使用肿瘤和配对正常样本进行HLA分型、HLA LOH和纯度分析；指纹未评估"
    if not prov.get("disease"):
        profile = str(prov.get("profile") or (generated.get("sample") or {}).get("profile") or "")
        stem = Path(profile).stem if profile else ""
        if stem:
            prov["disease"] = stem
    parallel = prov.get("parallel_rankings") if isinstance(prov.get("parallel_rankings"), Mapping) else {}
    if parallel.get("rules_version"):
        prov.setdefault("evidence_consensus", {})
        if isinstance(prov["evidence_consensus"], dict):
            prov["evidence_consensus"].setdefault("rules_version", parallel.get("rules_version"))
            prov["evidence_consensus"].setdefault("rules_name", parallel.get("rules_name"))
    prov["report_role"] = "production_evidence_consensus"
    return prov


def write_reports_from_final(
    final_dir: str | Path,
    *,
    event_top_n: int = 10,
    candidate_top_n: int = 50,
) -> dict[str, str]:
    root = Path(final_dir)
    production = root.parent
    scoring = root / "scoring"
    reports = root / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    manifest = _toml_load(production / "manifest/production.results.toml")
    generated = _toml_load(production / "run.production.generated.toml")
    provenance_path = root / "provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8")) if provenance_path.is_file() else {}
    materialize_hla_loh_layout(root, hla_loh=str((generated.get("inputs") or {}).get("hla_loh") or ""), manifest=manifest)
    provenance = enrich_report_provenance(root, provenance, manifest=manifest, generated=generated)
    events_path = scoring / "ranked_events.evidence_consensus.tsv"
    peptides_path = scoring / "ranked_peptides.evidence_consensus.tsv"
    if not events_path.is_file() or not peptides_path.is_file():
        raise FileNotFoundError("consensus ranking tables are required to rebuild the patient report")
    profile_name = str(
        (generated.get("sample") or {}).get("profile")
        or provenance.get("profile")
        or "sarcoma_rna_supported_v2_provisional"
    )
    try:
        profile = load_profile(profile_name)
    except Exception:
        profile = {"_profile_name": Path(profile_name).stem, "rules_version": str((provenance.get("parallel_rankings") or {}).get("rules_version") or "")}
    appm_rows = read_tsv(root / "appm/appm_summary.tsv") if (root / "appm/appm_summary.tsv").is_file() else []
    validation_rows = read_tsv(scoring / "validation_plan.tsv") if (scoring / "validation_plan.tsv").is_file() else []
    bundle = load_report_bundle(
        profile=profile,
        events=read_tsv(events_path),
        peptides=read_tsv(peptides_path),
        appm_summary=appm_rows[0] if appm_rows else {},
        validation_rows=validation_rows,
        outdir=root,
        provenance=provenance,
        sample_id=str(provenance.get("sample_id") or (generated.get("sample") or {}).get("id") or ""),
        entry_mode=str(provenance.get("entry_mode") or (generated.get("inputs") or {}).get("entry_mode") or ""),
    )
    patient_path = reports / "evidence_report.patient.html"
    technical_path = reports / "evidence_report.technical.html"
    make_patient_report(patient_path, bundle, event_top_n=event_top_n, candidate_top_n=candidate_top_n)
    make_technical_report(technical_path, bundle)
    legacy_path = reports / "evidence_report.html"
    legacy_path.write_text(technical_path.read_text(encoding="utf-8"), encoding="utf-8")
    write_json(provenance_path, provenance)
    return {
        "evidence_report_patient": str(patient_path),
        "evidence_report_technical": str(technical_path),
        "evidence_report": str(legacy_path),
        "provenance": str(provenance_path),
    }
