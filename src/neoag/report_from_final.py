"""Rebuild patient/technical reports from an existing production final/ directory."""

from __future__ import annotations

import json
import shutil
import statistics
from pathlib import Path
from typing import Any, Mapping

from .agent_skills.purity_cnv_review import collect_tool_results, consensus as purity_consensus
from .cohort_rules import discover_matching_cohort_contract, load_cohort_rule_contract, validate_cohort_rule_pair
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


def _purity_qc_status(tool: str, result_path: Path) -> str:
    """Return tool QC without confusing evidence provenance with QC."""
    if not result_path.exists():
        return "RESULT_PATH_MISSING"
    root = result_path if result_path.is_dir() else result_path.parent
    key = tool.upper()
    if key == "PURPLE":
        for path in sorted(root.rglob("*.purple.qc")):
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
                parts = line.split("\t", 1)
                if len(parts) == 2 and parts[0].strip().lower() == "qcstatus":
                    return parts[1].strip() or "ASSESSED"
    if key == "FACETS":
        for path in sorted(root.rglob("facets_omni2p5_summary.tsv")):
            for row in read_tsv(path):
                if str(row.get("metric") or "").strip().lower() == "status":
                    return str(row.get("value") or "ASSESSED").upper()
    return "ASSESSED"


def _canonical_purity_tool(value: Any) -> str:
    """Normalize a tool label without constraining which purity tools may exist."""
    name = str(value or "").strip()
    lowered = name.lower()
    for prefix in ("purity_", "cnv_"):
        if lowered.startswith(prefix):
            name = name[len(prefix):]
            lowered = name.lower()
    for suffix in ("_result", "_results", "_output", "_outputs"):
        if lowered.endswith(suffix):
            name = name[:-len(suffix)]
            lowered = name.lower()
            break
    aliases = {"facets": "FACETS", "purple": "PURPLE", "sequenza": "SEQUENZA", "ascat": "ASCAT"}
    return aliases.get(lowered, name.upper())


def _declared_purity_results(manifest: Mapping[str, Any]) -> dict[str, Path]:
    """Discover purity/CNV result stages from manifest metadata.

    Stage and output names are intentionally discovered instead of enumerated so
    reports remain portable when a case has a different subset of tools.
    """
    declared: dict[str, Path] = {}
    stages = manifest.get("stages") or {}
    if not isinstance(stages, Mapping):
        return declared
    for stage_name, stage in stages.items():
        if not isinstance(stage, Mapping):
            continue
        stage_key = str(stage_name).lower()
        outputs = stage.get("outputs") or {}
        if not isinstance(outputs, Mapping):
            continue
        is_purity_stage = stage_key.startswith("purity_") or stage_key.startswith("cnv_")
        for output_name, value in outputs.items():
            output_key = str(output_name).lower()
            is_tool_result = output_key.endswith(("_result", "_results", "_output", "_outputs"))
            if not value or not is_purity_stage or not is_tool_result:
                continue
            label = output_name if is_tool_result else stage_name
            tool = _canonical_purity_tool(label)
            if tool and tool not in {"CONSENSUS", "PURITY", "CNV", "TOOL_SUMMARY"}:
                declared[tool] = Path(str(value).replace("{outdir}", str((manifest.get("run") or {}).get("outdir") or "{outdir}")))
    return declared


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
    declared = _declared_purity_results(manifest)
    evidence = dict(manifest.get("evidence") or {})
    external_consensus_dirs: list[Path] = []
    for key in ("purity", "cnv"):
        value = str(evidence.get(key) or "").replace("{outdir}", str(production))
        if value:
            path = Path(value)
            external_consensus_dirs.append(path if path.is_dir() else path.parent)
    consensus_dirs = [
        *external_consensus_dirs,
        production / "evidence" / "purity_cnv",
        production / "purity" / "consensus",
        final_dir / "evidence" / "purity_cnv",
        final_dir / "purity" / "consensus",
    ]
    tool_map: dict[str, dict[str, str]] = {}
    seen_summaries: set[Path] = set()
    for consensus_dir in consensus_dirs:
        tool_summary = consensus_dir / "purity_cnv_tool_summary.tsv"
        try:
            summary_key = tool_summary.resolve()
        except OSError:
            summary_key = tool_summary
        if summary_key in seen_summaries or not tool_summary.is_file():
            continue
        seen_summaries.add(summary_key)
        summary_rows = read_tsv(tool_summary)
        for row in summary_rows:
            tool = _canonical_purity_tool(row.get("tool") or row.get("source_tool"))
            if not tool or str(row.get("status") or "").upper() == "MISSING":
                continue
            parsed = {
                "tool": tool,
                "purity": str(row.get("purity") or ""),
                "ploidy": str(row.get("ploidy") or ""),
                "status": str(row.get("status") or "ASSESSED"),
                "note": str(row.get("notes") or row.get("parse_method") or "from purity_cnv_tool_summary.tsv"),
            }
            existing = tool_map.get(tool)
            if not existing or (parsed["purity"] and not existing.get("purity")):
                tool_map[tool] = parsed
    search: list[Path] = list(declared.values())
    provenance_tools = provenance.get("tools") or {}
    if isinstance(provenance_tools, Mapping):
        for record in provenance_tools.values():
            if isinstance(record, Mapping) and record.get("file"):
                search.append(Path(str(record["file"])))
    search = [path for path in search if path.exists()]
    rows = collect_tool_results(search, sample_id=None) if search else []
    for row in rows:
        if str(row.get("status") or "").upper() == "MISSING":
            continue
        tool = _canonical_purity_tool(row.get("tool"))
        source_path = declared.get(tool) or Path(str(row.get("source_file") or ""))
        parsed = {
            "tool": tool,
            "purity": "" if row.get("purity") in {None, ""} else f"{row['purity']}",
            "ploidy": "" if row.get("ploidy") in {None, ""} else f"{row['ploidy']}",
            "status": _purity_qc_status(tool, source_path),
            "note": "已从工具原始结果解析纯度/倍性；用于多工具交叉核对。",
        }
        existing = tool_map.get(tool)
        if not existing or parsed["purity"] or not existing.get("purity"):
            tool_map[tool] = parsed
    preferred_order = ("FACETS", "PURPLE", "SEQUENZA", "ASCAT")
    tools = [tool_map[key] for key in preferred_order if key in tool_map]
    tools.extend(tool_map[key] for key in sorted(tool_map) if key not in preferred_order)
    present = set(tool_map)
    for tool, path in declared.items():
        if tool in present:
            continue
        status = "RESULT_PATH_MISSING" if not path.exists() else "NO_VALID_ESTIMATE"
        missing_record = {
            "tool": tool,
            "purity": "",
            "ploidy": "",
            "status": status,
            "note": (
                "结果路径不存在，未完成评估。"
                if status == "RESULT_PATH_MISSING"
                else "工具结果目录存在，但未形成可用的纯度/倍性估计。"
            ),
        }
        tools.append(missing_record)
    if not tools:
        return [], {}
    assessed_rows = [
        {
            "tool": row["tool"],
            "status": "FOUND",
            "purity": float(row["purity"]),
            "ploidy": float(row["ploidy"]) if row.get("ploidy") else None,
        }
        for row in tools if row.get("purity")
    ]
    cons = purity_consensus(assessed_rows)
    selected = next((row for row in tools if row.get("purity")), tools[0])
    consensus_tool_names = set((cons.get("tool_values") or {}).keys())
    consensus_ploidies = [
        float(row["ploidy"])
        for row in tools
        if row.get("ploidy") and (not consensus_tool_names or row["tool"] in consensus_tool_names)
    ]
    ploidy = f"{statistics.median(consensus_ploidies):g}" if consensus_ploidies else ""
    values = [float(row["purity"]) for row in tools if row.get("purity")]
    value_text = "、".join(f"{row['tool']}={float(row['purity']):.4f}" for row in tools if row.get("purity"))
    range_text = f"{min(values):.4f}-{max(values):.4f}" if values else "未形成"
    status = str(cons.get("status") or ("MULTI_TOOL_REVIEW" if len(values) > 1 else "SINGLE_TOOL_NO_CROSSCHECK"))
    if status == "CONCORDANT":
        basis = f"{value_text}；范围 {range_text}，多工具结果基本一致，工作值采用中位数。"
    elif status == "CONCORDANT_WITH_OUTLIER":
        excluded = str(cons.get("excluded_outliers") or "异常工具")
        basis = f"{value_text}；范围 {range_text}，{excluded} 与多数工具明显冲突。保留全部结果，工作值采用其余工具中位数并标记低置信度。"
    elif status in {"MODERATE_DISCORDANCE", "STRONG_DISCORDANCE"}:
        degree = "中等" if status == "MODERATE_DISCORDANCE" else "明显"
        basis = f"{value_text}；范围 {range_text}，工具间存在{degree}差异。工作值采用中位数并标记低置信度，需结合BAF、深度和CNV拟合图审阅。"
    else:
        basis = f"{value_text or '未获得有效纯度值'}；未形成充分的多工具数值共识。"
    consensus = {
        "recommended_purity": str(cons.get("recommended_purity") or selected.get("purity") or ""),
        "recommended_ploidy": ploidy,
        "selected_tool": "多工具中位数（排除离群值）" if status == "CONCORDANT_WITH_OUTLIER" else ("多工具中位数" if len(values) > 1 else selected.get("tool") or ""),
        "status": status,
        "basis": basis,
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
        ("reference_proteome", "reference_proteome"),
        ("gencode_gtf", "gencode_gtf"),
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
    profile = str(
        prov.get("profile")
        or (generated.get("sample") or {}).get("profile")
        or (manifest.get("run") or {}).get("profile")
        or ""
    )
    stem = Path(profile).stem if profile else ""
    if stem and not prov.get("analysis_profile"):
        prov["analysis_profile"] = stem
    parallel = prov.get("parallel_rankings") if isinstance(prov.get("parallel_rankings"), Mapping) else {}
    if parallel.get("rules_version"):
        prov.setdefault("evidence_consensus", {})
        if isinstance(prov["evidence_consensus"], dict):
            prov["evidence_consensus"].setdefault("rules_version", parallel.get("rules_version"))
            prov["evidence_consensus"].setdefault("rules_name", parallel.get("rules_name"))
    run_metadata = dict(manifest.get("run") or {})
    clinical_context = run_metadata.get("clinical_context")
    if isinstance(clinical_context, Mapping):
        # A cohort contract selects a knowledge file, but never supplies a patient diagnosis.
        prov["clinical_context"] = {**dict(prov.get("clinical_context") or {}), **dict(clinical_context)}
    clinical_context_source = str(run_metadata.get("clinical_context_source") or "").strip()
    if clinical_context_source:
        prov["clinical_context_source"] = clinical_context_source
    evidence_metadata = dict(manifest.get("evidence") or {})
    ranking_profile = str(run_metadata.get("profile") or (generated.get("sample") or {}).get("profile") or "")
    consensus_rules = str(
        evidence_metadata.get("evidence_consensus_rules")
        or (generated.get("inputs") or {}).get("evidence_consensus_rules")
        or parallel.get("rules")
        or ""
    )
    explicit_contract = str(run_metadata.get("cohort_rule_set") or "")
    contract = None
    contract_errors: list[str] = []
    try:
        if explicit_contract:
            contract = load_cohort_rule_contract(explicit_contract)
            contract_errors = validate_cohort_rule_pair(
                contract,
                ranking_profile=ranking_profile,
                evidence_consensus_rules=consensus_rules,
            )
        elif ranking_profile and consensus_rules:
            contract = discover_matching_cohort_contract(ranking_profile, consensus_rules)
    except (OSError, ValueError) as exc:
        contract_errors = [str(exc)]
    if contract:
        hash_fields = (
            "cohort_rule_set_sha256", "ranking_profile_sha256",
            "evidence_consensus_rules_sha256",
        )
        expected_hashes = {
            "cohort_rule_set_sha256": contract["contract_sha256"],
            "ranking_profile_sha256": contract["ranking_profile_sha256"],
            "evidence_consensus_rules_sha256": contract["evidence_consensus_rules_sha256"],
        }
        for key in hash_fields:
            recorded = str(run_metadata.get(key) or "")
            if recorded and recorded != expected_hashes[key]:
                contract_errors.append(f"{key} mismatch")
        comparability = (
            "MISMATCH"
            if contract_errors
            else ("LOCKED_COMPARABLE" if explicit_contract else "INFERRED_NOT_LOCKED")
        )
        prov["cohort_rule_contract"] = {
            "id": contract["id"],
            "version": contract["version"],
            "status": contract["status"],
            "disease": contract["disease"],
            "path": contract["path"],
            "sha256": contract["contract_sha256"],
            "ranking_profile": contract["ranking_profile"],
            "ranking_profile_sha256": contract["ranking_profile_sha256"],
            "evidence_consensus_rules": contract["evidence_consensus_rules"],
            "evidence_consensus_rules_sha256": contract["evidence_consensus_rules_sha256"],
            "report_contract_version": contract["report_contract_version"],
            "release_audit_policy": contract["release_audit_policy"],
            "comparability_status": comparability,
            "errors": contract_errors,
            "disease_knowledge_file": contract["disease_knowledge_file"],
        }
        if contract["disease_knowledge_file"]:
            prov.setdefault("disease_knowledge_file", contract["disease_knowledge_file"])
            prov["molecular_knowledge_selection"] = {
                "disease": contract["disease"],
                "source": "cohort_rule_contract",
                "knowledge_file": contract["disease_knowledge_file"],
            }
    else:
        prov["cohort_rule_contract"] = {
            "comparability_status": "UNASSESSED_NO_MATCHING_CONTRACT",
            "errors": contract_errors,
        }
    prov["report_role"] = "production_evidence_consensus"
    return prov


def write_reports_from_final(
    final_dir: str | Path,
    *,
    event_top_n: int = 20,
    candidate_top_n: int = 100,
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
    consensus_rule_path = Path(str((provenance.get("parallel_rankings") or {}).get("rules") or ""))
    consensus_rules = _toml_load(consensus_rule_path)
    manual_review = (
        consensus_rules.get("manual_review")
        or (provenance.get("evidence_consensus") or {}).get("manual_review")
        or {}
    )
    if manual_review:
        profile = dict(profile)
        profile["manual_review"] = manual_review
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
