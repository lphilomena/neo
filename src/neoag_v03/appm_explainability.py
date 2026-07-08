"""APPM v0.4.2 explainability helpers.

This module adds confidence calls and MHC-I submodule scoring on top of APPM
2.0 sidecars. It is deliberately rule based: the output is mechanism evidence
for computational triage, not a clinical resistance diagnosis.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Mapping

from .utils import clamp, read_tsv, to_float, write_tsv

APPM_SUBMODULE_SCORE_FIELDS = [
    "sample_id", "parent_module", "submodule", "score", "status", "driver_defects",
    "defect_severity", "appm_call_confidence", "appm_call_confidence_score",
    "confidence_reason", "action_hint", "reason",
]

CONFIDENCE_FIELDS = [
    "appm_call_confidence", "appm_call_confidence_score", "confidence_reason",
    "critical_missing_evidence", "evidence_conflict_impact",
]

# Genes are intentionally duplicated from appm_v2 to avoid circular imports.
MHC_I_CORE = {"B2M", "HLA-A", "HLA-B", "HLA-C"}
MHC_I_PROCESSING = {"TAP1", "TAP2", "TAPBP", "ERAP1", "ERAP2", "PSMB8", "PSMB9", "PSMB10", "CALR", "CANX", "PDIA3"}
MHC_I_REGULATION = {"NLRC5"}
MHC_II_CORE = {"CIITA", "RFX5", "RFXANK", "RFXAP", "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1"}
IFNG_SIGNALING = {"IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1", "STAT2", "IRF1", "IRF9", "PTPN2", "SOCS1"}


def _conf_label(score: float) -> str:
    if score >= 0.75:
        return "high"
    if score >= 0.50:
        return "medium"
    if score >= 0.25:
        return "low"
    return "insufficient"


def _severity_from_score(score: float) -> str:
    if score <= 0.10:
        return "lethal"
    if score <= 0.35:
        return "strong"
    if score <= 0.70:
        return "moderate"
    if score < 0.95:
        return "weak"
    return "none"


def _status_from_score(prefix: str, score: float, assessed: bool = True) -> str:
    if not assessed:
        return f"{prefix}_UNASSESSED"
    if score <= 0.25:
        return f"{prefix}_DEFECTIVE"
    if score < 0.85:
        return f"{prefix}_CAUTION"
    return f"{prefix}_INTACT"


def _row_gene_map(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {r.get("gene", "").upper(): r for r in rows if r.get("gene")}


def _is_bial(row: Mapping[str, str] | None) -> bool:
    return bool(row and row.get("biallelic_status") == "BIALLELIC_LOSS")


def _is_caution(row: Mapping[str, str] | None) -> bool:
    return bool(row and row.get("functional_status") in {"caution", "defective"})


def _has_data(row: Mapping[str, str] | None) -> bool:
    if not row:
        return False
    return any(str(row.get(k, "")).strip() not in {"", "unknown", "not_assessed", "unassessed", "NO_EVIDENCE", "intact", "0.0000"} for k in [
        "mutation_consequence", "copy_number_status", "loh_status", "expression_tpm", "expression_status", "biallelic_status",
    ])


def _input_status(input_rows: list[dict[str, str]], input_type: str) -> str:
    for r in input_rows:
        if r.get("input_type") == input_type:
            return r.get("input_status", "")
    return "not_provided"


def _assessed(input_rows: list[dict[str, str]], input_type: str) -> bool:
    return _input_status(input_rows, input_type) == "assessed"


def _conflict_genes(conflict_rows: list[dict[str, str]]) -> set[str]:
    return {r.get("gene", "") for r in conflict_rows if r.get("evidence_conflict_flag") == "yes" and r.get("gene") != "NA"}


def _critical_missing(input_rows: list[dict[str, str]], module: str) -> list[str]:
    missing: list[str] = []
    if module in {"MHC-I", "MHC_I_CORE", "MHC_I_PROCESSING", "MHC_I_REGULATION", "MHC_I_HLA_LOH"}:
        if not _assessed(input_rows, "cnv_loh"):
            missing.append("MISSING_CNV")
        if not _assessed(input_rows, "rna_expression"):
            missing.append("MISSING_RNA")
        if not _assessed(input_rows, "hla_loh"):
            missing.append("HLA_LOH_NOT_ASSESSED")
    if module in {"MHC-II", "MHC_II_CORE", "MHC_II_REGULATION"}:
        if not _assessed(input_rows, "rna_expression"):
            missing.append("MISSING_RNA")
        if not _assessed(input_rows, "cnv_loh"):
            missing.append("MISSING_CNV")
    if module in {"IFNG-JAK-STAT", "IFNG_SIGNALING"}:
        if not _assessed(input_rows, "rna_expression"):
            missing.append("MISSING_RNA")
        if not _assessed(input_rows, "cnv_loh"):
            missing.append("MISSING_CNV")
    return missing


def _defect_strength(driver_defects: str, score: float, rows: Iterable[Mapping[str, str]]) -> tuple[float, list[str]]:
    reasons: list[str] = []
    drivers = (driver_defects or "").lower()
    genes = list(rows)
    bial = [r.get("gene", "") for r in genes if r.get("biallelic_status") == "BIALLELIC_LOSS"]
    if "b2m" in drivers or "B2M" in bial:
        return 1.0, ["B2M_LOSS_HIGH_CONFIDENCE"]
    if any(g in {"TAP1", "TAP2"} for g in bial):
        return 0.85, ["TAP_PROCESSING_DEFECT_WITH_CNV_SUPPORT"]
    if "hla_loh" in drivers or "hla" in drivers:
        return 0.75, ["HLA_LOH_SIGNAL"]
    if any(g in {"JAK1", "JAK2", "IFNGR1", "IFNGR2", "STAT1", "CIITA", "RFX5", "RFXANK", "RFXAP"} for g in bial):
        return 0.75, ["BIALLELIC_PATHWAY_DEFECT"]
    if "nlrc5" in drivers:
        return 0.45, ["REGULATORY_CAUTION"]
    if score < 0.75:
        return 0.55, ["MODERATE_DEFECT_SIGNAL"]
    return 0.50, ["NO_MAJOR_SIGNAL"]


def _validation_score(input_rows: list[dict[str, str]]) -> tuple[float, list[str]]:
    if _assessed(input_rows, "hla_ligandome"):
        return 1.0, ["LIGANDOME_SUPPORTED_OR_ASSESSED"]
    if _assessed(input_rows, "proteomics") or _assessed(input_rows, "phosphoproteomics"):
        return 0.75, ["PROTEOMICS_SUPPORTED_OR_ASSESSED"]
    return 0.40, ["COMPUTATIONAL_PROXY_ONLY"]


def _concordance_score(rows: Iterable[Mapping[str, str]], conflict_genes: set[str]) -> tuple[float, list[str]]:
    rows = list(rows)
    if any(r.get("gene") in conflict_genes for r in rows):
        return 0.25, ["CONFLICTING_DNA_RNA_EVIDENCE"]
    if any(r.get("biallelic_status") == "BIALLELIC_LOSS" and r.get("expression_status") in {"low", "missing_from_expression_matrix"} for r in rows):
        return 1.0, ["DNA_CNV_RNA_CONCORDANT"]
    if any(r.get("biallelic_status") == "BIALLELIC_LOSS" for r in rows):
        return 0.80, ["DNA_CNV_SUPPORTED"]
    if any(r.get("expression_status") == "low" and r.get("mutation_status") == "intact" and r.get("copy_number_status") in {"unknown", "not_assessed", ""} for r in rows):
        return 0.45, ["RNA_ONLY_LOW_EXPRESSION"]
    return 0.70, ["NO_MAJOR_CONFLICT"]


def _module_confidence(
    *,
    sample_id: str,
    module: str,
    score: float,
    status: str,
    driver_defects: str,
    gene_rows: list[dict[str, str]],
    input_rows: list[dict[str, str]],
    conflict_rows: list[dict[str, str]],
) -> dict[str, str]:
    completeness_score = 0.0
    # appm_evidence_completeness.tsv is passed as an input row surrogate in the caller; input rows alone are sufficient fallback.
    assessed_count = sum(1 for t in ["mutation_vep", "cnv_loh", "rna_expression", "hla_loh", "proteomics", "hla_ligandome"] if _assessed(input_rows, t))
    completeness_score = min(1.0, assessed_count / 4.0)
    defect_score, defect_reasons = _defect_strength(driver_defects, score, gene_rows)
    validation, validation_reasons = _validation_score(input_rows)
    concordance, concordance_reasons = _concordance_score(gene_rows, _conflict_genes(conflict_rows))
    missing = _critical_missing(input_rows, module)
    conflict_penalty = 0.20 if any("CONFLICTING" in r for r in concordance_reasons) else 0.0
    missing_penalty = min(0.25, 0.08 * len(missing))
    conf = 0.35 * completeness_score + 0.30 * defect_score + 0.20 * validation + 0.15 * concordance - conflict_penalty - missing_penalty
    # RNA-only low expression should not become a high-confidence defect claim.
    if any("RNA_ONLY_LOW_EXPRESSION" in r for r in concordance_reasons):
        conf = min(conf, 0.64)
    if not gene_rows:
        conf = min(conf, 0.20)
    conf = clamp(conf)
    reasons = defect_reasons + validation_reasons + concordance_reasons + (missing or ["NO_CRITICAL_MISSING_EVIDENCE"])
    return {
        "appm_call_confidence": _conf_label(conf),
        "appm_call_confidence_score": f"{conf:.4f}",
        "confidence_reason": ";".join(dict.fromkeys(reasons)),
        "critical_missing_evidence": ";".join(missing) if missing else "none",
        "evidence_conflict_impact": "penalized" if conflict_penalty else "none",
    }


def _submodule_score(sample_id: str, parent: str, submodule: str, genes: set[str], gene_map: Mapping[str, Mapping[str, str]], lost_hla: list[str], input_rows: list[dict[str, str]], conflict_rows: list[dict[str, str]]) -> dict[str, str]:
    rows = [dict(gene_map[g]) for g in sorted(genes) if g in gene_map]
    score = 1.0
    reasons: list[str] = []
    action = "none"
    if submodule == "MHC_I_CORE":
        if _is_bial(gene_map.get("B2M")):
            score = min(score, 0.05); reasons.append("B2M_biallelic_loss"); action = "hard_cap_mhc_i"
        hla_losses = [g for g in ["HLA-A", "HLA-B", "HLA-C"] if _is_bial(gene_map.get(g))]
        if hla_losses:
            score = min(score, 0.40 if len(hla_losses) >= 2 else 0.70); reasons.append("hla_locus_core_loss:" + ",".join(hla_losses)); action = action or "cap_affected_hla"
        if any(_is_caution(gene_map.get(g)) for g in ["HLA-A", "HLA-B", "HLA-C", "B2M"]):
            score = min(score, 0.75); reasons.append("core_gene_caution")
    elif submodule == "MHC_I_PROCESSING":
        if _is_bial(gene_map.get("TAP1")) or _is_bial(gene_map.get("TAP2")):
            score = min(score, 0.25); reasons.append("TAP1_TAP2_biallelic_defect"); action = "processing_cap"
        if any(_is_caution(gene_map.get(g)) for g in genes):
            score = min(score, 0.75); reasons.append("processing_gene_caution")
    elif submodule == "MHC_I_REGULATION":
        if _is_bial(gene_map.get("NLRC5")):
            score = min(score, 0.40); reasons.append("NLRC5_biallelic_loss"); action = "regulation_cap"
        elif _is_caution(gene_map.get("NLRC5")):
            score = min(score, 0.65); reasons.append("NLRC5_caution")
    elif submodule == "MHC_I_HLA_LOH":
        if lost_hla:
            score = min(score, 0.70 if len(lost_hla) == 1 else 0.45); reasons.append("hla_allele_loh:" + ",".join(lost_hla)); action = "peptide_level_hla_loh_policy"
        elif not _assessed(input_rows, "hla_loh"):
            reasons.append("HLA_LOH_NOT_ASSESSED"); action = "review_missing_hla_loh"
    elif submodule == "MHC_II_CORE":
        if any(_is_bial(gene_map.get(g)) for g in ["CIITA", "RFX5", "RFXANK", "RFXAP"]):
            score = min(score, 0.15); reasons.append("CIITA_RFX_biallelic_defect"); action = "mhc_ii_cap"
        elif any(_is_caution(gene_map.get(g)) for g in genes):
            score = min(score, 0.65); reasons.append("mhc_ii_core_caution")
    elif submodule == "IFNG_SIGNALING":
        if any(_is_bial(gene_map.get(g)) for g in ["JAK1", "JAK2", "IFNGR1", "IFNGR2", "STAT1"]):
            score = min(score, 0.20); reasons.append("JAK_STAT_IFNGR_biallelic_defect"); action = "ifng_response_cap"
        elif any(_is_caution(gene_map.get(g)) for g in genes):
            score = min(score, 0.65); reasons.append("ifng_signaling_caution")
    if not reasons:
        reasons.append("no_major_signal")
    status_prefix = submodule
    status = _status_from_score(status_prefix, score, assessed=bool(rows or submodule == "MHC_I_HLA_LOH"))
    conf = _module_confidence(
        sample_id=sample_id,
        module=submodule,
        score=score,
        status=status,
        driver_defects=";".join(reasons),
        gene_rows=rows,
        input_rows=input_rows,
        conflict_rows=conflict_rows,
    )
    return {
        "sample_id": sample_id,
        "parent_module": parent,
        "submodule": submodule,
        "score": f"{score:.4f}",
        "status": status,
        "driver_defects": ";".join(reasons),
        "defect_severity": _severity_from_score(score),
        **{k: conf[k] for k in ["appm_call_confidence", "appm_call_confidence_score", "confidence_reason"]},
        "action_hint": action or "none",
        "reason": ";".join(reasons),
    }


def _augment_rows_with_confidence(path: Path, key_name: str, gene_map: Mapping[str, Mapping[str, str]], input_rows: list[dict[str, str]], conflict_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = read_tsv(path) if path.exists() else []
    out: list[dict[str, str]] = []
    for row in rows:
        module = row.get(key_name, row.get("module", row.get("pathway", "")))
        score = to_float(row.get("score", row.get("pathway_score", 1.0)), 1.0)
        status = row.get("status", row.get("pathway_status", ""))
        drivers = row.get("driver_defects", "")
        if module == "MHC-I":
            genes = MHC_I_CORE | MHC_I_PROCESSING | MHC_I_REGULATION
        elif module == "MHC-II":
            genes = MHC_II_CORE
        elif module == "IFNG-JAK-STAT":
            genes = IFNG_SIGNALING
        else:
            genes = set(gene_map)
        conf = _module_confidence(
            sample_id=row.get("sample_id", ""), module=module, score=score, status=status,
            driver_defects=drivers, gene_rows=[dict(gene_map[g]) for g in genes if g in gene_map],
            input_rows=input_rows, conflict_rows=conflict_rows,
        )
        out.append({**row, **conf})
    if out:
        # Preserve original field order and append confidence fields.
        fields = list(rows[0].keys()) + [f for f in CONFIDENCE_FIELDS if f not in rows[0]]
        write_tsv(path, out, fields)
    return out


def enhance_appm_outputs_v042(outdir: str | Path, *, raw_peptides: str | Path | None = None, profile: Mapping[str, Any] | None = None) -> dict[str, str]:
    """Post-process APPM 2.0 sidecars with v0.4.2 P1 explainability.

    The function is idempotent and safe on partial outputs: missing sidecars are
    left untouched. It writes appm_submodule_scores.tsv and appends confidence
    fields to summary/module/pathway/modifier sidecars when present.
    """
    out = Path(outdir)
    gene_path = out / "appm_gene_status.tsv"
    if not gene_path.exists():
        return {}
    gene_rows = read_tsv(gene_path)
    gene_map = _row_gene_map(gene_rows)
    input_rows = read_tsv(out / "appm_input_status.tsv") if (out / "appm_input_status.tsv").exists() else []
    conflict_rows = read_tsv(out / "appm_conflicts.tsv") if (out / "appm_conflicts.tsv").exists() else []
    summary_rows = read_tsv(out / "appm_summary.tsv") if (out / "appm_summary.tsv").exists() else []
    sample_id = gene_rows[0].get("sample_id", "SAMPLE") if gene_rows else (summary_rows[0].get("sample_id", "SAMPLE") if summary_rows else "SAMPLE")
    lost_hla = []
    if summary_rows:
        lost_hla = [x for x in summary_rows[0].get("hla_loh_alleles", "").split(",") if x]

    sub_rows = [
        _submodule_score(sample_id, "MHC-I", "MHC_I_CORE", MHC_I_CORE, gene_map, lost_hla, input_rows, conflict_rows),
        _submodule_score(sample_id, "MHC-I", "MHC_I_PROCESSING", MHC_I_PROCESSING, gene_map, lost_hla, input_rows, conflict_rows),
        _submodule_score(sample_id, "MHC-I", "MHC_I_REGULATION", MHC_I_REGULATION, gene_map, lost_hla, input_rows, conflict_rows),
        _submodule_score(sample_id, "MHC-I", "MHC_I_HLA_LOH", set(), gene_map, lost_hla, input_rows, conflict_rows),
        _submodule_score(sample_id, "MHC-II", "MHC_II_CORE", MHC_II_CORE, gene_map, lost_hla, input_rows, conflict_rows),
        _submodule_score(sample_id, "IFNG-JAK-STAT", "IFNG_SIGNALING", IFNG_SIGNALING, gene_map, lost_hla, input_rows, conflict_rows),
    ]
    sub_path = out / "appm_submodule_scores.tsv"
    write_tsv(sub_path, sub_rows, APPM_SUBMODULE_SCORE_FIELDS)

    module_rows = _augment_rows_with_confidence(out / "appm_module_scores.tsv", "module", gene_map, input_rows, conflict_rows)
    pathway_rows = _augment_rows_with_confidence(out / "appm_pathway_status.tsv", "pathway", gene_map, input_rows, conflict_rows)

    if summary_rows:
        # Overall confidence is the minimum confidence among core modules, with reason aggregation.
        conf_scores = [to_float(r.get("appm_call_confidence_score"), 0.0) for r in module_rows if r.get("module") in {"MHC-I", "MHC-II", "IFNG-JAK-STAT"}]
        score = min(conf_scores) if conf_scores else 0.0
        reasons = []
        for r in module_rows:
            if r.get("confidence_reason"):
                reasons.extend(r["confidence_reason"].split(";"))
        summary_rows[0].update({
            "appm_call_confidence": _conf_label(score),
            "appm_call_confidence_score": f"{score:.4f}",
            "confidence_reason": ";".join(dict.fromkeys(reasons)) or "INSUFFICIENT_INPUTS",
            "critical_missing_evidence": ";".join(dict.fromkeys(x for r in module_rows for x in r.get("critical_missing_evidence", "").split(";") if x and x != "none")) or "none",
            "evidence_conflict_impact": "penalized" if any(r.get("evidence_conflict_impact") == "penalized" for r in module_rows) else "none",
        })
        fields = list(summary_rows[0].keys())
        write_tsv(out / "appm_summary.tsv", summary_rows, fields)

    # Add confidence context to peptide modifiers without changing their core logic.
    mod_path = out / "appm_peptide_modifiers.tsv"
    if mod_path.exists():
        mods = read_tsv(mod_path)
        mhc_i_conf = next((r for r in module_rows if r.get("module") == "MHC-I"), {})
        mhc_ii_conf = next((r for r in module_rows if r.get("module") == "MHC-II"), {})
        ifng_conf = next((r for r in module_rows if r.get("module") == "IFNG-JAK-STAT"), {})
        out_mods = []
        for m in mods:
            conf_src = mhc_ii_conf if str(m.get("mhc_class", "I")).upper() in {"II", "MHC-II", "CLASSII"} else mhc_i_conf
            m = dict(m)
            m.update({
                "appm_call_confidence": conf_src.get("appm_call_confidence", "insufficient"),
                "appm_call_confidence_score": conf_src.get("appm_call_confidence_score", "0.0000"),
                "confidence_reason": conf_src.get("confidence_reason", "INSUFFICIENT_INPUTS"),
                "ifng_call_confidence": ifng_conf.get("appm_call_confidence", "insufficient"),
            })
            out_mods.append(m)
        fields = list(mods[0].keys()) + [f for f in ["appm_call_confidence", "appm_call_confidence_score", "confidence_reason", "ifng_call_confidence"] if f not in mods[0]] if mods else []
        if out_mods:
            write_tsv(mod_path, out_mods, fields)

    return {
        "appm_submodule_scores": str(sub_path),
        "appm_module_scores": str(out / "appm_module_scores.tsv"),
        "appm_pathway_status": str(out / "appm_pathway_status.tsv"),
        "appm_summary": str(out / "appm_summary.tsv"),
    }
