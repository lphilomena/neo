"""System-level benchmark utilities for APPM/escape/scoring behavior.

These benchmarks are release checks, not clinical validation. They exercise
synthetic perturbations, sensitivity grids, and optional ligandome/MS validation
inputs while keeping missing external data explicitly marked as pending.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .benchmark_improve import auprc, auroc
from .config import load_profile
from .scoring_v03 import apply_priority_cap_value, compute_peptide_efficacy, priority
from .utils import read_tsv, to_float, write_tsv

SYNTHETIC_FIELDS = [
    "scenario", "perturbation", "benchmark_status", "claim_scope",
    "mhc_class", "appm_multiplier", "escape_multiplier", "ccf_multiplier",
    "priority_cap", "efficacy_score", "final_priority", "expected_behavior",
    "observed_behavior", "pass_fail", "notes",
]

SENSITIVITY_FIELDS = [
    "analysis", "parameter", "value", "peptide_id", "baseline_rank", "new_rank",
    "baseline_score", "new_score", "score_delta", "baseline_priority", "new_priority",
    "top3_jaccard", "rank_changed", "interpretation",
]

LIGANDOME_FIELDS = [
    "dataset", "benchmark_status", "claim_scope", "input_path", "n", "n_presented",
    "predictor", "auroc", "auprc", "mean_score_presented", "mean_score_not_presented",
    "required_columns", "notes",
]

APPM_MS_STRATIFIED_FIELDS = [
    "dataset", "stratum_type", "stratum", "n", "n_presented", "presentation_rate",
    "mean_score", "mean_score_presented", "mean_score_not_presented", "auroc", "auprc",
    "claim_scope", "benchmark_status", "notes",
]

APPM_MULTIPLIER_DELTA_FIELDS = [
    "dataset", "metric", "without_appm", "with_appm", "delta", "n", "benchmark_status", "interpretation",
]

HLA_LIGAND_DETECTION_FIELDS = [
    "dataset", "hla_allele", "hla_lost", "mhc_class", "n_candidates", "n_ms_detected",
    "detection_rate", "appm_status", "escape_status", "interpretation",
]

SUMMARY_FIELDS = ["component", "status", "path", "n_rows", "notes"]


def _base_event() -> dict[str, str]:
    return {
        "event_id": "BENCH_EVENT", "sample_id": "BENCH_SAMPLE", "event_type": "SNV",
        "mutation_source": "synthetic", "gene": "BENCH", "event_confidence": "0.95",
        "event_expression": "30", "driver_relevance": "0.7", "tumor_vaf": "0.45",
        "clonality": "0.9", "persistence": "0.9", "tumor_specificity": "0.9",
        "safety_status": "PASS",
    }


def _base_peptide(peptide_id: str = "BENCH_PEP", mhc_class: str = "I") -> dict[str, str]:
    return {
        "peptide_id": peptide_id, "event_id": "BENCH_EVENT", "sample_id": "BENCH_SAMPLE",
        "event_type": "SNV", "mutation_source": "synthetic", "gene": "BENCH",
        "peptide": "SLYNTVATL", "hla_allele": "HLA-A*02:01" if mhc_class == "I" else "HLA-DRB1*04:01",
        "mhc_class": mhc_class, "source_tool": "synthetic", "binding_rank": "0.2",
        "el_rank": "0.2", "presentation_score": "0.9", "immunogenicity_score": "0.8",
        "wildtype_binding_rank": "99", "self_similarity_score": "0.05",
        "normal_hla_ligand_overlap": "no", "safety_status": "PASS", "safety_multiplier": "1.0",
    }


def _base_presentation(peptide_id: str = "BENCH_PEP") -> dict[str, str]:
    return {
        "peptide_id": peptide_id, "netmhcpan_ba_rank": "0.2", "netmhcpan_el_rank": "0.2",
        "mhcflurry_affinity_percentile": "0.2", "mhcflurry_presentation_score": "0.9",
        "binding_evidence_score": "0.98", "presentation_evidence_score": "0.95",
        "presentation_evidence_grade": "A", "immunogenicity_composite_score": "0.8",
        "immunogenicity_source": "synthetic",
    }


def _score_scenario(
    profile: Mapping[str, Any], *, scenario: str, perturbation: str, mhc_class: str = "I",
    appm: float = 1.0, escape: float = 1.0, ccf: float = 1.0, cap: str = "",
    expected: str = "", notes: str = "",
) -> dict[str, str]:
    event = _base_event()
    event["clonality_multiplier"] = f"{ccf:.4f}"
    pep = _base_peptide(mhc_class=mhc_class)
    pep["escape_multiplier"] = f"{escape:.4f}"
    if cap:
        pep["priority_cap"] = cap
    scored = compute_peptide_efficacy(pep, event, _base_presentation(), profile, appm=appm, ccf=ccf)
    score = to_float(scored.get("efficacy_score"), 0.0)
    raw_priority = priority("PASS", score)
    final_priority = apply_priority_cap_value(raw_priority, cap)
    observed_parts = ["score_preserved"]
    if appm < 1.0 or escape < 1.0 or ccf < 1.0:
        observed_parts = ["score_reduced"]
    if cap:
        observed_parts.append(f"priority_at_or_below_cap_{cap}")
    if final_priority != raw_priority:
        observed_parts.append(f"priority_capped_to_{final_priority}")
    observed = ";".join(observed_parts)
    return {
        "scenario": scenario, "perturbation": perturbation, "benchmark_status": "completed",
        "claim_scope": "synthetic_system_behavior_check", "mhc_class": mhc_class,
        "appm_multiplier": f"{appm:.4f}", "escape_multiplier": f"{escape:.4f}",
        "ccf_multiplier": f"{ccf:.4f}", "priority_cap": cap, "efficacy_score": f"{score:.4f}",
        "final_priority": final_priority, "expected_behavior": expected, "observed_behavior": observed,
        "pass_fail": "PASS" if not expected or expected in observed else "REVIEW", "notes": notes,
    }


def run_synthetic_perturbation(*, outdir: str | Path, profile_name: str = "default") -> Path:
    profile = load_profile(profile_name)
    rows = [
        _score_scenario(profile, scenario="baseline", perturbation="none", expected="score_preserved"),
        _score_scenario(profile, scenario="b2m_biallelic_loss", perturbation="global_mhc_i_loss", appm=0.0, escape=0.0, cap="D", expected="priority_at_or_below_cap_D", notes="MHC-I peptides should not remain high priority under global B2M/APM loss."),
        _score_scenario(profile, scenario="hla_allele_loh", perturbation="restricting_hla_lost", appm=1.0, escape=0.2, cap="C", expected="priority_at_or_below_cap_C", notes="Allele-specific loss is handled through escape/cap rather than double-counted in APPM."),
        _score_scenario(profile, scenario="jak_ifng_defect", perturbation="ifng_jak_stat_defect", appm=0.75, escape=0.7, cap="B_CAUTION", expected="priority_at_or_below_cap_B_CAUTION"),
        _score_scenario(profile, scenario="ciita_mhc_ii_defect", perturbation="mhc_ii_defect", mhc_class="II", appm=0.25, escape=0.8, cap="C", expected="priority_at_or_below_cap_C"),
        _score_scenario(profile, scenario="appm_unassessed", perturbation="missing_appm_inputs", appm=1.0, escape=1.0, expected="score_preserved", notes="Missing APPM inputs require review but should not create an artificial penalty."),
    ]
    path = Path(outdir) / "synthetic_perturbation.tsv"
    write_tsv(path, rows, SYNTHETIC_FIELDS)
    return path


def _rank_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(rows, key=lambda r: to_float(r.get("score"), 0.0), reverse=True)


def _top_jaccard(a: list[str], b: list[str]) -> float:
    sa, sb = set(a), set(b)
    if not sa and not sb:
        return 1.0
    return len(sa & sb) / len(sa | sb)


def run_sensitivity_analysis(*, outdir: str | Path, profile_name: str = "default") -> Path:
    profile = load_profile(profile_name)
    base_specs = [("P1", 1.00, 1.00, 1.00, ""), ("P2", 0.85, 1.00, 0.90, ""), ("P3", 1.00, 0.80, 0.80, ""), ("P4", 0.65, 1.00, 0.70, "B_CAUTION")]

    def scored_set(parameter: str, value: float | str | None = None) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        for pid, appm, escape, ccf, cap in base_specs:
            a, e, c, pc = appm, escape, ccf, cap
            if parameter == "appm_multiplier":
                a = float(value)
            elif parameter == "escape_multiplier":
                e = float(value)
            elif parameter == "ccf_multiplier":
                c = float(value)
            elif parameter == "priority_cap":
                pc = str(value or "")
            pep = _base_peptide(peptide_id=pid)
            pep["escape_multiplier"] = f"{e:.4f}"
            if pc:
                pep["priority_cap"] = pc
            event = _base_event()
            event["clonality_multiplier"] = f"{c:.4f}"
            score = to_float(compute_peptide_efficacy(pep, event, _base_presentation(pid), profile, appm=a, ccf=c).get("efficacy_score"), 0.0)
            rows.append({"peptide_id": pid, "score": score, "priority": apply_priority_cap_value(priority("PASS", score), pc)})
        return _rank_rows(rows)

    baseline = scored_set("baseline")
    base_rank = {r["peptide_id"]: i + 1 for i, r in enumerate(baseline)}
    base_score = {r["peptide_id"]: r["score"] for r in baseline}
    base_pri = {r["peptide_id"]: r["priority"] for r in baseline}
    base_top3 = [r["peptide_id"] for r in baseline[:3]]
    grid: list[tuple[str, float | str]] = [("appm_multiplier", 1.0), ("appm_multiplier", 0.75), ("appm_multiplier", 0.5), ("escape_multiplier", 1.0), ("escape_multiplier", 0.7), ("escape_multiplier", 0.3), ("ccf_multiplier", 1.0), ("ccf_multiplier", 0.75), ("ccf_multiplier", 0.5), ("priority_cap", ""), ("priority_cap", "B_CAUTION"), ("priority_cap", "C")]
    rows: list[dict[str, str]] = []
    for parameter, value in grid:
        new = scored_set(parameter, value)
        new_rank = {r["peptide_id"]: i + 1 for i, r in enumerate(new)}
        new_score = {r["peptide_id"]: r["score"] for r in new}
        new_pri = {r["peptide_id"]: r["priority"] for r in new}
        jacc = _top_jaccard(base_top3, [r["peptide_id"] for r in new[:3]])
        for pid in base_rank:
            rows.append({
                "analysis": "one_parameter_sensitivity", "parameter": parameter, "value": str(value),
                "peptide_id": pid, "baseline_rank": str(base_rank[pid]), "new_rank": str(new_rank[pid]),
                "baseline_score": f"{base_score[pid]:.4f}", "new_score": f"{new_score[pid]:.4f}",
                "score_delta": f"{new_score[pid] - base_score[pid]:.4f}", "baseline_priority": base_pri[pid],
                "new_priority": new_pri[pid], "top3_jaccard": f"{jacc:.4f}",
                "rank_changed": "yes" if base_rank[pid] != new_rank[pid] else "no",
                "interpretation": "rank_or_priority_sensitive" if base_rank[pid] != new_rank[pid] or base_pri[pid] != new_pri[pid] else "stable_under_this_setting",
            })
    path = Path(outdir) / "sensitivity_analysis.tsv"
    write_tsv(path, rows, SENSITIVITY_FIELDS)
    return path


def _label_from_row(row: Mapping[str, Any]) -> int | None:
    for key in ["observed_label", "presented", "ms_detected", "ligandome_detected", "label"]:
        if key not in row:
            continue
        val = str(row.get(key, "")).strip().lower()
        if val in {"1", "yes", "true", "present", "presented", "detected", "positive"}:
            return 1
        if val in {"0", "no", "false", "absent", "not_presented", "undetected", "negative"}:
            return 0
    return None


def _score_from_row(row: Mapping[str, Any]) -> float:
    for key in ["efficacy_score", "presentation_evidence_score", "presentation_score", "binding_evidence_score", "score"]:
        if str(row.get(key, "")).strip() != "":
            return to_float(row.get(key), 0.0)
    return 0.0


def run_ligandome_ms_validation(*, outdir: str | Path, ligandome_ms: str | Path | None = None) -> Path:
    out = Path(outdir) / "ligandome_ms_validation.tsv"
    required = "peptide,hla_allele,observed_label_or_presented,score_or_prediction"
    if not ligandome_ms:
        write_tsv(out, [{"dataset": "ligandome_ms_external", "benchmark_status": "external_required", "claim_scope": "pending_external_ligandome_ms_validation", "input_path": "", "n": "0", "n_presented": "0", "predictor": "presentation_or_efficacy_score", "auroc": "", "auprc": "", "mean_score_presented": "", "mean_score_not_presented": "", "required_columns": required, "notes": "Provide a peptide-HLA ligandome/MS table to compute external validation metrics."}], LIGANDOME_FIELDS)
        return out
    path = Path(ligandome_ms)
    pairs = [(label, _score_from_row(r)) for r in read_tsv(path) if (label := _label_from_row(r)) is not None]
    if not pairs:
        write_tsv(out, [{"dataset": path.stem, "benchmark_status": "pending_labels_missing", "claim_scope": "pending_external_ligandome_ms_validation", "input_path": str(path), "n": "0", "n_presented": "0", "predictor": "presentation_or_efficacy_score", "auroc": "", "auprc": "", "mean_score_presented": "", "mean_score_not_presented": "", "required_columns": required, "notes": "Input was readable but no binary ligandome/MS labels were found."}], LIGANDOME_FIELDS)
        return out
    labels = [x for x, _ in pairs]
    scores = [x for _, x in pairs]
    pos = [s for y, s in pairs if y == 1]
    neg = [s for y, s in pairs if y == 0]
    status = "completed" if sum(labels) and sum(1 for y in labels if y == 0) else "pending_two_class_labels"
    write_tsv(out, [{
        "dataset": path.stem, "benchmark_status": status,
        "claim_scope": "external_ligandome_ms_validation" if status == "completed" else "pending_external_ligandome_ms_validation",
        "input_path": str(path), "n": str(len(pairs)), "n_presented": str(sum(labels)),
        "predictor": "presentation_or_efficacy_score",
        "auroc": f"{auroc(labels, scores):.4f}" if status == "completed" else "nan",
        "auprc": f"{auprc(labels, scores):.4f}" if status == "completed" else "nan",
        "mean_score_presented": f"{(sum(pos) / len(pos)) if pos else float('nan'):.4f}",
        "mean_score_not_presented": f"{(sum(neg) / len(neg)) if neg else float('nan'):.4f}",
        "required_columns": required, "notes": "Metrics are dataset-level technical validation, not patient outcome validation.",
    }], LIGANDOME_FIELDS)
    return out


def _peptide_key(row: Mapping[str, str]) -> tuple[str, str]:
    return (str(row.get("peptide", "")).strip(), str(row.get("hla_allele", row.get("allele", ""))).strip())


def _labelled_ligand_rows(ligandome_ms: str | Path | None) -> tuple[list[dict[str, str]], str]:
    if not ligandome_ms:
        return [], "external_required"
    p = Path(ligandome_ms)
    if not p.exists():
        return [], "external_required"
    rows = read_tsv(p)
    labelled = []
    for r in rows:
        lab = _label_from_row(r)
        if lab is None:
            continue
        rr = dict(r)
        rr["_label"] = str(lab)
        rr["_score"] = str(_score_from_row(r))
        labelled.append(rr)
    if not labelled:
        return [], "pending_labels_missing"
    if not any(r["_label"] == "1" for r in labelled) or not any(r["_label"] == "0" for r in labelled):
        return labelled, "pending_two_class_labels"
    return labelled, "completed"


def _join_ranked_to_ligandome(ranked_peptides: str | Path | None, labelled: list[dict[str, str]]) -> list[dict[str, str]]:
    if not ranked_peptides or not Path(ranked_peptides).exists():
        return labelled
    ranked = read_tsv(ranked_peptides)
    by_pid = {r.get("peptide_id", ""): r for r in ranked if r.get("peptide_id")}
    by_pair = {_peptide_key(r): r for r in ranked if r.get("peptide")}
    out=[]
    for r in labelled:
        rr=dict(r)
        match = by_pid.get(r.get("peptide_id", "")) or by_pair.get(_peptide_key(r)) or {}
        rr.update({f"ranked_{k}": v for k, v in match.items()})
        if "_score" not in rr or rr.get("_score") in {"", "0.0"}:
            rr["_score"] = match.get("efficacy_score", match.get("presentation_evidence_score", rr.get("_score", "0")))
        out.append(rr)
    return out


def run_appm_ms_stratified_validation(*, outdir: str | Path, ligandome_ms: str | Path | None = None, ranked_peptides: str | Path | None = None, appm_summary: str | Path | None = None, appm_module_scores: str | Path | None = None, appm_submodule_scores: str | Path | None = None, peptide_appm_flags: str | Path | None = None, peptide_escape_flags: str | Path | None = None) -> Path:
    out = Path(outdir) / "appm_ms_stratified_validation.tsv"
    labelled, status = _labelled_ligand_rows(ligandome_ms)
    dataset = Path(ligandome_ms).stem if ligandome_ms else "ligandome_ms_external"
    if status != "completed":
        write_tsv(out, [{"dataset": dataset, "stratum_type": "APPM", "stratum": "external_required" if status == "external_required" else status, "n": "0", "n_presented": "0", "presentation_rate": "", "mean_score": "", "mean_score_presented": "", "mean_score_not_presented": "", "auroc": "", "auprc": "", "claim_scope": "pending_appm_stratified_ligandome_ms_validation", "benchmark_status": status, "notes": "Provide binary ligandome/MS labels plus ranked peptides to stratify by APPM state."}], APPM_MS_STRATIFIED_FIELDS)
        return out
    joined = _join_ranked_to_ligandome(ranked_peptides, labelled)
    appm_status = "UNASSESSED"
    if appm_summary and Path(appm_summary).exists():
        rows = read_tsv(appm_summary)
        if rows:
            appm_status = rows[0].get("mhc_i_integrity_status", rows[0].get("appm_overall_status", "UNASSESSED"))
    strata: dict[tuple[str, str], list[dict[str, str]]] = {("APPM", appm_status): joined}
    if peptide_escape_flags and Path(peptide_escape_flags).exists():
        esc = {r.get("peptide_id", ""): r for r in read_tsv(peptide_escape_flags) if r.get("peptide_id")}
        for r in joined:
            pid = r.get("peptide_id", r.get("ranked_peptide_id", ""))
            st = esc.get(pid, {}).get("escape_status", "ESCAPE_UNMAPPED")
            strata.setdefault(("escape_status", st), []).append(r)
    rows_out=[]
    for (stype, stratum), rows in strata.items():
        pairs=[(int(r["_label"]), to_float(r.get("_score", r.get("ranked_efficacy_score", "0")), 0.0)) for r in rows]
        labels=[x for x,_ in pairs]; scores=[s for _,s in pairs]
        pos=[s for y,s in pairs if y==1]; neg=[s for y,s in pairs if y==0]
        complete = bool(pos and neg)
        rows_out.append({
            "dataset": dataset, "stratum_type": stype, "stratum": stratum,
            "n": str(len(rows)), "n_presented": str(sum(labels)),
            "presentation_rate": f"{(sum(labels)/len(rows)) if rows else 0:.4f}",
            "mean_score": f"{(sum(scores)/len(scores)) if scores else 0:.4f}",
            "mean_score_presented": f"{(sum(pos)/len(pos)) if pos else float('nan'):.4f}",
            "mean_score_not_presented": f"{(sum(neg)/len(neg)) if neg else float('nan'):.4f}",
            "auroc": f"{auroc(labels, scores):.4f}" if complete else "nan",
            "auprc": f"{auprc(labels, scores):.4f}" if complete else "nan",
            "claim_scope": "appm_stratified_ligandome_ms_technical_validation",
            "benchmark_status": "completed" if complete else "pending_two_class_labels",
            "notes": "MS-positive is strong presentation evidence; MS-negative is not evidence of absence.",
        })
    write_tsv(out, rows_out, APPM_MS_STRATIFIED_FIELDS)
    return out


def run_appm_multiplier_delta(*, outdir: str | Path, ligandome_ms: str | Path | None = None, ranked_peptides: str | Path | None = None) -> Path:
    out = Path(outdir) / "appm_multiplier_delta.tsv"
    labelled, status = _labelled_ligand_rows(ligandome_ms)
    dataset = Path(ligandome_ms).stem if ligandome_ms else "ligandome_ms_external"
    if status != "completed" or not ranked_peptides or not Path(ranked_peptides).exists():
        write_tsv(out, [{"dataset": dataset, "metric": "AUPRC", "without_appm": "", "with_appm": "", "delta": "", "n": "0", "benchmark_status": status if status != "completed" else "pending_ranked_peptides", "interpretation": "Provide ligandome/MS labels and ranked_peptides to compute APPM multiplier delta."}], APPM_MULTIPLIER_DELTA_FIELDS)
        return out
    joined = _join_ranked_to_ligandome(ranked_peptides, labelled)
    labels=[]; with_scores=[]; without_scores=[]
    for r in joined:
        y=int(r["_label"]); labels.append(y)
        eff=to_float(r.get("_score", r.get("ranked_efficacy_score", "0")), 0.0)
        appm=to_float(r.get("ranked_appm_multiplier", "1.0"), 1.0) or 1.0
        with_scores.append(eff)
        without_scores.append(min(1.0, eff / appm))
    metrics=[("AUROC", auroc(labels, without_scores), auroc(labels, with_scores)), ("AUPRC", auprc(labels, without_scores), auprc(labels, with_scores))]
    rows=[{"dataset": dataset, "metric": m, "without_appm": f"{a:.4f}", "with_appm": f"{b:.4f}", "delta": f"{(b-a):.4f}", "n": str(len(labels)), "benchmark_status": "completed", "interpretation": "Positive delta means APPM multiplier improved this technical label ranking; no clinical outcome claim."} for m,a,b in metrics]
    write_tsv(out, rows, APPM_MULTIPLIER_DELTA_FIELDS)
    return out


def run_hla_ligand_detection_by_appm(*, outdir: str | Path, ligandome_ms: str | Path | None = None, ranked_peptides: str | Path | None = None, peptide_escape_flags: str | Path | None = None, appm_summary: str | Path | None = None) -> Path:
    out = Path(outdir) / "hla_ligand_detection_by_appm.tsv"
    labelled, status = _labelled_ligand_rows(ligandome_ms)
    dataset = Path(ligandome_ms).stem if ligandome_ms else "ligandome_ms_external"
    if status == "external_required":
        write_tsv(out, [{"dataset": dataset, "hla_allele": "", "hla_lost": "", "mhc_class": "", "n_candidates": "0", "n_ms_detected": "0", "detection_rate": "", "appm_status": "", "escape_status": "", "interpretation": "external_required"}], HLA_LIGAND_DETECTION_FIELDS)
        return out
    joined = _join_ranked_to_ligandome(ranked_peptides, labelled)
    esc = {r.get("peptide_id", ""): r for r in read_tsv(peptide_escape_flags)} if peptide_escape_flags and Path(peptide_escape_flags).exists() else {}
    appm_status="UNASSESSED"
    if appm_summary and Path(appm_summary).exists() and read_tsv(appm_summary):
        appm_status=read_tsv(appm_summary)[0].get("mhc_i_integrity_status", "UNASSESSED")
    groups: dict[tuple[str,str,str,str], list[dict[str,str]]] = {}
    for r in joined:
        pid=r.get("peptide_id", r.get("ranked_peptide_id", "")); e=esc.get(pid,{})
        key=(r.get("hla_allele", r.get("ranked_hla_allele", "")), e.get("restricting_hla_lost", "unmapped"), r.get("mhc_class", r.get("ranked_mhc_class", "I")), e.get("escape_status", "unmapped"))
        groups.setdefault(key,[]).append(r)
    rows=[]
    for (hla,lost,mhc,estatus), rs in groups.items():
        n=len(rs); det=sum(1 for r in rs if r.get("_label")=="1")
        rows.append({"dataset": dataset, "hla_allele": hla, "hla_lost": lost, "mhc_class": mhc, "n_candidates": str(n), "n_ms_detected": str(det), "detection_rate": f"{(det/n) if n else 0:.4f}", "appm_status": appm_status, "escape_status": estatus, "interpretation": "technical_ligand_detection_by_hla_and_escape_status"})
    write_tsv(out, rows or [{"dataset": dataset, "hla_allele": "", "hla_lost": "", "mhc_class": "", "n_candidates": "0", "n_ms_detected": "0", "detection_rate": "", "appm_status": appm_status, "escape_status": "", "interpretation": status}], HLA_LIGAND_DETECTION_FIELDS)
    return out


def _write_report(outdir: Path, summary_rows: list[dict[str, str]]) -> Path:
    report = outdir / "benchmark_system_report.md"
    lines = ["# System benchmark report", "", "This report summarizes release-level benchmark checks for APPM/escape/scoring behavior.", "These outputs do not establish clinical efficacy; they document synthetic behavior, optional ligandome/MS technical validation, and parameter sensitivity.", "", "| Component | Status | Rows | Output |", "|---|---|---:|---|"]
    for r in summary_rows:
        lines.append(f"| {r['component']} | {r['status']} | {r['n_rows']} | {r['path']} |")
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


def run_system_benchmark(*, outdir: str | Path, profile_name: str = "default", ligandome_ms: str | Path | None = None, mode: str = "all", ranked_peptides: str | Path | None = None, appm_summary: str | Path | None = None, appm_module_scores: str | Path | None = None, appm_submodule_scores: str | Path | None = None, peptide_appm_flags: str | Path | None = None, peptide_escape_flags: str | Path | None = None) -> dict[str, Any]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    if mode in {"all", "synthetic"}:
        paths["synthetic_perturbation"] = run_synthetic_perturbation(outdir=out, profile_name=profile_name)
    if mode in {"all", "sensitivity"}:
        paths["sensitivity_analysis"] = run_sensitivity_analysis(outdir=out, profile_name=profile_name)
    if mode in {"all", "ligandome-ms", "ligandome_ms"}:
        paths["ligandome_ms_validation"] = run_ligandome_ms_validation(outdir=out, ligandome_ms=ligandome_ms)
        paths["appm_ms_stratified_validation"] = run_appm_ms_stratified_validation(outdir=out, ligandome_ms=ligandome_ms, ranked_peptides=ranked_peptides, appm_summary=appm_summary, appm_module_scores=appm_module_scores, appm_submodule_scores=appm_submodule_scores, peptide_appm_flags=peptide_appm_flags, peptide_escape_flags=peptide_escape_flags)
        paths["appm_multiplier_delta"] = run_appm_multiplier_delta(outdir=out, ligandome_ms=ligandome_ms, ranked_peptides=ranked_peptides)
        paths["hla_ligand_detection_by_appm"] = run_hla_ligand_detection_by_appm(outdir=out, ligandome_ms=ligandome_ms, ranked_peptides=ranked_peptides, peptide_escape_flags=peptide_escape_flags, appm_summary=appm_summary)
    summary_rows: list[dict[str, str]] = []
    for key, path in paths.items():
        rows = read_tsv(path)
        status = rows[0].get("benchmark_status", "completed") if rows else "empty"
        if key == "sensitivity_analysis":
            status = "completed" if rows else "empty"
        summary_rows.append({"component": key, "status": status, "path": str(path), "n_rows": str(len(rows)), "notes": ""})
    summary_path = out / "benchmark_system_summary.tsv"
    write_tsv(summary_path, summary_rows, SUMMARY_FIELDS)
    report = _write_report(out, summary_rows)
    return {"summary_tsv": summary_path, "report_md": report, **paths}
