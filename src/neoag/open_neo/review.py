from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any

from neoag.controlled_execution.io_utils import load_limited_yaml, markdown_table, read_tsv, write_json, write_tsv
from neoag.skill_taxonomy.review_skills import run_experiment_design, run_patient_report, run_ranking_compare, run_technical_report

from .contracts import MacroResult, MacroStep
from .execution_adapters import discover_result_artifacts
from .state import RunLayout, audit, new_run_id, safe_identifier, update_case_state


def _get(row: dict[str, str], *keys: str, default: str = "") -> str:
    for key in keys:
        value = row.get(key)
        if value not in {None, ""}:
            return str(value)
    return default


def _priority(grade: str, manual: str) -> str:
    grade = grade.upper()
    if grade == "R1":
        return "HIGH"
    if grade == "R2":
        return "MEDIUM_HIGH"
    if grade == "R3":
        return "EVIDENCE_COMPLETION_FIRST"
    if manual.lower() in {"yes", "true", "1"}:
        return "MANUAL_REVIEW_ONLY"
    return "HOLD"


def _validation(event_type: str, next_steps: str) -> str:
    text = f"{event_type} {next_steps}".lower()
    if "fusion" in text:
        return "RT-PCR/Sanger + breakpoint long peptide/minigene"
    if any(x in text for x in ["splice", "junction", "exon"]):
        return "targeted RNA + abnormal-junction long peptide/minigene"
    if any(x in text for x in ["frameshift", "novel tail"]):
        return "novel-tail long peptide/minigene"
    if any(x in text for x in ["sv", "structural"]):
        return "DNA breakpoint + RNA transcript + long peptide/minigene"
    return "MT/WT paired short peptide ELISpot/tetramer"


def _review_rows(event_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for row in event_rows:
        grade = _get(row, "best_evidence_grade", "evidence_grade", default="R4")
        event_type = _get(row, "event_type", "biological_event_track", "evidence_track")
        manual = _get(row, "manual_review_required", default="no")
        next_steps = _get(row, "recommended_next_steps")
        out.append({
            "pipeline_event_rank": _get(row, "event_evidence_rank"),
            "event_group_id": _get(row, "event_group_id", "event_id"),
            "event_id": _get(row, "event_id"),
            "gene": _get(row, "gene"),
            "event_type": event_type,
            "pipeline_r_grade": grade,
            "pareto_front": _get(row, "best_pareto_front"),
            "representative_peptide": _get(row, "representative_1_peptide", "best_peptide"),
            "representative_hla": _get(row, "representative_1_hla_allele", "best_hla_allele"),
            "representative_peptide_id": _get(row, "representative_1_peptide_id", "best_peptide_id"),
            "manual_review_required": manual,
            "review_status": _priority(grade, manual),
            "experiment_priority": _priority(grade, manual),
            "recommended_validation": _validation(event_type, next_steps),
            "pipeline_recommended_next_steps": next_steps,
            "consensus_action": _get(row, "consensus_action"),
            "review_reason": _get(row, "event_consensus_trace"),
            "member_event_count": _get(row, "member_event_count"),
            "peptide_count": _get(row, "peptide_count"),
        })
    return out


def _first_batch(rows: list[dict[str, str]], top_n: int) -> list[dict[str, str]]:
    # Deterministic heuristic, not a mathematical vaccine-set optimizer.
    priority_order = {"HIGH": 0, "MEDIUM_HIGH": 1, "EVIDENCE_COMPLETION_FIRST": 2, "MANUAL_REVIEW_ONLY": 3, "HOLD": 4}
    ranked = sorted(rows, key=lambda r: (priority_order.get(r["experiment_priority"], 9), int(r.get("pipeline_event_rank") or 10**9), r.get("event_group_id", "")))
    chosen: list[dict[str, str]] = []
    genes: set[str] = set()
    types: dict[str, int] = {}
    hlas: dict[str, int] = {}
    deferred: list[dict[str, str]] = []
    for row in ranked:
        if row["experiment_priority"] not in {"HIGH", "MEDIUM_HIGH"}:
            continue
        gene = row.get("gene", "")
        event_type = row.get("event_type", "OTHER")
        hla = row.get("representative_hla", "")
        diversity_ok = gene not in genes and types.get(event_type, 0) < max(2, top_n // 3) and (not hla or hlas.get(hla, 0) < max(2, top_n // 3))
        if diversity_ok:
            chosen.append(row)
            genes.add(gene)
            types[event_type] = types.get(event_type, 0) + 1
            if hla:
                hlas[hla] = hlas.get(hla, 0) + 1
        else:
            deferred.append(row)
        if len(chosen) >= top_n:
            break
    for row in deferred:
        if len(chosen) >= top_n:
            break
        if row.get("event_group_id") not in {x.get("event_group_id") for x in chosen}:
            chosen.append(row)
    for i, row in enumerate(chosen, 1):
        row["first_batch_rank"] = str(i)
        row["selection_method"] = "deterministic_event_dedup_diversity_heuristic_not_optimized"
    return chosen


def _read_context(path: str | Path | None) -> dict[str, Any]:
    if not path:
        return {}
    p = Path(path)
    if not p.is_file():
        return {}
    try:
        return load_limited_yaml(p)
    except Exception:
        return {}


def _context_summary(context: dict[str, Any]) -> str:
    preferred = ("disease", "diagnosis", "sample_type", "analysis_goal", "therapy_context", "notes")
    rows = []
    for key in preferred:
        value = context.get(key)
        if isinstance(value, (str, int, float, bool)) and str(value).strip():
            rows.append(f"- {key}: {value}")
    return "\n".join(rows) if rows else "No structured clinical context was supplied."


def _write_docx(path: Path, title: str, sections: list[tuple[str, str]], rows: list[dict[str, str]]) -> bool:
    try:
        from docx import Document
        from docx.shared import Pt
    except Exception:
        return False
    doc = Document()
    doc.add_heading(title, 0)
    doc.add_paragraph("研究性计算筛选结果；不构成临床诊断、治疗建议或疗效承诺。")
    for heading, body in sections:
        doc.add_heading(heading, level=1)
        doc.add_paragraph(body)
    doc.add_heading("事件级候选与实验优先级", level=1)
    cols = ["first_batch_rank", "gene", "event_type", "pipeline_r_grade", "representative_peptide", "representative_hla", "experiment_priority", "recommended_validation"]
    table = doc.add_table(rows=1, cols=len(cols))
    table.style = "Table Grid"
    for i, col in enumerate(cols):
        table.rows[0].cells[i].text = col
    for row in rows:
        cells = table.add_row().cells
        for i, col in enumerate(cols):
            cells[i].text = str(row.get(col, ""))
    styles = doc.styles
    if "Normal" in styles:
        styles["Normal"].font.size = Pt(10)
    doc.save(path)
    return True


def _write_reports(layout: RunLayout, context: dict[str, Any], review_rows: list[dict[str, str]], first_batch: list[dict[str, str]], artifacts: dict[str, str]) -> dict[str, str]:
    patient_sections = [
        ("样本与分析背景", _context_summary(context)),
        ("分析范围", "本报告基于事件级 Evidence consensus 结果，对候选事件进行去重审阅，并给出下一步实验验证优先级。"),
        ("结果边界", "R1–R4 是计算证据等级，不等同于已验证新抗原。缺失证据标记为未评估或部分证据，不能解释为阴性。"),
        ("候选选择", f"本次从事件级候选中形成 {len(first_batch)} 个第一批实验候选。选择同时考虑 R 等级、事件去重、HLA 和事件类型多样性；该集合是确定性启发式清单，不是已经验证的最优治疗组合。"),
        ("下一步", "错义突变优先采用 MT/WT 成对短肽；移码、融合和异常剪接优先采用 targeted RNA、long peptide 或 minigene，以验证异常序列能否自然加工和呈递。"),
    ]
    patient_md = ["# Open-Neo 患者沟通版审阅报告", "", "> 仅供研究使用。候选新抗原需要进一步实验验证，不能直接用于临床治疗决策。", ""]
    for heading, body in patient_sections:
        patient_md += [f"## {heading}", "", body, ""]
    patient_md += ["## 第一批实验候选", "", markdown_table(first_batch, columns=["first_batch_rank", "gene", "event_type", "pipeline_r_grade", "representative_peptide", "representative_hla", "experiment_priority", "recommended_validation"], max_rows=30)]
    patient_path = layout.reports / "patient_report.md"
    patient_path.write_text("\n".join(patient_md) + "\n", encoding="utf-8")
    patient_html = layout.reports / "patient_report.html"
    patient_html.write_text("<html><body><pre>" + html.escape("\n".join(patient_md)) + "</pre></body></html>", encoding="utf-8")

    technical_md = [
        "# Open-Neo technical review report",
        "",
        "## Evidence boundary",
        "Pipeline ranks remain unchanged. This review writes independent experiment-priority outputs and does not mutate the Skill2 result directory.",
        "",
        "## Result artifacts",
        "",
        markdown_table([{"artifact": k, "path": v} for k, v in sorted(artifacts.items())], max_rows=50),
        "",
        "## Event-level review",
        "",
        markdown_table(review_rows, max_rows=50),
    ]
    technical_path = layout.reports / "technical_report.md"
    technical_path.write_text("\n".join(technical_md) + "\n", encoding="utf-8")
    technical_html = layout.reports / "technical_report.html"
    technical_html.write_text("<html><body><pre>" + html.escape("\n".join(technical_md)) + "</pre></body></html>", encoding="utf-8")
    outputs = {
        "patient_report_md": str(patient_path),
        "patient_report_html": str(patient_html),
        "technical_report_md": str(technical_path),
        "technical_report_html": str(technical_html),
    }
    patient_docx = layout.reports / "patient_report.docx"
    if _write_docx(patient_docx, "Open-Neo 新抗原筛选报告（患者沟通版）", patient_sections, first_batch):
        outputs["patient_report_docx"] = str(patient_docx)
    return outputs


def run_review(args: dict[str, Any]) -> dict[str, Any]:
    result_dir = Path(args.get("result_dir") or args.get("input") or "").resolve()
    case_id = safe_identifier(str(args.get("case_id") or result_dir.name or "REVIEW"))
    layout = RunLayout.create(args.get("outdir") or f"work/open-neo-review/{case_id}")
    top_n = int(args.get("top_n") or 12)
    result = MacroResult("open-neo-review", case_id, new_run_id(case_id, "review"), "review")
    audit(layout, "open_neo_review.start", "START", result_dir=str(result_dir), top_n=top_n)

    if layout.root.resolve() == result_dir:
        result.blocking_issues.append("REPORT_BOUNDARY_VIOLATION")
        result.steps.append(MacroStep("00", "report-output-boundary", "BLOCKED", "Review outdir must differ from the source result directory", failure_code="REPORT_BOUNDARY_VIOLATION"))
        result.finish("BLOCKED").write(layout.skill_result)
        return result.to_dict()

    artifacts = discover_result_artifacts(result_dir)
    required = ["consensus_events", "consensus_peptides"]
    missing = [k for k in required if not artifacts.get(k)]
    result.steps.append(MacroStep("01", "result-integrity-check", "PASS" if not missing else "BLOCKED", detail="missing=" + ",".join(missing), outputs=artifacts))
    if missing:
        result.blocking_issues.append("CONSENSUS_RANKING_MISSING")
        result.finish("BLOCKED").write(layout.skill_result)
        return result.to_dict()

    _, event_rows = read_tsv(artifacts["consensus_events"])
    review_rows = _review_rows(event_rows)
    if not review_rows:
        result.blocking_issues.append("EVENT_MAPPING_FAILED")
        result.finish("BLOCKED").write(layout.skill_result)
        return result.to_dict()
    write_tsv(layout.review / "candidate_review.tsv", review_rows)
    first_batch = _first_batch([dict(x) for x in review_rows], top_n)
    write_tsv(layout.review / "first_batch_experiment_set.tsv", first_batch)
    evidence_completion = [x for x in review_rows if x["experiment_priority"] == "EVIDENCE_COMPLETION_FIRST"]
    write_tsv(layout.review / "evidence_completion_queue.tsv", evidence_completion)
    manual = [x for x in review_rows if x["manual_review_required"].lower() in {"yes", "true", "1"}]
    write_tsv(layout.review / "manual_review_candidates.tsv", manual)
    result.steps.append(MacroStep("02", "event-level-review-and-first-batch", "PASS", detail=f"events={len(review_rows)}; first_batch={len(first_batch)}; evidence_completion={len(evidence_completion)}", outputs={"candidate_review": str(layout.review / "candidate_review.tsv"), "first_batch": str(layout.review / "first_batch_experiment_set.tsv"), "evidence_completion_queue": str(layout.review / "evidence_completion_queue.tsv")}))

    # Reuse the internal experiment-design Skill as the implementation layer.
    exp = run_experiment_design({
        "outdir": str(layout.review / "experiment_design"),
        "ranked_events": artifacts["consensus_events"],
        "ranked_peptides": artifacts["consensus_peptides"],
        "top_n": top_n,
        "therapy_context": args.get("therapy_context") or "research",
    })
    result.steps.append(MacroStep("03", "experiment-design", exp.get("status", "PARTIAL"), outputs=exp.get("outputs", {})))

    comparison_outputs: dict[str, str] = {}
    if artifacts.get("weighted_baseline") and artifacts.get("consensus_peptides"):
        cmp = run_ranking_compare({
            "outdir": str(layout.review / "ranking_compare"),
            "left": artifacts["weighted_baseline"],
            "left_name": "weighted_baseline",
            "right": artifacts["consensus_peptides"],
            "right_name": "evidence_consensus",
        })
        comparison_outputs = cmp.get("outputs", {})
        result.steps.append(MacroStep("04", "weighted-vs-consensus-review", cmp.get("status", "PARTIAL"), outputs=comparison_outputs))

    context = _read_context(args.get("clinical_context"))
    report_outputs = _write_reports(layout, context, review_rows, first_batch, artifacts)
    patient = run_patient_report({
        "outdir": str(layout.reports / "production_patient"),
        "evidence_report": artifacts.get("evidence_report", ""),
    })
    technical = run_technical_report({
        "outdir": str(layout.reports / "production_technical"),
        "result_dir_or_summary": str(result_dir),
    })
    production_reports = {
        **{f"production_patient_{k}": v for k, v in patient.get("outputs", {}).items()},
        **{f"production_technical_{k}": v for k, v in technical.get("outputs", {}).items()},
    }
    result.steps.append(MacroStep("05", "reports", "PASS", outputs={**report_outputs, **production_reports}))
    result.outputs.update({
        "candidate_review": str(layout.review / "candidate_review.tsv"),
        "first_batch_experiment_set": str(layout.review / "first_batch_experiment_set.tsv"),
        "evidence_completion_queue": str(layout.review / "evidence_completion_queue.tsv"),
        "manual_review_candidates": str(layout.review / "manual_review_candidates.tsv"),
        **report_outputs,
        **production_reports,
        **{f"comparison_{k}": v for k, v in comparison_outputs.items()},
        **{f"pipeline_{k}": v for k, v in artifacts.items()},
    })
    result.warnings.append("first_batch_experiment_set uses a deterministic diversity heuristic and is not a validated vaccine-set optimizer")
    write_json(layout.run_manifest, {
        "schema_version": "open-neo-review-manifest-v1",
        "run_id": result.run_id,
        "case_id": case_id,
        "source_result_dir": str(result_dir),
        "source_artifacts": artifacts,
        "review_outputs": result.outputs,
        "top_n": top_n,
        "status": "PASS_WITH_WARNINGS",
    })
    result.outputs["review_manifest"] = str(layout.run_manifest)
    update_case_state(layout, case_id=case_id, current_intent="review", status="PASS_WITH_WARNINGS", source_result_dir=str(result_dir), outputs=result.outputs)
    audit(layout, "open_neo_review.finish", "PASS_WITH_WARNINGS", candidates=len(review_rows), first_batch=len(first_batch))
    result.finish("PASS_WITH_WARNINGS").write(layout.skill_result)
    return result.to_dict()
