from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from neoag.open_neo.html_render import markdown_to_html

from .io import ensure_dir, markdown_table, read_table, row_get, safe_float, write_json, write_tsv


def run_ranking_compare(args: dict[str, Any]) -> dict[str, Any]:
    from neoag.agent_skills.ranking_compare import main as ranking_compare_main
    outdir = ensure_dir(args["outdir"])
    left = args.get("left") or args.get("netmhcpan42")
    right = args.get("right") or args.get("recommendation")
    left_name = args.get("left_name") or ("netmhcpan42" if args.get("netmhcpan42") else "left")
    right_name = args.get("right_name") or ("recommendation" if args.get("recommendation") else "right")
    if not left or not right:
        res = {"status": "FAIL", "skill": "neoag-ranking-compare", "failure_reason": "MISSING_INPUT:left,right"}
        write_json(outdir / "skill_result.json", res)
        return res
    rc = ranking_compare_main([
        "--left", str(left), "--left-name", str(left_name),
        "--right", str(right), "--right-name", str(right_name),
        "--outdir", str(outdir),
    ])
    status = "PASS" if rc == 0 else "FAIL"
    outputs = {
        "report": str(outdir / "ranking_compare_report.md"),
        "topn_overlap": str(outdir / "topn_overlap.tsv"),
        "candidate_rank_changes": str(outdir / "candidate_rank_changes.tsv"),
        "high_rank_hard_fail": str(outdir / "high_rank_hard_fail.tsv"),
        "top_composition": str(outdir / "top_composition.tsv"),
        "evidence_qc_summary": str(outdir / "evidence_qc_summary.tsv"),
        "manual_review_candidates": str(outdir / "manual_review_candidates.tsv"),
        "summary_json": str(outdir / "ranking_comparison_summary.json"),
    }
    res = {"status": status, "skill": "neoag-ranking-compare", "outputs": outputs, "summary": f"Compared {left_name} vs {right_name}"}
    write_json(outdir / "skill_result.json", res)
    return res


def _validation_route(row: dict[str, str]) -> str:
    priority = row_get(row, ["experiment_priority"], "").upper()
    if priority == "FUSION_CONFIRMATION_FIRST":
        return "fusion_rt_pcr_sanger_before_peptide_or_minigene"
    if priority == "TARGETED_RNA_FIRST":
        return "targeted_rna_before_peptide_assay"
    if priority == "PHASING_FIRST":
        return "read_backed_phasing_before_peptide_design"
    if priority == "SAFETY_REVIEW_FIRST":
        return "normal_reference_and_safety_review_before_assay"
    if priority == "PRESENTATION_REVIEW_FIRST":
        return "independent_presentation_prediction_before_assay"
    # Explicit biological identity outranks annotation consequence. In
    # particular, a VCF SNV/InDel annotated with splice_junction still needs
    # DNA-variant validation rather than a splice-junction workflow.
    explicit = (row_get(row, ["event_kind", "event_type", "source_event_type", "mutation_source"], "") or "").lower()
    if any(x in explicit for x in ["snv", "missense", "substitution"]):
        return "short_peptide_plus_wt_control" if row_get(row, ["wildtype_peptide", "wt_peptide"], "") else "short_peptide_elispot_with_wt_control_if_available"
    if any(x in explicit for x in ["indel", "frameshift", "insertion", "deletion"]):
        return "novel_tail_long_peptide_or_minigene" if "frame" in explicit else "short_peptide_plus_wt_control"
    if "fusion" in explicit:
        return "fusion_rt_pcr_sanger_then_long_peptide_or_minigene"
    if any(x in explicit for x in ["splice", "junction"]):
        return "targeted_rna_then_junction_long_peptide_or_minigene"
    src = (row_get(row, ["source_type", "peptide_consequence", "consequence"], "") or "").lower()
    if "fusion" in src:
        return "fusion_rt_pcr_sanger_then_long_peptide_or_minigene"
    if any(x in src for x in ["splice", "junction"]):
        return "targeted_rna_then_junction_long_peptide_or_minigene"
    if "frameshift" in src:
        return "novel_tail_long_peptide_or_minigene"
    if any(x in src for x in ["dna_sv", "structural", " bnd", " sv"]):
        return "dna_breakpoint_plus_rna_transcript_then_minigene"
    if row_get(row, ["wildtype_peptide", "wt_peptide"], ""):
        return "short_peptide_plus_wt_control"
    return "short_peptide_elispot_with_wt_control_if_available"


def _event_representatives(events: list[dict[str, str]], top_n: int) -> list[dict[str, str]]:
    candidates: list[dict[str, str]] = []
    for index in (1, 2):
        for event in events:
            prefix = f"representative_{index}_"
            peptide_id = row_get(event, [f"{prefix}peptide_id"], "")
            if not peptide_id:
                if index == 1:
                    peptide_id = row_get(event, ["best_peptide_id"], "")
                else:
                    continue
            candidate = dict(event)
            candidate.update({
                "peptide_id": peptide_id,
                "peptide": row_get(event, [f"{prefix}peptide", "best_peptide"], ""),
                "hla_allele": row_get(event, [f"{prefix}hla_allele", "best_hla_allele"], ""),
                "event_id": row_get(event, [f"{prefix}event_id", "event_id"], ""),
                "evidence_rank": row_get(event, [f"{prefix}evidence_rank", "best_peptide_evidence_rank"], ""),
                "evidence_grade": row_get(event, [f"{prefix}evidence_grade", "best_evidence_grade"], ""),
                "pareto_front": row_get(event, [f"{prefix}pareto_front", "best_pareto_front"], ""),
                "representative_index": str(index),
                "candidate_source": "ranked_events_representative",
            })
            candidates.append(candidate)
            if len(candidates) >= top_n:
                return candidates
    return candidates


def run_experiment_design(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    top_n = int(args.get("top_n") or 20)
    candidate_review = Path(args.get("candidate_review") or "")
    first_batch = Path(args.get("first_batch") or "")
    ranked_peptides = Path(args.get("ranked_peptides") or args.get("input") or "")
    ranked_events_value = args.get("ranked_events")
    if not ranked_events_value and ranked_peptides.is_file():
        for name in ("ranked_events.evidence_consensus.tsv", "ranked_events.tsv"):
            candidate = ranked_peptides.with_name(name)
            if candidate.is_file():
                ranked_events_value = str(candidate)
                break
    if first_batch.is_file():
        _, rows = read_table(first_batch)
        input_source = "first_batch_event_review"
    elif candidate_review.is_file():
        _, candidate_rows = read_table(candidate_review)
        rows = candidate_rows[:top_n]
        input_source = "candidate_review"
    elif ranked_events_value:
        ranked_events = Path(ranked_events_value)
        _, event_rows = read_table(ranked_events)
        rows = _event_representatives(event_rows, top_n)
        input_source = "ranked_events"
    else:
        _, peptide_rows = read_table(ranked_peptides)
        rows = peptide_rows[:top_n]
        input_source = "ranked_peptides_fallback"
    candidates = []
    short = []
    longp = []
    mini = []
    targeted = []
    manual = []
    for i, row in enumerate(rows, 1):
        route = _validation_route(row)
        rec = {
            "rank": i,
            "pipeline_event_rank": row_get(row, ["pipeline_event_rank", "event_evidence_rank"], ""),
            "event_group_id": row_get(row, ["event_group_id"], ""),
            "phase_group_id": row_get(row, ["phase_group_id"], ""),
            "redundancy_group": row_get(row, ["redundancy_group"], ""),
            "event_id": row_get(row, ["event_id"], ""),
            "peptide_id": row_get(row, ["peptide_id", "representative_1_peptide_id"], ""),
            "gene": row_get(row, ["gene"], ""),
            "peptide": row_get(row, ["peptide", "representative_1_peptide"], ""),
            "wildtype_peptide": row_get(row, ["wildtype_peptide", "wt_peptide"], ""),
            "hla_allele": row_get(row, ["hla_allele", "representative_1_hla", "representative_1_hla_allele"], ""),
            "event_type": row_get(row, ["event_type", "source_type", "peptide_consequence"], ""),
            "event_kind": row_get(row, ["event_kind"], ""),
            "pipeline_r_grade": row_get(row, ["pipeline_r_grade", "evidence_grade", "best_evidence_grade"], ""),
            "review_status": row_get(row, ["review_status"], ""),
            "experiment_priority": row_get(row, ["experiment_priority"], ""),
            "review_reason": row_get(row, ["review_reason"], ""),
            "recommended_validation": route,
            "reason": "event-level review triage; requires wet-lab validation",
        }
        candidates.append(rec)
        if rec["experiment_priority"] == "MANUAL_REVIEW_ONLY":
            manual.append({**rec, "manual_review_action": "mechanism review only; do not auto-promote"})
        if route.startswith("short"):
            short.append({**rec, "pairing_requirement": "mutant and wild-type peptide must be tested together"})
        if any(token in route for token in ("long_peptide", "novel_tail")):
            longp.append({**rec, "long_peptide_design": "cover junction/novel region with 25-30 aa peptide; overlap 10-15 aa if long tail"})
        if "minigene" in route:
            mini.append({**rec, "minigene_design": "include junction/novel sequence with 45-90 nt flanks when possible"})
        if any(token in route for token in ("targeted_rna", "rt_pcr", "rna_transcript")):
            targeted.append({**rec, "targeted_rna": "confirm RNA junction/alt expression before immunoassay"})
    write_tsv(outdir / "experiment_candidates.tsv", candidates)
    write_tsv(outdir / "short_peptide_pool.tsv", short)
    write_tsv(outdir / "long_peptide_design.tsv", longp)
    write_tsv(outdir / "minigene_design.tsv", mini)
    write_tsv(outdir / "targeted_rna_validation_plan.tsv", targeted)
    write_tsv(outdir / "manual_review_candidates.tsv", manual)
    (outdir / "targeted_rna_validation_plan.md").write_text("# Targeted RNA validation plan\n\n" + markdown_table(targeted, max_rows=top_n) + "\nBoundary: this is an experimental validation plan, not a treatment recommendation.\n", encoding="utf-8")
    outputs = {
        "experiment_candidates": str(outdir / "experiment_candidates.tsv"),
        "short_peptide_pool": str(outdir / "short_peptide_pool.tsv"),
        "long_peptide_design": str(outdir / "long_peptide_design.tsv"),
        "minigene_design": str(outdir / "minigene_design.tsv"),
        "targeted_rna_validation_plan": str(outdir / "targeted_rna_validation_plan.tsv"),
        "targeted_rna_validation_plan_md": str(outdir / "targeted_rna_validation_plan.md"),
        "manual_review_candidates": str(outdir / "manual_review_candidates.tsv"),
    }
    res = {"status": "PASS", "skill": "neoag-experiment-design", "summary": f"Designed validation routes for top {len(candidates)} event-prioritized representatives", "input_source": input_source, "outputs": outputs}
    write_json(outdir / "skill_result.json", res)
    return res


def run_patient_report(args: dict[str, Any]) -> dict[str, Any]:
    from neoag.agent_skills.patient_report import main as patient_report_main
    outdir = ensure_dir(args["outdir"])
    argv = ["--outdir", str(outdir)]
    recommendation = args.get("recommendation") or args.get("ranked_peptides_or_summary")
    if recommendation:
        argv += ["--recommendation", str(recommendation)]
    if args.get("evidence_report"):
        argv += ["--evidence-report", str(args["evidence_report"])]
    for key, flag in (("ranking_compare_report", "--ranking-compare-report"), ("appm_review", "--appm-review"), ("ccf_review", "--ccf-review")):
        if args.get(key):
            argv += [flag, str(args[key])]
    rc = patient_report_main(argv)
    status = "PASS" if rc == 0 else "FAIL"
    res = {"status": status, "skill": "neoag-patient-report", "summary": "Generated patient-facing draft report", "outputs": {"patient_report_md": str(outdir / "patient_report.md")}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_technical_report(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    result_dir = Path(args.get("result_dir_or_summary") or args.get("input") or ".")
    files = sorted([p for p in result_dir.rglob("*") if p.is_file()])[:200] if result_dir.exists() and result_dir.is_dir() else []
    rows = []
    for path in files:
        digest = hashlib.sha256(path.read_bytes()).hexdigest() if path.stat().st_size < 50 * 1024 * 1024 else "not_computed_large_file"
        rows.append({"relative_path": str(path.relative_to(result_dir)), "size_bytes": path.stat().st_size, "sha256": digest})
    write_tsv(outdir / "technical_report_files.tsv", rows)
    manifest_summary = "No pipeline manifest supplied."
    manifest_path = Path(args.get("pipeline_manifest") or "")
    if manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest_summary = markdown_table([{"field": key, "value": json.dumps(value, ensure_ascii=False) if isinstance(value, (dict, list)) else value} for key, value in manifest.items()], max_rows=50)
        except (OSError, json.JSONDecodeError) as exc:
            manifest_summary = f"Manifest parse failed: {exc}"
    md = ["# NeoAg Technical Report Draft", "", "## Evidence boundary", "This report is for technical review. Candidate neoantigens are computational triage outputs and require experimental validation.", "", "## Pipeline manifest / provenance", manifest_summary, "", "## Result files and hashes", markdown_table(rows, max_rows=200)]
    (outdir / "technical_report.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    (outdir / "technical_report.html").write_text(
        markdown_to_html("\n".join(md), title="NeoAg Technical Report Draft"),
        encoding="utf-8",
    )
    res = {"status": "PASS", "skill": "neoag-technical-report", "summary": f"Generated technical report over {len(rows)} files", "outputs": {"technical_report": str(outdir / "technical_report.md")}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_concept_explainer(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    concept = str(args.get("concept") or args.get("input") or "").strip().lower()
    audience = str(args.get("audience") or "patient")
    explanations = {
        "appm": "APPM（抗原加工与呈递系统）指肿瘤细胞把内部异常蛋白片段加工后通过 HLA 分子展示到细胞表面的分子系统。APPM 完整并不等同于新抗原一定有效，APPM 缺陷也需要 DNA、RNA、HLA LOH 和蛋白层证据综合判断。",
        "ccf": "CCF（Cancer Cell Fraction）是估计有多少比例癌细胞携带某个突变或事件。它不是 VAF；低纯度样本中 CCF 置信度会下降，只能作为排序辅助证据。",
        "hla loh": "HLA LOH 指肿瘤丢失某个 HLA 等位基因。若候选肽段依赖的 restricting HLA 已丢失，该候选通常应强降权或排除。未检测到 LOH 不等于绝对没有 LOH。",
        "minigene": "minigene 验证是把突变、融合或异常剪接片段构建成小型表达载体，让细胞内源性加工并呈递，适合 fusion、splice、frameshift 等不能只靠短肽证明的候选。",
        "elispot": "ELISpot 是检测 T 细胞受到候选抗原刺激后是否释放 IFN-γ 等细胞因子的实验，用于验证候选肽段是否能诱导免疫反应。",
    }
    text = explanations.get(concept, "未找到预设术语。建议在 technical report 中人工补充，避免编造概念解释。")
    md = f"# Concept explanation: {concept or 'NA'}\n\nAudience: {audience}\n\n{text}\n"
    (outdir / "concept_explanation.md").write_text(md, encoding="utf-8")
    res = {"status": "PASS", "skill": "neoag-concept-explainer", "summary": f"Generated explanation for {concept}", "outputs": {"concept_explanation": str(outdir / "concept_explanation.md")}}
    write_json(outdir / "skill_result.json", res)
    return res
