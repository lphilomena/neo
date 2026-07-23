from __future__ import annotations

from pathlib import Path
from typing import Any

from .io import ensure_dir, markdown_table, read_table, row_get, safe_float, write_json, write_tsv


def run_ranking_compare(args: dict[str, Any]) -> dict[str, Any]:
    from neoag.agent_skills.ranking_compare import main as ranking_compare_main
    outdir = ensure_dir(args["outdir"])
    rc = ranking_compare_main(["--recommendation", str(args["recommendation"]), "--netmhcpan42", str(args["netmhcpan42"]), "--outdir", str(outdir)])
    status = "PASS" if rc == 0 else "FAIL"
    res = {"status": status, "skill": "neoag-ranking-compare", "outputs": {"report": str(outdir / "ranking_compare_report.md")}, "summary": "Generated ranking comparison report"}
    write_json(outdir / "skill_result.json", res)
    return res


def _validation_route(row: dict[str, str]) -> str:
    src = (row_get(row, ["source_type", "peptide_consequence", "event_type"], "") or "").lower()
    if any(x in src for x in ["fusion", "splice", "junction", "frameshift", "sv"]):
        return "long_peptide_or_minigene_plus_targeted_rna"
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
    ranked_peptides = Path(args.get("ranked_peptides") or args.get("input") or "")
    ranked_events_value = args.get("ranked_events")
    if not ranked_events_value and ranked_peptides.is_file():
        for name in ("ranked_events.evidence_consensus.tsv", "ranked_events.tsv"):
            candidate = ranked_peptides.with_name(name)
            if candidate.is_file():
                ranked_events_value = str(candidate)
                break
    if ranked_events_value:
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
    for i, row in enumerate(rows, 1):
        route = _validation_route(row)
        rec = {"rank": i, "event_evidence_rank": row_get(row, ["event_evidence_rank"], ""), "event_group_id": row_get(row, ["event_group_id"], ""), "representative_index": row_get(row, ["representative_index"], ""), "peptide_id": row_get(row, ["peptide_id"], ""), "gene": row_get(row, ["gene"], ""), "peptide": row_get(row, ["peptide"], ""), "hla_allele": row_get(row, ["hla_allele"], ""), "event_type": row_get(row, ["event_type", "source_type", "peptide_consequence"], ""), "evidence_grade": row_get(row, ["evidence_grade", "best_evidence_grade"], ""), "final_priority": row_get(row, ["final_priority"], ""), "efficacy_score": row_get(row, ["efficacy_score"], ""), "recommended_validation": route, "reason": "event-deduplicated computed triage; requires wet-lab validation"}
        candidates.append(rec)
        if route.startswith("short"):
            short.append(rec)
        else:
            longp.append({**rec, "long_peptide_design": "cover junction/novel region with 25-30 aa peptide; overlap 10-15 aa if long tail"})
            mini.append({**rec, "minigene_design": "include junction/novel sequence with 45-90 nt flanks when possible"})
            targeted.append({**rec, "targeted_rna": "confirm RNA junction/alt expression before immunoassay"})
    write_tsv(outdir / "experiment_candidates.tsv", candidates)
    write_tsv(outdir / "short_peptide_pool.tsv", short)
    write_tsv(outdir / "long_peptide_design.tsv", longp)
    write_tsv(outdir / "minigene_design.tsv", mini)
    (outdir / "targeted_rna_validation_plan.md").write_text("# Targeted RNA validation plan\n\n" + markdown_table(targeted, max_rows=top_n) + "\nBoundary: this is an experimental validation plan, not a treatment recommendation.\n", encoding="utf-8")
    res = {"status": "PASS", "skill": "neoag-experiment-design", "summary": f"Designed validation routes for top {len(candidates)} event-prioritized representatives", "input_source": input_source, "outputs": {"experiment_candidates": str(outdir / "experiment_candidates.tsv")}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_patient_report(args: dict[str, Any]) -> dict[str, Any]:
    from neoag.agent_skills.patient_report import main as patient_report_main
    outdir = ensure_dir(args["outdir"])
    argv = ["--outdir", str(outdir)]
    if args.get("ranked_peptides_or_summary"):
        argv += ["--ranked-peptides", str(args["ranked_peptides_or_summary"])]
    if args.get("evidence_report"):
        argv += ["--evidence-report", str(args["evidence_report"])]
    rc = patient_report_main(argv)
    status = "PASS" if rc == 0 else "FAIL"
    res = {"status": status, "skill": "neoag-patient-report", "summary": "Generated patient-facing draft report", "outputs": {"patient_report_md": str(outdir / "patient_report.md")}}
    write_json(outdir / "skill_result.json", res)
    return res


def run_technical_report(args: dict[str, Any]) -> dict[str, Any]:
    outdir = ensure_dir(args["outdir"])
    result_dir = Path(args.get("result_dir_or_summary") or args.get("input") or ".")
    files = sorted([p for p in result_dir.rglob("*") if p.is_file()])[:200] if result_dir.exists() and result_dir.is_dir() else []
    rows = [{"relative_path": str(p.relative_to(result_dir)), "size_bytes": p.stat().st_size} for p in files]
    write_tsv(outdir / "technical_report_files.tsv", rows)
    md = ["# NeoAg Technical Report Draft", "", "## Evidence boundary", "This report is for technical review. Candidate neoantigens are computational triage outputs and require experimental validation.", "", "## Result files", markdown_table(rows, max_rows=50)]
    (outdir / "technical_report.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    (outdir / "technical_report.html").write_text("<html><body><pre>" + "\n".join(md).replace("&", "&amp;").replace("<", "&lt;") + "</pre></body></html>", encoding="utf-8")
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
