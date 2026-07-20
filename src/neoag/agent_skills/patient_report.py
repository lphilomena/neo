from __future__ import annotations

import argparse
from pathlib import Path
from .common import count_by, ensure_dir, markdown_table, read_tsv, strip_html_text
from .appm_review import parse_evidence_html


def load_top_candidates(path: str | None, n: int = 20) -> tuple[list[dict[str, str]], dict[str, int], dict[str, int]]:
    if not path or not Path(path).exists():
        return [], {}, {}
    _, rows = read_tsv(path)
    p_counts = dict(count_by(rows, "final_priority"))
    event_counts = dict(count_by(rows, "event_type"))
    top = []
    for r in rows[:n]:
        top.append({"priority": r.get("final_priority", ""), "gene": r.get("gene", ""), "peptide": r.get("peptide", ""), "hla": r.get("hla_allele", ""), "type": r.get("event_type", ""), "use": r.get("recommended_use", "")[:120]})
    return top, p_counts, event_counts


def write_docx_if_available(path: Path, title: str, sections: list[tuple[str, str]]) -> bool:
    try:
        from docx import Document  # type: ignore
        from docx.shared import Pt  # type: ignore
    except Exception:
        return False
    doc = Document()
    doc.add_heading(title, 0)
    for heading, body in sections:
        doc.add_heading(heading, level=1)
        for para in body.split("\n"):
            if para.strip():
                doc.add_paragraph(para.strip())
    doc.save(path)
    return True


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Generate patient-facing report from Project B evidence outputs")
    ap.add_argument("--recommendation")
    ap.add_argument("--netmhcpan42")
    ap.add_argument("--evidence-report")
    ap.add_argument("--ranking-compare-report")
    ap.add_argument("--appm-review")
    ap.add_argument("--ccf-review")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--title", default="肿瘤新抗原预测分析及检测结果报告｜患者沟通版")
    args = ap.parse_args(argv)
    outdir = ensure_dir(args.outdir)
    html_summary = parse_evidence_html(args.evidence_report)
    top, priority_counts, event_counts = load_top_candidates(args.recommendation, 20)
    n_rows = sum(priority_counts.values())
    unique_summary = ""
    if args.recommendation and Path(args.recommendation).exists():
        _, rows = read_tsv(args.recommendation)
        unique_peptides = len({r.get("peptide", "") for r in rows})
        unique_events = len({r.get("event_id", "") for r in rows})
        unique_summary = f"本次候选 peptide-HLA 组合共 {len(rows):,} 条，独特肽段约 {unique_peptides:,} 条，候选事件约 {unique_events:,} 个。"

    priority_text = "；".join([f"{k}={v}" for k, v in sorted(priority_counts.items())]) if priority_counts else "未提供 final_priority 统计"
    appm_bits = []
    for k in ["MHC-I status", "MHC-II status", "IFNG response status", "appm_call_confidence", "evidence_completeness", "immune_escape_risk"]:
        if k in html_summary:
            appm_bits.append(f"{k}: {html_summary[k]}")
    appm_text = "；".join(appm_bits) if appm_bits else "未提供完整 APPM evidence report。"
    top_md = markdown_table(top, ["priority", "gene", "peptide", "hla", "type", "use"], max_rows=20)

    sections = []
    sections.append(("重要说明", "本报告是基于测序数据和计算模型形成的研究性分析，用于候选靶点筛选和后续实验设计参考。报告中的候选新抗原不等同于已经验证的新抗原，也不等同于已经确定的治疗方案。最终是否用于个体化疫苗或 T 细胞相关治疗，需要结合 RNA 表达复核、HLA LOH、APPM、ELISpot、minigene 或 long peptide 等验证，以及临床医生综合判断。"))
    sections.append(("本次综合结果", f"{unique_summary}\n综合推荐分级为：{priority_text}。A/B/C/D 级表示实验验证优先级，不是临床疗效等级。"))
    sections.append(("APPM 与免疫逃逸", f"{appm_text}\n这表示当前没有发现足以直接推翻 MHC-I 候选的强免疫逃逸证据；但如果 evidence completeness 为 PARTIAL，应理解为证据仍不完整，而不是完全证明抗原呈递系统正常。"))
    sections.append(("候选肽段分层", "B 级候选可作为优先验证对象；C 级中 RNA 支持较强或生物学意义明确者也可进入验证池；C_CAUTION 候选需要 WT peptide、normal proteome/ligandome 或安全性复核；D 级一般不优先推进，但 KRAS、TP53、EWSR1::WT1 等机制重要事件仍可人工审阅。"))
    sections.append(("Top 候选示例", top_md))
    sections.append(("后续验证建议", "SNV 错义候选可先做 mutant short peptide + WT peptide 对照的 ELISpot 或 T 细胞激活检测；frameshift、fusion、splice/exon junction 类候选更建议 targeted RNA、long peptide 或 minigene 验证；第一批实验建议选择 10–20 个候选，而不是一次性推进所有非 D 候选。"))

    md_lines = [f"# {args.title}", ""]
    for h, b in sections:
        md_lines += [f"## {h}", b, ""]
    (outdir / "patient_report.md").write_text("\n".join(md_lines), encoding="utf-8")
    html = "<html><head><meta charset='utf-8'><title>Patient Report</title></head><body>" + "\n".join(
        [f"<h1>{args.title}</h1>"] + [f"<h2>{h}</h2><pre style='white-space:pre-wrap'>{b}</pre>" for h, b in sections]
    ) + "</body></html>"
    (outdir / "patient_report.html").write_text(html, encoding="utf-8")
    wrote = write_docx_if_available(outdir / "patient_report.docx", args.title, sections)
    (outdir / "patient_report_outputs.txt").write_text("\n".join([str(outdir / "patient_report.md"), str(outdir / "patient_report.html"), str(outdir / "patient_report.docx") if wrote else "DOCX not generated: python-docx unavailable"]) + "\n", encoding="utf-8")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
