from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Iterable


def read_tsv(path: Path, limit: int | None = None) -> list[dict[str, str]]:
    if not path.exists():
        return []
    rows: list[dict[str, str]] = []
    with path.open('r', encoding='utf-8', errors='replace', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            rows.append({str(k): (v or '') for k, v in row.items()})
            if limit and len(rows) >= limit:
                break
    return rows


def count_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open('r', encoding='utf-8', errors='replace') as fh:
        n = sum(1 for _ in fh)
    return max(0, n - 1)


def find_latest_run(candidates: Iterable[Path]) -> Path | None:
    found: list[Path] = []
    for base in candidates:
        if not base:
            continue
        if (base / 'scoring/ranked_peptides.v03.tsv').exists():
            found.append(base)
        found.extend(p for p in base.glob('**/run-full') if (p / 'scoring/ranked_peptides.v03.tsv').exists())
    if not found:
        return None
    return max(found, key=lambda p: (p / 'scoring/ranked_peptides.v03.tsv').stat().st_mtime)


def safe_float(v: str) -> float | None:
    try:
        if v == '' or v.upper() == 'NA':
            return None
        return float(v)
    except Exception:
        return None


def fmt_candidate(row: dict[str, str], idx: int) -> str:
    bits = [
        f"{idx}. {row.get('gene','')}",
        row.get('peptide',''),
        row.get('hla_allele',''),
        row.get('event_type',''),
        f"priority={row.get('final_priority') or row.get('priority') or 'NA'}",
        f"presentation={row.get('presentation_evidence_grade') or row.get('presentation_grade') or 'NA'}",
        f"efficacy={row.get('efficacy_score') or 'NA'}",
        f"immunology={row.get('immunology_composite_score') or 'NA'}",
    ]
    return '- ' + ' | '.join(x for x in bits if x)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description='Inspect latest NeoAg result directory and write a result summary')
    ap.add_argument('--result-dir')
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--project-root', default='.')
    args = ap.parse_args(argv)

    project_root = Path(args.project_root).resolve()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    parent_outdir = outdir.parent
    candidates = []
    if args.result_dir:
        candidates.append(Path(args.result_dir))
    candidates += [parent_outdir, project_root / 'results/llm_agent_web', project_root / 'results']
    run_dir = find_latest_run(candidates)
    if not run_dir:
        msg = '# 结果分析\n\n未找到可分析的 sliding-window 结果目录。请提供 result-dir，或先运行 neoag sliding-window。\n'
        (outdir / 'result_inspection.md').write_text(msg, encoding='utf-8')
        (outdir / 'result_inspection.json').write_text(json.dumps({'status': 'NO_RESULT'}, indent=2) + '\n', encoding='utf-8')
        print(msg)
        return 0

    ranked_peptides = run_dir / 'scoring/ranked_peptides.v03.tsv'
    ranked_events = run_dir / 'scoring/ranked_events.v03.tsv'
    validation_plan = run_dir / 'scoring/validation_plan.v03.tsv'
    report = run_dir / 'reports/evidence_report.v03.html'
    patient_report = run_dir / 'reports/evidence_report.patient.html'
    provenance = run_dir / 'provenance.v03.json'

    top_peptides = read_tsv(ranked_peptides, limit=20)
    top_validation = read_tsv(validation_plan, limit=20)
    event_rows = read_tsv(ranked_events, limit=100000)
    peptide_rows_for_counts = read_tsv(ranked_peptides, limit=100000)

    priority_counts = Counter((r.get('final_priority') or r.get('priority') or 'NA') for r in peptide_rows_for_counts)
    grade_counts = Counter((r.get('presentation_evidence_grade') or r.get('presentation_grade') or 'NA') for r in peptide_rows_for_counts)
    event_counts = Counter((r.get('event_type') or 'NA') for r in event_rows)
    safety_counts = Counter((r.get('safety_status') or 'NA') for r in peptide_rows_for_counts)

    summary = {
        'status': 'PASS',
        'run_dir': str(run_dir),
        'ranked_peptides': str(ranked_peptides),
        'ranked_events': str(ranked_events),
        'validation_plan': str(validation_plan),
        'evidence_report': str(report) if report.exists() else '',
        'patient_report': str(patient_report) if patient_report.exists() else '',
        'provenance': str(provenance) if provenance.exists() else '',
        'n_ranked_peptides': count_rows(ranked_peptides),
        'n_ranked_events': count_rows(ranked_events),
        'n_validation_candidates': count_rows(validation_plan),
        'priority_counts': dict(priority_counts),
        'presentation_grade_counts': dict(grade_counts),
        'event_type_counts': dict(event_counts),
        'safety_status_counts': dict(safety_counts),
    }

    md = [
        '# NeoAg sliding-window 结果分析',
        '',
        f"结果目录：`{run_dir}`",
        '',
        '## 总体情况',
        f"- ranked peptides: {summary['n_ranked_peptides']}",
        f"- ranked events: {summary['n_ranked_events']}",
        f"- validation candidates: {summary['n_validation_candidates']}",
        f"- evidence report: `{report}`" if report.exists() else '- evidence report: 未找到',
        '',
        '## 分布概览',
        '- priority: ' + ', '.join(f'{k}={v}' for k, v in sorted(priority_counts.items())) if priority_counts else '- priority: NA',
        '- presentation grade: ' + ', '.join(f'{k}={v}' for k, v in sorted(grade_counts.items())) if grade_counts else '- presentation grade: NA',
        '- event type: ' + ', '.join(f'{k}={v}' for k, v in sorted(event_counts.items())) if event_counts else '- event type: NA',
        '- safety status: ' + ', '.join(f'{k}={v}' for k, v in sorted(safety_counts.items())) if safety_counts else '- safety status: NA',
        '',
        '## Top ranked peptides',
    ]
    if top_peptides:
        for i, row in enumerate(top_peptides[:10], 1):
            md.append(fmt_candidate(row, i))
    else:
        md.append('- 未读取到 ranked peptide 记录')

    md += ['', '## 验证建议 Top candidates']
    if top_validation:
        for i, row in enumerate(top_validation[:10], 1):
            md.append('- ' + ' | '.join(x for x in [
                f"{i}. {row.get('gene','')}",
                row.get('peptide',''),
                row.get('hla_allele',''),
                f"priority={row.get('priority','NA')}",
                row.get('recommended_assay',''),
                row.get('validation_notes',''),
            ] if x))
    else:
        md.append('- 未读取到 validation_plan 记录')

    md += [
        '',
        '## 解释要点',
        '- 这批结果已经完成 sliding-window 主流程；可优先查看 ranked_peptides、ranked_events、validation_plan 和 evidence_report。',
        '- 推荐优先级综合了呈递预测、免疫原性、APPM/免疫逃逸、安全性和 CCF/持久性等证据，不等同于单一 NetMHCpan 排名。',
        '- 当前结果属于计算筛选候选，需要后续实验验证；不要把候选直接表述为已确认新抗原。',
        '',
        '## 关键文件',
        f"- ranked_peptides: `{ranked_peptides}`",
        f"- ranked_events: `{ranked_events}`",
        f"- validation_plan: `{validation_plan}`",
    ]
    if report.exists():
        md.append(f"- evidence_report: `{report}`")
    if patient_report.exists():
        md.append(f"- patient_report: `{patient_report}`")
    if provenance.exists():
        md.append(f"- provenance: `{provenance}`")

    text = '\n'.join(md) + '\n'
    (outdir / 'result_inspection.md').write_text(text, encoding='utf-8')
    (outdir / 'result_inspection.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(text)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
