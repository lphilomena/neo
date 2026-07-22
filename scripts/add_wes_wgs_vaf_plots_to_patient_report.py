#!/usr/bin/env python3
"""Add WES and WGS protein-altering SNV/InDel VAF plots to a patient report."""

from __future__ import annotations

import argparse
import csv
import html
import math
import shutil
import statistics
from pathlib import Path


START = "<!-- WES_WGS_VAF_PLOTS_START -->"
END = "<!-- WES_WGS_VAF_PLOTS_END -->"
BINS = (0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, math.inf)
LABELS = ("<1%", "1-2%", "2-5%", "5-10%", "10-20%", "20-50%", ">=50%")
COLORS = {"SNV": "#2878b5", "InDel": "#e07a2d"}


def read_variants(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            try:
                af = float(row["tumor_af"])
            except (KeyError, TypeError, ValueError):
                continue
            if math.isfinite(af) and 0 <= af <= 1:
                rows.append({"af": af, "variant_type": row.get("variant_type", "Unknown")})
    return rows


def quantile(values: list[float], q: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    pos = (len(ordered) - 1) * q
    low = math.floor(pos)
    high = math.ceil(pos)
    if low == high:
        return ordered[low]
    return ordered[low] * (high - pos) + ordered[high] * (pos - low)


def summarize(rows: list[dict[str, object]]) -> dict[str, object]:
    values = [float(row["af"]) for row in rows]
    return {
        "n": len(values),
        "snv": sum(row["variant_type"] == "SNV" for row in rows),
        "indel": sum(row["variant_type"] == "InDel" for row in rows),
        "median": statistics.median(values) if values else math.nan,
        "q1": quantile(values, 0.25),
        "q3": quantile(values, 0.75),
        "min": min(values) if values else math.nan,
        "max": max(values) if values else math.nan,
    }


def bin_counts(rows: list[dict[str, object]]) -> dict[str, list[int]]:
    counts = {"SNV": [0] * len(LABELS), "InDel": [0] * len(LABELS)}
    for row in rows:
        kind = str(row["variant_type"])
        if kind not in counts:
            continue
        af = float(row["af"])
        for idx, (lower, upper) in enumerate(zip(BINS[:-1], BINS[1:])):
            if lower <= af < upper:
                counts[kind][idx] += 1
                break
    return counts


def write_svg(path: Path, label: str, rows: list[dict[str, object]]) -> dict[str, object]:
    stats = summarize(rows)
    counts = bin_counts(rows)
    width, height = 760, 430
    left, right, top, bottom = 70, 25, 78, 72
    plot_w, plot_h = width - left - right, height - top - bottom
    totals = [counts["SNV"][i] + counts["InDel"][i] for i in range(len(LABELS))]
    ymax = max(totals + [1])
    ymax = max(5, int(math.ceil(ymax / 5.0) * 5))
    bar_w = plot_w / len(LABELS) * 0.62
    gap = plot_w / len(LABELS)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" style="width:100%;height:auto;display:block">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="{left}" y="28" font-family="Arial,sans-serif" font-size="20" font-weight="700" fill="#17324d">{html.escape(label)} 蛋白改变型 SNV/InDel VAF 分布</text>',
        f'<text x="{left}" y="52" font-family="Arial,sans-serif" font-size="13" fill="#4b5563">n={stats["n"]}；SNV={stats["snv"]}；InDel={stats["indel"]}；中位 VAF={stats["median"]:.1%}；IQR={stats["q1"]:.1%}-{stats["q3"]:.1%}</text>',
    ]
    for tick in range(0, 6):
        value = ymax * tick / 5
        y = top + plot_h - plot_h * tick / 5
        parts.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left + plot_w}" y2="{y:.1f}" stroke="#e5e7eb"/>')
        parts.append(f'<text x="{left - 10}" y="{y + 4:.1f}" text-anchor="end" font-family="Arial,sans-serif" font-size="11" fill="#6b7280">{value:.0f}</text>')
    for idx, bin_label in enumerate(LABELS):
        x = left + gap * (idx + 0.5) - bar_w / 2
        current_y = top + plot_h
        for kind in ("SNV", "InDel"):
            count = counts[kind][idx]
            h = plot_h * count / ymax
            current_y -= h
            parts.append(f'<rect x="{x:.1f}" y="{current_y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{COLORS[kind]}"/>')
        if totals[idx]:
            parts.append(f'<text x="{x + bar_w / 2:.1f}" y="{current_y - 6:.1f}" text-anchor="middle" font-family="Arial,sans-serif" font-size="11" fill="#374151">{totals[idx]}</text>')
        parts.append(f'<text x="{x + bar_w / 2:.1f}" y="{top + plot_h + 22}" text-anchor="middle" font-family="Arial,sans-serif" font-size="11" fill="#374151">{html.escape(bin_label)}</text>')
    parts.extend([
        f'<text x="{left + plot_w / 2}" y="{height - 14}" text-anchor="middle" font-family="Arial,sans-serif" font-size="12" fill="#374151">肿瘤变异等位基因频率（VAF）</text>',
        f'<text x="18" y="{top + plot_h / 2}" transform="rotate(-90 18 {top + plot_h / 2})" text-anchor="middle" font-family="Arial,sans-serif" font-size="12" fill="#374151">变异数</text>',
        f'<rect x="{width - 165}" y="20" width="13" height="13" fill="{COLORS["SNV"]}"/><text x="{width - 145}" y="31" font-family="Arial,sans-serif" font-size="12">SNV</text>',
        f'<rect x="{width - 95}" y="20" width="13" height="13" fill="{COLORS["InDel"]}"/><text x="{width - 75}" y="31" font-family="Arial,sans-serif" font-size="12">InDel</text>',
        '</svg>',
    ])
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")
    return stats


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--wes-tsv", type=Path, required=True)
    parser.add_argument("--wgs-tsv", type=Path, required=True)
    parser.add_argument("--patient-report", type=Path, required=True)
    parser.add_argument("--assets-dir", type=Path, required=True)
    parser.add_argument("--backup", type=Path)
    args = parser.parse_args()

    args.assets_dir.mkdir(parents=True, exist_ok=True)
    wes_rows = read_variants(args.wes_tsv)
    wgs_rows = read_variants(args.wgs_tsv)
    if not wes_rows or not wgs_rows:
        raise SystemExit("Both WES and WGS tables must contain valid tumor_af values")

    wes_svg = args.assets_dir / "wes_protein_altering_snv_indel_vaf.svg"
    wgs_svg = args.assets_dir / "wgs_protein_altering_snv_indel_vaf.svg"
    wes_stats = write_svg(wes_svg, "今年 WES", wes_rows)
    wgs_stats = write_svg(wgs_svg, "今年 WGS", wgs_rows)

    summary_path = args.assets_dir / "vaf_distribution_summary.tsv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["dataset", "n", "snv", "indel", "median_vaf", "q1_vaf", "q3_vaf", "min_vaf", "max_vaf"])
        for dataset, stats in (("current_WES", wes_stats), ("current_WGS", wgs_stats)):
            writer.writerow([dataset, stats["n"], stats["snv"], stats["indel"], *(f'{stats[k]:.6g}' for k in ("median", "q1", "q3", "min", "max"))])

    report = args.patient_report.read_text(encoding="utf-8")
    if args.backup:
        args.backup.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(args.patient_report, args.backup)
    wes_svg_inline = wes_svg.read_text(encoding="utf-8").strip()
    wgs_svg_inline = wgs_svg.read_text(encoding="utf-8").strip()
    block = f"""{START}
<h3>今年 WES 与今年 WGS 的 VAF 分布</h3>
<div style='display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:14px;margin:12px 0 18px'>
<figure style='margin:0;border:1px solid #ddd'>{wes_svg_inline}<figcaption class='small' style='padding:6px'>今年 WES：n={wes_stats['n']}，中位 VAF {wes_stats['median']:.1%}（IQR {wes_stats['q1']:.1%}-{wes_stats['q3']:.1%}）。</figcaption></figure>
<figure style='margin:0;border:1px solid #ddd'>{wgs_svg_inline}<figcaption class='small' style='padding:6px'>今年 WGS：n={wgs_stats['n']}，中位 VAF {wgs_stats['median']:.1%}（IQR {wgs_stats['q1']:.1%}-{wgs_stats['q3']:.1%}）。</figcaption></figure>
</div>
<p class='small'>图中仅包含本节严格口径的蛋白改变型 SNV/InDel。WES 与 WGS 的测序深度、文库和检出阈值不同，因此 VAF 分布差异不能单独解释为肿瘤克隆演化。</p>
{END}"""
    if START in report and END in report:
        before, rest = report.split(START, 1)
        _, after = rest.split(END, 1)
        report = before + block + after
    else:
        anchor = "<h3>差异位点回查结果</h3>"
        if anchor not in report:
            raise SystemExit(f"Report anchor not found: {anchor}")
        report = report.replace(anchor, block + "\n" + anchor, 1)
    args.patient_report.write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
