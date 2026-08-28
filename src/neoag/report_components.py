"""Shared HTML report components."""

from __future__ import annotations

from collections import Counter
import html
from typing import Any, Mapping

DEFAULT_REPORT_CSS = """
<style>
body{font-family:Arial,sans-serif;margin:32px;color:#222;line-height:1.45;max-width:1100px}
h1,h2,h3{color:#17324d}.section{margin-top:28px}
table{border-collapse:collapse;width:100%;margin:12px 0 24px}th,td{border:1px solid #ddd;padding:7px;font-size:12px;vertical-align:top}th{background:#f3f6f9;position:sticky;top:0}
.badge{padding:3px 7px;border-radius:8px;font-size:12px;display:inline-block}.PASS{background:#d6f5d6}.CAUTION{background:#fff1b8}.FAIL{background:#ffd6d6}.UNASSESSED{background:#eee;color:#555}
.card{border:1px solid #ddd;border-radius:10px;padding:14px;margin:12px 0;box-shadow:0 1px 4px #eee}
.small{color:#555;font-size:13px}.mono{font-family:Menlo,Consolas,monospace;font-size:11px;word-break:break-all}
.warn{background:#fff7e6;border-left:4px solid #e6a700;padding:12px;margin:14px 0}
.info{background:#f0f7ff;border-left:4px solid #3b82f6;padding:12px;margin:14px 0}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px}
.metric{border:1px solid #ddd;border-radius:8px;padding:10px;background:#fafafa}
.metric b.value{font-size:1.45rem;color:#17324d}
ul.compact{margin:8px 0 8px 20px}
.patient h1{font-size:1.6rem}.patient .lead{font-size:1.05rem;color:#333}
@media print{body{margin:18px;max-width:none}.card{box-shadow:none}th{position:static}.no-print{display:none}}
</style>
"""


def esc(value: Any) -> str:
    return html.escape(str(value if value is not None else ""))


def badge(text: Any) -> str:
    raw = str(text or "")
    upper = raw.upper()
    cls = "UNASSESSED"
    if any(token in upper for token in ["PASS", "INTACT", "HIGH", "A", "B"]):
        cls = "PASS"
    if any(token in upper for token in ["CAUTION", "REVIEW", "MEDIUM", "LOW", "C"]):
        cls = "CAUTION"
    if any(token in upper for token in ["DEFECT", "REJECT", "FAIL", "LOST", "GLOBAL", "D"]):
        cls = "FAIL"
    if any(token in upper for token in ["UNASSESSED", "INSUFFICIENT", "MISSING", "INCONCLUSIVE", "UNKNOWN"]):
        cls = "UNASSESSED"
    return f"<span class='badge {cls}'>{esc(raw)}</span>"


def table(rows: list[Mapping[str, Any]], headers: list[str], *, max_rows: int | None = None) -> str:
    view = rows[:max_rows] if max_rows else rows
    out = ["<table><tr>" + "".join(f"<th>{esc(header)}</th>" for header in headers) + "</tr>"]
    for row in view:
        out.append("<tr>" + "".join(f"<td>{esc(row.get(header, ''))}</td>" for header in headers) + "</tr>")
    if max_rows and len(rows) > max_rows:
        out.append(f"<tr><td colspan='{len(headers)}' class='small'>Showing {max_rows} of {len(rows)} rows.</td></tr>")
    out.append("</table>")
    return "\n".join(out)


def count_values(rows: list[Mapping[str, Any]], field: str) -> Counter[str]:
    counts: Counter[str] = Counter()
    for row in rows:
        counts[str(row.get(field) or "UNASSESSED")] += 1
    return counts


def count_matching(rows: list[Mapping[str, Any]], field: str, tokens: tuple[str, ...]) -> int:
    total = 0
    for row in rows:
        value = str(row.get(field) or "").upper()
        if any(token in value for token in tokens):
            total += 1
    return total


def metric_card(label: str, value: Any, detail: Any = "") -> str:
    detail_html = f"<div class='small'>{esc(detail)}</div>" if detail else ""
    return f"<div class='metric'><div>{esc(label)}</div><b class='value'>{esc(value)}</b>{detail_html}</div>"
