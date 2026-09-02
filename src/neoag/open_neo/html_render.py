from __future__ import annotations

import html
import re


_TABLE_SEPARATOR = re.compile(r"^:?-{3,}:?$")


def _inline(text: str) -> str:
    escaped = html.escape(text, quote=True)
    escaped = re.sub(r"`([^`]+)`", r"<code>\1</code>", escaped)
    escaped = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", escaped)
    return escaped


def _table_cells(line: str) -> list[str] | None:
    stripped = line.strip()
    if not stripped.startswith("|") or not stripped.endswith("|"):
        return None
    return [cell.strip() for cell in stripped[1:-1].split("|")]


def _is_table_separator(cells: list[str] | None) -> bool:
    return bool(cells) and all(_TABLE_SEPARATOR.fullmatch(cell.replace(" ", "")) for cell in cells)


def markdown_to_html(markdown: str, *, title: str = "Open-Neo report", lang: str = "zh-CN") -> str:
    """Render the Markdown subset emitted by Open-Neo reports as standalone HTML."""
    lines = markdown.splitlines()
    body: list[str] = []
    index = 0
    while index < len(lines):
        line = lines[index]
        stripped = line.strip()
        if not stripped:
            index += 1
            continue

        heading = re.match(r"^(#{1,6})\s+(.+)$", stripped)
        if heading:
            level = len(heading.group(1))
            body.append(f"<h{level}>{_inline(heading.group(2))}</h{level}>")
            index += 1
            continue

        header_cells = _table_cells(line)
        separator_cells = _table_cells(lines[index + 1]) if index + 1 < len(lines) else None
        if header_cells is not None and _is_table_separator(separator_cells):
            body.append('<div class="table-wrap"><table><thead><tr>')
            body.extend(f"<th>{_inline(cell)}</th>" for cell in header_cells)
            body.append("</tr></thead><tbody>")
            index += 2
            while index < len(lines):
                row_cells = _table_cells(lines[index])
                if row_cells is None:
                    break
                padded = (row_cells + [""] * len(header_cells))[: len(header_cells)]
                body.append("<tr>" + "".join(f"<td>{_inline(cell)}</td>" for cell in padded) + "</tr>")
                index += 1
            body.append("</tbody></table></div>")
            continue

        if stripped.startswith("- "):
            body.append("<ul>")
            while index < len(lines) and lines[index].strip().startswith("- "):
                body.append(f"<li>{_inline(lines[index].strip()[2:])}</li>")
                index += 1
            body.append("</ul>")
            continue

        paragraph = [stripped]
        index += 1
        while index < len(lines):
            candidate = lines[index].strip()
            if not candidate or candidate.startswith("#") or candidate.startswith("- "):
                break
            next_cells = _table_cells(lines[index + 1]) if index + 1 < len(lines) else None
            if _table_cells(lines[index]) is not None and _is_table_separator(next_cells):
                break
            paragraph.append(candidate)
            index += 1
        body.append("<p>" + "<br>".join(_inline(item) for item in paragraph) + "</p>")

    css = """
body{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","PingFang SC","Microsoft YaHei",Arial,sans-serif;color:#202124;line-height:1.6;max-width:1480px;margin:0 auto;padding:32px;background:#fff}
h1{font-size:30px;border-bottom:2px solid #1f6f8b;padding-bottom:12px}h2{font-size:22px;margin-top:32px;color:#174c5f}h3{font-size:18px}p{margin:10px 0 16px}.table-wrap{overflow-x:auto;margin:16px 0 24px;border:1px solid #d9e1e5}table{border-collapse:collapse;width:100%;font-size:13px}th,td{border:1px solid #d9e1e5;padding:8px 10px;text-align:left;vertical-align:top;white-space:normal;overflow-wrap:anywhere}th{background:#e8f2f5;font-weight:650;position:sticky;top:0}tbody tr:nth-child(even){background:#f7f9fa}code{background:#f1f3f4;padding:1px 4px;border-radius:3px}ul{padding-left:24px}@media print{body{max-width:none;padding:12px}.table-wrap{overflow:visible}th{position:static}}
""".strip()
    return (
        "<!doctype html>\n"
        f'<html lang="{html.escape(lang, quote=True)}"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        f"<title>{html.escape(title)}</title><style>{css}</style></head><body>"
        + "\n".join(body)
        + "</body></html>\n"
    )
