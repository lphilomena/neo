from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def read_tsv(path: str | Path, limit: int | None = None) -> tuple[list[str], list[dict[str, str]]]:
    p = Path(path)
    with p.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        rows: list[dict[str, str]] = []
        for i, row in enumerate(reader):
            if limit is not None and i >= limit:
                break
            rows.append({k: (v if v is not None else "") for k, v in row.items()})
    return header, rows


def write_tsv(path: str | Path, rows: Sequence[Mapping[str, Any]], fieldnames: Sequence[str] | None = None) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        keys: list[str] = []
        for row in rows:
            for k in row.keys():
                if k not in keys:
                    keys.append(k)
        fieldnames = keys
    with p.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(fieldnames), delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def write_json(path: str | Path, data: Any) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def read_json(path: str | Path) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def safe_float(value: Any, default: float | None = None) -> float | None:
    try:
        if value is None:
            return default
        s = str(value).strip()
        if s == "" or s.lower() in {"na", "nan", "none", "null"}:
            return default
        return float(s)
    except Exception:
        return default


def safe_int(value: Any, default: int = 0) -> int:
    f = safe_float(value, None)
    if f is None:
        return default
    return int(f)


def count_by(rows: Sequence[Mapping[str, Any]], field: str) -> Counter[str]:
    c: Counter[str] = Counter()
    for row in rows:
        val = str(row.get(field, "") or "").strip() or "NA"
        c[val] += 1
    return c


def find_files(root: str | Path, patterns: Sequence[str]) -> list[Path]:
    base = Path(root)
    if not base.exists():
        return []
    out: list[Path] = []
    for pat in patterns:
        out.extend(base.rglob(pat))
    seen: set[str] = set()
    uniq: list[Path] = []
    for p in out:
        key = str(p.resolve())
        if key not in seen and p.is_file():
            seen.add(key)
            uniq.append(p)
    return sorted(uniq)


def executable_status(exe: str) -> dict[str, str]:
    path = shutil.which(exe)
    return {
        "tool": exe,
        "status": "OK" if path else "MISSING",
        "path": path or "",
        "message": "found on PATH" if path else "not found on PATH",
    }


def run_command(cmd: Sequence[str], cwd: str | Path | None = None, timeout: int = 60) -> dict[str, Any]:
    try:
        proc = subprocess.run(list(cmd), cwd=str(cwd) if cwd else None, text=True, capture_output=True, timeout=timeout)
        return {
            "cmd": " ".join(cmd),
            "returncode": proc.returncode,
            "stdout": proc.stdout[-5000:],
            "stderr": proc.stderr[-5000:],
            "ok": proc.returncode == 0,
        }
    except Exception as e:
        return {"cmd": " ".join(cmd), "returncode": 999, "stdout": "", "stderr": str(e), "ok": False}


def markdown_table(rows: Sequence[Mapping[str, Any]], columns: Sequence[str] | None = None, max_rows: int = 30) -> str:
    if not rows:
        return "\n（无记录）\n"
    if columns is None:
        columns = list(rows[0].keys())
    def cell(x: Any) -> str:
        s = str(x if x is not None else "")
        return s.replace("|", "\\|").replace("\n", " ")
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows[:max_rows]:
        lines.append("| " + " | ".join(cell(row.get(c, "")) for c in columns) + " |")
    if len(rows) > max_rows:
        lines.append(f"\n仅显示前 {max_rows} 行，共 {len(rows)} 行。")
    return "\n".join(lines) + "\n"


def strip_html_text(html: str) -> str:
    text = re.sub(r"<script.*?</script>", " ", html, flags=re.S | re.I)
    text = re.sub(r"<style.*?</style>", " ", text, flags=re.S | re.I)
    text = re.sub(r"<[^>]+>", " ", text)
    text = re.sub(r"&nbsp;", " ", text)
    text = re.sub(r"&gt;", ">", text)
    text = re.sub(r"&lt;", "<", text)
    text = re.sub(r"&amp;", "&", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def percentile(sorted_vals: list[float], pct: float) -> float | None:
    if not sorted_vals:
        return None
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    k = (len(sorted_vals) - 1) * pct
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return sorted_vals[int(k)]
    return sorted_vals[f] * (c - k) + sorted_vals[c] * (k - f)
