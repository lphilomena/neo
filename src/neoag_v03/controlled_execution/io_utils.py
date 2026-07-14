from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import sys
import time
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def to_jsonable(obj: Any) -> Any:
    if is_dataclass(obj):
        return asdict(obj)
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, dict):
        return {str(k): to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [to_jsonable(v) for v in obj]
    return obj


def write_json(path: str | Path, data: Any) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(to_jsonable(data), ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return p


def read_json(path: str | Path) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_tsv(path: str | Path, rows: Sequence[Mapping[str, Any]], fieldnames: Sequence[str] | None = None) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    if fieldnames is None:
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


def append_jsonl(path: str | Path, item: Any) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(to_jsonable(item), ensure_ascii=False, sort_keys=True) + "\n")
    return p


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def file_size(path: str | Path) -> int:
    return Path(path).stat().st_size


def markdown_table(rows: Sequence[Mapping[str, Any]], columns: Sequence[str] | None = None, max_rows: int = 50) -> str:
    if not rows:
        return "\n（无记录）\n"
    if columns is None:
        columns = list(rows[0].keys())
    def cell(v: Any) -> str:
        s = "" if v is None else str(v)
        return s.replace("|", "\\|").replace("\n", " ")
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows[:max_rows]:
        lines.append("| " + " | ".join(cell(row.get(c, "")) for c in columns) + " |")
    if len(rows) > max_rows:
        lines.append(f"\n... truncated {len(rows) - max_rows} rows ...")
    return "\n".join(lines)


def safe_rel(path: str | Path, root: str | Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(Path(root).resolve()))
    except Exception:
        return str(path)


def detect_private_path(text: str) -> bool:
    patterns = [
        r"/mnt/[A-Za-z0-9_\-./]+",
        r"/home/[A-Za-z0-9_\-./]+",
        r"/root/[A-Za-z0-9_\-./]+",
        r"[A-Z]:\\\\",
    ]
    return any(re.search(p, text) for p in patterns)


def load_limited_yaml(path: str | Path) -> dict[str, Any]:
    """Dependency-light YAML reader for simple manifest files.

    Supports JSON directly. For YAML it supports nested dictionaries via
    indentation, scalar values, and top-level lists written as '- item'. This is
    intentionally small; if PyYAML is installed it is used automatically.
    """
    p = Path(path)
    text = p.read_text(encoding="utf-8")
    if p.suffix.lower() == ".json":
        return json.loads(text)
    try:
        import yaml  # type: ignore
        loaded = yaml.safe_load(text)
        return loaded or {}
    except Exception:
        pass
    root: dict[str, Any] = {}
    stack: list[tuple[int, dict[str, Any]]] = [(-1, root)]
    for raw in text.splitlines():
        line = raw.rstrip()
        if not line or line.strip().startswith("#"):
            continue
        indent = len(line) - len(line.lstrip(" "))
        stripped = line.strip()
        while stack and indent <= stack[-1][0]:
            stack.pop()
        parent = stack[-1][1]
        if stripped.startswith("- "):
            # minimal support: append to special _list on parent
            parent.setdefault("_list", []).append(_parse_scalar(stripped[2:].strip()))
            continue
        if ":" not in stripped:
            continue
        key, val = stripped.split(":", 1)
        key = key.strip()
        val = val.strip()
        if val == "":
            child: dict[str, Any] = {}
            parent[key] = child
            stack.append((indent, child))
        else:
            parent[key] = _parse_scalar(val)
    return root


def _parse_scalar(s: str) -> Any:
    if s in {"null", "None", "~"}:
        return None
    if s.lower() == "true":
        return True
    if s.lower() == "false":
        return False
    if (s.startswith('"') and s.endswith('"')) or (s.startswith("'") and s.endswith("'")):
        return s[1:-1]
    try:
        if "." in s:
            return float(s)
        return int(s)
    except Exception:
        return s
