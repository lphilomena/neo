"""Shared adapter helpers."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Iterable, Mapping

from neoag.utils import MISSING, first, open_text_maybe_gz, to_float


def clean(value: Any) -> str:
    text = str(value if value is not None else "").strip()
    return "" if text in MISSING else text


def as_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(clean(value)))
    except Exception:
        return default


def as_float_text(value: Any) -> str:
    text = clean(value)
    if not text:
        return ""
    try:
        return f"{float(text):.12g}"
    except Exception:
        return ""


def truth_text(value: Any) -> str:
    return "true" if clean(value).casefold() in {"1", "true", "yes", "y", "pass", "present"} else "false"


def split_tokens(value: Any) -> list[str]:
    result: list[str] = []
    for token in re.split(r"[;,]", clean(value)):
        item = token.strip()
        if item and item not in result:
            result.append(item)
    return result


def join_tokens(values: Iterable[Any]) -> str:
    tokens: set[str] = set()
    for value in values:
        if isinstance(value, (list, tuple, set, frozenset)):
            tokens.update(clean(x) for x in value if clean(x))
        else:
            tokens.update(split_tokens(value))
    return ";".join(sorted(tokens))


def row_hash(row: Mapping[str, Any]) -> str:
    payload = json.dumps({str(k): clean(v) for k, v in sorted(row.items())}, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def source_record_id(tool: str, path: str | Path, row_number: int, row: Mapping[str, Any]) -> str:
    return f"{tool}|{Path(path).name}|{row_number}|{row_hash(row)[:16]}"


def parse_attributes(value: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in value.strip().strip(";").split(";"):
        if not item.strip():
            continue
        if "=" in item:
            key, val = item.split("=", 1)
        elif " " in item:
            key, val = item.split(" ", 1)
            val = val.strip().strip('"')
        else:
            continue
        attrs[key.strip()] = val.strip().strip('"')
    return attrs


def read_delimited(path: str | Path) -> list[dict[str, str]]:
    """Read a headered tab/comma-delimited file, including gzip content."""
    p = Path(path)
    with open_text_maybe_gz(p) as fh:
        sample = fh.read(8192)
        fh.seek(0)
        delimiter = "\t" if sample.count("\t") >= sample.count(",") else ","
        reader = csv.DictReader(fh, delimiter=delimiter)
        return [{str(k): clean(v) for k, v in row.items() if k is not None} for row in reader]


def get(row: Mapping[str, Any], *aliases: str, default: str = "") -> str:
    return first(row, list(aliases), default=default)


def infer_mhc_class(allele: str) -> str:
    upper = clean(allele).upper()
    if any(marker in upper for marker in ("DPA", "DPB", "DQA", "DQB", "DRA", "DRB")):
        return "II"
    return "I" if upper else ""
