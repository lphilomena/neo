from __future__ import annotations

import csv
import gzip
import json
import re
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def open_text(path: str | Path):
    p = Path(path)
    if str(p).endswith(".gz"):
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return p.open("r", encoding="utf-8-sig", errors="replace")


def detect_delimiter(path: str | Path) -> str:
    with open_text(path) as fh:
        sample = "".join([fh.readline() for _ in range(5)])
    if sample.count("\t") >= sample.count(","):
        return "\t"
    return ","


def read_table(path: str | Path, delimiter: str | None = None, limit: int | None = None) -> tuple[list[str], list[dict[str, str]]]:
    delim = delimiter or detect_delimiter(path)
    with open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        header = reader.fieldnames or []
        rows: list[dict[str, str]] = []
        for i, row in enumerate(reader):
            if limit is not None and i >= limit:
                break
            rows.append({str(k): (v if v is not None else "") for k, v in row.items() if k is not None})
    return header, rows


def write_tsv(path: str | Path, rows: Sequence[Mapping[str, Any]], fieldnames: Sequence[str] | None = None) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        keys: list[str] = []
        for row in rows:
            for key in row.keys():
                if key not in keys:
                    keys.append(str(key))
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


def safe_float(value: Any, default: float | None = None) -> float | None:
    try:
        s = str(value if value is not None else "").strip()
        if not s or s.lower() in {"na", "nan", "none", "null", "."}:
            return default
        return float(s)
    except Exception:
        return default


def normalize_hla(value: str) -> str:
    s = str(value or "").strip().replace("_", "*")
    s = s.replace("HLA-", "HLA-").replace("hla-", "HLA-")
    if not s:
        return ""
    s = s.replace(" ", "")
    # Common variants: A*02:01, HLA_A_02_01, HLA-A0201, A0201
    if s.startswith("HLA") and "*" not in s and ":" not in s:
        m = re.match(r"HLA[-_]?([A-Z0-9]+)[-_]?(\d{2})(\d{2})", s)
        if m:
            return f"HLA-{m.group(1)}*{m.group(2)}:{m.group(3)}"
    if re.match(r"^[A-Z0-9]+\*\d{2}:\d{2}$", s):
        return "HLA-" + s
    if re.match(r"^[A-Z0-9]+\d{4}$", s):
        return f"HLA-{s[:-4]}*{s[-4:-2]}:{s[-2:]}"
    if s.startswith("HLA-"):
        return s.replace("_", ":")
    return s


def row_get(row: Mapping[str, Any], candidates: Sequence[str], default: str = "") -> str:
    lower = {str(k).lower(): k for k in row.keys()}
    for c in candidates:
        if c in row:
            return str(row.get(c, "") or "")
        lc = c.lower()
        if lc in lower:
            return str(row.get(lower[lc], "") or "")
    return default


def parse_vcf_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in info.split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        else:
            out[part] = "true"
    return out


def parse_csq_header(path: str | Path) -> list[str]:
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                m = re.search(r"Format: ([^\"]+)", line)
                if m:
                    return [x.strip() for x in m.group(1).split("|")]
            if line.startswith("#CHROM"):
                break
    return []


def read_vcf_records(path: str | Path, limit: int | None = None) -> tuple[list[str], list[dict[str, str]]]:
    csq_fields = parse_csq_header(path)
    records: list[dict[str, str]] = []
    sample_names: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                sample_names = parts[9:]
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt, qual, flt, info = parts[:8]
            info_map = parse_vcf_info(info)
            alts = [a for a in alt.split(",") if a]
            csq_entries: list[dict[str, str]] = []
            csq = info_map.get("CSQ", "")
            if csq and csq_fields:
                for entry in csq.split(","):
                    vals = entry.split("|")
                    csq_entries.append({k: vals[i] if i < len(vals) else "" for i, k in enumerate(csq_fields)})
            for alt_i, one_alt in enumerate(alts or [alt]):
                rec = {
                    "chrom": chrom, "pos": pos, "id": vid, "ref": ref, "alt": one_alt, "qual": qual, "filter": flt, "info": info,
                    "svtype": info_map.get("SVTYPE", ""), "gene": "", "transcript": "", "consequence": "", "protein_change": "", "allele_index": str(alt_i + 1),
                }
                csq_map: dict[str, str] = {}
                if csq_entries:
                    csq_map = next((x for x in csq_entries if x.get("Allele") == one_alt), csq_entries[min(alt_i, len(csq_entries) - 1)])
                    rec["gene"] = csq_map.get("SYMBOL") or csq_map.get("Gene") or csq_map.get("Feature") or ""
                    rec["transcript"] = csq_map.get("Feature", "")
                    rec["consequence"] = csq_map.get("Consequence", "")
                    rec["protein_change"] = csq_map.get("HGVSp") or csq_map.get("Protein_position") or ""
                records.append(rec)
                if limit is not None and len(records) >= limit:
                    return sample_names, records
    return sample_names, records


def read_fasta_sequences(path: str | Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name:
        seqs[name] = "".join(chunks).upper()
    return seqs


def markdown_table(rows: Sequence[Mapping[str, Any]], columns: Sequence[str] | None = None, max_rows: int = 20) -> str:
    if not rows:
        return "（无记录）\n"
    if columns is None:
        columns = list(rows[0].keys())
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows[:max_rows]:
        lines.append("| " + " | ".join(str(row.get(c, "")).replace("|", "\\|").replace("\n", " ") for c in columns) + " |")
    if len(rows) > max_rows:
        lines.append(f"\n仅显示前 {max_rows} 行，共 {len(rows)} 行。")
    return "\n".join(lines) + "\n"
