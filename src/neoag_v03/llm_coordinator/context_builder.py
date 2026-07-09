from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from neoag_v03.agents.skill_router import find_named_files
from .schemas import FileRecord


def infer_file_kind(path: Path) -> str:
    n = path.name.lower()
    if "ranked_peptides.recommendation" in n:
        return "ranked_peptides_recommendation"
    if "ranked_peptides.netmhcpan42" in n:
        return "ranked_peptides_netmhcpan42"
    if n.startswith("evidence_report") and n.endswith(".html"):
        return "evidence_report"
    if "spechla" in n or "lohhla" in n or "hla_loh" in n:
        return "hla_loh"
    if "facets" in n or "purple" in n or "purity" in n:
        return "purity_or_cnv"
    if "appm" in n:
        return "appm"
    if "ccf" in n:
        return "ccf"
    if n.endswith(".docx"):
        return "docx"
    if n.endswith(".pptx"):
        return "pptx"
    if n.endswith(".tsv"):
        return "tsv"
    if n.endswith(".vcf") or n.endswith(".vcf.gz"):
        return "vcf"
    return "unknown"


def build_file_records(files: list[str] | None = None, result_dir: str | None = None) -> list[FileRecord]:
    paths: list[Path] = []
    for f in files or []:
        paths.append(Path(f))
    if result_dir:
        base = Path(result_dir)
        if base.exists():
            paths.extend([p for p in base.rglob("*") if p.is_file()])
    seen: set[str] = set()
    records: list[FileRecord] = []
    for p in paths:
        key = str(p)
        if key in seen:
            continue
        seen.add(key)
        exists = p.exists()
        size = p.stat().st_size if exists else None
        records.append(FileRecord(path=str(p), name=p.name, kind=infer_file_kind(p), exists=exists, size_bytes=size))
    return records


def build_context(message: str, files: list[str] | None = None, result_dir: str | None = None, project_root: str = ".") -> dict[str, Any]:
    records = build_file_records(files, result_dir)
    named = find_named_files(result_dir, files)
    return {
        "message": message,
        "project_root": str(Path(project_root).resolve()),
        "result_dir": result_dir,
        "available_files": [r.to_dict() for r in records[:200]],
        "file_count": len(records),
        "named_files": named,
        "cwd": os.getcwd(),
    }
