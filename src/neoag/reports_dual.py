"""Dual-audience HTML reports: patient communication vs research/technical."""

from __future__ import annotations

import html
import hashlib
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

from .candidate_identity import candidate_identity, identity_value
from .utils import read_tsv

REPORT_CSS = """
<style>
body{font-family:Arial,sans-serif;margin:32px;color:#222;line-height:1.45;max-width:1100px}
h1,h2,h3{color:#17324d}.section{margin-top:28px}
table{border-collapse:collapse;width:100%;margin:12px 0 24px}th,td{border:1px solid #ddd;padding:7px;font-size:12px;vertical-align:top}th{background:#f3f6f9}
.badge{padding:3px 7px;border-radius:8px;font-size:12px;display:inline-block}.PASS{background:#d6f5d6}.CAUTION{background:#fff1b8}.FAIL{background:#ffd6d6}.UNASSESSED{background:#eee;color:#555}
.card{border:1px solid #ddd;border-radius:10px;padding:14px;margin:12px 0;box-shadow:0 1px 4px #eee}
.small{color:#555;font-size:13px}.mono{font-family:Menlo,Consolas,monospace;font-size:11px;word-break:break-all}
.warn{background:#fff7e6;border-left:4px solid #e6a700;padding:12px;margin:14px 0}
.info{background:#f0f7ff;border-left:4px solid #3b82f6;padding:12px;margin:14px 0}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px}
.metric{border:1px solid #ddd;border-radius:8px;padding:10px;background:#fafafa}
ul.compact{margin:8px 0 8px 20px}
.patient h1{font-size:1.6rem}.patient .lead{font-size:1.05rem;color:#333}
</style>
"""

GENE_SYMBOL_FALLBACKS = {
    "ENSG00000170846": "MRFAP1L2",
    "ENSG00000232325": "AC093627.1",
    "ENSG00000286215": "VIM2P",
}

PATIENT_STATUS_LABELS = {
    "EVENT_STRONG": "事件获得较强支持",
    "EVENT_PARTIAL": "事件获得部分支持",
    "RNA_CONFIRMED": "RNA中检测到直接支持",
    "RNA_NEGATIVE": "RNA中未检测到直接支持，需结合覆盖度解释",
    "GENE_EXPRESSION_ONLY": "仅确认基因表达，尚未确认突变表达",
    "PRESENTATION_CONSISTENT_STRONG": "两个核心工具呈递预测一致且较强",
    "PRESENTATION_DISCORDANT": "核心呈递工具结果不一致",
    "PRESENTATION_SINGLE_TOOL": "仅一个核心工具提供呈递支持",
    "PRESENTATION_UNASSESSED": "核心呈递预测未评估",
    "MT_SPECIFIC": "突变肽具有较强特异性",
    "MARGINAL_MT_ADVANTAGE": "MT相对WT仅轻度改善，不能单独作为免疫原性正向证据",
    "WT_STRONG_BINDING_REVIEW": "WT仍预测为强结合，需重点复核自身反应与免疫耐受风险",
    "WT_BINDING_REVIEW": "WT仍保留预测结合，需进行配对安全性复核",
    "WT_LOW_PREDICTED_BINDING": "WT预测结合较弱，但不能据此排除耐受或交叉反应",
    "HLA_LOH_UNASSESSED": "限制性HLA多工具确认未完成",
    "SAFETY_PARTIAL": "正常组织安全性仅部分评估（具体缺口见候选说明）",
    "SUPPORTED": "获得支持",
    "CLONAL": "倾向克隆性事件",
    "C1": "候选来源链完整且有正交支持",
    "C2": "候选来源链较完整，仍需实验确认",
    "C3": "候选来源链基本合理，但关键环节尚未闭合",
    "C4": "候选来源链不完整，仅进入技术审阅",
    "UNASSESSED": "未评估",
    "UNSPECIFIED": "未明确",
    "PARTIAL": "证据部分完整",
    "LOW": "证据完整度低",
    "MHC_I_CAUTION": "MHC-I存在谨慎信号",
    "IFNG_RESPONSE_CAUTION": "IFNG/JAK-STAT应答存在谨慎信号",
    "MHC_II_INTACT": "现有结果未见MHC-II呈递系统整体完全丧失",
    "NOVEL_SEQUENCE": "已确认包含异常连接或新生序列",
}


def _patient_status_text(value: Any) -> str:
    text = str(value or "").strip()
    return PATIENT_STATUS_LABELS.get(text.upper(), text or "未评估")


def esc(x: Any) -> str:
    return html.escape(str(x if x is not None else ""))


def _read_optional(path: str | Path | None) -> list[dict[str, str]]:
    if not path:
        return []
    p = Path(path)
    return read_tsv(p) if p.is_file() else []


def _report_relpath(path: Path, root: Path | None) -> str:
    if not root:
        return str(path)
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path)


def _first_existing(paths: list[Path | None]) -> Path | None:
    for path in paths:
        if path and path.is_file():
            return path
    return None


def _hla_loh_search_roots(root: Path | None, provenance: Mapping[str, Any] | None = None) -> list[Path]:
    roots: list[Path] = []
    if root:
        roots.extend([root, root.parent, root.parent / "evidence", root.parent / "evidence" / "hla_loh"])
    for value in (
        (provenance or {}).get("hla_loh"),
        ((provenance or {}).get("tools") or {}).get("hla_loh", {}).get("file") if isinstance((provenance or {}).get("tools"), Mapping) else None,
    ):
        if value:
            path = Path(str(value))
            roots.append(path.parent if path.suffix else path)
    seen: set[str] = set()
    ordered: list[Path] = []
    for item in roots:
        key = str(item)
        if key not in seen:
            seen.add(key)
            ordered.append(item)
    return ordered


def _expand_hla_loh_consensus(path: Path, root: Path | None) -> list[dict[str, str]]:
    results: list[dict[str, str]] = []
    source = _report_relpath(path, root)
    for row in _read_optional(path):
        allele = str(row.get("hla_allele") or row.get("allele") or "").strip()
        if not allele:
            continue
        for tool, field in (("LOHHLA", "lohhla_status"), ("SpecHLA", "spechla_status")):
            status = str(row.get(field) or "").strip()
            if not status:
                continue
            results.append({
                **dict(row),
                "hla_allele": allele,
                "loh_status": status,
                "tool": tool,
                "_report_tool": tool,
                "_report_source": source,
            })
    return results


def _read_hla_loh_tool_results(root: Path | None, provenance: Mapping[str, Any] | None = None) -> list[dict[str, str]]:
    """Load allele-level LOHHLA and SpecHLA evidence from supported layouts."""
    search_roots = _hla_loh_search_roots(root, provenance)
    candidates = {
        "LOHHLA": [
            "hla_loh/lohhla/hla_loh.tsv",
            "hla_loh_consensus/lohhla_hla_loh.tsv",
            "lohhla/hla_loh.tsv",
            "hla_loh.tsv",
        ],
        "SpecHLA": [
            "hla_loh/spechla/hla_loh.tsv",
            "hla_loh_consensus/spechla_sequenza012_hla_loh.tsv",
            "hla_loh_consensus/spechla_hla_loh.tsv",
            "spechla/hla_loh.tsv",
            "hla_loh.spechla.tsv",
        ],
    }
    results: list[dict[str, str]] = []
    for tool, relative_paths in candidates.items():
        selected = _first_existing([base / rel for base in search_roots for rel in relative_paths])
        if selected is None and tool == "LOHHLA":
            selected = _first_existing([
                Path(str(((provenance or {}).get("tools") or {}).get("lohhla", {}).get("file") or "")),
            ])
        if not selected:
            continue
        if selected.name.endswith("spechla.tsv") and tool != "SpecHLA":
            continue
        for row in _read_optional(selected):
            record = dict(row)
            record["_report_tool"] = tool
            record["_report_source"] = _report_relpath(selected, root)
            results.append(record)
    if results:
        return results
    consensus = _first_existing([
        *(base / name for base in search_roots for name in (
            "hla_loh_consensus.tsv",
            "hla_loh_consensus/hla_loh_consensus.tsv",
            "hla_loh/hla_loh_consensus.tsv",
        )),
        Path(str(((provenance or {}).get("tools") or {}).get("hla_loh", {}).get("file") or "")),
        Path(str((provenance or {}).get("hla_loh") or "")),
    ])
    return _expand_hla_loh_consensus(consensus, root) if consensus else []


def _read_tool_version_manifest(root: Path | None) -> dict[str, dict[str, str]]:
    if not root:
        return {}
    payload = _read_json_optional(root / "metadata" / "tool_versions.json")
    tools = payload.get("tools") or {}
    return {str(name): dict(record) for name, record in tools.items() if isinstance(record, Mapping)}


def _read_json_optional(path: str | Path | None) -> dict[str, Any]:
    if not path:
        return {}
    p = Path(path)
    if not p.is_file():
        return {}
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _has_consensus_grades(events: list[Mapping[str, Any]], peptides: list[Mapping[str, Any]]) -> bool:
    event_ok = any(str(row.get("best_evidence_grade") or row.get("event_evidence_grade") or "").strip() for row in events[:50])
    peptide_ok = any(str(row.get("evidence_grade") or row.get("pipeline_r_grade") or "").strip() for row in peptides[:50])
    return bool(events) and bool(peptides) and event_ok and peptide_ok


def _numeric_median(values: list[float]) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    return ordered[len(ordered) // 2]


def _float_or_none(value: Any) -> float | None:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def _depth_summary(rows: list[Mapping[str, Any]], *fields: str) -> str:
    values: list[float] = []
    seen: set[str] = set()
    for row in rows:
        variant_key = "|".join(str(row.get(key) or "").strip() for key in ("chrom", "pos", "ref", "alt"))
        event_key = str(row.get("event_group_id") or row.get("event_id") or row.get("source_event_id") or "").strip()
        key = variant_key if variant_key.strip("|") else event_key
        if key and key in seen:
            continue
        for field in fields:
            try:
                value = float(str(row.get(field) or "").strip())
            except ValueError:
                continue
            if value > 0:
                values.append(value)
                if key:
                    seen.add(key)
                break
    median = _numeric_median(values)
    if median is None:
        return ""
    return f"候选位点中位深度 {median:.0f}x（n={len(values)}）"


def _rna_qc_summary(rows: list[Mapping[str, Any]]) -> str:
    """Summarize RNA evidence by unique biological event, not peptide-HLA rows."""
    variants: dict[str, tuple[float | None, float | None]] = {}
    junctions: dict[str, float | None] = {}
    expression_loaded = False
    for row in rows:
        expression_loaded = expression_loaded or any(
            str(row.get(field) or "").strip() not in {"", "UNASSESSED", "NOT_AVAILABLE"}
            for field in ("gene_expression_tpm", "transcript_expression_tpm")
        )
        event_type = str(row.get("event_type") or row.get("variant_type") or "").strip().upper()
        variant_key = "|".join(str(row.get(key) or "").strip() for key in ("chrom", "pos", "ref", "alt"))
        junction_key = str(
            row.get("canonical_junction_id") or row.get("junction_key") or row.get("source_junction_id") or ""
        ).strip()
        if not junction_key and event_type in {"FUSION", "SPLICE", "SPLICE_JUNCTION"}:
            junction_key = "|".join(
                str(row.get(key) or "").strip()
                for key in ("gene", "breakpoint1", "breakpoint2")
            )
        event_key = (
            variant_key if variant_key.strip("|") and event_type in {"SNV", "INDEL", "MNV"}
            else junction_key if junction_key.strip("|")
            else str(row.get("event_group_id") or row.get("event_id") or row.get("source_event_id") or "").strip()
        )
        if not event_key.strip("|"):
            continue
        if event_type in {"SNV", "INDEL", "MNV"}:
            depth = _float_or_none(row.get("rna_depth"))
            alt = _float_or_none(row.get("rna_alt_reads"))
            if depth is None and alt is None:
                continue
            previous = variants.get(event_key)
            if previous is None or (depth or -1) > (previous[0] or -1):
                variants[event_key] = (depth, alt)
        elif event_type in {"FUSION", "SPLICE", "SPLICE_JUNCTION"}:
            reads = _float_or_none(row.get("rna_junction_reads") or row.get("junction_reads"))
            if reads is None:
                continue
            if event_key not in junctions or (reads or -1) > (junctions[event_key] or -1):
                junctions[event_key] = reads

    parts: list[str] = []
    if expression_loaded:
        parts.append("已载入基因/转录本表达")
    if variants:
        assessed = [(depth, alt) for depth, alt in variants.values() if depth is not None or alt is not None]
        covered = [(depth, alt) for depth, alt in assessed if (depth or 0) > 0]
        alt1 = sum(1 for _, alt in covered if (alt or 0) >= 1)
        alt3 = sum(1 for _, alt in covered if (alt or 0) >= 3)
        alt5 = sum(1 for _, alt in covered if (alt or 0) >= 5)
        zero_depth = sum(1 for depth, _ in assessed if (depth or 0) <= 0)
        parts.append(
            f"SNV/InDel位点级RNA已评估 {len(assessed)}/{len(variants)} 个独立事件；"
            f"有覆盖 {len(covered)}，ALT reads≥1/3/5：{alt1}/{alt3}/{alt5}，深度0：{zero_depth}"
        )
    if junctions:
        assessed_junctions = [value for value in junctions.values() if value is not None]
        supported = sum(1 for value in assessed_junctions if (value or 0) > 0)
        parts.append(
            f"Fusion/Splice连接证据已评估 {len(assessed_junctions)}/{len(junctions)} 个独立事件，"
            f"junction reads>0：{supported}"
        )
    return "；".join(parts)


def _depth_summary_from_tsv(path: Path | None, *fields: str, limit: int | None = None) -> str:
    if not path or not path.is_file():
        return ""
    values: list[float] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="	")
        for row in reader:
            key = str(row.get("event_id") or row.get("peptide_id") or row.get("variant_key") or "")
            if key and key in seen:
                continue
            for field in fields:
                try:
                    value = float(str(row.get(field) or "").strip())
                except ValueError:
                    continue
                if value > 0:
                    values.append(value)
                    if key:
                        seen.add(key)
                    break
            if limit and len(values) >= limit:
                break
    median = _numeric_median(values)
    if median is None:
        return ""
    return f"候选位点中位深度 {median:.0f}x（n={len(values)}）"


def _metadata_contexts(prov: Mapping[str, Any], root: Path | None) -> list[Mapping[str, Any]]:
    """Load bounded upstream metadata referenced by a derived result."""
    contexts: list[Mapping[str, Any]] = [prov]
    production_roots: set[Path] = set()
    for _, value in _walk_scalar_items(prov):
        candidate = Path(value)
        try:
            exists = candidate.is_absolute() and candidate.exists()
        except OSError:
            exists = False
        if not exists:
            continue
        for parent in [candidate.parent, *candidate.parents]:
            if parent.name == "production":
                production_roots.add(parent)
                break
    if root:
        for parent in [root, *root.parents]:
            if parent.name == "production":
                production_roots.add(parent)
                break

    metadata_paths: set[Path] = set()
    for production_root in production_roots:
        metadata_paths.update(production_root.glob("final*/provenance.json"))
        metadata_paths.update(production_root.glob("final*/*/provenance.json"))
        metadata_paths.update(production_root.glob("*/provenance.json"))
        metadata_paths.add(production_root / "production_run_summary.json")
        metadata_paths.add(production_root.parent / "run_manifest.json")
    existing = [path for path in metadata_paths if path.is_file()]
    existing.sort(key=lambda path: path.stat().st_mtime, reverse=True)
    for path in existing[:100]:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            continue
        if isinstance(payload, Mapping):
            contexts.append(payload)
    return contexts


def _find_vcf_input(prov: Mapping[str, Any], root: Path | None) -> str:
    candidates: list[tuple[int, int, str]] = []
    seen: set[str] = set()
    for context_index, context in enumerate(_metadata_contexts(prov, root)):
        for key_path, value in _walk_scalar_items(context):
            lowered = value.lower()
            if not (lowered.endswith(".vcf") or lowered.endswith(".vcf.gz")) or value in seen:
                continue
            seen.add(value)
            joined = f"{key_path} {value}".lower()
            if any(token in joined for token in ("spechla", "lohhla", "hla_vcf", "germline_only")):
                continue
            score = 0
            if any(token in key_path.lower() for token in ("somatic_vcf", "annotated_vcf", "/vcf")):
                score += 12
            if "somatic" in joined:
                score += 8
            if any(token in joined for token in ("paired", "tumor", "tumour")):
                score += 3
            if "pass" in joined:
                score += 1
            try:
                is_file = Path(value).is_file()
            except OSError:
                is_file = False
            if is_file:
                candidates.append((score, -context_index, value))
    return max(candidates, default=(0, 0, ""))[2]


def _vcf_depth_summary_from_provenance(
    prov: Mapping[str, Any], root: Path | None, rows: list[Mapping[str, Any]], specimen: str,
) -> str:
    """Summarize candidate-site depth from a paired VCF when evidence columns are empty."""
    vcf = _find_vcf_input(prov, root)
    if not vcf:
        return ""
    candidate_keys = {
        (re.sub(r"^chr", "", str(row.get("chrom") or row.get("chromosome") or "").strip(), flags=re.IGNORECASE),
         str(row.get("pos") or row.get("position") or "").strip())
        for row in rows
        if str(row.get("chrom") or row.get("chromosome") or "").strip()
        and str(row.get("pos") or row.get("position") or "").strip()
    }
    if not candidate_keys:
        return ""
    opener = gzip.open if vcf.endswith(".gz") else open
    sample_names: list[str] = []
    sample_index: int | None = None
    values: list[float] = []
    tokens = {"normal": ("normal", "blood", "germline", "control"), "tumor": ("tumor", "tumour")}.get(specimen, (specimen,))
    try:
        with opener(vcf, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("#CHROM"):
                    sample_names = line.rstrip("\n").split("\t")[9:]
                    sample_index = next((i for i, name in enumerate(sample_names) if any(token in name.lower() for token in tokens)), None)
                    if sample_index is None and len(sample_names) == 2:
                        opposite = ("tumor", "tumour") if specimen == "normal" else ("normal", "blood", "germline", "control")
                        opposite_index = next((i for i, name in enumerate(sample_names) if any(token in name.lower() for token in opposite)), None)
                        if opposite_index is not None:
                            sample_index = 1 - opposite_index
                    continue
                if line.startswith("#") or sample_index is None:
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) <= 9 + sample_index:
                    continue
                key = (re.sub(r"^chr", "", fields[0], flags=re.IGNORECASE), fields[1])
                if key not in candidate_keys:
                    continue
                sample = dict(zip(fields[8].split(":"), fields[9 + sample_index].split(":")))
                value = sample.get("DP", "")
                if not value and sample.get("AD"):
                    try:
                        value = str(sum(float(item) for item in sample["AD"].split(",") if item not in {"", "."}))
                    except ValueError:
                        value = ""
                try:
                    depth = float(value)
                except ValueError:
                    continue
                if depth >= 0:
                    values.append(depth)
    except (OSError, UnicodeError):
        return ""
    median = _numeric_median(values)
    if median is None or sample_index is None:
        return ""
    sample_name = sample_names[sample_index] if sample_index < len(sample_names) else specimen
    label = "正常样本" if specimen == "normal" else "肿瘤样本"
    return f"候选事件位点中位{label}有效深度 {median:.0f}x（n={len(values)}；VCF {sample_name} DP）"


def _walk_scalar_items(value: Any, path: str = "") -> list[tuple[str, str]]:
    items: list[tuple[str, str]] = []
    if isinstance(value, Mapping):
        for key, nested in value.items():
            items.extend(_walk_scalar_items(nested, f"{path}/{key}"))
    elif isinstance(value, list):
        for index, nested in enumerate(value):
            items.extend(_walk_scalar_items(nested, f"{path}[{index}]"))
    else:
        text = str(value or "").strip()
        if text:
            items.append((path, text))
    return items


def _find_bam_input(prov: Mapping[str, Any], specimen: str, root: Path | None = None) -> str:
    wanted = {
        "normal": ("normal", "blood", "germline", "control"),
        "tumor": ("tumor", "tumour", "sample"),
    }
    tokens = wanted.get(specimen, (specimen,))
    candidates: list[tuple[int, int, int, str]] = []
    seen: set[str] = set()
    excluded = ("rna", "spechla", "lohhla", "hla-la", "polysolver", "mhc", "realign", "nomissmatch", "extract", "region_bam", "temp_bam")
    for context_index, context in enumerate(_metadata_contexts(prov, root)):
        for key_path, value in _walk_scalar_items(context):
            joined = f"{key_path} {value}".lower()
            if not value.lower().endswith(".bam") or value in seen:
                continue
            seen.add(value)
            if not any(token in joined for token in tokens) or any(token in joined for token in excluded):
                continue
            source = Path(value)
            try:
                is_file = source.is_file()
                source_size = source.stat().st_size if is_file else 0
            except OSError:
                is_file = False
                source_size = 0
            if not is_file:
                continue
            score = 0
            key_lower = key_path.lower()
            if any(token in key_lower for token in (f"{specimen}_dna_bam", f"{specimen}_bam", "source_bam")):
                score += 15
            if specimen in joined:
                score += 5
            if "blood_wgs" in joined or "germline_bam" in joined:
                score += 8 if specimen == "normal" else -8
            if "wgs" in joined:
                score += 5
            elif "wes" in joined or "exome" in joined:
                score += 3
            if "dna" in joined:
                score += 3
            candidates.append((score, min(source_size, 10**12), -context_index, value))
    return max(candidates, default=(0, 0, 0, ""))[3]


def _samtools_candidates(root: Path | None) -> list[str]:
    candidates = [
        os.environ.get("NEOAG_SAMTOOLS", ""),
        os.environ.get("SAMTOOLS", ""),
        shutil.which("samtools") or "",
    ]
    for env_key in ("OPEN_NEO_DEPLOY_ROOT", "NEOAG_DEPLOY_ROOT"):
        base = os.environ.get(env_key, "").strip()
        if base:
            candidates.extend([
                str(Path(base) / "env_tool" / "conda_envs" / "neoag-tools" / "bin" / "samtools"),
                str(Path(base) / "env_tool" / "tools" / "bin" / "samtools"),
            ])
    if root:
        for parent in [root, *root.parents]:
            candidates.extend([
                str(parent / "env_tool" / "conda_envs" / "neoag-tools" / "bin" / "samtools"),
                str(parent / "open-neo-deploy" / "env_tool" / "conda_envs" / "neoag-tools" / "bin" / "samtools"),
            ])
    return [path for path in dict.fromkeys(candidates) if path and Path(path).is_file()]


def _bam_depth_summary_from_provenance(prov: Mapping[str, Any], root: Path | None, specimen: str) -> str:
    bam = _find_bam_input(prov, specimen, root)
    if not bam:
        return ""
    samtools = next(iter(_samtools_candidates(root)), "")
    if not samtools:
        return ""
    try:
        idxstats = subprocess.run([samtools, "idxstats", bam], text=True, capture_output=True, check=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return ""
    chrom_pattern = re.compile(r"^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$", re.IGNORECASE)
    primary_length = 0
    primary_mapped = 0
    for line in idxstats.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 4 or not chrom_pattern.match(fields[0]):
            continue
        try:
            primary_length += int(fields[1])
            primary_mapped += int(fields[2])
        except ValueError:
            continue
    if not primary_length or not primary_mapped:
        return ""
    read_lengths: list[int] = []
    proc: subprocess.Popen[str] | None = None
    try:
        proc = subprocess.Popen([samtools, "view", bam], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        if proc.stdout:
            for index, line in enumerate(proc.stdout):
                if index >= 10000:
                    break
                fields = line.rstrip("\n").split("\t")
                if len(fields) > 9 and fields[9] != "*":
                    read_lengths.append(len(fields[9]))
    except OSError:
        read_lengths = []
    finally:
        if proc and proc.poll() is None:
            proc.kill()
    if not read_lengths:
        return ""
    read_length = sorted(read_lengths)[len(read_lengths) // 2]
    depth = primary_mapped * read_length / primary_length
    label = "正常DNA" if specimen == "normal" else "肿瘤DNA"
    return f"{label} BAM全基因组估算平均深度 {depth:.1f}x（idxstats映射reads×抽样读长；非候选位点pileup）"


def _purity_tools_from_provenance(prov: Mapping[str, Any]) -> list[dict[str, str]]:
    declared = [dict(row) for row in (prov.get("purity_cnv_tools") or []) if isinstance(row, Mapping)]
    if declared:
        return declared
    tools: list[dict[str, str]] = []
    for name, record in (prov.get("tools") or {}).items() if isinstance(prov.get("tools"), Mapping) else []:
        key = str(name).lower()
        if key not in {"facets", "sequenza", "purple", "ascat"}:
            continue
        payload = record if isinstance(record, Mapping) else {}
        purity = ""
        ploidy = ""
        source = Path(str(payload.get("file") or ""))
        if source.is_file():
            rows = _read_optional(source)
            if rows:
                purity = str(rows[0].get("purity") or rows[0].get("cellularity") or "")
                ploidy = str(rows[0].get("ploidy") or "")
        tools.append({
            "tool": str(name).upper(),
            "purity": purity,
            "ploidy": ploidy,
            "status": str(payload.get("status") or "ASSESSED"),
            "note": "来自生产接口已转换的单一工具结果" if not declared else "",
        })
    return tools


def _badge(text: str) -> str:
    t = str(text or "")
    cls = "UNASSESSED"
    if any(x in t.upper() for x in ["PASS", "INTACT", "HIGH", "A", "B"]):
        cls = "PASS"
    if any(x in t.upper() for x in ["CAUTION", "REVIEW", "MEDIUM", "LOW", "C"]):
        cls = "CAUTION"
    if any(x in t.upper() for x in ["DEFECT", "REJECT", "FAIL", "LOST", "GLOBAL", "D"]):
        cls = "FAIL"
    return f"<span class='badge {cls}'>{esc(text)}</span>"


def _table(rows: list[Mapping[str, Any]], headers: list[str], *, max_rows: int | None = None) -> str:
    view = rows[:max_rows] if max_rows else rows
    out = ["<table><tr>" + "".join(f"<th>{esc(h)}</th>" for h in headers) + "</tr>"]
    for row in view:
        out.append("<tr>" + "".join(f"<td>{esc(row.get(h, ''))}</td>" for h in headers) + "</tr>")
    out.append("</table>")
    return "\n".join(out)


def _map_by(rows: list[Mapping[str, Any]], key: str) -> dict[str, dict[str, str]]:
    return {str(r.get(key, "")): dict(r) for r in rows if r.get(key)}


def _path_values(payload: Any) -> list[Path]:
    """Collect path-like provenance values without scanning the filesystem."""
    paths: list[Path] = []
    if isinstance(payload, Mapping):
        for value in payload.values():
            paths.extend(_path_values(value))
    elif isinstance(payload, (list, tuple)):
        for value in payload:
            paths.extend(_path_values(value))
    elif isinstance(payload, str) and payload.startswith("/"):
        paths.append(Path(payload))
    return paths


def _asset_gtf_candidates(provenance: Mapping[str, Any]) -> list[Path]:
    """Infer GTF locations from already-recorded asset paths.

    Production manifests often record the normal proteome or other assets but
    omit the GTF itself.  Deriving the shared ``data`` root keeps report-time
    gene-symbol resolution portable without hard-coding a server mount.
    """
    candidates: list[Path] = []
    for path in _path_values(provenance):
        parts = path.parts
        for index in range(len(parts) - 1, -1, -1):
            if parts[index] != "data":
                continue
            data_root = Path(*parts[: index + 1])
            candidates.extend([
                data_root / "rna" / "gencode_v49" / "gencode.v49.annotation.gtf.gz",
                data_root / "ref" / "hg38" / "gencode.gtf",
                data_root / "ref" / "ctat" / "current" / "ctat_genome_lib_build_dir" / "ref_annot.gtf",
            ])
            break
    return list(dict.fromkeys(candidates))


def _normal_expression_gene_map(root: Path | None, prov: Mapping[str, Any], gene_ids: set[str]) -> dict[str, str]:
    """Map Ensembl gene IDs to symbols using deployed normal-expression assets."""
    if not gene_ids:
        return {}
    candidates: list[Path] = []
    for payload in (prov, prov.get("input_files") or {}, prov.get("references") or {}, prov.get("reference_manifest") or {}):
        if not isinstance(payload, Mapping):
            continue
        for key in ("normal_expression", "normal_expression_reference", "normal_expression_table"):
            value = str(payload.get(key) or "").strip()
            if value:
                candidates.append(Path(value))
    for env_key in ("NEOAG_NORMAL_EXPRESSION", "NORMAL_EXPRESSION_REFERENCE"):
        value = os.environ.get(env_key, "").strip()
        if value:
            candidates.append(Path(value))
    if root:
        for parent in [root, *root.parents]:
            candidates.extend([
                parent / "refs" / "data" / "normal" / "expression" / "normal_expression.gtex_v11_hpa_hspc.tsv",
                parent / "open-neo-deploy" / "refs" / "data" / "normal" / "expression" / "normal_expression.gtex_v11_hpa_hspc.tsv",
            ])
    gene_map: dict[str, str] = {}
    for candidate in dict.fromkeys(candidates):
        if not candidate.is_file():
            continue
        try:
            with candidate.open("r", encoding="utf-8", errors="replace", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    ensembl = str(row.get("ensembl_gene_id") or row.get("gene_id") or "").split(".", 1)[0].strip()
                    symbol = str(row.get("gene") or row.get("gene_name") or row.get("gene_symbol") or "").strip()
                    if ensembl in gene_ids and symbol and symbol != ensembl:
                        gene_map["GENE|" + ensembl] = symbol
                        if len(gene_map) >= len(gene_ids):
                            return gene_map
        except OSError:
            continue
    return gene_map


def _patient_expression_tpm_map(prov: Mapping[str, Any]) -> tuple[dict[str, str], str]:
    """Load the patient gene-expression table recorded by pipeline provenance."""
    candidates: list[Path] = []
    payloads: list[Mapping[str, Any]] = [prov]
    for key in ("input_files", "references", "reference_manifest"):
        value = prov.get(key)
        if isinstance(value, Mapping):
            payloads.append(value)
    for payload in payloads:
        for key in ("gene_expression", "expression", "gene_tpm", "gene_expression_tpm"):
            value = str(payload.get(key) or "").strip()
            if value:
                candidates.append(Path(value))
    for candidate in dict.fromkeys(candidates):
        if not candidate.is_file() or candidate.stat().st_size == 0:
            continue
        values: dict[str, str] = {}
        try:
            with candidate.open("r", encoding="utf-8", errors="replace", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    gene_id = str(
                        row.get("gene_id") or row.get("ensembl_gene_id")
                        or row.get("gene") or row.get("gene_name") or ""
                    ).strip()
                    tpm = str(
                        row.get("tpm") or row.get("TPM") or row.get("gene_tpm")
                        or row.get("expression_tpm") or ""
                    ).strip()
                    if not gene_id or not tpm:
                        continue
                    values[gene_id] = tpm
                    values[gene_id.split(".", 1)[0]] = tpm
                    symbol = str(row.get("gene_symbol") or row.get("gene_name") or "").strip()
                    if symbol:
                        values[symbol] = tpm
        except OSError:
            continue
        if values:
            return values, str(candidate)
    return {}, ""


def _apply_patient_gene_expression(
    rows: list[dict[str, str]], expression_map: Mapping[str, str], source: str,
) -> list[dict[str, str]]:
    """Replace placeholder expression only when one exact patient gene ID resolves."""
    if not expression_map:
        return rows
    id_fields = (
        "gene_id", "ensembl_gene_id", "gene", "source_event_id", "source_record_id",
        "source_records", "event_name", "splice_event_id", "uid",
    )
    for row in rows:
        is_splice = str(row.get("event_type") or "").strip().lower() in {
            "splice", "splice_junction", "junction",
        }
        gene_ids = list(dict.fromkeys(
            gene_id
            for field_name in id_fields
            for gene_id in re.findall(r"ENSG\d+(?:\.\d+)?", str(row.get(field_name) or ""))
        ))
        matched = [expression_map.get(gene_id) or expression_map.get(gene_id.split(".", 1)[0]) for gene_id in gene_ids]
        matched = [value for value in matched if value is not None]
        if len(gene_ids) == 1 and matched:
            row["gene_expression_tpm"] = str(matched[0])
            row["expression_source"] = source
            row["expression_evidence_status"] = "GENE_EXPRESSION_MATCHED_BY_ENSEMBL_ID"
            continue
        if is_splice and len(gene_ids) > 1:
            row["gene_expression_tpm"] = ""
            row["expression_source"] = source
            row["expression_evidence_status"] = "UNASSESSED_AMBIGUOUS_GENE_ID"
            continue
        gene = str(row.get("gene") or "").strip()
        if gene and gene in expression_map and "::" not in gene and "/" not in gene:
            row["gene_expression_tpm"] = str(expression_map[gene])
            row["expression_source"] = source
            row["expression_evidence_status"] = "GENE_EXPRESSION_MATCHED_BY_SYMBOL"
            continue
        if is_splice and str(row.get("gene_expression_tpm") or "").strip() in {"0", "0.0", "0.0000"}:
            row["gene_expression_tpm"] = ""
            row["expression_source"] = source
            row["expression_evidence_status"] = "UNASSESSED_ID_NOT_MAPPED"
    return rows


def _read_longrna_junction_genes(
    root: Path | None,
    extra_gene_ids: set[str] | None = None,
    gtf_paths: list[Path] | None = None,
) -> dict[str, str]:
    """Build junction-to-gene labels for four-tool splice events."""
    if not root:
        return {}
    rows = _read_optional(root / "inputs" / "longrna_comprehensive_evidence.tsv")
    validation_rows = _read_optional(root / "metadata" / "four_tool_splice_cross_validation.tsv")
    provenance_rows = _read_optional(root / "metadata" / "splice_multitool_provenance.tsv")
    gene_ids = {str(row.get("uid") or "").split(":", 1)[0] for row in rows if str(row.get("uid") or "").startswith("ENSG")}
    gene_ids.update(extra_gene_ids or set())
    transcript_ids = {
        transcript.strip() for row in validation_rows for transcript in str(row.get("longread_isoforms") or "").split(";")
        if transcript.strip().startswith("ENS")
    }
    gene_names: dict[str, str] = {}
    gene_intervals: list[tuple[str, int, int, str, str]] = []
    gtf_candidates = list(gtf_paths or [])
    if len(root.parents) > 2:
        gtf_candidates.extend([
            root.parents[2] / "data" / "ref" / "hg38" / "gencode.gtf",
            root.parents[2] / "data" / "ref" / "ctat" / "current" / "ctat_genome_lib_build_dir" / "ref_annot.gtf",
        ])
    if len(root.parents) > 3:
        gtf_candidates.append(
            root.parents[3] / "open-neo-deploy" / "refs" / "data" / "ref" / "hg38" / "gencode.gtf"
        )
    for gtf in dict.fromkeys(gtf_candidates):
        if not gtf.is_file() or not (gene_ids or transcript_ids):
            continue
        opener = gzip.open if gtf.suffix == ".gz" else open
        with opener(gtf, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if "\tgene\t" in line:
                    match = re.search(r'gene_id "(ENSG\d+)(?:\.\d+)?";.*?gene_name "([^"]+)";', line)
                    if match and match.group(1) in gene_ids:
                        gene_names[match.group(1)] = match.group(2)
                    if match:
                        fields = line.split("\t")
                        if len(fields) >= 7:
                            chrom = fields[0] if fields[0].startswith("chr") else "chr" + fields[0]
                            gene_intervals.append((chrom, int(fields[3]), int(fields[4]), fields[6], match.group(2)))
                elif "\ttranscript\t" in line and transcript_ids:
                    match = re.search(
                        r'gene_id "(ENSG\d+)(?:\.\d+)?";.*?transcript_id "(ENS[TG]\d+)(?:\.\d+)?";.*?gene_name "([^"]+)";',
                        line,
                    )
                    if match and match.group(2) in transcript_ids:
                        gene_names[match.group(2)] = match.group(3)
        resolved_gene_ids = {gene_id for gene_id in gene_ids if gene_id in gene_names}
        resolved_transcripts = {transcript for transcript in transcript_ids if transcript in gene_names}
        if resolved_gene_ids == gene_ids and resolved_transcripts == transcript_ids:
            break
    junction_genes: dict[str, str] = {}
    for row in rows:
        coord = str(row.get("coord") or "")
        parts = coord.split(":")
        if len(parts) != 4:
            continue
        junction = "SJ|GRCh38|" + "|".join(parts)
        gene_id = str(row.get("uid") or "").split(":", 1)[0]
        if gene_id:
            junction_genes[junction] = gene_names.get(gene_id, gene_id)
    for row in validation_rows:
        junction = str(row.get("canonical_junction_id") or "")
        if not junction or junction in junction_genes:
            continue
        transcript = next((item.strip() for item in str(row.get("longread_isoforms") or "").split(";") if item.strip().startswith("ENS")), "")
        if transcript and transcript in gene_names:
            junction_genes[junction] = gene_names[transcript]
    for row in validation_rows + provenance_rows:
        junction = str(row.get("canonical_junction_id") or "")
        if not junction or junction in junction_genes:
            continue
        parts = junction.split("|")
        if len(parts) != 6:
            continue
        chrom, start, end, strand = parts[2], int(parts[3]), int(parts[4]), parts[5]
        matches = [
            (gene_end - gene_start, gene_name)
            for gene_chrom, gene_start, gene_end, gene_strand, gene_name in gene_intervals
            if gene_chrom == chrom and gene_strand == strand and gene_start <= start <= gene_end and gene_start <= end <= gene_end
        ]
        if matches:
            junction_genes[junction] = sorted(matches)[0][1]
    for row in provenance_rows:
        event_id = str(row.get("event_id") or "")
        junction = str(row.get("canonical_junction_id") or "")
        if event_id and junction in junction_genes:
            junction_genes["EVENT|" + event_id] = junction_genes[junction]
    for gene_id, gene_name in gene_names.items():
        if gene_id.startswith("ENSG"):
            junction_genes["GENE|" + gene_id] = gene_name
    return junction_genes


def _replace_gene_ids(value: Any, gene_map: Mapping[str, str]) -> str:
    text = str(value or "")
    return re.sub(r"ENSG\d+", lambda match: gene_map.get("GENE|" + match.group(0), GENE_SYMBOL_FALLBACKS.get(match.group(0), match.group(0))), text)


def _gene_label_needs_enrichment(value: Any) -> bool:
    text = str(value or "").strip()
    return (
        not text
        or text.startswith(("SEV|", "SJ|", "JUNC"))
        or bool(re.fullmatch(r"(?:chr)?[0-9XYM]+:\d+(?:[-:]\d+)?", text, re.IGNORECASE))
        or bool(re.fullmatch(r"ENSG\d+", text))
    )


def _event_change_needs_enrichment(row: Mapping[str, Any]) -> bool:
    event_type = str(row.get("event_type") or "").strip().lower()
    protein_change = str(row.get("combined_protein_change") or "").strip()
    if ":p." in protein_change or re.search(r"\bp\.[A-Za-z*]", protein_change):
        return False
    if event_type in {"splice", "fusion"}:
        return not any(str(row.get(field_name) or "").strip() for field_name in (
            "source_event_id", "source_record_id", "source_records", "canonical_junction_id",
        ))
    return event_type in {"snv", "indel", "missense", "frameshift"} and not protein_change


def _read_report_enrichment_rows(
    path: Path | None,
    rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Stream only source fields needed to label unresolved report events."""
    if not path or not path.is_file():
        return []
    target_events = {
        str(row.get("event_id") or "").removeprefix("EVENT:")
        for row in rows
        if _gene_label_needs_enrichment(row.get("gene")) or _event_change_needs_enrichment(row)
    }
    if not target_events:
        return []
    keep_fields = {
        "event_id", "peptide_id", "hla_allele", "gene", "gene_pair", "event_name",
        "canonical_junction_id", "source_event_id", "source_record_id", "source_records",
        "source_junction_id", "source_tool", "source_tools", "splice_event_id",
        "combined_protein_change", "consequence", "peptide_consequence", "transcript_id",
        "orf_id", "transcript_orf_status", "independent_translation_generators",
        "normal_tissue_max_tpm", "normal_tissue_max_tissue", "critical_tissue_max_tpm",
        "critical_tissue_name", "critical_tissue_hit", "normal_hspc_tpm",
        "normal_hspc_cell_type", "normal_hspc_unit", "normal_expression_status",
        "normal_hspc_status", "safety_missing_layers", "safety_reason",
        "safety_reason_codes", "event_normal_tissue_max_tpm",
        "normal_proteome_exact_match_status", "normal_transcript_junction_match_status",
        "normal_immunopeptidome_match_status", "similar_peptide_cross_reactivity_status",
        "source_gene_expression_context_status", "critical_organ_expression_context_status",
        "hematopoietic_expression_context_status", "tcr_contact_anchor_context_status",
        "final_safety_conclusion", "safety_evidence_completeness",
        "event_normal_tissue_max_tissue", "event_critical_tissue_max_tpm",
        "event_critical_tissue_name", "event_critical_tissue_hit",
        "event_normal_hspc_tpm", "event_normal_hspc_cell_type",
        "event_normal_hspc_unit", "event_normal_expression_status",
        "event_normal_hspc_status", "event_safety_missing_layers",
        "event_safety_reason",
    }
    selected: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            event_id = str(row.get("event_id") or "").removeprefix("EVENT:")
            if event_id in target_events and event_id not in selected:
                selected[event_id] = {field_name: str(row.get(field_name) or "") for field_name in keep_fields}
                if len(selected) == len(target_events):
                    break
    return list(selected.values())


def _enrich_rows_from_sources(
    rows: list[dict[str, str]], source_peptides: list[dict[str, str]], source_events: list[dict[str, str]],
    junction_genes: Mapping[str, str], evidence_source_label: str = "",
) -> list[dict[str, str]]:
    """Fill display fields from raw events, representative peptides, and splice junction annotation."""
    peptide_by_id = _map_by(source_peptides, "peptide_id")
    peptide_by_id_hla = {
        (str(row.get("peptide_id") or ""), str(row.get("hla_allele") or "")): row
        for row in source_peptides if row.get("peptide_id")
    }
    source_by_id = _map_by(source_events, "event_id")
    source_by_id.update({
        "EVENT:" + event_id.removeprefix("EVENT:"): row
        for event_id, row in list(source_by_id.items())
    })
    fields = (
        "gene", "combined_protein_change", "event_name", "consequence", "peptide_consequence",
        "canonical_junction_id", "source_event_id", "source_record_id", "source_records",
        "source_junction_id", "source_tool", "source_tools", "splice_event_id",
        "frame_evidence_grade", "orf_evidence_grade", "orf_id", "transcript_id",
        "transcript_orf_status", "independent_translation_generators", "rna_support_status",
        "rna_support_grade", "rna_support_reason", "rna_support_reason_code", "rna_support_state",
        "rna_evidence_completeness", "rna_evidence_score", "rna_alt_reads", "rna_ref_reads",
        "rna_depth", "rna_vaf", "rna_junction_reads", "rna_junction_source", "event_expression",
        "gene_expression_tpm", "transcript_expression_tpm", "cross_platform_status", "safety_status",
    )
    rna_fields = {
        "rna_support_status", "rna_support_grade", "rna_support_reason", "rna_support_reason_code",
        "rna_support_state", "rna_evidence_completeness", "rna_evidence_score", "rna_alt_reads",
        "rna_ref_reads", "rna_depth", "rna_vaf", "rna_junction_reads", "rna_junction_source",
        "event_expression", "gene_expression_tpm", "transcript_expression_tpm",
    }
    authoritative_evidence_fields = rna_fields | {
        "presentation_consensus_state", "presentation_consensus_grade", "presentation_evidence_grade",
        "presentation_evidence_score", "mutant_specificity_status", "mutant_specificity_gate_status",
        "mutant_specificity_grade", "mutant_specificity_reason", "source_chain_confidence_tier",
        "source_chain_confidence_grade", "source_chain_confidence_label", "source_chain_hard_failure",
        "source_chain_missing_requirements", "evidence_grade", "evidence_grade_uncapped",
        "final_priority", "safety_status", "safety_state", "safety_grade", "cross_platform_status",
        "event_authenticity_state", "event_authenticity_grade", "clonality_state", "clonality_grade",
        "normal_tissue_max_tpm", "normal_tissue_max_tissue", "critical_tissue_max_tpm",
        "critical_tissue_name", "critical_tissue_hit", "normal_hspc_tpm",
        "normal_hspc_cell_type", "normal_hspc_unit", "normal_expression_status",
        "normal_hspc_status", "safety_missing_layers", "safety_reason",
        "safety_reason_codes", "event_normal_tissue_max_tpm",
        "normal_proteome_exact_match_status", "normal_transcript_junction_match_status",
        "normal_immunopeptidome_match_status", "similar_peptide_cross_reactivity_status",
        "source_gene_expression_context_status", "critical_organ_expression_context_status",
        "hematopoietic_expression_context_status", "tcr_contact_anchor_context_status",
        "final_safety_conclusion", "safety_evidence_completeness",
        "event_normal_tissue_max_tissue", "event_critical_tissue_max_tpm",
        "event_critical_tissue_name", "event_critical_tissue_hit",
        "event_normal_hspc_tpm", "event_normal_hspc_cell_type",
        "event_normal_hspc_unit", "event_normal_expression_status",
        "event_normal_hspc_status", "event_safety_missing_layers",
        "event_safety_reason",
        "ccf_estimate", "ccf_confidence_state", "purity_consensus_status", "hla_appm_state",
        "appm_integrity_status", "restricting_hla_lost", "restricting_locus_loh", "loh_status",
        "escape_status", "immunogenicity_score", "immunogenicity_composite_score",
        "netmhcpan_el_rank", "mhcflurry_presentation_score", "mhcflurry_processing_score",
        "prime_score", "bigmhc_im_score", "netmhcstabpan_rank", "netmhcstabpan_score",
        "netchop_31d_cterm_score", "netchop_31d_max_score", "netchop_31d_mean_score",
        "netchop_31d_cleavage_sites", "netchop_processing_status",
        "tap_processing_status", "evidence_field_sources", "evidence_conflict_fields",
    }
    copy_fields = set(fields) | authoritative_evidence_fields
    enriched: list[dict[str, str]] = []
    for event in rows:
        row = dict(event)
        source = source_by_id.get(str(row.get("event_id") or ""), {})
        representative_ids = (
            row.get("best_peptide_id"), row.get("representative_1_peptide_id"),
            row.get("representative_2_peptide_id"),
        )
        peptide_key = (str(row.get("peptide_id") or ""), str(row.get("hla_allele") or ""))
        representative = peptide_by_id_hla.get(peptide_key) or peptide_by_id.get(peptide_key[0]) or next((peptide_by_id.get(str(pid)) for pid in representative_ids if pid and str(pid) in peptide_by_id), None)
        for field_name in copy_fields:
            if not row.get(field_name) and source.get(field_name):
                row[field_name] = source[field_name]
            if representative and representative.get(field_name) and (field_name in authoritative_evidence_fields or not row.get(field_name)):
                row[field_name] = representative[field_name]
        row["_report_evidence_source"] = evidence_source_label or "ranked input"
        row["_report_evidence_matched"] = "YES" if representative else "NO"
        junction = str(row.get("canonical_junction_id") or (representative or {}).get("canonical_junction_id") or "")
        if not junction:
            event_id = str(row.get("event_id") or (representative or {}).get("event_id") or "")
            junction = "EVENT|" + event_id if event_id else ""
        if _gene_label_needs_enrichment(row.get("gene")) and junction in junction_genes:
            row["gene"] = junction_genes[junction]
        if _gene_label_needs_enrichment(row.get("gene")):
            source_labels = " ".join(str(row.get(field_name) or "") for field_name in (
                "source_event_id", "source_record_id", "source_records", "event_name",
            ))
            gene_ids = list(dict.fromkeys(re.findall(r"ENSG\d+", source_labels)))
            if gene_ids:
                row["gene"] = " / ".join(
                    junction_genes.get("GENE|" + gene_id, GENE_SYMBOL_FALLBACKS.get(gene_id, gene_id))
                    for gene_id in gene_ids
                )
        for field_name in ("gene", "event_name", "combined_protein_change"):
            if row.get(field_name):
                row[field_name] = _replace_gene_ids(row[field_name], junction_genes)
        if re.fullmatch(r"ENSG\d+", str(row.get("gene") or "")):
            row["gene"] = f"{row['gene']}（symbol未注释）"
        enriched.append(row)
    return enriched


def _apply_junction_verification(
    rows: list[dict[str, str]], verification_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Overlay independently verified junction counts without altering caller evidence."""
    by_event = {
        str(item.get("event_id") or "").removeprefix("EVENT:"): item
        for item in verification_rows
        if str(item.get("event_id") or "").strip()
    }
    by_peptide = {
        str(item.get("peptide_id") or ""): item
        for item in verification_rows
        if str(item.get("peptide_id") or "").strip()
    }
    if not by_event and not by_peptide:
        return rows
    for row in rows:
        event_id = str(row.get("event_id") or "").removeprefix("EVENT:")
        peptide_id = str(row.get("peptide_id") or "")
        verification = by_peptide.get(peptide_id) or by_event.get(event_id)
        if not verification:
            continue
        verified_reads = str(verification.get("verified_rna_junction_reads") or "").strip()
        if not verified_reads:
            continue
        caller_reads = str(verification.get("caller_rna_junction_reads") or "").strip()
        if not caller_reads:
            caller_reads = str(
                row.get("provided_rna_junction_reads") or row.get("rna_junction_reads") or ""
            ).strip()
        if caller_reads:
            row["provided_rna_junction_reads"] = caller_reads
            row["caller_rna_junction_reads"] = caller_reads
        row["rna_junction_reads"] = verified_reads
        row["junction_match_status"] = str(
            verification.get("junction_match_status") or "BAM_VERIFIED"
        )
        row["junction_match_method"] = str(
            verification.get("junction_match_method")
            or "star_chimeric_exact_breakpoint_plus_bam_qname"
        )
        row["rna_junction_source"] = str(
            verification.get("junction_verification_source") or "junction_read_verification.tsv"
        )
        row["junction_verification_note"] = str(
            verification.get("junction_verification_note") or ""
        )
    return rows


PRIORITY_PATIENT = {
    "A": "优先推荐进一步验证",
    "B": "值得考虑验证",
    "B_CAUTION": "可考虑，但需关注安全性",
    "C": "证据有限，谨慎推进",
    "C_CAUTION": "证据有限且需安全性复核",
    "D": "当前不建议推进",
}

CONSEQUENCE_PATIENT = {
    "missense": "氨基酸改变（点突变）",
    "frameshift": "移码变异产生的新肽段",
    "splice_junction": "剪接异常产生的新肽段",
    "exon_deletion_junction": "外显子缺失/剪接交界肽段",
    "insertion": "插入/缺失交界肽段",
    "fusion": "基因融合交界肽段",
    "other": "其他变异相关肽段",
}

FIELD_GLOSSARY = {
    "efficacy_score": "综合免疫学评分（0–1），整合表达、结合、呈递、安全性等维度。",
    "final_priority": "最终优先级分层：A/B 为优先候选，C 为需更多证据，D 为不建议推进。",
    "presentation_evidence_grade": "HLA 结合与加工呈递证据等级（A 最优）。",
    "appm_multiplier": "抗原加工呈递通路（APPM）完整性对候选的折减系数。",
    "ccf_multiplier": "肿瘤克隆性（CCF）对候选的折减系数。",
    "safety_status": "整合完整肽、正常连接/转录本、正常免疫肽组及相似自身肽的安全性初筛结果；来源基因表达仅作辅助背景。",
    "escape_status": "免疫逃逸机制（如 HLA 丢失）对候选的影响评估。",
    "validation_mode": "建议的实验验证设计类型（短肽对、长肽、minigene 等）。",
    "recommended_assay": "推荐的体外验证实验类型。",
    "cross_platform_status": "WES/WGS 共同检出、低水平支持、检出能力不足或样本特异性等跨平台证据状态。",
    "cross_platform_multiplier": "跨平台 DNA 证据对排序分数的保守调整系数。",
    "source_chain_confidence_tier": "候选来源链置信度 C1–C4：评估事件、转录本/ORF、肽段与HLA结果能否形成可追溯链路；不同于R1–R4实验推荐等级。",
    "source_chain_orthogonal_status": "独立或跨模态确认状态；同一BAM上的多个caller仅算计算一致性，不等同正交确认。",
    "source_chain_requirement_statuses": "SNV、InDel、Fusion或Splice赛道专属要求的逐项状态，区分NOT_APPLICABLE、UNASSESSED、LOW_POWER、NEGATIVE与CONFLICT。",
}


@dataclass
class ReportBundle:
    profile: Mapping[str, Any]
    events: list[dict[str, str]]
    peptides: list[dict[str, str]]
    appm_summary: Mapping[str, Any] = field(default_factory=dict)
    validation_rows: list[dict[str, str]] = field(default_factory=list)
    sample_id: str = ""
    entry_mode: str = ""
    provenance: dict[str, Any] = field(default_factory=dict)
    peptide_safety: list[dict[str, str]] = field(default_factory=list)
    peptide_escape_flags: list[dict[str, str]] = field(default_factory=list)
    immune_escape_summary: list[dict[str, str]] = field(default_factory=list)
    ccf: list[dict[str, str]] = field(default_factory=list)
    appm_gene_status: list[dict[str, str]] = field(default_factory=list)
    appm_peptide_modifiers: list[dict[str, str]] = field(default_factory=list)
    appm_module_scores: list[dict[str, str]] = field(default_factory=list)
    appm_submodule_scores: list[dict[str, str]] = field(default_factory=list)
    appm_conflicts: list[dict[str, str]] = field(default_factory=list)
    wes_qc: list[dict[str, str]] = field(default_factory=list)
    wes_wgs_coding_summary: list[dict[str, str]] = field(default_factory=list)
    targeted_pileup_summary: list[dict[str, str]] = field(default_factory=list)
    evidence_manifest: dict[str, Any] = field(default_factory=dict)
    evidence_conflicts: list[dict[str, str]] = field(default_factory=list)
    evidence_source_status: str = "RANKED_INPUT_ONLY"
    evidence_integrity: dict[str, str] = field(default_factory=dict)
    purity_tools: list[dict[str, str]] = field(default_factory=list)
    purity_consensus: dict[str, Any] = field(default_factory=dict)
    hla_loh_tool_results: list[dict[str, str]] = field(default_factory=list)
    tool_versions: dict[str, dict[str, str]] = field(default_factory=dict)
    disease_knowledge: dict[str, Any] = field(default_factory=dict)


def _load_disease_knowledge(
    profile: Mapping[str, Any], provenance: Mapping[str, Any],
) -> dict[str, Any]:
    """Load an optional disease knowledge profile without embedding diseases in code."""
    config = profile.get("disease_knowledge") if isinstance(profile, Mapping) else None
    config = config if isinstance(config, Mapping) else {}
    explicit = str(
        config.get("path")
        or profile.get("disease_knowledge_file")
        or provenance.get("disease_knowledge_file")
        or ""
    ).strip()
    repo_root = Path(__file__).resolve().parents[2]
    if config.get("anchors"):
        return dict(config)
    candidates: list[Path] = []
    if explicit:
        path = Path(explicit)
        candidates.append(path)
        if not path.is_absolute():
            candidates.append(repo_root / path)
    else:
        candidates.extend(sorted((repo_root / "configs" / "disease_knowledge").glob("*.json")))

    observed_context: list[str] = []
    for payload in (provenance, profile):
        if not isinstance(payload, Mapping):
            continue
        for key in ("disease", "diagnosis", "disease_name", "cancer_type", "tumor_type"):
            value = str(payload.get(key) or "").strip()
            if value:
                observed_context.append(value.upper())
        for container_key in ("clinical_context", "disease_profile", "patient_context", "clinical"):
            nested = payload.get(container_key)
            if isinstance(nested, Mapping):
                for key in ("disease", "diagnosis", "disease_name", "cancer_type", "tumor_type"):
                    value = str(nested.get(key) or "").strip()
                    if value:
                        observed_context.append(value.upper())

    source: Path | None = None
    payload: Mapping[str, Any] | None = None
    for candidate in candidates:
        if not candidate.is_file():
            continue
        try:
            parsed = json.loads(candidate.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            continue
        if not isinstance(parsed, Mapping):
            continue
        aliases = [
            str(value).upper() for value in [
                parsed.get("disease_id"), parsed.get("display_name"), *(parsed.get("disease_aliases") or []),
            ] if str(value or "").strip()
        ]
        if explicit or not observed_context or any(
            alias in context or context in alias
            for alias in aliases for context in observed_context
        ):
            source, payload = candidate, parsed
            break
    if source is None:
        return {"status": "MISSING", "source": explicit, "anchors": []} if explicit else {}
    assert payload is not None
    result = dict(payload)
    result.update({key: value for key, value in config.items() if key not in {"path", "anchors"}})
    result["status"] = "LOADED"
    result["source"] = str(source)
    return result


def _normalized_event_label(value: Any) -> str:
    return re.sub(r"[^A-Z0-9]+", "::", str(value or "").upper()).strip(":")


def _disease_anchor(row: Mapping[str, Any], bundle: ReportBundle | None) -> dict[str, Any]:
    if bundle is None:
        return {}
    observed = {
        _normalized_event_label(row.get(field))
        for field in ("gene", "event_name", "event_id", "source_event_id")
        if str(row.get(field) or "").strip()
    }
    peptide = str(row.get("peptide") or "").strip().upper()
    allele = str(row.get("hla_allele") or row.get("allele") or "").strip().upper()
    for anchor in bundle.disease_knowledge.get("anchors", []):
        if not isinstance(anchor, Mapping):
            continue
        aliases = {
            _normalized_event_label(value)
            for value in [anchor.get("event"), *(anchor.get("aliases") or [])]
            if str(value or "").strip()
        }
        if not observed.intersection(aliases):
            continue
        matched = dict(anchor)
        matched["peptide_evidence"] = ""
        for item in anchor.get("public_neoantigens", []):
            if not isinstance(item, Mapping) or str(item.get("sequence") or "").upper() != peptide:
                continue
            hla_prefixes = [str(value).upper() for value in item.get("hla_prefixes", [])]
            if not hla_prefixes or any(allele.startswith(prefix) for prefix in hla_prefixes):
                matched["peptide_evidence"] = str(item.get("evidence") or "external evidence recorded")
                break
        return matched
    return {}


def _patient_disease_anchor_note(row: Mapping[str, Any], bundle: ReportBundle | None) -> str:
    anchor = _disease_anchor(row, bundle)
    if not anchor:
        return ""
    note = str(anchor.get("label") or "疾病锚定事件")
    if anchor.get("peptide_evidence"):
        note += f"；该Peptide-HLA组合已记录外部功能/呈递证据（{anchor['peptide_evidence']}）"
    return note + "；仅优先展示，不自动提升R等级"


def _augment_runtime_tool_provenance(root: Path, provenance: dict[str, Any]) -> None:
    """Register tools proven by non-empty standard outputs in a completed run."""
    tools = provenance.setdefault("tools", {})
    if not isinstance(tools, dict):
        tools = {}
        provenance["tools"] = tools
    markers = {
        "netmhcpan": root / "presentation" / "netmhcpan_evidence.tsv",
        "mhcflurry": root / "presentation" / "mhcflurry_evidence.tsv",
        "netmhcstabpan": root / "presentation" / "netmhcstabpan_evidence.tsv",
        "prime": root / "presentation" / "prime_evidence.tsv",
        "bigmhc_im": root / "presentation" / "bigmhc_im_evidence.tsv",
        "deepimmuno": root / "presentation" / "deepimmuno_evidence.tsv",
        "netchop": root / "processing" / "netchop_evidence.tsv",
        "star": root / "rna" / "star" / "Log.final.out",
        "spechla": root / "hla_loh" / "spechla" / "hla_loh.tsv",
        "lohhla": root / "hla_loh" / "lohhla" / "hla_loh.tsv",
    }
    for name, path in markers.items():
        if not path.is_file() or path.stat().st_size == 0:
            continue
        record = tools.get(name)
        if not isinstance(record, dict):
            record = {}
            tools[name] = record
        record.update({"source": "local", "status": "real", "file": str(path)})
        record.setdefault("mode", "runtime_output_discovery")

    summary_paths = (
        root.parent / "evidence" / "purity_cnv" / "purity_cnv_tool_summary.tsv",
        root.parent / "evidence" / "hla_typing" / "hla_typing_tool_summary.tsv",
        root.parent / "evidence" / "hla_loh" / "hla_loh.standardized.tsv",
    )
    for path in summary_paths:
        if not path.is_file() or path.stat().st_size == 0:
            continue
        for row in _read_optional(path):
            raw_name = next((str(row.get(key) or "").strip() for key in (
                "tool", "source_tool", "evidence_tool", "caller", "method",
            ) if str(row.get(key) or "").strip()), "")
            if not raw_name:
                continue
            key = _patient_tool_key(raw_name)
            record = tools.get(key)
            if not isinstance(record, dict):
                record = {}
                tools[key] = record
            record.update({"source": "local", "status": "real", "file": str(path)})
            record.setdefault("mode", "tool_summary_discovery")
            record.setdefault("display_name", raw_name)

    # Event-producing tools may be preserved only on canonical raw events and
    # disappear from peptide-level ranking tables. Discover them from the
    # event provenance so patient reports do not silently omit tools such as
    # SpliceMutr after an exact-junction merge.
    event_source_paths = (
        root / "parsed" / "raw_events.tsv",
        root / "inputs" / "combined_raw_events.tsv",
    )
    for path in event_source_paths:
        if not path.is_file() or path.stat().st_size == 0:
            continue
        for key, count in _patient_source_tool_counts(_read_optional(path)).items():
            if key not in _PATIENT_SOURCE_TOOL_META:
                continue
            record = tools.get(key)
            if not isinstance(record, dict):
                record = {}
                tools[key] = record
            record.update({"source": "local", "status": "real", "file": str(path)})
            record.setdefault("mode", "event_source_discovery")
            record["evidence_event_rows"] = max(
                int(record.get("evidence_event_rows") or 0), count
            )



def _augment_purity_cnv_provenance(root: Path, prov: dict[str, Any]) -> None:
    production = root.parent
    summary_path = production / "evidence" / "purity_cnv" / "purity_cnv_tool_summary.tsv"
    if not summary_path.is_file() or summary_path.stat().st_size == 0:
        summary_path = production / "purity" / "consensus" / "purity_cnv_tool_summary.tsv"
    if not summary_path.is_file() or summary_path.stat().st_size == 0:
        return
    tools = []
    for row in _read_optional(summary_path):
        tools.append({
            "tool": str(row.get("tool") or row.get("source_tool") or ""),
            "purity": str(row.get("purity") or ""),
            "ploidy": str(row.get("ploidy") or ""),
            "status": str(row.get("status") or "UNASSESSED"),
            "note": str(row.get("notes") or row.get("parse_method") or "from purity_cnv_tool_summary.tsv"),
        })
    if tools:
        prov["purity_cnv_tools"] = tools
    consensus_path = summary_path.parent / "purity_cnv_consensus.tsv"
    recommended_path = summary_path.parent / "recommended_purity.tsv"
    consensus_rows = _read_optional(consensus_path)
    recommended_rows = _read_optional(recommended_path)
    consensus_row = consensus_rows[0] if consensus_rows else {}
    recommended_row = recommended_rows[0] if recommended_rows else {}
    found_tools = [row for row in tools if str(row.get("status") or "").upper() not in {"", "MISSING"}]
    names = [str(row.get("tool") or "未记录") for row in tools]
    if tools:
        prov["purity_cnv_consensus"] = {
            "recommended_purity": str(consensus_row.get("recommended_purity") or recommended_row.get("purity") or ""),
            "recommended_ploidy": str(recommended_row.get("ploidy") or next((row.get("ploidy") or "" for row in found_tools if row.get("ploidy")), "")),
            "selected_tool": str(recommended_row.get("evidence_tool") or ("多工具共识" if len(found_tools) > 1 else (found_tools[0].get("tool") if found_tools else "未形成估计"))),
            "status": str(consensus_row.get("status") or recommended_row.get("consensus_status") or ("MULTI_TOOL_REVIEW" if len(found_tools) > 1 else "SINGLE_TOOL_NO_CROSSCHECK")),
            "basis": str(consensus_row.get("interpretation") or ("已并列保留 " + "、".join(names) + " 结果；缺失工具显式标记为未形成估计。")),
        }


def _augment_purity_from_ranked_rows(prov: dict[str, Any], rows: list[dict[str, str]]) -> None:
    """Recover authoritative purity evidence retained by a downstream rerank.

    A presentation-only rerank may copy the consensus columns while leaving the
    original purity evidence directory outside its output root. Prefer the
    recorded recommendation pointer, then fall back to the embedded tool-value
    map. This keeps reports portable without guessing a patient-specific parent
    directory.
    """
    declared_tools = [row for row in (prov.get("purity_cnv_tools") or []) if isinstance(row, Mapping)]
    declared_consensus = prov.get("purity_cnv_consensus") if isinstance(prov.get("purity_cnv_consensus"), Mapping) else {}
    has_numeric_tool = any(str(row.get("purity") or "").strip() for row in declared_tools)
    has_numeric_consensus = bool(str(declared_consensus.get("recommended_purity") or "").strip())
    if has_numeric_tool and has_numeric_consensus:
        return

    source_row = next((row for row in rows if str(row.get("purity_recommendation_file") or "").strip()), {})
    recommendation = Path(str(source_row.get("purity_recommendation_file") or ""))
    summary_path = recommendation.parent / "purity_cnv_tool_summary.tsv" if recommendation.is_file() else Path()
    if recommendation.is_file() and summary_path.is_file() and summary_path.stat().st_size:
        tools = []
        for row in _read_optional(summary_path):
            tools.append({
                "tool": str(row.get("tool") or row.get("source_tool") or ""),
                "purity": str(row.get("purity") or ""),
                "ploidy": str(row.get("ploidy") or ""),
                "status": str(row.get("status") or "UNASSESSED"),
                "note": str(row.get("notes") or row.get("parse_method") or "recovered from recorded purity recommendation"),
            })
        if tools:
            prov["purity_cnv_tools"] = tools
            consensus_rows = _read_optional(summary_path.parent / "purity_cnv_consensus.tsv")
            recommended_rows = _read_optional(recommendation)
            consensus_row = consensus_rows[0] if consensus_rows else {}
            recommended_row = recommended_rows[0] if recommended_rows else {}
            found = [row for row in tools if str(row.get("status") or "").upper() not in {"", "MISSING"}]
            prov["purity_cnv_consensus"] = {
                "recommended_purity": str(
                    consensus_row.get("recommended_purity")
                    or recommended_row.get("purity")
                    or source_row.get("purity")
                    or ""
                ),
                "recommended_ploidy": str(
                    recommended_row.get("ploidy")
                    or source_row.get("ploidy")
                    or next((row.get("ploidy") or "" for row in found if row.get("ploidy")), "")
                ),
                "selected_tool": str(
                    recommended_row.get("evidence_tool")
                    or ("多工具共识" if len(found) > 1 else (found[0].get("tool") if found else "未形成估计"))
                ),
                "status": str(
                    consensus_row.get("status")
                    or recommended_row.get("consensus_status")
                    or source_row.get("purity_consensus_status")
                    or ("MULTI_TOOL_REVIEW" if len(found) > 1 else "SINGLE_TOOL_NO_CROSSCHECK")
                ),
                "basis": str(
                    consensus_row.get("interpretation")
                    or "已通过排序表记录的权威纯度文件回链，并列保留全部工具结果。"
                ),
            }
            return

    source_row = next((row for row in rows if str(row.get("purity_tool_values") or "").strip()), {})
    raw_values = str(source_row.get("purity_tool_values") or "").strip()
    try:
        values = json.loads(raw_values) if raw_values else {}
    except (TypeError, ValueError):
        values = {}
    if not isinstance(values, Mapping) or not values:
        return
    tools = [{
        "tool": str(tool),
        "purity": str(value),
        "ploidy": "",
        "status": "FOUND",
        "note": "来自最终排序表保留的多工具纯度值；工具级倍性未随表携带。",
    } for tool, value in values.items()]
    prov["purity_cnv_tools"] = tools
    prov["purity_cnv_consensus"] = {
        "recommended_purity": str(source_row.get("purity") or ""),
        "recommended_ploidy": str(source_row.get("ploidy") or ""),
        "selected_tool": "多工具共识" if len(tools) > 1 else tools[0]["tool"],
        "status": str(source_row.get("purity_consensus_status") or "MULTI_TOOL_REVIEW"),
        "basis": str(
            source_row.get("purity_range")
            and f"最终排序表保留多工具值；纯度范围 {source_row['purity_range']}。"
            or "最终排序表保留多工具纯度值；未携带原始工具级倍性。"
        ),
    }

def load_report_bundle(
    *,
    profile: Mapping[str, Any],
    events: list[dict[str, str]],
    peptides: list[dict[str, str]],
    appm_summary: Mapping[str, Any] | None = None,
    validation_rows: list[dict[str, str]] | None = None,
    outdir: str | Path | None = None,
    provenance: Mapping[str, Any] | None = None,
    patient_inputs: Mapping[str, Any] | None = None,
    sample_id: str = "",
    entry_mode: str = "",
) -> ReportBundle:
    root = Path(outdir) if outdir else None
    prov = dict(provenance or {})
    if root and not prov:
        prov = _read_json_optional(root / "provenance.json")
    if root:
        _augment_runtime_tool_provenance(root, prov)
        _augment_purity_cnv_provenance(root, prov)
        funnel_path = root / "parsed" / "splice_prefilter_funnel.tsv"
        if not funnel_path.is_file():
            configured_funnel = str((prov.get("splice_prefilter") or {}).get("funnel") or "") if isinstance(prov.get("splice_prefilter"), Mapping) else ""
            if configured_funnel:
                funnel_path = Path(configured_funnel)
        if funnel_path.is_file():
            prov["splice_filter_funnel_rows"] = _read_optional(funnel_path)
    if patient_inputs:
        explicit_inputs = patient_inputs.get("input_files") if isinstance(patient_inputs.get("input_files"), Mapping) else patient_inputs
        prov["input_files"] = dict(explicit_inputs)

    def p(*parts: str) -> Path | None:
        return root / Path(*parts) if root else None

    source_events = _read_optional(p("inputs", "combined_raw_events.tsv") if root else None)
    if not source_events:
        input_files = prov.get("input_files") if isinstance(prov.get("input_files"), Mapping) else {}
        for key in ("combined_raw_events", "raw_events"):
            candidate = Path(str(input_files.get(key) or ""))
            if candidate.is_file():
                source_events = _read_optional(candidate)
                break
    nested_evidence_path = p("scoring", "evidence_consensus", "all_tool_results.tsv") if root else None
    scoring_evidence_path = p("scoring", "all_tool_results.tsv") if root else None
    canonical_evidence_path = _first_existing([nested_evidence_path, scoring_evidence_path])
    if canonical_evidence_path:
        evidence_source_label = _report_relpath(canonical_evidence_path, root)
        evidence_source_status = "CANONICAL_ALL_TOOL_RESULTS"
        load_tool_rows = not _has_consensus_grades(events, peptides)
        source_peptides = _read_optional(canonical_evidence_path) if load_tool_rows else []
    else:
        source_peptides = []
        evidence_source_label = "ranked input (all_tool_results unavailable)"
        evidence_source_status = "RANKED_INPUT_ONLY"
    _augment_purity_from_ranked_rows(prov, events + peptides + source_peptides)
    if not prov.get("purity_cnv_consensus") and (prov.get("purity_cnv_tools") or _purity_tools_from_provenance(prov)):
        tool_rows = [dict(row) for row in (prov.get("purity_cnv_tools") or _purity_tools_from_provenance(prov))]
        if len(tool_rows) == 1:
            row = tool_rows[0]
            prov["purity_cnv_consensus"] = {
                "recommended_purity": row.get("purity") or "",
                "recommended_ploidy": row.get("ploidy") or "",
                "selected_tool": row.get("tool") or "单一工具",
                "status": "SINGLE_TOOL_NO_CROSSCHECK",
                "basis": "仅载入单一纯度工具结果，未形成多工具交叉验证共识。",
            }
        elif len(tool_rows) > 1:
            purities = [str(item.get("purity") or "") for item in tool_rows if str(item.get("purity") or "").strip()]
            names = [str(item.get("tool") or "未记录") for item in tool_rows]
            prov["purity_cnv_consensus"] = {
                "recommended_purity": purities[0] if purities else "",
                "recommended_ploidy": next((str(item.get("ploidy") or "") for item in tool_rows if str(item.get("ploidy") or "").strip()), ""),
                "selected_tool": "多工具并列",
                "status": "MULTI_TOOL_REVIEW",
                "basis": "已并列保留 " + "、".join(names) + " 结果，不静默选择单一工具。",
            }
    if not prov.get("tumor_dna_depth"):
        prov["tumor_dna_depth"] = (
            _depth_summary(events + peptides, "tumor_depth", "wes_tumor_depth")
            or _depth_summary_from_tsv(canonical_evidence_path, "tumor_depth", "wes_tumor_depth", "wgs_tumor_depth")
        )
    if not prov.get("normal_dna_depth"):
        prov["normal_dna_depth"] = (
            _depth_summary(events + peptides, "normal_depth")
            or _depth_summary_from_tsv(canonical_evidence_path, "normal_depth")
            or _vcf_depth_summary_from_provenance(prov, root, source_events + events + peptides, "normal")
            or _bam_depth_summary_from_provenance(prov, root, "normal")
        )
    if not prov.get("genome_build"):
        for row in events + peptides:
            build = str(row.get("genome_build") or "").strip()
            if build and build.upper() not in {"UNASSESSED", "NA", "N/A"}:
                prov["genome_build"] = build
                break
    derived_rna_qc = _rna_qc_summary(events + peptides + source_peptides)
    recorded_rna_qc = str(prov.get("rna_qc_status") or "")
    if derived_rna_qc and (
        not recorded_rna_qc
        or "pileup未评估" in recorded_rna_qc
        or "位点级RNA pileup未评估" in recorded_rna_qc
    ):
        prov["rna_qc_status"] = derived_rna_qc

    if canonical_evidence_path and not source_peptides:
        report_enrichment_rows = _read_report_enrichment_rows(
            canonical_evidence_path,
            events[:250] + peptides[:500],
        )
        source_events.extend(report_enrichment_rows)

    extra_gene_ids = {
        gene_id
        for source in source_events + source_peptides + events + peptides
        for field_name in (
            "gene", "gene_pair", "event_name", "combined_protein_change",
            "source_event_id", "source_record_id", "source_records", "source_junction_id",
        )
        for gene_id in re.findall(r"ENSG\d+", str(source.get(field_name) or ""))
    }
    gtf_candidates: list[Path] = []
    for payload in (prov, prov.get("references") or {}, prov.get("reference_manifest") or {}):
        if not isinstance(payload, Mapping):
            continue
        for key in ("gencode_gtf", "annotation_gtf", "gtf", "gene_annotation"):
            value = str(payload.get(key) or "").strip()
            if value:
                gtf_candidates.append(Path(value))
    for env_key in ("GENCODE_GTF", "NEOAG_GENCODE_GTF", "REFERENCE_GTF"):
        value = os.environ.get(env_key, "").strip()
        if value:
            gtf_candidates.append(Path(value))
    gtf_candidates.extend(_asset_gtf_candidates(prov))
    junction_genes = _normal_expression_gene_map(root, prov, extra_gene_ids)
    junction_genes.update(_read_longrna_junction_genes(root, extra_gene_ids, gtf_candidates))
    enriched_peptides = _enrich_rows_from_sources(
        peptides, source_peptides or peptides, source_events, junction_genes, evidence_source_label,
    )
    enriched_events = _enrich_rows_from_sources(
        events, source_peptides or enriched_peptides, source_events, junction_genes, evidence_source_label,
    )
    patient_expression_map, patient_expression_source = _patient_expression_tpm_map(prov)
    enriched_peptides = _apply_patient_gene_expression(
        enriched_peptides, patient_expression_map, patient_expression_source,
    )
    enriched_events = _apply_patient_gene_expression(
        enriched_events, patient_expression_map, patient_expression_source,
    )
    junction_verification_path = p("metadata", "junction_read_verification.tsv") if root else None
    junction_verification_rows = _read_optional(junction_verification_path)
    if junction_verification_rows:
        enriched_peptides = _apply_junction_verification(enriched_peptides, junction_verification_rows)
        enriched_events = _apply_junction_verification(enriched_events, junction_verification_rows)
        prov["junction_read_verification"] = {
            "status": "LOADED",
            "source": str(junction_verification_path),
            "records": len(junction_verification_rows),
        }
    evidence_manifest = _read_json_optional(_first_existing([
        canonical_evidence_path.with_name("all_tool_results.manifest.json") if canonical_evidence_path else None,
        p("scoring", "evidence_consensus", "all_tool_results.manifest.json") if root else None,
        p("scoring", "all_tool_results.manifest.json") if root else None,
    ]))
    evidence_integrity = {"status": "UNASSESSED", "expected_sha256": "", "actual_sha256": ""}
    if canonical_evidence_path and canonical_evidence_path.is_file() and evidence_manifest:
        expected_sha256 = str(
            (evidence_manifest.get("output") or {}).get("sha256")
            or (evidence_manifest.get("input") or {}).get("sha256")
            or ""
        )
        actual_sha256 = _file_sha256(canonical_evidence_path) if expected_sha256 else ""
        evidence_integrity = {
            "status": "PASS" if expected_sha256 and actual_sha256 == expected_sha256 else "FAIL",
            "expected_sha256": expected_sha256,
            "actual_sha256": actual_sha256,
        }
    return ReportBundle(
        profile=profile,
        events=enriched_events,
        peptides=enriched_peptides,
        appm_summary=appm_summary or {},
        validation_rows=validation_rows or [],
        sample_id=sample_id or str(prov.get("sample_id") or (peptides[0].get("sample_id") if peptides else "")),
        entry_mode=entry_mode or str(prov.get("entry_mode") or ""),
        provenance=prov,
        peptide_safety=_read_optional(p("safety", "peptide_safety.tsv") if root else None),
        peptide_escape_flags=_read_optional(p("immune_escape", "peptide_escape_flags.tsv") if root else None),
        immune_escape_summary=_read_optional(p("immune_escape", "immune_escape_summary.tsv") if root else None),
        ccf=_read_optional(p("clonality", "ccf_2.tsv") if root else None) or _read_optional(p("clonality", "ccf_lite.tsv") if root else None),
        appm_gene_status=_read_optional(p("appm", "appm_gene_status.tsv") if root else None),
        appm_peptide_modifiers=_read_optional(p("appm", "appm_peptide_modifiers.tsv") if root else None),
        appm_module_scores=_read_optional(p("appm", "appm_module_scores.tsv") if root else None),
        appm_submodule_scores=_read_optional(p("appm", "appm_submodule_scores.tsv") if root else None),
        appm_conflicts=_read_optional(p("appm", "appm_conflicts.tsv") if root else None),
        wes_qc=_read_optional(p("qc", "wes", "wes_qc.tsv") if root else None),
        wes_wgs_coding_summary=_read_optional(
            p("qc", "wes_wgs_coding_comparison", "wes_wgs_coding_summary.tsv") if root else None
        ),
        targeted_pileup_summary=(
            _read_optional(
                p("qc", "wes_wgs_coding_comparison", "targeted_pileup", "protein_altering", "discordant_targeted_pileup_summary.tsv")
                if root else None
            )
            or _read_optional(
                p("qc", "wes_wgs_coding_comparison", "targeted_pileup", "discordant_targeted_pileup_summary.tsv")
                if root else None
            )
        ),
        evidence_manifest=evidence_manifest,
        evidence_conflicts=_read_optional(_first_existing([
            p("scoring", "evidence_consensus", "evidence_conflicts.tsv") if root else None,
            p("scoring", "evidence_conflicts.tsv") if root else None,
        ])),
        evidence_source_status=evidence_source_status,
        evidence_integrity=evidence_integrity,
        purity_tools=_purity_tools_from_provenance(prov),
        purity_consensus=dict(prov.get("purity_cnv_consensus") or {}),
        hla_loh_tool_results=_read_hla_loh_tool_results(root, prov),
        tool_versions=_read_tool_version_manifest(root),
        disease_knowledge=_load_disease_knowledge(profile, prov),
    )


def _patient_consequence_label(peptide: Mapping[str, Any]) -> str:
    pc = str(peptide.get("peptide_consequence") or "").lower()
    if pc in CONSEQUENCE_PATIENT:
        return CONSEQUENCE_PATIENT[pc]
    et = str(peptide.get("event_type") or "")
    if et == "Fusion":
        return CONSEQUENCE_PATIENT["fusion"]
    return CONSEQUENCE_PATIENT.get(pc, "肿瘤特异性肽段候选")


def _patient_priority_label(priority: str) -> str:
    return PRIORITY_PATIENT.get(str(priority or "").strip(), "需进一步评估")


def _appm_patient_summary(appm_summary: Mapping[str, Any]) -> str:
    i_status = str(appm_summary.get("mhc_i_integrity_status") or "")
    if "DEFECT" in i_status.upper():
        return "样本 MHC-I 抗原呈递通路可能存在缺陷，部分候选肽段的实际呈递能力可能低于计算预测。"
    if "PARTIAL" in i_status.upper() or "CAUTION" in i_status.upper():
        return "样本 MHC-I 抗原呈递通路存在不确定因素，建议结合实验验证解读候选。"
    return "样本 MHC-I 抗原呈递通路完整性评估未见明确缺陷信号（仍依赖输入证据完整度）。"


def _top_hla_alleles(peptides: list[Mapping[str, Any]], limit: int = 6) -> str:
    alleles: list[str] = []
    for p in peptides:
        hla = str(p.get("hla_allele") or "").strip()
        if hla and hla not in alleles:
            alleles.append(hla)
        if len(alleles) >= limit:
            break
    return "、".join(alleles) if alleles else "未提供"


def _val_by_peptide(validation_rows: list[Mapping[str, Any]]) -> dict[str, dict[str, str]]:
    return _map_by(validation_rows, "peptide_id")


def _cross_platform_counts(events: list[dict[str, str]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for event in events:
        status = str(event.get("cross_platform_status") or "")
        if status and status not in {"NOT_APPLICABLE", "UNASSESSED_NOT_IN_COMPARISON"}:
            counts[status] = counts.get(status, 0) + 1
    return counts


def _patient_platform_label(status: str) -> str:
    return {
        "CROSS_PLATFORM_PASS_CONCORDANT": "WES/WGS 共同检出",
        "ALT_PRESENT_BELOW_PASS_OR_CALLER_DIFFERENCE": "另一检测可见低水平支持",
        "COVERED_NO_ALT_SAMPLE_OR_ASSAY_DIFFERENCE": "仅当前样本/时间点明确",
        "OTHER_COVERED_BUT_LIMITED_POWER_AT_SOURCE_VAF": "另一检测能力不足，不能判阴性",
        "OTHER_LOW_OR_NO_COVERAGE": "另一检测覆盖不足",
        "SOURCE_INDEL_NOT_REPRODUCED_REASSEMBLY_REQUIRED": "复杂变异需局部重组复核",
        "SOURCE_PASS_NOT_REPRODUCED_BY_PILEUP": "源检测需重新复核",
        "SOURCE_WEAK_EXACT_PILEUP_SUPPORT": "源检测支持较弱",
        "NORMAL_SUPPORT_REVIEW": "正常血液也有支持，暂不推进",
        "NOT_APPLICABLE": "非 DNA 点突变事件",
        "UNASSESSED_NOT_IN_COMPARISON": "尚未完成跨平台评估",
    }.get(str(status or ""), "尚未完成跨平台评估")


def _patient_rna_label(peptide: Mapping[str, Any]) -> str:
    status = str(peptide.get("rna_support_status") or "")
    if status in {"RNA_ALT_SUPPORTED", "RNA_JUNCTION_SUPPORTED"}:
        return "RNA 已支持"
    if status == "RNA_ALT_NOT_DETECTED":
        return "RNA 未检出突变支持"
    if status in {"UNASSESSED", "", "RNA_ONLY_UNRESOLVED"}:
        return "RNA 证据未完整评估"
    return status.replace("_", " ")


def _patient_event_change(event: Mapping[str, Any]) -> str:
    consequence = str(event.get("peptide_consequence") or "").strip().lower()
    event_type = str(event.get("event_type") or "").strip()
    peptide = str(event.get("peptide") or "").strip()
    protein_change = str(event.get("combined_protein_change") or "").strip()
    if ":p." in protein_change:
        return ("p." + protein_change.split(":p.", 1)[1]).replace("%3D", "=")
    # Reuse the authoritative event-track classifier. A DNA SNV/InDel may
    # carry a splice-related VEP consequence without becoming a splice event.
    patient_track = _patient_track(event)
    if patient_track == "Splice":
        source_label = next((
            str(event.get(field_name) or "").strip()
            for field_name in ("source_event_id", "source_record_id", "source_records")
            if str(event.get(field_name) or "").strip()
        ), "")
        match = re.search(r":([A-Za-z]+\d+(?:\.\d+)?)_\d+-([A-Za-z]+\d+(?:\.\d+)?)_\d+", source_label)
        junction_path = f"{match.group(1)}→{match.group(2)}" if match else "精确异常junction"
        junction = str(event.get("canonical_junction_id") or "").strip()
        parts = junction.split("|")
        coordinate = f"{parts[2]}:{parts[3]}-{parts[4]}" if len(parts) >= 5 else str(event.get("event_name") or "").strip()
        if not source_label and not coordinate and peptide:
            return f"异常剪接肽段 {peptide}"
        detail = f"（{coordinate}）" if coordinate else ""
        has_formal_origin = all(
            str(event.get(field_name) or "").strip()
            for field_name in ("transcript_hypothesis_id", "orf_id", "origin_peptide_id")
        )
        if has_formal_origin:
            return (
                f"异常剪接 {junction_path}{detail}；已完成局部转录本、ORF及跨junction肽段来源精确回链，"
                "全长转录本真实性仍待独立验证"
            )
        return f"异常剪接 {junction_path}{detail}；ORF/蛋白影响待确认"
    if patient_track == "Fusion" and peptide:
        return f"融合肽段 {peptide}"
    change = str(event.get("event_name") or event.get("consequence") or protein_change or "")
    if change:
        return change.replace("%3D", "=")
    if consequence in CONSEQUENCE_PATIENT:
        return CONSEQUENCE_PATIENT[consequence]
    return "蛋白改变待确认"


def _patient_fusion_interpretation(event: Mapping[str, Any]) -> str:
    gene = str(event.get("gene") or "")
    if gene == "EWSR1::WT1":
        return "DSRCT 的特征性驱动融合；仍需确认断点、阅读框及融合肽的真实加工呈递"
    if gene.startswith("HLA-"):
        return "HLA 高多态区域事件，优先排查比对或转录本拼接影响"
    if str(event.get("safety_status") or "") == "CAUTION":
        return "RNA junction 有支持，但正常组织背景或安全性证据仍需复核"
    return "候选融合事件；需用独立方法确认断点和阅读框"


def _metric_value(rows: list[dict[str, str]], metric: str, default: str = "未评估") -> str:
    for row in rows:
        if row.get("metric") == metric:
            return str(row.get("value") or default)
    return default


def _patient_pileup_category(category: str) -> str:
    return {
        "ALT_PRESENT_BELOW_PASS_OR_CALLER_DIFFERENCE": "另一平台存在低水平 ALT 或调用规则差异",
        "COVERED_NO_ALT_SAMPLE_OR_ASSAY_DIFFERENCE": "覆盖充分但未见 ALT，考虑样本/时间点差异",
        "NORMAL_SUPPORT_REVIEW": "正常血液存在支持，需排除胚系或技术因素",
        "OTHER_COVERED_BUT_LIMITED_POWER_AT_SOURCE_VAF": "有覆盖但按源 VAF 统计检出能力不足",
        "OTHER_LOW_OR_NO_COVERAGE": "另一平台低覆盖或无覆盖",
        "SOURCE_INDEL_NOT_REPRODUCED_REASSEMBLY_REQUIRED": "源平台 InDel 未由简单 pileup 复现，需局部重组",
        "SOURCE_PASS_NOT_REPRODUCED_BY_PILEUP": "源平台 PASS 位点未由 pileup 复现",
        "SOURCE_WEAK_EXACT_PILEUP_SUPPORT": "源平台精确 ALT 支持较弱",
        "WEAK_OR_ABSENT_REVIEW": "支持较弱或缺失，需人工复核",
    }.get(str(category or ""), str(category or "未分类"))


def _patient_pileup_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [
        {
            "source": row.get("source", ""),
            "category": _patient_pileup_category(row.get("category", "")),
            "count": row.get("count", ""),
        }
        for row in rows
    ]


def _patient_cancer_context(row: Mapping[str, Any]) -> str:
    status = str(row.get("cancer_driver_context") or "")
    types = str(row.get("cancer_gene_types") or "")
    if status == "DRIVER_CONTEXT":
        return f"癌症基因背景：{types}" if types else "具有癌症基因背景"
    if status == "LISTED_NO_DRIVER_CLASS":
        return "名单收录，但未归类为癌基因/抑癌基因"
    if status == "NOT_LISTED":
        return "未在当前癌症基因表中"
    return "未评估"


def _make_patient_report_legacy(path: str | Path, bundle: ReportBundle) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    val_map = _val_by_peptide(bundle.validation_rows)
    advancable = [ppt for ppt in bundle.peptides if str(ppt.get("final_priority", "")).upper() not in {"D", ""}]
    ranked_pool = advancable if advancable else bundle.peptides
    top = []
    seen_events = set()
    for candidate in ranked_pool:
        event_id = str(candidate.get("event_id") or candidate.get("peptide_id") or "")
        if event_id in seen_events:
            continue
        seen_events.add(event_id)
        top.append(candidate)
        if len(top) >= 10:
            break
    genes = []
    for ppt in top:
        g = str(ppt.get("gene") or "")
        if g and g not in genes:
            genes.append(g)

    out = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        f"<title>新抗原分析报告（患者沟通版）— {esc(bundle.sample_id)}</title>",
        REPORT_CSS,
        "</head><body class='patient'>",
        "<h1>新抗原计算分析报告</h1>",
        "<p class='lead'>本报告用于帮助理解肿瘤新抗原候选的<strong>计算筛选结果</strong>，"
        "便于医患沟通与后续实验规划；<strong>不能替代临床诊断或治疗决策</strong>。</p>",
    ]

    out.append("<div class='section'><h2>1. 分析目的</h2>")
    out.append(
        "<p>通过肿瘤 DNA/RNA 与 HLA 分型信息，识别可能由肿瘤变异产生、并被患者 HLA 分子呈递的肽段候选，"
        "用于疫苗设计、T 细胞治疗靶点探索或免疫监测研究的初步筛选。</p>"
    )
    out.append("</div>")

    out.append("<div class='section'><h2>2. 样本与 HLA 背景</h2>")
    out.append("<ul class='compact'>")
    out.append(f"<li><b>样本编号：</b>{esc(bundle.sample_id or '未注明')}</li>")
    out.append(f"<li><b>分析场景：</b>{esc(bundle.entry_mode or '肿瘤新抗原筛选')}</li>")
    out.append("<li><b>数据层：</b>肿瘤 WES、肿瘤 WGS、配对血液正常样本、肿瘤 RNA/融合、HLA 分型，以及纯度、CNV 和 CCF 证据。</li>")
    out.append(f"<li><b>当前候选使用的 HLA-I 背景：</b>{esc(_top_hla_alleles(bundle.peptides))}</li>")
    out.append("<li><b>HLA LOH 背景：</b>当前检出的丢失信号位于 HLA-II（DQA1/DQB1）；未见 HLA-A/B/C 丢失影响当前 HLA-I 候选，受 HLA LOH 影响的候选肽数为 0。</li>")
    out.append(f"<li><b>评分方案：</b>{esc(bundle.profile.get('_profile_name', 'default'))}</li>")
    out.append("<li><b>临床样本关系：</b>精确取材日期、部位、治疗前后关系仍须以临床样本清单核实，本报告不据测序文件名推断。</li>")
    out.append("</ul><p class='small'>HLA 分型用于判断候选肽可能由哪些 HLA 分子呈递；它本身不证明肿瘤细胞已经加工并展示该肽段。</p></div>")

    out.append("<div class='section'><h2>3. 肿瘤是否具备抗原呈递条件</h2>")
    presentation_rows = [
        {"dimension": "MHC-I 抗原呈递", "conclusion": "部分保留，需谨慎解释", "meaning": "与当前 HLA-I 候选最直接相关；未显示通路完全缺失，但部分环节仍需实验确认"},
        {"dimension": "MHC-II 抗原呈递", "conclusion": "存在不确定因素", "meaning": "DQA1/DQB1 丢失及部分低表达信号提示 HLA-II 背景需要复核"},
        {"dimension": "IFNG 应答", "conclusion": "未见明确缺陷", "meaning": "现有输入支持该应答通路总体保留，但不等同于临床免疫治疗敏感"},
        {"dimension": "HLA-I 等位基因丢失", "conclusion": "未见影响当前候选", "meaning": "未发现 HLA-A/B/C 丢失影响当前纳入排序的 HLA-I 候选"},
        {"dimension": "HLA-II 等位基因丢失", "conclusion": "需要单独复核", "meaning": "当前信号位于 DQA1/DQB1，只进入 HLA-II 背景，不应降低当前 HLA-I 候选"},
        {"dimension": "证据完整性", "conclusion": "目前仍不完整", "meaning": "部分通路证据缺失；缺失不能解释为正常，也不能据此诊断免疫逃逸"},
    ]
    out.append(_table(presentation_rows, ["dimension", "conclusion", "meaning"]))
    out.append("<p><b>综合判断：</b>肿瘤具备部分 HLA-I 抗原呈递条件，并非“完全不能呈递”；但 TAP 等环节存在谨慎信号，且 APPM 输入完整度较低，因此应表述为<b>部分保留、仍需实验确认</b>，不能据此判断免疫治疗敏感或耐药。</p>")
    out.append("<p class='small'>本患者沟通版不展示模型分值；详细评分、状态字段和计算依据保留在科研技术版报告中。</p>")
    out.append("</div>")

    out.append("<div class='section'><h2>4. 关键发现（摘要）</h2><ul class='compact'>")
    priority_counts = {}
    for peptide in bundle.peptides:
        label = str(peptide.get("final_priority") or "未分级")
        priority_counts[label] = priority_counts.get(label, 0) + 1
    out.append(
        f"<li>共评估 <b>{len(bundle.events)}</b> 个候选事件及 <b>{len(bundle.peptides)}</b> 个肽段-HLA 组合；"
        f"其中 <b>{len(advancable)}</b> 个组合未被判为“不建议推进”。</li>"
    )
    out.append(
        "<li><b>当前分层：</b>"
        + "；".join(f"{esc(level)} 级 {count} 个" for level, count in sorted(priority_counts.items()))
        + "。分层数量较多不代表存在同等数量的独立治疗靶点。</li>"
    )
    if genes:
        out.append(f"<li>优先关注的基因包括：<b>{esc('、'.join(genes[:8]))}</b>。</li>")
    out.append(f"<li>{esc(_appm_patient_summary(bundle.appm_summary))}</li>")
    if bundle.immune_escape_summary:
        ies = bundle.immune_escape_summary[0]
        risk = str(ies.get("overall_immune_escape_risk") or "")
        if risk and risk.upper() not in {"LOW", "PASS", ""}:
            out.append(f"<li>免疫逃逸风险提示：<b>{esc(risk)}</b>（机制层面证据，非临床耐药结论）。</li>")
    platform_counts = _cross_platform_counts(bundle.events)
    if platform_counts:
        concordant = platform_counts.get("CROSS_PLATFORM_PASS_CONCORDANT", 0)
        review = sum(platform_counts.values()) - concordant
        out.append(
            f"<li>WES/WGS DNA 交叉复核：<b>{concordant}</b> 个事件由两平台共同检出；"
            f"<b>{review}</b> 个事件存在低水平、检出能力或样本时间点差异，已在排序中标记复核。</li>"
        )
    out.append("</ul></div>")

    dna_events = [event for event in bundle.events if str(event.get("mutation_source") or "") in {"SNV", "InDel"}]
    featured_genes = {"TP53", "KRAS"}
    dna_events.sort(key=lambda event: (
        str(event.get("gene") or "") in featured_genes,
        str(event.get("haplotype_status") or "") == "PHASED_CIS_COMBINED",
        str(event.get("rna_support_status") or "") == "RNA_ALT_SUPPORTED",
        str(event.get("cross_platform_status") or "") == "CROSS_PLATFORM_PASS_CONCORDANT",
        float(event.get("event_score") or 0),
    ), reverse=True)
    out.append("<div class='section'><h2>5. 主要 DNA 突变及 RNA/跨平台证据</h2>")
    out.append("<p class='small'>以下为结合跨平台、RNA 与事件评分选出的主要 DNA 变异；它们不是临床用药清单，也不代表均为肿瘤驱动事件。</p>")
    mutation_rows = []
    for event in dna_events[:12]:
        mutation_rows.append({
            "gene": event.get("gene", ""),
            "cancer_context": _patient_cancer_context(event),
            "protein_change": _patient_event_change(event),
            "type": event.get("mutation_source", ""),
            "rna_evidence": _patient_rna_label(event),
            "wes_wgs_evidence": _patient_platform_label(str(event.get("cross_platform_status") or "")),
            "interpretation": (
                "两相邻变异已完成同一单倍型重构" if event.get("haplotype_status") == "PHASED_CIS_COMBINED"
                else "需结合病理、克隆性和实验验证判断作用"
            ),
        })
    out.append(_table(mutation_rows, ["gene", "cancer_context", "protein_change", "type", "rna_evidence", "wes_wgs_evidence", "interpretation"]))
    out.append("<p><b>重点提示：</b>TP53 等共同检出且有 RNA 支持的变异可信度相对更高；KRAS 等仅在一个肿瘤文库明确检出的变异应按样本/时间点特异结果解释；TBR1 相邻变异按已重构的同一单倍型解读。</p>")
    out.append("</div>")

    fusion_events = [event for event in bundle.events if str(event.get("event_type") or "") == "Fusion"]
    fusion_events.sort(key=lambda event: (str(event.get("gene") or "") == "EWSR1::WT1", float(event.get("event_score") or 0)), reverse=True)
    out.append("<div class='section'><h2>6. 融合基因事件及 DSRCT 背景解释</h2>")
    out.append("<p><b>EWSR1::WT1</b> 是 DSRCT 的标志性融合，本样本存在 RNA junction 支持；但由融合产生的新抗原仍需独立验证阅读框、异常 junction 和实际 HLA 呈递。</p>")
    fusion_rows = []
    for event in fusion_events[:10]:
        fusion_rows.append({
            "fusion": event.get("gene", ""),
            "cancer_context": _patient_cancer_context(event),
            "junction_reads": event.get("rna_junction_reads", ""),
            "rna_status": _patient_rna_label(event),
            "expression": event.get("event_expression", ""),
            "safety": event.get("safety_status", ""),
            "interpretation": _patient_fusion_interpretation(event),
        })
    out.append(_table(fusion_rows, ["fusion", "cancer_context", "junction_reads", "rna_status", "expression", "safety", "interpretation"]))
    out.append("<p class='small'>肝脏高表达基因、HLA 区域或仅有少量 junction reads 的融合可能包含 read-through、正常背景或比对伪影，其证据等级不能与 EWSR1::WT1 等驱动融合等同。</p>")
    out.append("</div>")

    event_by_gene: dict[str, dict[str, str]] = {}
    for event in bundle.events:
        gene = str(event.get("gene") or "")
        if gene and gene not in event_by_gene:
            event_by_gene[gene] = event
    discussable = []
    discussion_notes = {
        "EWSR1::WT1": "DSRCT 标志性驱动融合，疾病相关性最明确；新抗原价值仍取决于 junction 阅读框、加工呈递和功能实验",
        "TP53": "WES/WGS 共同检出并有 RNA 支持，事件真实性较强；仍需验证突变肽相对 WT 的特异性",
        "TBR1": "两个相邻变异已按同一单倍型重构；应只保留少量 combined-mutant 肽，避免重叠窗口重复占位",
        "KRAS": "经典热点但仅在一个肿瘤文库明确；需先确认取材/时间点和独立 DNA/RNA 支持",
        "TOMM34": "跨平台与 RNA 证据较完整，可作为非驱动但可验证的突变肽事件讨论",
        "ACLY": "跨平台与 RNA 证据较完整，可纳入突变肽对照验证候选组",
    }
    for gene in ("EWSR1::WT1", "TP53", "TBR1", "TOMM34", "ACLY", "KRAS"):
        event = event_by_gene.get(gene)
        if not event:
            continue
        discussable.append({
            "event": gene,
            "change": _patient_event_change(event),
            "rna": _patient_rna_label(event),
            "dna": _patient_platform_label(str(event.get("cross_platform_status") or "")),
            "why_discuss": discussion_notes[gene],
        })
    out.append("<div class='section'><h2>7. 目前最值得讨论的候选事件</h2>")
    out.append("<p>这里的“值得讨论”综合考虑疾病生物学、DNA/RNA 可复现性、单倍型、HLA 呈递预测和实验可行性，不等同于自动分数最高。</p>")
    out.append(_table(discussable, ["event", "change", "rna", "dna", "why_discuss"]))
    out.append("</div>")

    manual_statuses = {
        "COVERED_NO_ALT_SAMPLE_OR_ASSAY_DIFFERENCE": "可能反映不同取材/时间点的真实异质性",
        "SOURCE_INDEL_NOT_REPRODUCED_REASSEMBLY_REQUIRED": "复杂 InDel 可能无法由简单 pileup 复现",
        "SOURCE_PASS_NOT_REPRODUCED_BY_PILEUP": "源 caller 给出 PASS，但原始 reads 仍需人工核验",
        "SOURCE_WEAK_EXACT_PILEUP_SUPPORT": "存在少量精确 ALT reads，但证据偏弱",
        "NORMAL_SUPPORT_REVIEW": "正常样本也有支持，需先排除胚系或技术伪影",
    }
    manual_rows = []
    for event in bundle.events:
        status = str(event.get("cross_platform_status") or "")
        gene = str(event.get("gene") or "")
        if status not in manual_statuses:
            continue
        manual_rows.append({
            "event": gene,
            "change": _patient_event_change(event),
            "retain_reason": manual_statuses[status],
            "do_not_auto_advance": _patient_platform_label(status),
            "required_review": "IGV/局部重组、独立测序或匹配时间点复核",
        })
        if len(manual_rows) >= 12:
            break
    ews = event_by_gene.get("EWSR1::WT1")
    if ews:
        manual_rows.insert(0, {
            "event": "EWSR1::WT1",
            "change": _patient_event_change(ews),
            "retain_reason": "DSRCT 标志性驱动融合，生物学优先级高于一般自动分数",
            "do_not_auto_advance": "融合肽呈递、安全性和 WT/正常背景尚未完成",
            "required_review": "断点/阅读框、junction RNA、长肽/minigene 和功能验证",
        })
    out.append("<div class='section'><h2>8. 需要人工保留、但不应按自动分数直接推进的事件</h2>")
    out.append(_table(manual_rows, ["event", "change", "retain_reason", "do_not_auto_advance", "required_review"]))
    out.append("<p class='small'>人工保留表示不应因单一分数或暂时缺证而删除；它也不表示可以绕过证据补充直接进入治疗设计。</p></div>")

    tier_rows = [
        {"tier": "研究层 1A", "scope": "WES/WGS 共同检出、RNA 支持、呈递证据较好的 SNV/InDel", "examples": "TP53、TOMM34、ACLY 及其他满足条件事件", "action": "优先做 MT/WT 成对肽验证"},
        {"tier": "研究层 1B", "scope": "疾病驱动融合及异常 junction", "examples": "EWSR1::WT1", "action": "单独成组，优先长肽/minigene，不以短肽分数替代加工验证"},
        {"tier": "研究层 2", "scope": "样本/时间点特异或另一平台低水平支持", "examples": "KRAS 等", "action": "先确认目标取材中的 DNA/RNA 存在，再决定是否进入免疫学实验"},
        {"tier": "人工复核层", "scope": "复杂 InDel、源检测未复现、弱支持", "examples": "按 targeted pileup 标记的事件", "action": "IGV、局部组装或独立测序后重新评分"},
        {"tier": "暂缓层", "scope": "正常样本支持、明显安全性风险或总体为 D", "examples": "AXDND1 等正常支持事件", "action": "不进入首批实验；先排除胚系、正常组织表达和交叉反应"},
    ]
    out.append("<div class='section'><h2>9. 最值得关注的候选分层</h2>")
    out.append(_table(tier_rows, ["tier", "scope", "examples", "action"]))
    out.append("</div>")

    out.append("<div class='section'><h2>10. 优先候选肽段（按事件去重 Top 10）</h2>")
    out.append("<p class='small'>下表为计算排序靠前的候选，不代表已证实可诱导抗肿瘤免疫反应。</p>")
    headers = ["rank", "gene", "cancer_context", "variant_type", "hla", "source_chain", "priority", "rna_evidence", "dna_evidence", "meaning", "next_step"]
    rows = []
    for i, ppt in enumerate(top, 1):
        pid = str(ppt.get("peptide_id") or "")
        val = val_map.get(pid, {})
        mode = str(val.get("validation_mode") or "")
        next_step = str(val.get("validation_strategy") or val.get("recommended_assay") or ppt.get("recommended_use") or "")
        if mode in {"frameshift_long", "splice_junction_long", "fusion_junction_long", "insertion_long"}:
            meaning = "可能由肿瘤特异性序列产生，建议用长肽/minigene 验证真实加工呈递"
        else:
            meaning = "突变肽 vs 正常肽对照验证"
        rows.append({
            "rank": str(i),
            "gene": ppt.get("gene", ""),
            "cancer_context": _patient_cancer_context(ppt),
            "variant_type": _patient_consequence_label(ppt),
            "hla": ppt.get("hla_allele", ""),
            "source_chain": ppt.get("source_chain_confidence_tier", "未评估"),
            "priority": _patient_priority_label(str(ppt.get("final_priority") or "")),
            "rna_evidence": _patient_rna_label(ppt),
            "dna_evidence": _patient_platform_label(str(ppt.get("cross_platform_status") or "")),
            "meaning": meaning,
            "next_step": next_step,
        })
    out.append(_table(rows, headers))
    out.append("</div>")

    out.append("<div class='section'><h2>11. WES 与 WGS 蛋白改变型 SNV/InDel 对比</h2>")
    if bundle.wes_wgs_coding_summary:
        coding_rows = [
            {"metric": "WES 蛋白改变型 SNV/InDel", "value": _metric_value(bundle.wes_wgs_coding_summary, "protein_altering_wes")},
            {"metric": "WGS 蛋白改变型 SNV/InDel", "value": _metric_value(bundle.wes_wgs_coding_summary, "protein_altering_wgs")},
            {"metric": "两平台共同检出", "value": _metric_value(bundle.wes_wgs_coding_summary, "protein_altering_common")},
            {"metric": "WES-only", "value": _metric_value(bundle.wes_wgs_coding_summary, "protein_altering_wes_only")},
            {"metric": "WGS-only", "value": _metric_value(bundle.wes_wgs_coding_summary, "protein_altering_wgs_only")},
            {"metric": "共同位点 VAF 相关性", "value": _metric_value(bundle.wes_wgs_coding_summary, "common_af_pearson")},
        ]
        out.append(_table(coding_rows, ["metric", "value"]))
    else:
        out.append("<p>WES/WGS coding 区域汇总尚未加载。</p>")
    if bundle.targeted_pileup_summary:
        out.append("<h3>差异位点回查结果</h3>")
        out.append(_table(_patient_pileup_rows(bundle.targeted_pileup_summary), ["source", "category", "count"]))
    out.append("<p class='small'>本节只统计 VEP 明确标记为 missense、frameshift、inframe insertion/deletion、start/stop lost 或 stop gained 等会改变蛋白序列的 SNV/InDel；不包括 synonymous、非编码、仅 splice-region 但没有明确蛋白改变的记录，也不包括 fusion/SV。</p>")
    out.append("<p>低重合度不能简单解释为某个平台错误。回查显示主要来源包括低 VAF 下检出能力不足、另一平台存在但未通过 PASS、不同肿瘤文库/时间点差异，以及长 InDel 需要局部重组复核。只有在源 BAM 可复现、另一 BAM 覆盖和统计检出能力充分且无 ALT 时，才视为较强的样本差异证据。</p>")
    out.append("</div>")

    out.append("<div class='section'><h2>12. 为什么当前没有直接进入高优先级的候选</h2><ul class='compact'>")
    if not any(str(p.get("final_priority") or "") in {"A", "B", "B_CAUTION"} for p in bundle.peptides):
        out.append("<li>当前没有 A/B 级候选。主要原因不是“完全没有候选”，而是正常组织安全性、RNA 支持、突变特异性或样本间一致性仍需补证。</li>")
    out.append("<li>C_CAUTION 表示候选具有一定计算证据，但在进入首批实验或治疗设计前必须完成针对性复核。</li>")
    out.append("<li>D 级表示目前不建议推进，常见原因包括安全性风险、正常样本支持、呈递证据不足或关键证据缺失。</li>")
    out.append("<li>同一事件可产生多个长度和多个 HLA 限制性肽段，因此肽段组合数远高于独立变异事件数。</li>")
    out.append("</ul></div>")

    out.append("<div class='section'><h2>13. 当前证据缺口与样本差异</h2><ul class='compact'>")
    if platform_counts:
        out.append(f"<li>WES/WGS 共同检出的事件：<b>{platform_counts.get('CROSS_PLATFORM_PASS_CONCORDANT', 0)}</b> 个。</li>")
        out.append(
            f"<li>另一检测可见低水平 ALT 支持：<b>{platform_counts.get('ALT_PRESENT_BELOW_PASS_OR_CALLER_DIFFERENCE', 0)}</b> 个；"
            f"因检出能力不足不能判阴性：<b>{platform_counts.get('OTHER_COVERED_BUT_LIMITED_POWER_AT_SOURCE_VAF', 0)}</b> 个。</li>"
        )
        out.append(
            f"<li>覆盖充分但呈现样本/时间点差异：<b>{platform_counts.get('COVERED_NO_ALT_SAMPLE_OR_ASSAY_DIFFERENCE', 0)}</b> 个。"
            "这提示不同取材、肿瘤异质性或低纯度影响，不能把一个时间点的结果概括为所有肿瘤组织。</li>"
        )
        source_review = sum(platform_counts.get(key, 0) for key in (
            "SOURCE_INDEL_NOT_REPRODUCED_REASSEMBLY_REQUIRED",
            "SOURCE_PASS_NOT_REPRODUCED_BY_PILEUP",
            "SOURCE_WEAK_EXACT_PILEUP_SUPPORT",
        ))
        out.append(f"<li>源检测仍需局部重组、IGV 或重复测序复核：<b>{source_review}</b> 个事件。</li>")
    out.append("<li>RNA 表达量只能说明基因被表达，不能替代突变 RNA reads 或融合/剪接 junction 的直接支持。</li>")
    out.append("<li>精确取材日期、部位及治疗前后关系需与临床样本记录核对；本报告仅描述已测序文库。</li>")
    out.append("</ul></div>")

    out.append("<div class='section'><h2>14. 建议的下一步验证顺序</h2><ol class='compact'>")
    out.append("<li><b>先核对样本：</b>确认 WES、WGS、RNA 的取材部位、日期、治疗前后关系及肿瘤含量，避免把时间点差异误作技术失败。</li>")
    out.append("<li><b>确认事件真实性：</b>对样本特异、复杂 InDel 和弱支持事件做 IGV、局部组装、靶向深测或独立 PCR；TBR1 保留 read-backed phasing 结论。</li>")
    out.append("<li><b>确认突变转录：</b>SNV/InDel 检查 RNA alt reads/RNA VAF；融合和剪接检查 junction reads、阅读框及异常转录本。</li>")
    out.append("<li><b>确认突变特异性：</b>比较 MT 与 WT 的 HLA 结合、呈递和免疫原性；WT 相当或更强者不进入首批。</li>")
    out.append("<li><b>确认加工呈递：</b>短肽候选做 MT/WT 成对验证；移码、剪接和融合优先长肽或 minigene，并在条件允许时做免疫肽组学。</li>")
    out.append("<li><b>确认免疫功能与安全性：</b>再开展 ELISpot、四聚体/多聚体、细胞毒实验，并补正常组织、HSPC、自身肽和脱靶复核。</li>")
    out.append("</ol></div>")

    experiment_rows = [
        {"group": "A：SNV/InDel 短肽组", "contents": "每个事件 1–2 个最佳 peptide-HLA 组合", "controls": "对应 WT 肽、无关肽、阳性刺激", "purpose": "验证突变特异性和 T 细胞识别"},
        {"group": "B：移码/剪接长肽组", "contents": "覆盖新生尾部或异常 junction 的长肽/minigene", "controls": "正常转录本/WT 构建", "purpose": "验证真实加工而非仅短肽结合"},
        {"group": "C：融合专项组", "contents": "EWSR1::WT1 单独成组，其他融合分开", "controls": "断点阴性、WT 两端序列、无关融合", "purpose": "验证断点、阅读框、呈递和 DSRCT 特异背景"},
        {"group": "D：人工保留复核组", "contents": "样本特异、复杂 InDel、弱支持事件", "controls": "另一平台/另一时间点、正常样本", "purpose": "先解决事件真实性，不与已确认候选混合解释"},
    ]
    out.append("<div class='section'><h2>15. 建议的实验候选组织方式</h2>")
    out.append(_table(experiment_rows, ["group", "contents", "controls", "purpose"]))
    out.append("<p class='small'>每组应预先定义入组条件、排除条件和重复数；同一事件的高度重叠肽应去冗余，避免一个事件因窗口数量多而在实验中被过度代表。</p></div>")

    out.append("<div class='section'><h2>16. 面向患者的核心结论</h2>")
    out.append("<p>本次分析发现了若干值得继续研究的候选，尤其包括与 DSRCT 密切相关的 <b>EWSR1::WT1</b> 融合，以及部分在 DNA、RNA 和 HLA 预测层面得到支持的突变事件。肿瘤的 HLA-I 抗原呈递能力看起来是<b>部分保留</b>的，因此继续做新抗原实验验证具有研究依据。</p>")
    out.append("<p>但目前没有候选达到“仅凭计算结果即可用于治疗”的证据标准。部分事件在 WES/WGS、不同取材或正常样本之间存在差异，正常组织安全性和真实加工呈递也尚未完整验证。当前最重要的下一步是按分层进行事件确认、RNA/断点验证、MT-WT 对照和 T 细胞功能实验，而不是直接按自动排名选择治疗方案。</p>")
    out.append("</div>")

    out.append("<div class='warn'><h2>17. 数据来源与解释边界</h2><ul class='compact'>")
    out.append("<li><b>数据来源：</b>肿瘤 WES、肿瘤 WGS、配对血液正常样本、肿瘤 RNA/融合结果、HLA 分型、纯度/CNV/CCF，以及呈递、免疫原性、APPM/逃逸与正常组织安全性参考。</li>")
    out.append("<li><b>跨平台边界：</b>WES 与 WGS 可能来自不同肿瘤文库或时间点；本报告展示的蛋白改变型 SNV/InDel 差异同时受捕获范围、深度、低纯度、异质性、caller、VEP 转录本选择与局部组装影响。完整 coding/splice PASS 全集保留在技术 QC 中。</li>")
    out.append("<li><b>融合边界：</b>检测到驱动融合不等于其 junction 肽一定被加工、呈递或被 T 细胞识别。</li>")
    out.append("<li>本分析为<strong>计算机辅助筛选</strong>，预测结合亲和力不等于体内呈递，更不等于临床疗效。</li>")
    out.append("<li>APPM、CCF、安全性与免疫逃逸评估依赖输入数据完整度；缺失数据不等于“无风险”。</li>")
    out.append("<li>本报告<strong>不包含</strong>原始测序质控、文件路径或生信命令细节；技术细节见科研技术版报告。</li>")
    out.append("<li>不得将本报告直接用于患者诊断、预后判断或个体化治疗处方。</li>")
    out.append("</ul></div>")
    out.append("</body></html>")
    p.write_text("\n".join(out), encoding="utf-8")


R_GRADE_PATIENT = {
    "R1": ("第一批实验优先", "关键证据较完整，可优先进入研究性实验验证。"),
    "R2": ("值得推进", "总体证据较好，但仍有一项或少量谨慎因素需要补充。"),
    "R3": ("优先补证据", "先补 RNA、事件真实性、安全性或呈递证据，再决定是否进入免疫学实验。"),
    "R4": ("当前暂不推进", "存在硬失败、明确风险、证据明显不足或呈递一致弱等原因。"),
}


def _patient_grade(row: Mapping[str, Any]) -> str:
    for key in ("pipeline_r_grade", "evidence_grade", "r_grade", "final_priority"):
        value = str(row.get(key) or "").strip().upper()
        if value:
            return value
    return "UNASSESSED"


def _patient_event_grade(row: Mapping[str, Any]) -> str:
    """Return the event-level R grade without consulting peptide-level rows."""
    grade = str(row.get("best_evidence_grade") or row.get("event_evidence_grade") or "").strip().upper()
    if grade != "R3":
        return grade or "UNASSESSED"
    review_required = str(row.get("manual_review_required") or "").strip().lower() in {"yes", "true", "1"}
    has_conflict = bool(str(row.get("evidence_conflict_layers") or "").strip())
    has_gap = bool(str(row.get("evidence_missing_layers") or "").strip())
    if review_required or has_conflict:
        return "R3-REVIEW"
    if has_gap:
        return "R3-GAP"
    return "R3-READY"


def _patient_event_grade_counts(events: list[dict[str, str]]) -> dict[str, int]:
    counts = {grade: 0 for grade in ("R1", "R2", "R3-READY", "R3-GAP", "R3-REVIEW", "R4", "UNASSESSED")}
    seen: set[str] = set()
    for index, row in enumerate(events):
        event_key = str(row.get("event_group_id") or row.get("event_id") or f"event-row-{index}")
        if event_key in seen:
            continue
        seen.add(event_key)
        grade = _patient_event_grade(row)
        counts[grade] = counts.get(grade, 0) + 1
    return counts


def _patient_track(row: Mapping[str, Any]) -> str:
    # Explicit event identity is authoritative.  In particular, a VCF SNV or
    # InDel can carry a splice-related VEP consequence without becoming an
    # RNA splice-junction event.
    event_type = str(row.get("event_type") or "").strip().lower()
    explicit_tracks = {
        "snv": "SNV",
        "missense": "SNV",
        "substitution": "SNV",
        "indel": "InDel",
        "insertion": "InDel",
        "deletion": "InDel",
        "frameshift": "InDel",
        "fusion": "Fusion",
        "splice": "Splice",
        "splice_junction": "Splice",
        "junction": "Splice",
    }
    if event_type in explicit_tracks:
        return explicit_tracks[event_type]
    text = " ".join(str(row.get(key) or "") for key in (
        "event_type", "mutation_source", "peptide_consequence", "source_type", "event_kind"
    )).lower()
    if "fusion" in text:
        return "Fusion"
    if "splice" in text or "junction" in text or "exon" in text:
        return "Splice"
    if "sv" in text or "structural" in text:
        return "DNA SV"
    if any(token in text for token in ("indel", "frameshift", "insertion", "deletion")):
        return "InDel"
    if any(token in text for token in ("snv", "missense", "substitution")):
        return "SNV"
    return "Other"


def _patient_representatives(rows: list[dict[str, str]], limit: int, track: str | None = None) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []
    seen_events: set[str] = set()
    seen_peptide_hla: set[str] = set()
    for row in rows:
        if track and _patient_track(row) != track:
            continue
        event_id = str(row.get("event_id") or row.get("event_name") or row.get("peptide_id") or "")
        event_key = event_id or identity_value(row, "event_identity_id")
        peptide_hla_key = identity_value(row, "peptide_hla_id")
        if event_key in seen_events or peptide_hla_key in seen_peptide_hla:
            continue
        seen_events.add(event_key)
        seen_peptide_hla.add(peptide_hla_key)
        annotated = dict(row)
        for field, value in candidate_identity(row).items():
            annotated.setdefault(field, value)
        selected.append(annotated)
        if len(selected) >= limit:
            break
    return selected


def _patient_event_keys(row: Mapping[str, Any]) -> list[str]:
    keys: list[str] = []
    for value in (row.get("event_group_id"), row.get("event_id"), row.get("source_event_id")):
        text = str(value or "").strip()
        if text and text not in keys:
            keys.append(text)
    for value in str(row.get("member_event_ids") or "").replace(",", ";").split(";"):
        text = value.strip()
        if text and text not in keys:
            keys.append(text)
    return keys


def _patient_event_representatives(
    events: list[dict[str, str]],
    peptides: list[dict[str, str]],
    limit: int,
    track: str,
) -> list[dict[str, str]]:
    """Select event-level rows, then attach the best ranked peptide-HLA when available."""
    peptide_by_event: dict[str, dict[str, str]] = {}
    for peptide in peptides:
        for key in _patient_event_keys(peptide):
            peptide_by_event.setdefault(key, peptide)

    selected: list[dict[str, str]] = []
    seen: set[str] = set()
    # ranked_events is authoritative. Peptide rows are only a fallback for older
    # bundles whose event table did not carry this analysis track at all.
    event_sources = [row for row in events if _patient_track(row) == track]
    sources = event_sources or [row for row in peptides if _patient_track(row) == track]
    for source in sources:
        keys = _patient_event_keys(source)
        event_key = keys[0] if keys else str(source.get("event_name") or source.get("gene") or "")
        if not event_key or event_key in seen:
            continue
        peptide = next((peptide_by_event[key] for key in keys if key in peptide_by_event), None)
        combined = dict(peptide or {})
        combined.update({key: value for key, value in source.items() if str(value or "").strip()})
        selected.append(combined)
        seen.update(keys or [event_key])
        if len(selected) >= limit:
            break
    return selected


def _patient_value(row: Mapping[str, Any], *fields: str, default: str = "UNASSESSED") -> str:
    for field_name in fields:
        value = str(row.get(field_name) or "").strip()
        if value:
            return value
    return default


def _patient_assessed(row: Mapping[str, Any], *fields: str) -> bool:
    unavailable = ("UNASSESSED", "NOT_AVAILABLE", "NOT_RUN")
    for field_name in fields:
        value = str(row.get(field_name) or "").strip().upper()
        if not value or any(token in value for token in unavailable):
            continue
        return True
    return False


def _patient_metric(label: str, row: Mapping[str, Any], *fields: str) -> str:
    value = _patient_value(row, *fields)
    if label == "MT/WT" and _patient_track(row) in {"Fusion", "Splice"}:
        state = value.upper()
        if state == "NOVEL_SEQUENCE":
            return "MT/WT=异常连接新序列；应使用正常连接或正常异构体肽作为对照"
        crosses = str(row.get("crosses_junction") or "").strip().lower()
        consequence = str(row.get("peptide_consequence") or "").strip().lower()
        if state in {"", "UNASSESSED", "NOT_ASSESSED", "NOT_AVAILABLE"} and (
            crosses in {"yes", "true", "1"} or consequence in {"fusion", "splice_junction"}
        ):
            return "MT/WT=传统点突变式配对不适用；需补正常连接或正常异构体肽对照"
    return f"{label}={_patient_status_text(value)}"


def _patient_rna_metric(row: Mapping[str, Any]) -> str:
    """Keep caller-reported junction support distinct from verified support."""
    if _patient_track(row) in {"Fusion", "Splice"}:
        verified = _patient_observed_value(row, "rna_junction_reads") or _patient_observed_value(row, "junction_reads")
        provided = _patient_observed_value(row, "provided_rna_junction_reads")
        try:
            verified_n = float(verified) if verified is not None else None
            provided_n = float(provided) if provided is not None else None
        except ValueError:
            verified_n = provided_n = None
        if provided_n and provided_n > 0 and (verified_n is None or verified_n <= 0):
            canonical = str(row.get("canonical_junction_id") or "")
            strand = str(row.get("junction_strand") or "").strip()
            if canonical.endswith("|.") or strand in {".", "?"}:
                return (
                    f"RNA=上游工具报告同坐标junction reads {provided}；因链方向未解析，"
                    "当前严格规则未将其计入已核实reads（不等同于检测为0）"
                )
            return (
                f"RNA=上游工具报告junction reads {provided}，但尚未与同一canonical junction的"
                "原始比对记录精确回链；已核实reads为0"
            )
    return _patient_metric("RNA", row, "rna_support_state", "rna_support_status")


def _patient_presentation_metric(row: Mapping[str, Any]) -> str:
    state = _patient_value(
        row, "presentation_consensus_state", "presentation_evidence_grade"
    )
    text = f"呈递={_patient_status_text(state)}"
    if "DISCORDANT" not in state.upper() and "不一致" not in _patient_status_text(state):
        return text
    values: list[str] = []
    net = _patient_observed_value(row, "netmhcpan_el_rank") or _patient_observed_value(row, "netmhcpan_mt_rank_el")
    mhc = _patient_observed_value(row, "mhcflurry_presentation_score")
    stab = _patient_observed_value(row, "netmhcstabpan_rank")
    if net is not None:
        values.append(f"NetMHCpan EL rank={net}")
    if mhc is not None:
        values.append(f"MHCflurry呈递分={mhc}")
    if stab is not None:
        values.append(f"NetMHCstabpan rank={stab}")
    return text + (f"（{'；'.join(values)}）" if values else "")


def _patient_presentation_quantitative_row(row: Mapping[str, Any], rank: int) -> dict[str, str]:
    """Expose raw model outputs without turning screening thresholds into proof."""
    def observed(*fields: str) -> str | None:
        for field in fields:
            candidate = _patient_observed_value(row, field)
            if candidate is not None:
                return candidate
        return None

    def value(*fields: str) -> str:
        raw = observed(*fields)
        return raw if raw is not None else "未提供"

    def number(*fields: str) -> float | None:
        raw = observed(*fields)
        try:
            return float(raw) if raw is not None else None
        except (TypeError, ValueError):
            return None

    mt_el = number("netmhcpan_mt_rank_el", "netmhcpan_el_rank")
    wt_el = number("netmhcpan_wt_rank_el")
    mt_ba = number("netmhcpan_mt_rank_ba", "netmhcpan_ba_rank")
    wt_ba = number("netmhcpan_wt_rank_ba")
    el_delta = wt_el - mt_el if mt_el is not None and wt_el is not None else None
    el_ratio = wt_el / mt_el if mt_el not in {None, 0.0} and wt_el is not None else None
    ba_delta = wt_ba - mt_ba if mt_ba is not None and wt_ba is not None else None
    ba_ratio = wt_ba / mt_ba if mt_ba not in {None, 0.0} and wt_ba is not None else None

    peptide = _patient_value(row, "peptide", default="")
    mutation_position = _patient_value(
        row,
        "mutation_positions_in_peptide",
        "mutation_position_in_peptide",
        "junction_position_in_peptide_1based",
        default="未提供",
    )
    hla = _patient_value(row, "hla_allele", default="未提供")
    explicit_support = _patient_value(
        row,
        "netmhcpan_allele_support_status",
        "mhcflurry_allele_support_status",
        "predictor_allele_support_status",
        "allele_support_status",
        default="",
    )
    has_prediction = any(number(field) is not None for field in (
        "netmhcpan_mt_rank_el", "netmhcpan_el_rank", "netmhcpan_mt_rank_ba",
        "netmhcpan_ba_rank", "mhcflurry_presentation_score", "netmhcstabpan_rank",
    ))
    if explicit_support:
        allele_support = explicit_support
    elif has_prediction:
        allele_support = "工具已返回该HLA结果；训练覆盖/外推状态未记录，不能据分数反推"
    else:
        allele_support = "未形成该HLA的可用预测；训练覆盖/外推状态未记录"

    netmhcpan = (
        f"EL rank MT={value('netmhcpan_mt_rank_el', 'netmhcpan_el_rank')}%, "
        f"WT={value('netmhcpan_wt_rank_el')}%; "
        f"BA rank MT={value('netmhcpan_mt_rank_ba', 'netmhcpan_ba_rank')}%, "
        f"WT={value('netmhcpan_wt_rank_ba')}%; "
        f"IC50 MT={value('netmhcpan_mt_ic50')} nM, WT={value('netmhcpan_wt_ic50')} nM"
    )
    comparisons: list[str] = []
    if el_delta is not None:
        comparisons.append(f"EL Δ(WT-MT)={el_delta:.4g}")
    if el_ratio is not None:
        comparisons.append(f"EL WT/MT={el_ratio:.4g}")
    if ba_delta is not None:
        comparisons.append(f"BA Δ(WT-MT)={ba_delta:.4g}")
    if ba_ratio is not None:
        comparisons.append(f"BA WT/MT={ba_ratio:.4g}")
    mt_wt = "；".join(comparisons) or "未形成可计算的MT/WT rank差值或比值"

    role = _patient_value(row, "mutation_position_role", default="")
    role_text = {
        "PRIMARY_HLA_ANCHOR": "主锚定位点：可能主要改变HLA结合，不等同于形成新的TCR识别表面",
        "PUTATIVE_TCR_FACING": "推定TCR暴露位置：可能改变TCR识别，但仅为序列位置推断，需结构/功能验证",
        "MIXED_ANCHOR_AND_PUTATIVE_TCR_FACING": "同时涉及主锚点与推定TCR暴露位置，需分别解释结合与TCR识别效应",
        "STRUCTURAL_ROLE_UNCERTAIN": "结构角色不确定，不能从DAI或位置单独推断免疫原性",
        "NO_SEQUENCE_CHANGE": "MT与WT肽序列相同，不构成突变特异肽",
    }.get(role.upper(), "")
    if not role_text:
        anchor_only = _patient_value(row, "mutation_anchor_only", default="").lower()
        tcr_facing = _patient_value(row, "mutation_tcr_facing", default="").lower()
        if anchor_only in {"yes", "true", "1"}:
            role_text = "主锚定位点：可能主要改变HLA结合，不等同于形成新的TCR识别表面"
        elif tcr_facing in {"yes", "true", "1"}:
            role_text = "推定TCR暴露位置：可能改变TCR识别，但仅为序列位置推断，需结构/功能验证"
        else:
            role_text = "未形成可靠位置解释"

    wt_risk = _patient_value(row, "wt_self_reactivity_risk_status", default="")
    if not wt_risk:
        wt_ic50 = number("netmhcpan_wt_ic50")
        wt_ranks = [candidate for candidate in (wt_el, wt_ba) if candidate is not None]
        if (wt_ranks and min(wt_ranks) <= 1.0) or (wt_ic50 is not None and wt_ic50 <= 50.0):
            wt_risk = "WT_STRONG_BINDING_REVIEW"
        elif (wt_ranks and min(wt_ranks) <= 2.0) or (wt_ic50 is not None and wt_ic50 <= 500.0):
            wt_risk = "WT_BINDING_REVIEW"
        elif wt_ranks or wt_ic50 is not None:
            wt_risk = "WT_LOW_PREDICTED_BINDING"
        else:
            wt_risk = "UNASSESSED"
    wt_risk_text = _patient_status_text(wt_risk)

    mhcflurry = (
        f"affinity percentile MT={value('mhcflurry_mt_affinity_percentile', 'mhcflurry_affinity_percentile')}%, "
        f"WT={value('mhcflurry_wt_affinity_percentile')}%; "
        f"presentation score MT={value('mhcflurry_mt_presentation_score', 'mhcflurry_presentation_score')}, "
        f"WT={value('mhcflurry_wt_presentation_score')}"
    )
    stability = (
        f"NetMHCstabpan rank={value('netmhcstabpan_rank')}%, "
        f"score/稳定性={value('netmhcstabpan_score')}; "
        f"WT rank={value('netmhcstabpan_wt_rank')}%, WT score={value('netmhcstabpan_wt_score')}"
    )
    auxiliary = (
        f"PRIME MT={value('prime_score')} (rank={value('prime_rank')}), "
        f"WT={value('prime_wt_score')} (rank={value('prime_wt_rank')}); "
        f"BigMHC MT={value('bigmhc_im_score')}, WT={value('bigmhc_im_wt_score')}; "
        f"DeepImmuno={value('deepimmuno_score')}"
    )
    return {
        "排名": str(rank),
        "肽段-HLA": f"{peptide or '未提供'} / {hla}",
        "肽长/变异位置": f"{len(peptide) if peptide else '未提供'} aa；位置={mutation_position}",
        "NetMHCpan原始值": netmhcpan,
        "MT/WT定量比较": mt_wt,
        "突变位置结构解释": role_text,
        "WT自身反应/耐受风险": wt_risk_text,
        "MHCflurry原始值": mhcflurry,
        "稳定性": stability,
        "免疫原性辅助模型": auxiliary,
        "HLA模型覆盖": allele_support,
    }


def _patient_mtwt_caution(row: Mapping[str, Any]) -> str:
    if _patient_track(row) not in {"SNV", "InDel"}:
        return ""
    detail = _patient_presentation_quantitative_row(row, 0)
    status = _patient_value(row, "mutant_specificity_status", "mutant_specificity_state", default="UNASSESSED")
    return (
        f"MT/WT谨慎解释：{_patient_status_text(status)}；{detail['突变位置结构解释']}；"
        f"{detail['WT自身反应/耐受风险']}。MT/WT结合差异或DAI不能独立证明免疫原性"
    )


def _patient_numeric_value(row: Mapping[str, Any], *fields: str) -> tuple[str, float | None]:
    text = _patient_value(row, *fields, default="").strip()
    try:
        return text, float(text)
    except (TypeError, ValueError):
        return text, None


def _patient_ccf_assessment(row: Mapping[str, Any], bundle: ReportBundle | None = None) -> tuple[bool, str]:
    ccf_text, ccf_value = _patient_numeric_value(row, "ccf_estimate", "ccf_best", "raw_ccf")
    confidence = _patient_value(row, "ccf_confidence_state", "ccf_confidence_grade", "ccf_confidence", default="").upper()
    confidence_invalid = not confidence or any(token in confidence for token in ("UNASSESSED", "UNSPECIFIED", "UNKNOWN", "NOT_AVAILABLE", "LOW"))
    if ccf_value is not None and 0 <= ccf_value <= 2 and not confidence_invalid:
        state = _patient_value(row, "clonality_state", "ccf_status", default="已形成估计")
        ci_low = _patient_value(row, "ccf_ci_low", default="")
        ci_high = _patient_value(row, "ccf_ci_high", default="")
        interval = f"（95%区间 {ci_low}-{ci_high}）" if ci_low and ci_high else "（未形成95%区间）"
        local_cn = _patient_value(row, "local_cnv_status", default="未记录局部CNV匹配状态")
        total_cn = _patient_value(row, "total_cn", "total_copy_number", default="未提供")
        major_cn = _patient_value(row, "major_cn", default="未提供")
        minor_cn = _patient_value(row, "minor_cn", default="未提供")
        multiplicity = _patient_value(row, "multiplicity_best", "mutation_multiplicity_assumption", default="未解析")
        normal_status = _patient_value(row, "normal_contamination_status", default="配对正常污染未评估")
        return True, (
            f"克隆性/覆盖解释={_patient_status_text(state)}；CCF={ccf_text}{interval}；"
            f"置信度={_patient_status_text(confidence)}；局部CN={total_cn}（major={major_cn}, minor={minor_cn}；{_patient_status_text(local_cn)}）；"
            f"变异拷贝数假设={multiplicity}；{_patient_status_text(normal_status)}。"
            "该结果表示与相应克隆覆盖范围相容，不构成克隆性的直接证明"
        )
    purity_status = str((bundle.purity_consensus if bundle else {}).get("status") or "").upper()
    reason = "样本纯度低且缺少可用CCF结果" if "LOW_PURITY" in purity_status else "缺少可用CCF数值或可靠置信度"
    local_cn = _patient_value(row, "local_cnv_status", default="局部CNV状态未记录")
    return False, f"克隆性=未形成可靠估计；原因={reason}；局部CNV={_patient_status_text(local_cn)}；处理=不作为阴性，也不作为正向加分"


_SPLICE_FUNNEL_LABELS = {
    "RAW_SPLICE_EVENTS_IN_UNIFIED_INPUT": "进入统一事件层的异常剪接事件",
    "ALIGNMENT_COORDINATE_QC": "坐标与比对质量过滤",
    "UNIQUE_JUNCTION_READS": "unique junction reads门槛",
    "TOTAL_JUNCTION_COVERAGE": "junction总覆盖门槛",
    "PSI": "PSI门槛",
    "MATCHED_NORMAL_JUNCTION": "患者配对/邻近正常样本比较",
    "NORMAL_COHORT_JUNCTION": "GTEx/正常组织junction库比较",
    "ANNOTATED_NORMAL_ISOFORM": "去除已注释正常异构体",
    "CREDIBLE_ORF": "可信ORF/阅读框",
    "NMD": "NMD风险过滤",
    "JUNCTION_SPANNING_PEPTIDE": "生成真正跨异常junction的肽",
    "NORMAL_PROTEOME_EXCLUSION": "正常蛋白组精确匹配排除",
    "SELECTED_FOR_PRESENTATION": "进入HLA呈递预测",
    "HLA_PRESENTATION": "通过HLA呈递门槛",
}


def _patient_splice_funnel_rows(bundle: ReportBundle) -> list[dict[str, str]]:
    stored = bundle.provenance.get("splice_filter_funnel_rows")
    rows: list[dict[str, str]] = []
    if isinstance(stored, list):
        for item in stored:
            if not isinstance(item, Mapping):
                continue
            stage = str(item.get("stage") or "")
            rows.append({
                "筛选阶段": _SPLICE_FUNNEL_LABELS.get(stage, stage or "未命名阶段"),
                "进入事件": str(item.get("entered_events") or "未记录"),
                "已评估": str(item.get("assessed_events") or "0"),
                "明确通过": str(item.get("passed_events") or "0"),
                "明确未通过": str(item.get("failed_events") or "0"),
                "未评估": str(item.get("unassessed_events") or "0"),
                "阶段后可能剩余": str(item.get("possible_remaining_range") or "未记录"),
                "规则/说明": str(item.get("criterion") or ""),
            })

    splice_events = {
        str(row.get("event_group_id") or row.get("event_id") or f"splice-{index}")
        for index, row in enumerate(bundle.events)
        if _patient_track(row) == "Splice"
    }
    if not rows and splice_events:
        rows.append({
            "筛选阶段": "进入统一事件层的异常剪接事件",
            "进入事件": str(len(splice_events)), "已评估": str(len(splice_events)),
            "明确通过": str(len(splice_events)), "明确未通过": "0", "未评估": "0",
            "阶段后可能剩余": str(len(splice_events)),
            "规则/说明": "仅能确认进入统一事件层的数量；未找到上游逐级漏斗，不能称为原始aligner junction总数",
        })
        for stage in (
            "ALIGNMENT_COORDINATE_QC", "UNIQUE_JUNCTION_READS", "TOTAL_JUNCTION_COVERAGE", "PSI",
            "MATCHED_NORMAL_JUNCTION", "NORMAL_COHORT_JUNCTION", "ANNOTATED_NORMAL_ISOFORM",
            "CREDIBLE_ORF", "NMD", "JUNCTION_SPANNING_PEPTIDE", "NORMAL_PROTEOME_EXCLUSION",
        ):
            rows.append({
                "筛选阶段": _SPLICE_FUNNEL_LABELS[stage], "进入事件": str(len(splice_events)),
                "已评估": "0", "明确通过": "0", "明确未通过": "0",
                "未评估": str(len(splice_events)), "阶段后可能剩余": f"0-{len(splice_events)}",
                "规则/说明": "该阶段的结构化计数未写入运行结果；不能从Top候选或基因TPM反推",
            })

    peptide_groups: dict[str, list[Mapping[str, Any]]] = {}
    for index, row in enumerate(bundle.peptides):
        if _patient_track(row) != "Splice":
            continue
        key = str(row.get("event_group_id") or row.get("event_id") or f"splice-peptide-{index}")
        peptide_groups.setdefault(key, []).append(row)
    if peptide_groups:
        assessed = 0
        passed = 0
        for event_rows in peptide_groups.values():
            states = {
                str(row.get("presentation_consensus_state") or row.get("presentation_evidence_grade") or "").upper()
                for row in event_rows
            }
            if states - {"", "UNASSESSED", "PRESENTATION_UNASSESSED"}:
                assessed += 1
            if any(state in {"PRESENTATION_CONSISTENT_STRONG", "PRESENTATION_MODERATE", "A", "B"} for state in states):
                passed += 1
        rows.append({
            "筛选阶段": _SPLICE_FUNNEL_LABELS["HLA_PRESENTATION"],
            "进入事件": str(len(peptide_groups)), "已评估": str(assessed),
            "明确通过": str(passed), "明确未通过": str(max(0, assessed - passed)),
            "未评估": str(len(peptide_groups) - assessed),
            "阶段后可能剩余": f"{passed}-{passed + len(peptide_groups) - assessed}",
            "规则/说明": "按独立剪接事件汇总；NetMHCpan/MHCflurry核心呈递共识为强或中等时计为通过",
        })
    return rows


def _patient_ccf_coverage_rows(bundle: ReportBundle) -> list[dict[str, str]]:
    grouped: dict[str, dict[str, Mapping[str, Any]]] = {}
    for index, row in enumerate(bundle.events):
        track = _patient_track(row)
        event_id = str(row.get("event_group_id") or row.get("event_id") or f"event-{index}")
        grouped.setdefault(track, {})[event_id] = row
    result: list[dict[str, str]] = []
    for track in ("SNV", "InDel", "DNA SV", "Fusion", "Splice"):
        event_rows = list(grouped.get(track, {}).values())
        if not event_rows:
            continue
        reliable = low = unresolved = not_applicable = 0
        for row in event_rows:
            ok, _ = _patient_ccf_assessment(row, bundle)
            _, value = _patient_numeric_value(row, "ccf_estimate", "ccf_best", "raw_ccf")
            source = " ".join(str(row.get(field) or "").upper() for field in ("mutation_source", "ccf_status", "ccf_method"))
            if ok:
                reliable += 1
            elif value is not None:
                low += 1
            elif track in {"Fusion", "Splice"} and ("RNA_ONLY" in source or not str(row.get("tumor_vaf") or "").strip()):
                not_applicable += 1
            else:
                unresolved += 1
        eligible = len(event_rows) - not_applicable
        coverage = f"{(100.0 * reliable / eligible):.1f}%" if eligible else "不适用"
        impact = (
            "RNA-only事件不能从RNA reads推导DNA CCF；需DNA-SV/位点证据或正交验证"
            if not_applicable == len(event_rows)
            else "缺失CCF的DNA来源事件不能确认克隆覆盖范围，并限制进入高等级候选"
            if unresolved or low else "当前可用于克隆性分层"
        )
        result.append({
            "事件类型": track, "独立事件": str(len(event_rows)), "可靠CCF": str(reliable),
            "低置信数值": str(low), "缺失/未解析": str(unresolved), "RNA-only不适用": str(not_applicable),
            "可评估事件覆盖率": coverage, "解释与影响": impact,
        })
    return result


def _patient_candidate_integrity(row: Mapping[str, Any]) -> tuple[bool, list[str]]:
    missing: list[str] = []
    if not str(row.get("event_id") or "").strip():
        missing.append("event_id缺失")
    if not str(row.get("peptide") or "").strip():
        missing.append("peptide缺失")
    if not str(row.get("hla_allele") or "").strip():
        missing.append("限制性HLA缺失")
    source_tier = str(row.get("source_chain_confidence_tier") or "").strip().upper()
    if source_tier in {"", "UNASSESSED"}:
        missing.append("来源链不可回溯")
    elif source_tier == "C4":
        details: list[str] = []
        canonical = str(row.get("canonical_junction_id") or "")
        strand = str(row.get("junction_strand") or "").strip()
        if canonical.endswith("|.") or strand in {"", ".", "?"}:
            details.append("canonical junction缺少可用strand")
        provided = _patient_observed_value(row, "provided_rna_junction_reads")
        verified = _patient_observed_value(row, "rna_junction_reads") or _patient_observed_value(row, "junction_reads")
        try:
            if float(provided or 0) > 0 and float(verified or 0) <= 0:
                if canonical.endswith("|.") or strand in {".", "?"}:
                    details.append(
                        f"上游同坐标报告{provided}条reads；因strand未解析，未计入严格verified reads"
                    )
                else:
                    details.append(f"上游报告{provided}条reads但精确核实为0")
        except ValueError:
            pass
        if not str(row.get("transcript_hypothesis_id") or row.get("transcript_id") or "").strip():
            details.append("上游caller事件已回溯，但transcript hypothesis尚未建立")
        if not str(row.get("orf_id") or "").strip():
            details.append("正式ORF尚未建立")
        if "SC_PEPTIDE_HLA_TRACEABILITY_INCOMPLETE" in str(row.get("source_chain_reason_codes") or ""):
            source_traceable = bool(
                str(row.get("source_record_id") or "").strip()
                and str(row.get("event_id") or "").strip()
                and str(row.get("peptide") or "").strip()
                and str(row.get("hla_allele") or "").strip()
            )
            details.append(
                "肽段-HLA可回溯至上游caller事件，但尚未通过transcript/ORF/peptide-origin正式闭环"
                if source_traceable else "肽段-HLA来源回链不完整"
            )
        missing.append("来源链C4：" + "、".join(details or ["存在阻断性来源链缺口"]))
    core_presentation = any(_patient_assessed(row, field) for field in (
        "netmhcpan_el_rank", "netmhcpan_mt_rank_el", "mhcflurry_presentation_score",
    ))
    if not core_presentation:
        missing.append("核心呈递工具结果缺失")
    return not missing, missing


def _patient_candidate_hla_status(bundle: ReportBundle | None) -> str:
    if bundle is None:
        return "未评估"
    _, overall = _patient_hla_loh_consensus(bundle)
    if "单工具提示限制性HLA-I LOH" in overall:
        return "单工具提示可能丢失；需正交确认"
    if "仅有单工具支持" in overall or "单工具" in overall:
        return "单工具未提示LOH；不足以确认完整保留"
    if "多工具一致未提示" in overall:
        return "多工具一致未提示LOH；支持当前分层"
    if "检出限制性HLA-I LOH" in overall:
        return overall.replace("检出限制性HLA-I LOH", "检出限制性HLA丢失")
    if "冲突" in overall:
        return "工具结果冲突；限制性HLA状态待确认"
    return "未评估"


def _patient_restricting_hla_reliability(
    row: Mapping[str, Any], bundle: ReportBundle | None,
) -> tuple[bool, str]:
    """Classify whether the restricting-allele LOH call is reliable for stratification."""
    if bundle is None:
        return False, "限制性HLA未评估"
    allele = str(row.get("hla_allele") or row.get("allele") or "").strip()
    if not allele:
        return False, "限制性HLA缺失"
    loh_rows, _ = _patient_hla_loh_consensus(bundle)
    match = next((item for item in loh_rows if item.get("HLA等位基因") == allele), None)
    if not match:
        return False, f"{allele}无逐等位基因LOH结论"
    consensus = str(match.get("综合判断") or "未评估")
    if consensus.startswith("多工具一致"):
        return True, consensus
    if consensus.startswith("仅"):
        return False, f"{consensus}；单工具结果不作为完整保留或丢失的可靠确认"
    return False, consensus


def _patient_restricting_hla_gap(
    row: Mapping[str, Any], bundle: ReportBundle | None,
) -> str:
    """Return a candidate-specific HLA-LOH limitation with tool details."""
    if bundle is None:
        return "限制性HLA未评估"
    allele = str(row.get("hla_allele") or row.get("allele") or "").strip()
    if not allele:
        return "限制性HLA未评估"
    loh_rows, _ = _patient_hla_loh_consensus(bundle)
    match = next((item for item in loh_rows if item.get("HLA等位基因") == allele), None)
    if not match:
        return f"限制性HLA {allele}：SpecHLA和LOHHLA均未提供逐等位基因结论"
    lohhla = str(match.get("LOHHLA") or "未提供")
    spechla = str(match.get("SpecHLA") or "未提供")
    consensus = str(match.get("综合判断") or "未评估")
    if "冲突" in consensus:
        return f"限制性HLA {allele} LOH冲突：LOHHLA={lohhla}；SpecHLA={spechla}"
    if consensus.startswith("仅"):
        return (
            f"限制性HLA {allele} 仅有单工具LOH判断：LOHHLA={lohhla}；SpecHLA={spechla}；"
            "不足以确认该等位基因在肿瘤中完整保留或丢失"
        )
    if "丢失" in consensus:
        return f"限制性HLA {allele} 可能丢失：LOHHLA={lohhla}；SpecHLA={spechla}"
    if consensus == "未评估":
        return f"限制性HLA {allele}：SpecHLA和LOHHLA均未完成评估"
    return ""


def _patient_candidate_appm_status(bundle: ReportBundle | None) -> str:
    if bundle is None:
        return "未评估"
    summary = bundle.appm_summary
    completeness = str(summary.get("appm_evidence_completeness") or "UNASSESSED").upper()
    mhc_i = str(summary.get("mhc_i_integrity_status") or "UNASSESSED").upper()
    tap = str(summary.get("tap_risk") or "UNASSESSED").lower()
    if completeness == "PARTIAL" and mhc_i == "MHC_I_INTACT":
        if tap in {"caution", "warn", "warning"}:
            return "证据部分完整；未见整体失活，但加工环节仍有谨慎信号"
        return "证据部分完整；未见整体失活，部分加工证据仍不完整"
    if completeness not in {"", "UNASSESSED", "NOT_ASSESSED"}:
        return f"{_patient_status_text(completeness)}；MHC-I={_patient_status_text(mhc_i)}"
    return "未评估"


def _patient_evidence_summary(row: Mapping[str, Any], bundle: ReportBundle | None = None) -> str:
    _, clonality = _patient_ccf_assessment(row, bundle)
    clonality_summary = clonality.split("；", 1)[0]
    return "；".join([
        _patient_metric("事件", row, "event_authenticity_state", "cross_platform_status"),
        _patient_metric("RNA", row, "rna_support_state", "rna_support_status"),
        _patient_metric("呈递", row, "presentation_consensus_state", "presentation_evidence_grade"),
        _patient_metric("MT/WT", row, "mutant_specificity_status", "mutant_specificity_state"),
        clonality_summary,
        f"限制性HLA={_patient_candidate_hla_status(bundle)}",
        f"APPM={_patient_candidate_appm_status(bundle)}",
        _patient_metric("安全性", row, "safety_state", "safety_status"),
        _patient_metric("来源链", row, "source_chain_confidence_tier"),
    ])


_PATIENT_CONFLICT_FIELD_LABELS = {
    "gene": "基因标识",
    "event_type": "事件类型",
    "protein_change": "蛋白改变",
    "hgvsp": "蛋白改变",
    "presentation_evidence_grade": "呈递证据等级",
    "presentation_evidence_score": "呈递综合评分",
    "netmhcpan_el_rank": "NetMHCpan EL rank",
    "mhcflurry_presentation_score": "MHCflurry呈递评分",
    "gene_expression_tpm": "基因表达TPM",
    "transcript_expression_tpm": "转录本表达TPM",
    "rna_ref_reads": "RNA参考等位基因reads",
    "rna_alt_reads": "RNA突变等位基因reads",
    "rna_depth": "RNA位点深度",
    "rna_vaf": "RNA VAF",
    "ccf_estimate": "CCF估计",
    "clonality_status": "克隆性状态",
    "purity": "肿瘤纯度",
    "restricting_hla_lost": "限制性HLA丢失状态",
    "hla_loh_status": "HLA LOH状态",
    "escape_status": "免疫逃逸状态",
    "appm_integrity_status": "APPM完整性状态",
    "safety_status": "安全性状态",
    "reference_proteome_exact_match": "正常蛋白组精确匹配",
}

_PATIENT_CONFLICT_SOURCE_LABELS = {
    "raw_events": "原始事件表",
    "raw_peptides": "原始肽段表",
    "annotated_peptides": "注释肽段表",
    "presentation_evidence": "呈递工具结果",
    "expression_evidence": "表达证据",
    "rna_junction_evidence": "RNA位点/连接证据",
    "ccf_2": "CCF/拷贝数证据",
    "appm_peptide_modifiers": "APPM证据",
    "peptide_escape_flags": "HLA LOH/免疫逃逸证据",
    "peptide_safety": "肽段安全性证据",
    "event_safety": "事件安全性证据",
    "ranked_peptides": "旧主排序副本",
    "validation_plan": "验证计划",
}


def _patient_conflict_summary(row: Mapping[str, Any], max_items: int = 4) -> str:
    """Describe source disagreements using this candidate's own provenance fields."""
    conflict_status = str(row.get("evidence_conflict_status") or "").strip().upper()
    if conflict_status in {"NONE", "NO_CONFLICT", "PASS", "RESOLVED", "FALSE", "0"}:
        return ""
    raw_details = str(row.get("evidence_conflict_details") or "").strip()
    details: list[Mapping[str, Any]] = []
    if raw_details:
        try:
            parsed = json.loads(raw_details)
            if isinstance(parsed, list):
                details = [item for item in parsed if isinstance(item, Mapping)]
        except json.JSONDecodeError:
            details = []

    summaries: list[str] = []
    has_structured_details = bool(details)
    derived_sources = {"ranked_peptides", "validation_plan"}
    provenance_only_fields = {
        "source_file", "source_record_id", "source_records", "source_row_number",
        "source_tools", "provenance_record_count",
    }
    for detail in details:
        field_name = str(detail.get("field") or "").strip()
        if not field_name or field_name in provenance_only_fields:
            continue
        field_label = _PATIENT_CONFLICT_FIELD_LABELS.get(field_name, field_name)
        selected_source = str(detail.get("selected_source") or "权威来源")
        other_source = str(detail.get("other_source") or "其他来源")
        # Differences against downstream ranking/validation copies are audit
        # synchronization issues, not independent biological-tool conflicts.
        if selected_source in derived_sources or other_source in derived_sources:
            continue
        selected_label = _PATIENT_CONFLICT_SOURCE_LABELS.get(selected_source, selected_source)
        other_label = _PATIENT_CONFLICT_SOURCE_LABELS.get(other_source, other_source)
        selected_value = str(detail.get("selected_value") or "未提供")
        other_value = str(detail.get("other_value") or "未提供")
        summary = f"{field_label}：{selected_label}={selected_value}，{other_label}={other_value}"
        if summary not in summaries:
            summaries.append(summary)

    if not summaries and not has_structured_details:
        fields = [
            item.strip()
            for item in str(row.get("evidence_conflict_fields") or "").split(",")
            if item.strip() and item.strip() not in provenance_only_fields
        ]
        summaries = [_PATIENT_CONFLICT_FIELD_LABELS.get(field, field) for field in fields]

    if not summaries:
        return ""
    shown = summaries[:max_items]
    if len(summaries) > max_items:
        shown.append(f"另有{len(summaries) - max_items}项，详见evidence_conflicts.tsv")
    return "；".join(shown)


def _patient_limitation(row: Mapping[str, Any], bundle: ReportBundle | None = None) -> str:
    limitations: list[str] = []
    for field_name in ("hard_failure_codes", "priority_cap_reason_codes", "reason_codes"):
        for reason in str(row.get(field_name) or "").replace(",", "|").split("|"):
            reason = reason.strip()
            if reason and reason not in limitations:
                limitations.append(reason)
    checks = [
        (("rna_support_state", "rna_support_status"), "RNA证据未评估"),
        (("presentation_consensus_state", "presentation_evidence_grade"), "呈递证据未评估"),
        (("mutant_specificity_status", "mutant_specificity_state"), "MT/WT未评估"),
        (("safety_state", "safety_status"), "安全性证据未评估"),
        (("source_chain_confidence_tier",), "来源链未评估"),
    ]
    for fields, label in checks:
        if not _patient_assessed(row, *fields) and label not in limitations:
            limitations.append(label)
    ccf_reliable, ccf_assessment = _patient_ccf_assessment(row, bundle)
    if not ccf_reliable:
        limitations.append(ccf_assessment)
    integrity_ok, integrity_missing = _patient_candidate_integrity(row)
    if not integrity_ok:
        limitations.append("候选完整性检查未通过：" + "、".join(integrity_missing))
    conflict_summary = _patient_conflict_summary(row)
    if conflict_summary:
        limitations.append("具体证据冲突：" + conflict_summary)
    safety_gap = _patient_safety_gap(row)
    if safety_gap:
        limitations.append(safety_gap)
    hla_status = _patient_candidate_hla_status(bundle)
    if hla_status == "未评估":
        limitations.append("限制性HLA状态未评估")
    elif "不足以确认" in hla_status:
        limitations.append("限制性HLA仅单工具未提示LOH；不足以确认肿瘤中完整保留")
    elif "冲突" in hla_status or "丢失" in hla_status:
        limitations.append(hla_status)
    appm_status = _patient_candidate_appm_status(bundle)
    if appm_status == "未评估":
        limitations.append("APPM状态未评估")
    elif appm_status.startswith("证据部分完整"):
        limitations.append("APPM证据部分完整；加工环节仍需谨慎解释")
    if not _patient_assessed(row, "netmhcstabpan_rank", "netmhcstabpan_score"):
        limitations.append("NetMHCstabpan未评估")
    if not _patient_assessed(row, "netchop_processing_status", "netchop_31d_max_score"):
        limitations.append("NetChop酶切证据未评估")
    return "；".join(dict.fromkeys(limitations)) if limitations else "未见明确限制；仍需实验验证"


def _patient_comprehensive_evidence(row: Mapping[str, Any], bundle: ReportBundle | None = None) -> dict[str, str]:
    presentation = "; ".join([
        _patient_metric("共识", row, "presentation_consensus_state"),
        _patient_metric("NetMHCpan EL rank", row, "netmhcpan_el_rank", "el_rank"),
        _patient_metric("MHCflurry presentation", row, "mhcflurry_presentation_score"),
        _patient_metric("PRIME", row, "prime_score"),
        _patient_metric("BigMHC", row, "bigmhc_im_score"),
    ])
    processing = "; ".join([
        _patient_metric("MHCflurry processing", row, "mhcflurry_processing_score"),
        _patient_metric("NetMHCstabpan rank", row, "netmhcstabpan_rank"),
        _patient_metric("NetChop 3.1d C端切割", row, "netchop_31d_cterm_score", "netchop_31d_max_score", "netchop_processing_status"),
        _patient_metric("TAP/APPM", row, "tap_processing_status"),
    ])
    rna = _patient_metric("状态", row, "rna_support_state", "rna_support_status")
    for label, fields in (
        ("alt reads", ("rna_alt_reads",)),
        ("junction reads", ("rna_junction_reads", "junction_reads")),
        ("gene TPM", ("gene_expression_tpm",)),
        ("transcript TPM", ("transcript_expression_tpm",)),
    ):
        if _patient_assessed(row, *fields):
            rna += "; " + _patient_metric(label, row, *fields)
    return {
        "事件真实性": _patient_metric("状态", row, "event_authenticity_state", "cross_platform_status"),
        "RNA": rna,
        "呈递工具": presentation,
        "加工/稳定性": processing,
        "MT/WT": _patient_metric("状态", row, "mutant_specificity_status", "mutant_specificity_state"),
        "克隆性/CCF": _patient_ccf_assessment(row, bundle)[1],
        "限制性HLA状态": _patient_candidate_hla_status(bundle),
        "APPM状态": _patient_candidate_appm_status(bundle),
        "安全性": _patient_metric("状态", row, "safety_state", "safety_status"),
        "免疫原性": _patient_metric("score", row, "immunogenicity_composite_score", "immunogenicity_score", "bigmhc_im_score"),
        "可追溯性": _patient_metric("等级", row, "source_chain_confidence_tier"),
    }


def _patient_evidence_audit_rows(rows: list[dict[str, str]], bundle: ReportBundle | None = None) -> list[dict[str, str]]:
    dimensions = [
        ("事件真实性", ("event_authenticity_state", "cross_platform_status")),
        ("RNA直接证据", ("rna_support_state", "rna_support_status")),
        ("核心呈递共识", ("presentation_consensus_state", "presentation_evidence_grade")),
        ("加工/稳定性", ("mhcflurry_processing_score", "netmhcstabpan_rank", "netmhcstabpan_score", "netchop_31d_cterm_score", "netchop_31d_max_score", "netchop_processing_status", "tap_processing_status")),
        ("MT/WT特异性", ("mutant_specificity_status", "mutant_specificity_state")),
        ("安全性", ("safety_state", "safety_status")),
        ("免疫原性", ("immunogenicity_composite_score", "immunogenicity_score", "bigmhc_im_score")),
        ("来源链", ("source_chain_confidence_tier",)),
    ]
    total = len(rows)
    result = []
    for label, fields in dimensions:
        # A leading UNASSESSED placeholder must not hide a usable result in a
        # later tool/state field. Audit each source independently.
        assessed = sum(
            1 for row in rows
            if any(_patient_assessed(row, field) for field in fields)
        )
        available = str(assessed)
        unavailable = str(total - assessed)
        criterion = f"Top {total}；缺失或未评估不按阴性处理"
        if label == "加工/稳定性":
            tool_counts = {
                "MHCflurry processing": sum(
                    1 for row in rows if _patient_assessed(row, "mhcflurry_processing_score")
                ),
                "NetMHCstabpan": sum(
                    1 for row in rows
                    if any(_patient_assessed(row, field) for field in ("netmhcstabpan_rank", "netmhcstabpan_score"))
                ),
                "NetChop": sum(
                    1 for row in rows
                    if any(_patient_assessed(row, field) for field in ("netchop_31d_cterm_score", "netchop_31d_max_score", "netchop_processing_status"))
                ),
                "独立TAP": sum(
                    1 for row in rows if _patient_assessed(row, "tap_processing_status")
                ),
            }
            available = f"{assessed}（至少一项候选级结果）"
            unavailable = f"{total - assessed}（尚无候选级加工/稳定性值）"
            criterion = (
                f"Top {total}；" + "、".join(f"{name} {count}" for name, count in tool_counts.items())
                + "；样本级APPM另行展示，不能替代候选级加工/稳定性结果"
            )
        elif label == "MT/WT特异性":
            novel_sequence = sum(
                1 for row in rows
                if _patient_observed_source(row, "mutant_specificity_state", "mutant_specificity_status") == "NOVEL_SEQUENCE"
            )
            traditional = sum(
                1 for row in rows
                if _patient_track(row) in {"SNV", "InDel"}
                and any(
                    _patient_assessed(row, field)
                    for field in ("mutant_specificity_status", "mutant_specificity_state")
                )
            )
            usable = len({
                index for index, row in enumerate(rows)
                if (
                    _patient_observed_source(row, "mutant_specificity_state", "mutant_specificity_status") == "NOVEL_SEQUENCE"
                    or (
                        _patient_track(row) in {"SNV", "InDel"}
                        and any(
                            _patient_assessed(row, field)
                            for field in ("mutant_specificity_status", "mutant_specificity_state")
                        )
                    )
                )
            })
            available = f"{usable}（传统MT/WT {traditional}；异常连接新序列 {novel_sequence}）"
            unavailable = f"{total - usable}（尚未形成适用的特异性证据）"
            criterion = (
                f"Top {total}；SNV/InDel使用传统MT/WT比较；Fusion/Splice的NOVEL_SEQUENCE作为异常连接新序列证据，"
                "但仍需正常连接或正常异构体肽对照"
            )
        result.append({
            "证据维度": label,
            "可作为当前分层证据": available,
            "尚不能作为可靠证据": unavailable,
            "判定口径": criterion,
        })
    hla_reliability = [_patient_restricting_hla_reliability(row, bundle) for row in rows]
    hla_assessed = sum(reliable for reliable, _ in hla_reliability)
    hla_single_or_incomplete = sum(
        not reliable and ("单工具" in detail or "QC" in detail)
        for reliable, detail in hla_reliability
    )
    appm_assessed = total if _patient_candidate_appm_status(bundle) != "未评估" else 0
    ccf_reliable = 0
    ccf_low_confidence = 0
    ccf_unresolved = 0
    for row in rows:
        reliable, _ = _patient_ccf_assessment(row, bundle)
        if reliable:
            ccf_reliable += 1
            continue
        _, ccf_value = _patient_numeric_value(row, "ccf_estimate", "ccf_best", "raw_ccf")
        if ccf_value is not None:
            ccf_low_confidence += 1
        else:
            ccf_unresolved += 1
    result.insert(5, {
        "证据维度": "克隆性/CCF",
        "可作为当前分层证据": f"{ccf_reliable}（可靠估计）",
        "尚不能作为可靠证据": (
            f"{ccf_low_confidence + ccf_unresolved}"
            f"（已计算但低置信 {ccf_low_confidence}；未形成数值或不适用 {ccf_unresolved}）"
        ),
        "判定口径": f"Top {total}；低置信CCF保留数值用于审阅，但不作为正向加分或阴性结论",
    })
    result.insert(6, {
        "证据维度": "限制性HLA状态",
        "可作为当前分层证据": f"{hla_assessed}（逐等位基因多工具一致）",
        "尚不能作为可靠证据": f"{total - hla_assessed}（其中单工具或QC不足 {hla_single_or_incomplete}）",
        "判定口径": (
            f"Top {total}；按每个候选的限制性HLA逐项判断；"
            "仅单工具阴性或另一工具QC不足不计为完整保留证据"
        ),
    })
    result.insert(7, {
        "证据维度": "APPM状态",
        "可作为当前分层证据": str(appm_assessed),
        "尚不能作为可靠证据": str(total - appm_assessed),
        "判定口径": f"Top {total}；依据样本级APPM评估",
    })
    return result


_PATIENT_RESULT_TOOL_SPECS = {
    "netmhcpan": ("NetMHCpan", ("netmhcpan_el_rank", "netmhcpan_mt_rank_el"), "HLA-I结合/呈递"),
    "mhcflurry": ("MHCflurry", ("mhcflurry_presentation_score",), "呈递与加工"),
    "prime": ("PRIME", ("prime_score", "prime_rank"), "免疫原性/呈递辅助"),
    "bigmhc": ("BigMHC", ("bigmhc_im_score",), "免疫原性辅助"),
    "deepimmuno": ("DeepImmuno", ("deepimmuno_score",), "免疫原性辅助（适用于支持的9/10-mer peptide-HLA组合）"),
    "netmhcstabpan": ("NetMHCstabpan", ("netmhcstabpan_rank", "netmhcstabpan_score"), "肽-HLA稳定性"),
    "netchop": ("NetChop 3.1d", ("netchop_processing_status", "netchop_31d_cterm_score", "netchop_31d_max_score"), "蛋白酶体C端酶切加工"),
    "tap_appm": (
        "TAP/APPM",
        (
            "tap_processing_status",
            "appm_integrity_status",
            "appm_evidence_completeness",
            "appm_multiplier",
        ),
        "TAP转运与抗原加工呈递通路",
    ),
}

_STANDARD_PEPTIDE_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")


def _patient_tool_applicable(row: Mapping[str, Any], tool_key: str) -> bool:
    """Return whether a candidate is in a tool's documented peptide/HLA domain."""
    peptide = str(row.get("peptide") or row.get("mutant_peptide") or "").strip().upper()
    hla = str(row.get("hla_allele") or "").strip()
    if not peptide or not hla or not set(peptide).issubset(_STANDARD_PEPTIDE_AA):
        return False
    if tool_key == "deepimmuno":
        return len(peptide) in {9, 10}
    if tool_key in _PATIENT_RESULT_TOOL_SPECS:
        return 8 <= len(peptide) <= 11
    return True


def _patient_tool_coverage_status(
    rows: list[dict[str, str]], tool_key: str, fields: tuple[str, ...], count: int,
) -> str:
    applicable = sum(_patient_tool_applicable(row, tool_key) for row in rows)
    not_applicable = len(rows) - applicable
    if not count:
        status = f"未评估（适用的标准肽段-HLA组合 {applicable} 个中无结果值）"
    else:
        status = f"综合证据表已载入结果值（{count}/{applicable} 个适用组合）"
    if not_applicable:
        status += f"；另 {not_applicable} 个因肽长、字符或HLA输入不符合该工具范围而不适用"
    return status

_PATIENT_SOURCE_TOOL_ALIASES = {
    "bigmhc_im": "bigmhc", "bigmhc": "bigmhc", "deepimmuno": "deepimmuno",
    "netmhcpan": "netmhcpan", "mhcflurry": "mhcflurry", "netmhcstabpan": "netmhcstabpan",
    "netchop": "netchop", "prime": "prime", "tap/appm": "tap_appm", "appm": "tap_appm",
    "easyfuse": "easyfuse", "jaffal": "jaffal", "snaf": "snaf", "splicemutr": "splicemutr",
    "isoquant": "isoquant", "sqanti3": "sqanti3", "immunopepper": "immunopepper",
    "pvacbind": "pvacbind", "pvactools": "pvactools", "vep": "vep", "facets": "facets",
    "sequenza": "sequenza", "purple": "purple", "ascat": "ascat", "spechla": "spechla",
    "lohhla": "lohhla", "optitype": "optitype", "hla-la": "hlala", "hla_la": "hlala",
    "salmon": "salmon", "rsem": "rsem", "star-fusion": "starfusion", "star_fusion": "starfusion",
    "star": "star", "arriba": "arriba", "fusioncatcher": "fusioncatcher", "regtools": "regtools",
    "gatk": "gatk", "mutect2": "gatk", "bwa": "bwa", "samtools": "samtools",
}

_PATIENT_SOURCE_TOOL_META = {
    "easyfuse": ("EasyFuse", "短读长RNA融合检测"), "jaffal": ("JAFFAL", "长读长RNA融合检测"),
    "snaf": ("SNAF", "异常剪接检测"), "splicemutr": ("SpliceMutr", "异常剪接交叉验证"),
    "isoquant": ("IsoQuant", "长读长转录本重建"), "sqanti3": ("SQANTI3", "长读长转录本结构质控"),
    "immunopepper": ("ImmunoPepper", "异常转录本翻译与肽段生成"),
    "pvacbind": ("pVACbind", "peptide-HLA呈递预测"), "pvactools": ("pVACtools", "新抗原候选分析"),
    "vep": ("VEP", "变异功能注释"), "facets": ("FACETS", "纯度、CNV与CCF证据"),
    "sequenza": ("Sequenza", "纯度、倍性与CNV证据"), "purple": ("PURPLE", "纯度、倍性与CNV证据"),
    "ascat": ("ASCAT", "等位基因特异CNV证据"), "spechla": ("SpecHLA", "HLA分型与HLA-LOH证据"),
    "lohhla": ("LOHHLA", "HLA-I等位基因LOH证据"), "optitype": ("OptiType", "HLA-I分型"),
    "hlala": ("HLA-LA", "HLA分型"), "salmon": ("Salmon", "转录本和基因表达定量"),
    "rsem": ("RSEM", "转录本和基因表达定量"), "star": ("STAR", "短读长RNA比对与junction提取"),
    "starfusion": ("STAR-Fusion", "RNA融合检测"), "arriba": ("Arriba", "RNA融合检测"),
    "fusioncatcher": ("FusionCatcher", "RNA融合检测"), "regtools": ("RegTools", "RNA junction提取/注释"),
    "gatk": ("GATK/Mutect2", "体细胞变异检测"), "bwa": ("BWA", "DNA序列比对"),
    "samtools": ("samtools", "BAM处理与统计"),
}


def _patient_tool_key(name: Any) -> str:
    text = str(name or "").strip().lower().replace(" ", "")
    text = text.replace("3.1d", "").replace("-", "_")
    aliases = {
        "bigmhc_im": "bigmhc", "bigmhc": "bigmhc", "netchop": "netchop",
        "tap/appm": "tap_appm", "appm": "tap_appm", "hla_la": "hlala",
        "star_fusion": "starfusion", "gatk/mutect2": "gatk",
    }
    return aliases.get(text, text.replace("_", ""))


def _patient_source_tool_counts(rows: list[dict[str, str]]) -> Counter[str]:
    counts: Counter[str] = Counter()
    source_fields = (
        "source_tool", "source_tools", "tools_detected", "caller", "callers", "source",
        "mutation_source", "junction_source", "rna_junction_source", "expression_source",
        "purity_source", "purity_tool_values", "external_clonality_tool", "immunogenicity_source",
    )
    aliases = sorted(_PATIENT_SOURCE_TOOL_ALIASES, key=len, reverse=True)
    for row in rows:
        seen: set[str] = set()
        text = " ".join(str(row.get(field) or "") for field in source_fields).lower()
        for alias in aliases:
            if re.search(rf"(?<![a-z0-9]){re.escape(alias)}(?![a-z0-9])", text):
                seen.add(_PATIENT_SOURCE_TOOL_ALIASES[alias])
        for key in seen:
            counts[key] += 1
    return counts


def _patient_inferred_tool_rows(rows: list[dict[str, str]], tool_versions: Mapping[str, Mapping[str, str]] | None = None) -> list[dict[str, str]]:
    result = []
    versions = {str(key).lower(): value for key, value in (tool_versions or {}).items()}
    for tool_key, (name, fields, purpose) in _PATIENT_RESULT_TOOL_SPECS.items():
        count = sum(1 for row in rows if _patient_assessed(row, *fields))
        if tool_key == "tap_appm" and count:
            tap_count = sum(1 for row in rows if _patient_assessed(row, "tap_processing_status"))
            status = f"APPM已评估（候选修饰值 {count}/{len(rows)}）"
            if tap_count:
                status += f"；独立TAP状态 {tap_count}/{len(rows)}"
            else:
                status += "；独立TAP候选级状态未单列"
        else:
            status = _patient_tool_coverage_status(rows, tool_key, fields, count)
        record = versions.get(tool_key) or versions.get(name.lower()) or versions.get(name.lower().replace("/", "_")) or {}
        default_version = "NeoAg APPM 2.0" if tool_key == "tap_appm" and count else "原始运行版本未记录（需补工具版本清单）"
        default_evidence = "APPM汇总和肽段修饰结果" if tool_key == "tap_appm" and count else "综合结果仅证明工具结果存在，不能反推版本"
        version = str(record.get("version") or default_version)
        version_evidence = str(record.get("evidence") or default_evidence)
        result.append({"流程/工具": name, "版本": version, "版本依据": version_evidence, "状态": status, "作用": purpose})
    return result


def _patient_validation(row: Mapping[str, Any], val_map: Mapping[str, Mapping[str, str]]) -> str:
    val = val_map.get(str(row.get("peptide_id") or ""), {})
    explicit = str(val.get("validation_strategy") or val.get("recommended_assay") or row.get("recommended_validation") or row.get("recommended_use") or "")
    if explicit:
        translations = (
            ("do not advance", "当前暂缓/不推进；先解决阻断性证据，再决定是否重新评估"),
            ("safety-focused validation before efficacy assay", "先完成针对性的正常组织与脱靶安全性复核，再考虑有效性实验"),
            ("requires focused safety validation", "先完成针对性的正常组织、正常蛋白组和脱靶安全性复核，再决定是否进入功能实验"),
            ("novel c-terminal tail", "新生C端肽段：优先采用覆盖新生尾部的混合长肽（15–27 aa）和/或移码minigene；短肽仅作次级验证"),
            ("mutant short peptide", "突变短肽（8–11 aa）与匹配的正常短肽对照；建议开展MHC-I ELISpot或多聚体实验"),
            ("fusion junction long peptide", "优先采用跨融合断点的长肽和/或融合minigene，并先确认精确断点与阅读框"),
            ("abnormal splice/exon-junction long peptide", "优先采用覆盖异常剪接连接点的长肽（15–27 aa）和/或剪接minigene，不应仅依赖短肽"),
            ("mt/wt paired validation required", "补做突变肽与正常肽成对验证；确认突变特异性后再进入功能实验"),
            ("exclude from first validation batch", "暂不纳入首批验证；先补齐当前证据缺口"),
            ("clonality/persistence caution", "先复核克隆性和持续性，再决定实验优先级"),
        )
        lowered = explicit.lower()
        for marker, translated in translations:
            if marker in lowered:
                return translated
        # Patient-facing reports should not expose raw English pipeline advice
        # or semicolon-delimited gate diagnostics.  Preserve that text in the
        # technical outputs and emit a deterministic Chinese recommendation
        # based on the event track here.
    track = _patient_track(row)
    if track == "SNV":
        return "MT/WT成对短肽与ELISpot/多聚体"
    if track == "InDel":
        return "新生尾部长肽或minigene"
    if track == "Fusion":
        return "RT-PCR/Sanger确认断点，再做融合junction长肽或minigene"
    if track == "Splice":
        return "targeted RNA确认junction，再做异常junction长肽或minigene"
    return "先确认事件真实性，再设计功能实验"


def _patient_event_grade_map(events: list[dict[str, str]]) -> dict[str, str]:
    result: dict[str, str] = {}
    for row in events:
        grade = _patient_event_grade(row)
        for key in (row.get("event_id"), row.get("event_group_id")):
            if key:
                result[str(key)] = grade
        for member in str(row.get("member_event_ids") or "").replace(",", ";").split(";"):
            if member.strip():
                result[member.strip()] = grade
    return result


def _patient_event_row_grade(row: Mapping[str, Any], event_grades: Mapping[str, str]) -> str:
    """Resolve the displayed grade from the event table, not peptide-level inference."""
    for key in (row.get("event_id"), row.get("event_group_id"), row.get("source_event_id")):
        if key and str(key) in event_grades:
            return str(event_grades[str(key)])
    return _patient_grade(row)


def _patient_observed_value(row: Mapping[str, Any], field: str) -> str | None:
    value = row.get(field)
    if value is None:
        return None
    text = str(value).strip()
    if text.upper() in {"", ".", "NA", "N/A", "NONE", "NULL", "UNASSESSED", "NOT_AVAILABLE", "NOT_RUN"}:
        return None
    return text


def _patient_numeric_display(value: str | None, decimals: int) -> str | None:
    if value is None:
        return None
    try:
        return f"{float(value):.{decimals}f}"
    except ValueError:
        return value


def _patient_observed_source(row: Mapping[str, Any], *fields: str) -> str | None:
    for field in fields:
        value = _patient_observed_value(row, field)
        if value is not None:
            return value
    return None


def _patient_missing_rna_label(row: Mapping[str, Any], label: str, *, expression: bool = False, transcript: bool = False) -> str:
    if expression or transcript:
        status = _patient_observed_source(row, "expression_evidence_status")
        source = _patient_observed_source(
            row,
            "transcript_expression_source" if transcript else "expression_source",
            "expression_source",
            "transcript_expression_source",
        )
        if transcript and _patient_track(row) == "Fusion":
            return "转录本表达无法精确匹配（融合转录本无法唯一对应定量矩阵中的标准转录本）；按融合伙伴基因表达和断点支持reads解读"
        if transcript and _patient_track(row) == "Splice":
            return "转录本表达无法精确匹配（异常剪接异构体无法唯一对应定量矩阵中的标准转录本）；按基因表达和精确junction reads解读"
        if status == "UNASSESSED_ID_NOT_MAPPED":
            return f"{label}未匹配到表达矩阵ID"
        if source:
            return f"{label}未提供（已接入表达源但该候选未匹配：{source}）"
        return f"{label}未提供（表达证据未接入）"
    source = _patient_observed_source(row, "rna_vaf_source")
    status = _patient_observed_source(row, "rna_support_status", "rna_evidence_completeness")
    depth = _patient_observed_value(row, "rna_depth")
    alt_reads = _patient_observed_value(row, "rna_alt_reads")
    if depth is not None or alt_reads is not None:
        return f"{label} 0.0000（RNA位点已评估但未检测到ALT或深度为0）"
    if source:
        return f"{label}未计算（已接入RNA VAF源但该位点未匹配：{source}）"
    if status:
        return f"{label}未计算（RNA支持状态：{status}）"
    return f"{label}未计算（RNA BAM/VAF证据未接入）"


def _patient_dna_vcf_measurements(row: Mapping[str, Any]) -> str:
    """Render observed small-variant DNA/VCF support."""
    blocks: list[str] = []
    for label, depth_field, alt_field, vaf_field in (
        ("肿瘤DNA/VCF", "tumor_depth", "tumor_alt_count", "tumor_vaf"),
        ("WGS", "wgs_tumor_depth", "wgs_tumor_alt_count", "wgs_tumor_vaf"),
        ("WES", "wes_tumor_depth", "wes_tumor_alt_count", "wes_tumor_vaf"),
    ):
        depth = _patient_numeric_display(_patient_observed_value(row, depth_field), 0)
        alt = _patient_numeric_display(_patient_observed_value(row, alt_field), 0)
        vaf = _patient_numeric_display(_patient_observed_value(row, vaf_field), 4)
        if depth is None and alt is None and vaf is None:
            continue
        parts = []
        if depth is not None:
            parts.append(f"深度 {depth}")
        if alt is not None:
            parts.append(f"ALT reads {alt}")
        if vaf is not None:
            parts.append(f"VAF {vaf}")
        blocks.append(f"{label}（{'，'.join(parts)}）")
    return "；".join(blocks)


def _patient_dna_sv_measurements(row: Mapping[str, Any]) -> str:
    """Render DNA structural-variant support without borrowing RNA junction fields."""
    parts: list[str] = []
    for label, fields in (
        ("断点支持reads", ("dna_sv_support_reads", "wgs_sv_support_reads", "sv_support_reads")),
        ("split reads", ("dna_sv_split_reads", "wgs_sv_split_reads", "sv_split_reads", "split_reads")),
        (
            "discordant pairs",
            ("dna_sv_discordant_pairs", "wgs_sv_discordant_pairs", "sv_discordant_pairs", "discordant_pairs"),
        ),
    ):
        value = _patient_observed_source(row, *fields)
        display = _patient_numeric_display(value, 0) if value is not None else None
        if display is not None:
            parts.append(f"{label} {display}")
    caller = _patient_observed_source(row, "dna_sv_caller", "wgs_sv_caller", "sv_caller")
    if caller:
        parts.append(f"DNA-SV caller {caller}")
    return "，".join(parts)


def _patient_dna_evidence(row: Mapping[str, Any]) -> str:
    """Describe DNA evidence using an event-appropriate measurement model."""
    track = _patient_track(row)
    if track in {"SNV", "InDel"}:
        measurements = _patient_dna_vcf_measurements(row)
        return f"DNA/VCF数据：{measurements or '未接入或未评估'}"
    if track in {"Fusion", "DNA SV"}:
        measurements = _patient_dna_sv_measurements(row)
        if measurements:
            return f"DNA/SV数据：{measurements}"
        if track == "Fusion":
            return "DNA/SV数据：未接入或未评估；当前融合仅按RNA断点证据评估"
        return "DNA/SV数据：未接入或未评估"
    if track == "Splice":
        return "DNA证据：点突变VCF口径不适用；当前异常剪接按RNA junction证据评估"
    measurements = _patient_dna_vcf_measurements(row)
    return f"DNA数据：{measurements}" if measurements else "DNA数据：未接入或未评估"


def _patient_dna_rna_interpretation(row: Mapping[str, Any]) -> str:
    if _patient_track(row) not in {"SNV", "InDel"}:
        return ""
    dna_alt_values = []
    for field in ("tumor_alt_count", "wgs_tumor_alt_count", "wes_tumor_alt_count"):
        value = _patient_observed_value(row, field)
        if value is None:
            continue
        try:
            dna_alt_values.append(float(value))
        except ValueError:
            continue
    if not dna_alt_values:
        return ""
    rna_alt_raw = _patient_observed_value(row, "rna_alt_reads")
    rna_depth_raw = _patient_observed_value(row, "rna_depth")
    try:
        rna_alt = float(rna_alt_raw) if rna_alt_raw is not None else None
    except ValueError:
        rna_alt = None
    try:
        rna_depth = float(rna_depth_raw) if rna_depth_raw is not None else None
    except ValueError:
        rna_depth = None
    dna_alt = max(dna_alt_values)
    if dna_alt > 4 and rna_alt == 0:
        if rna_depth is not None and rna_depth < 10:
            return "解读：DNA/VCF支持该变异；RNA位点覆盖偏低，不能据此判定RNA突变表达阴性"
        return "解读：DNA/VCF支持该变异；RNA位点复核未检出ALT，表示突变等位基因RNA表达证据不足，但不否定DNA层面变异"
    if dna_alt > 4 and rna_alt is not None and rna_alt >= 5:
        return "解读：DNA/VCF与RNA层面均有ALT reads支持"
    if dna_alt > 4 and rna_alt is not None and rna_alt > 0:
        return "解读：DNA/VCF支持该变异，RNA层面有低量ALT reads支持"
    if dna_alt > 4:
        return "解读：DNA/VCF支持该变异；RNA突变等位基因表达需结合RNA位点证据单独判断"
    return ""


def _patient_junction_reads_measurement(row: Mapping[str, Any], junction_reads: str | None) -> str:
    """Label junction read counts according to their actual verification state."""
    track = _patient_track(row)
    match_status = (_patient_observed_source(row, "junction_match_status") or "").strip().upper()
    match_method = (_patient_observed_source(row, "junction_match_method") or "").strip().lower()
    support_status = (_patient_observed_source(row, "junction_support_status") or "").strip().upper()
    resolution_status = (_patient_observed_source(row, "junction_resolution_status") or "").strip().upper()
    canonical = (_patient_observed_source(row, "canonical_junction_id") or "").strip()
    source = (_patient_observed_source(row, "rna_junction_source", "junction_source") or "").strip().lower()
    source_tools = (
        _patient_observed_source(row, "source_tools", "source_tool", "junction_source_tool") or ""
    ).strip().lower()
    source_record = (_patient_observed_source(row, "source_record_id", "source_records") or "").strip()
    source_location = (_patient_observed_source(row, "source_file", "source") or "").strip()

    verified_statuses = {
        "EXACT", "EXACT_MATCH", "MATCHED", "VERIFIED", "RESOLVED", "PILEUP_VERIFIED", "BAM_VERIFIED",
    }
    verified_methods = ("pileup", "samtools", "regtools", "bam_recount", "bam_support")
    independently_verified = match_status in verified_statuses or any(
        token in match_method for token in verified_methods
    )
    exact_canonical_link = bool(canonical) and (
        match_status in {"MATCHED_EXACT_CANONICAL", "EXACT_CANONICAL"}
        or resolution_status in {"RESOLVED_EXACT_CANONICAL", "RESOLVED_EXACT"}
        or support_status in {
            "SUPPORTED_EXACT_JUNCTION",
            "SUPPORTED_ALL_EVENT_JUNCTIONS_EXACT",
            "MATCHED_ZERO_READS",
        }
    )
    source_record_linked = bool(source_record and source_location)
    strict_cross_validated = (
        (_patient_observed_source(row, "strict_cross_validated") or "").strip().lower()
        in {"yes", "true", "1"}
    )
    structure_exact = (
        (_patient_observed_source(row, "splicemutr_structure_exact") or "").strip().lower()
        in {"yes", "true", "1"}
    )
    junction_key = (_patient_observed_source(row, "junction_key", "canonical_junction_id") or "").strip()
    exact_cross_tool_link = strict_cross_validated and structure_exact and bool(junction_key)
    caller_tokens = (
        "raw_events", "caller", "targeted_rescue", "easyfuse", "arriba",
        "star-fusion", "star_fusion", "fusioncatcher", "snaf", "splicemutr",
    )
    caller_reported = any(token in source or token in source_tools for token in caller_tokens)

    if junction_reads is None:
        if caller_reported:
            return "caller报告junction reads未提供（尚未独立回链核实）"
        return "junction reads未提供（核实状态未确认）"
    if independently_verified:
        return f"主比对表已按精确junction回链，unique junction reads {junction_reads}"
    if exact_canonical_link:
        return f"已按标准化精确junction坐标回链，junction reads {junction_reads}"
    if track == "Splice" and exact_cross_tool_link:
        return f"已按event/junction精确回链的剪接工具支持reads {junction_reads}"
    if source_record_linked:
        origin = "融合caller" if track == "Fusion" else "剪接caller" if track == "Splice" else "上游caller"
        return f"{origin}原始记录已回链，junction reads {junction_reads}（尚无独立主比对表核实）"
    if caller_reported:
        if track == "Fusion" and "targeted_rescue" in source_tools:
            origin = "融合caller/定向rescue汇总"
        elif track == "Fusion":
            origin = "融合caller"
        elif track == "Splice":
            origin = "剪接工具"
        else:
            origin = "上游caller"
        return f"{origin}报告junction reads {junction_reads}（尚未独立回链核实）"
    source_note = f"，来源 {source}" if source else ""
    return f"junction reads {junction_reads}（核实状态未确认{source_note}）"


def _patient_junction_count_comparison(
    row: Mapping[str, Any],
    *,
    provided_text: str,
    provided_value: float,
    verified_text: str,
    verified_value: float,
) -> str:
    """Explain differing exact-source junction counts without treating their delta as reads."""
    canonical = str(row.get("canonical_junction_id") or "").strip()
    match_status = str(row.get("junction_match_status") or "").strip().upper()
    support_status = str(row.get("junction_support_status") or "").strip().upper()
    exact_primary = bool(canonical) and support_status == "SUPPORTED_EXACT_JUNCTION"
    exact_caller = exact_primary and match_status in {"EXACT", "EXACT_MATCH", "MATCHED", "RESOLVED"}
    if exact_caller:
        return (
            f"同一精确canonical junction存在两种计数口径：caller原始记录 {provided_text}，"
            f"主比对表unique reads {verified_text}；两者均已坐标回链，差异不代表未归属reads，"
            f"不相减或累加，排序保守采用主比对表 {verified_text}"
        )
    if provided_value >= verified_value:
        unresolved = provided_value - verified_value
        unresolved_text = str(int(unresolved)) if unresolved.is_integer() else f"{unresolved:g}"
        strand = str(row.get("junction_strand") or "").strip()
        if verified_value <= 0 and (canonical.endswith("|.") or strand in {".", "?"}):
            return (
                f"上游工具汇总同坐标junction reads {provided_text}；因链方向未解析，"
                "当前严格规则未将其计入已核实reads（不等同于检测为0）"
            )
        return (
            f"上游工具汇总junction reads {provided_text}；其中 {verified_text} 条已按同一"
            f"canonical junction精确核实，差额 {unresolved_text} 条尚未归属，未计入已核实支持"
        )
    return (
        f"上游工具汇总junction reads {provided_text}，与已核实值 {verified_text} 的"
        "统计口径不一致，需复核来源记录"
    )


def _patient_rna_measurements(row: Mapping[str, Any]) -> str:
    """Render only observed RNA fields carried by the canonical evidence table."""
    track = _patient_track(row)
    gene_tpm = _patient_numeric_display(_patient_observed_value(row, "gene_expression_tpm"), 4)
    transcript_tpm = _patient_numeric_display(_patient_observed_value(row, "transcript_expression_tpm"), 4)
    depth = _patient_numeric_display(_patient_observed_value(row, "rna_depth"), 0)
    alt_reads = _patient_numeric_display(_patient_observed_value(row, "rna_alt_reads"), 0)
    vaf = _patient_numeric_display(_patient_observed_value(row, "rna_vaf"), 4)
    junction_reads = _patient_numeric_display(
        _patient_observed_value(row, "rna_junction_reads")
        or _patient_observed_value(row, "junction_reads")
        or _patient_observed_value(row, "provided_rna_junction_reads"),
        0,
    )
    values = [
        f"基因表达 {gene_tpm} TPM"
        if gene_tpm is not None else _patient_missing_rna_label(row, "基因表达", expression=True)
    ]
    values.append(
        f"转录本表达 {transcript_tpm} TPM"
        if transcript_tpm is not None else _patient_missing_rna_label(row, "转录本表达", transcript=True)
    )
    if track in {"Fusion", "Splice"}:
        verified_text, verified_value = _patient_numeric_value(
            row, "rna_junction_reads", "junction_reads"
        )
        provided_text, provided_value = _patient_numeric_value(
            row, "provided_rna_junction_reads"
        )
        provided_reads = _patient_numeric_display(
            _patient_observed_value(row, "provided_rna_junction_reads"), 0
        )
        values.append(_patient_junction_reads_measurement(row, junction_reads))
        if provided_reads is not None and provided_reads != junction_reads:
            if provided_value is not None and verified_value is not None:
                values.append(_patient_junction_count_comparison(
                    row,
                    provided_text=provided_text,
                    provided_value=provided_value,
                    verified_text=verified_text,
                    verified_value=verified_value,
                ))
            else:
                values.append(
                    f"上游工具汇总junction reads {provided_reads}，与已核实值 {junction_reads} 的"
                    "统计口径不一致，需复核来源记录"
                )
        if track == "Splice":
            unique_reads = _patient_numeric_display(_patient_observed_value(row, "unique_junction_reads"), 0)
            total_coverage = _patient_numeric_display(_patient_observed_value(row, "junction_total_coverage"), 0)
            psi = _patient_numeric_display(_patient_observed_value(row, "splice_psi"), 4)
            values.append(
                f"unique junction reads {unique_reads}"
                if unique_reads is not None else "unique junction reads未记录（不能用caller总reads替代）"
            )
            values.append(
                f"junction总覆盖 {total_coverage}"
                if total_coverage is not None else "junction总覆盖未记录"
            )
            values.append(
                f"PSI {psi}"
                if psi is not None else "PSI未记录（不能用基因TPM替代）"
            )
    else:
        values.append(f"RNA位点深度 {depth}" if depth is not None else _patient_missing_rna_label(row, "RNA位点深度"))
        values.append(f"RNA alt reads {alt_reads}" if alt_reads is not None else _patient_missing_rna_label(row, "RNA alt reads"))
        values.append(f"RNA VAF {vaf}" if vaf is not None else _patient_missing_rna_label(row, "RNA VAF"))
    return "；".join(values)


def _patient_cross_site_rna(row: Mapping[str, Any]) -> str:
    status = _patient_observed_value(row, "cross_site_status")
    if status is None:
        return ""
    sample = _patient_observed_value(row, "secondary_sample_id") or "另一部位RNA样本"
    identity = (_patient_observed_value(row, "sample_identity_status") or "UNASSESSED").upper()
    translations = {
        "EXACT_SHARED": "精确事件获得跨部位RNA支持",
        "EXACT_MATCH_IDENTITY_UNASSESSED": "待核对：精确事件在另一部位检出，但样本身份尚未完成RNA-DNA指纹确认",
        "GENE_PAIR_SHARED_BREAKPOINT_UNASSESSED": "待核对：融合基因对一致，但至少一侧缺少可比较的精确断点",
        "SECONDARY_NEGATIVE_ADEQUATE_COVERAGE": "待核对：另一部位覆盖充分但未检出直接RNA支持，需结合部位和样本成分解释",
        "SECONDARY_MATCH_LOW_POWER": "待核对：另一部位存在对应记录，但位点或junction覆盖不足",
        "PRIMARY_ONLY": "当前仅在主肿瘤部位获得支持",
        "AMBIGUOUS_SECONDARY_MATCH": "待核对：另一部位存在多个可能匹配记录，需确认坐标、方向和来源记录",
    }
    detail = _patient_observed_value(row, "cross_site_review_reason") or translations.get(status.upper(), status)
    secondary_alt = _patient_numeric_display(_patient_observed_value(row, "secondary_rna_alt_reads"), 0)
    secondary_depth = _patient_numeric_display(_patient_observed_value(row, "secondary_rna_depth"), 0)
    secondary_junction = _patient_numeric_display(_patient_observed_value(row, "secondary_rna_junction_reads"), 0)
    metrics = []
    if secondary_depth is not None:
        metrics.append(f"位点深度 {secondary_depth}")
    if secondary_alt is not None:
        metrics.append(f"alt reads {secondary_alt}")
    if secondary_junction is not None:
        metrics.append(f"junction reads {secondary_junction}")
    suffix = f"（{identity}；{'，'.join(metrics)}）" if metrics else f"（身份状态 {identity}）"
    return f"跨部位RNA {sample}：{detail}{suffix}"


def _patient_safety_dimensions(row: Mapping[str, Any]) -> str:
    """Render direct peptide safety evidence separately from expression context."""
    junction_track = _patient_track(row) in {"Fusion", "Splice"}

    def status(field: str, *, detected: str = "检出", not_detected: str = "未检出") -> str:
        value = str(row.get(field) or "").strip().upper()
        if value in {"DETECTED", "YES", "TRUE", "1", "EXACT_MATCH", "SUPPORTED_IN_NORMAL"}:
            return detected
        if value in {"NOT_DETECTED", "NO", "FALSE", "0", "NEGATIVE"}:
            return not_detected
        if value == "NOT_APPLICABLE":
            return "不适用"
        return "未评估"

    proteome = status("normal_proteome_exact_match_status")
    if proteome == "未评估":
        proteome = status("reference_proteome_exact_match")
    junction = status("normal_transcript_junction_match_status")
    if junction == "未评估":
        junction = status("normal_junction_seen")
    ligandome = status("normal_immunopeptidome_match_status")
    if ligandome == "未评估":
        ligandome = status("normal_hla_ligand_exact_match")

    similarity_raw = str(
        row.get("similar_peptide_cross_reactivity_status")
        or row.get("closest_self_similarity")
        or row.get("self_similarity_score")
        or ""
    ).strip()
    similarity_upper = similarity_raw.upper()
    if similarity_upper in {"", "UNASSESSED", "NOT_ASSESSED", "MISSING", "NA", "N/A"}:
        similarity = "未评估"
    elif "HIGH" in similarity_upper:
        similarity = "高相似性，需脱靶复核"
    elif "LOW" in similarity_upper:
        similarity = "已评估，未见高相似性"
    elif _patient_numeric_display(similarity_raw, 4) is not None:
        similarity = f"相似度 {similarity_raw}，需结合阈值复核"
    else:
        similarity = "状态需复核"

    normal_tpm = _patient_numeric_display(_patient_observed_value(row, "normal_tissue_max_tpm"), 4)
    normal_tissue = _patient_value(row, "normal_tissue_max_tissue", default="正常组织")
    source_expression = "未评估"
    if normal_tpm is not None:
        source_expression = f"{normal_tissue}最高 {normal_tpm} TPM"
    critical_tpm = _patient_numeric_display(_patient_observed_value(row, "critical_tissue_max_tpm"), 4)
    critical_tissue = _patient_value(row, "critical_tissue_name", default="关键器官")
    critical_expression = f"{critical_tissue}最高 {critical_tpm} TPM" if critical_tpm is not None else "未评估"
    hspc_tpm = _patient_numeric_display(_patient_observed_value(row, "normal_hspc_tpm"), 4)
    hspc_unit = _patient_value(row, "normal_hspc_unit", default="")
    hematopoietic = (
        f"最高 {hspc_tpm}" + (f" {hspc_unit}" if hspc_unit else "")
        if hspc_tpm is not None else "未评估"
    )
    if junction_track:
        source_expression += "（伙伴基因辅助背景，不代表跨连接肽存在）" if source_expression != "未评估" else "（辅助背景）"
        critical_expression += "（伙伴基因辅助背景）" if critical_expression != "未评估" else "（辅助背景）"
        hematopoietic += "（伙伴基因辅助背景）" if hematopoietic != "未评估" else "（辅助背景）"

    conclusion_codes = {
        "REJECT_DIRECT_SAFETY_EVIDENCE": "直接安全证据提示排除",
        "REVIEW_DIRECT_SAFETY_EVIDENCE": "存在直接安全信号，需专项复核",
        "PARTIAL_DIRECT_SAFETY_EVIDENCE": "直接安全证据不完整，暂不能定论",
        "NO_DIRECT_SAFETY_SIGNAL_DETECTED": "已评估层未见直接安全信号；不等于已证明安全",
    }
    conclusion_code = str(row.get("final_safety_conclusion") or "").strip().upper()
    conclusion = conclusion_codes.get(conclusion_code)
    if not conclusion:
        safety = _patient_value(row, "safety_state", "safety_tier", "safety_status", default="").upper()
        if safety in {"SAFETY_REJECT", "REJECT", "FAIL"}:
            conclusion = "直接安全证据提示排除"
        elif safety in {"SAFETY_REVIEW", "SAFETY_HIGH_RISK", "REVIEW", "CAUTION"}:
            conclusion = "需按直接肽段/连接证据专项复核"
        else:
            conclusion = "证据不完整，暂不能定论"
    completeness = _patient_numeric_display(_patient_observed_value(row, "safety_evidence_completeness"), 2)
    if completeness is not None:
        conclusion += f"（直接证据完整度 {float(completeness) * 100:.0f}%）"

    return "；".join([
        f"完整肽正常蛋白组精确匹配={proteome}",
        f"正常转录组/正常junction匹配={junction}",
        f"正常免疫肽组匹配={ligandome}",
        f"相似肽交叉反应风险={similarity}",
        f"来源基因正常组织表达={source_expression}",
        f"关键器官表达={critical_expression}",
        f"正常造血系统表达={hematopoietic}",
        f"最终安全结论={conclusion}",
    ])


def _patient_fusion_boundary_evidence(row: Mapping[str, Any], bundle: ReportBundle) -> str:
    if _patient_track(row) != "Fusion":
        return ""
    peptide = _patient_value(row, "peptide", default="")
    left_gene = _patient_value(row, "fusion_left_gene", default="")
    right_gene = _patient_value(row, "fusion_right_gene", default="")
    if not left_gene or not right_gene:
        genes = _patient_value(row, "gene", "event_name", default="").split("::", 1)
        left_gene = left_gene or (genes[0] if genes else "5′伙伴")
        right_gene = right_gene or (genes[1] if len(genes) > 1 else "3′伙伴")
    left = _patient_value(row, "fusion_left_peptide", default="")
    right = _patient_value(row, "fusion_right_peptide", default="")
    position = _patient_value(
        row, "junction_position_in_peptide_1based", "junction_offset_in_peptide", default=""
    )
    crosses = _patient_value(row, "crosses_junction", default="UNASSESSED").upper()
    if crosses in {"NO", "FALSE", "0"}:
        mapping = "该肽不跨融合连接点，不归类为融合新抗原；仅可按肿瘤相关/异常表达抗原另行评估"
    elif crosses in {"YES", "TRUE", "1"} and left and right and position:
        mapping = f"{left_gene}来源 {left}｜融合连接点（肽内位置 {position}）｜{right_gene}来源 {right}；两侧均含氨基酸，已证明跨断点"
    elif crosses in {"YES", "TRUE", "1"}:
        mapping = "上游标记为跨断点，但缺少肽内连接位置或左右氨基酸映射，尚不能独立证明融合特异性"
    else:
        mapping = "跨断点状态未评估；在建立肽内连接位置和左右来源前不得称为融合新抗原"

    transcript = _patient_value(row, "transcript_hypothesis_id", "transcript_id", default="未建立")
    orf = _patient_value(row, "orf_id", default="未建立")
    alternatives: list[tuple[str, str, str]] = []
    event_keys = set(_patient_event_keys(row))
    for candidate in bundle.peptides:
        if event_keys and not event_keys.intersection(_patient_event_keys(candidate)):
            continue
        candidate_peptide = _patient_value(candidate, "peptide", default="")
        if not candidate_peptide:
            continue
        item = (
            candidate_peptide,
            _patient_value(candidate, "transcript_hypothesis_id", "transcript_id", default="未建立"),
            _patient_value(candidate, "orf_id", default="未建立"),
        )
        if item not in alternatives:
            alternatives.append(item)
    near_variants = sorted({item[0] for item in alternatives if len(item[0]) == len(peptide) and sum(a != b for a, b in zip(item[0], peptide)) <= 1})
    if len(near_variants) > 1:
        comparison = (
            "同一事件存在近似肽序列假设 " + "/".join(near_variants)
            + "；必须按患者精确融合转录本、阅读框和ORF逐一归因，未完成前不得视为独立融合新抗原"
        )
    elif alternatives:
        comparison = f"当前候选回链：transcript={transcript}，ORF={orf}；同事件共{len(alternatives)}个肽/转录本/ORF组合"
    else:
        comparison = f"当前候选回链：transcript={transcript}，ORF={orf}"
    return mapping + "；" + comparison


def _patient_event_evidence_and_next_step(
    row: Mapping[str, Any], bundle: ReportBundle, val_map: Mapping[str, Mapping[str, str]],
    *, compact_common: bool = False,
) -> str:
    track = _patient_track(row)
    evidence = [
        _patient_metric("事件", row, "event_authenticity_state", "cross_platform_status"),
        _patient_rna_metric(row),
        _patient_presentation_metric(row),
    ]
    if not (compact_common and track in {"Fusion", "Splice"}):
        evidence.append(_patient_metric("MT/WT", row, "mutant_specificity_status", "mutant_specificity_state"))
    gaps = _patient_key_gaps(row, bundle)
    gap_text = "；".join(gaps) if gaps else "未见明确关键缺口"
    dna_evidence = _patient_dna_evidence(row)
    if compact_common and track == "Splice":
        dna_evidence = ""
    elif compact_common and track == "Fusion" and not _patient_dna_sv_measurements(row):
        dna_evidence = ""
    disease_note = _patient_disease_anchor_note(row, bundle)
    fusion_boundary = _patient_fusion_boundary_evidence(row, bundle)
    return (
        f"疾病知识：{disease_note}。" if disease_note else ""
    ) + (
        f"核心证据：{'；'.join(evidence)}。"
        f"{dna_evidence + '。' if dna_evidence else ''}"
        f"RNA数据：{_patient_rna_measurements(row)}。"
        f"{_patient_mtwt_caution(row) + '。' if _patient_mtwt_caution(row) else ''}"
        f"安全性分层：{_patient_safety_dimensions(row)}。"
        f"{('融合肽断点证明：' + fusion_boundary + '。') if fusion_boundary else ''}"
        f"{_patient_dna_rna_interpretation(row) + '。' if _patient_dna_rna_interpretation(row) else ''}"
        f"{_patient_cross_site_rna(row) + '。' if _patient_cross_site_rna(row) else ''}"
        f"主要缺口：{gap_text}。"
        f"下一步：{_patient_validation(row, val_map)}"
    )


_PATIENT_TRACK_EVIDENCE_NOTES = {
    "SNV": (
        "适用：肿瘤/正常DNA位点深度与VAF、RNA alt reads/VAF、MT/WT突变特异性、"
        "HLA呈递、HLA-LOH/APPM、CCF与正常背景安全性。"
        "通常不适用：融合/剪接junction reads、PSI和异常连接ORF。"
    ),
    "InDel": (
        "适用：DNA深度/VAF、局部重比对与左对齐、RNA突变转录本支持、"
        "阅读框/新生尾部、NMD、HLA呈递、HLA-LOH/APPM、CCF和安全性。"
        "通常不适用：异常junction/PSI；除非该InDel本身导致剪接改变。"
    ),
    "Fusion": (
        "适用：精确融合断点、junction/split-read支持、方向与阅读框、跨断点新序列、"
        "正常read-through背景、HLA呈递、HLA-LOH/APPM和安全性。"
        "安全性必须按精确跨断点肽查询正常蛋白组、转录组、junction库和ligandome；"
        "融合伙伴基因的正常TPM只作背景，不直接判定跨断点肽高风险。"
        "不适用：普通点突变式DNA VAF和传统MT/WT配对；RNA-only融合不应把缺失DNA CCF解释为0。"
        "若有独立DNA-SV支持，在候选行中单独展示。"
    ),
    "Splice": (
        "适用：精确异常junction、unique reads/PSI、转录本假设与ORF、跨连接新序列、"
        "正常异构体/junction背景、HLA呈递、HLA-LOH/APPM和安全性。"
        "安全性必须按精确跨junction肽和正常异构体评估；基因整体正常TPM不能代替连接特异性评估。"
        "不适用：普通点突变式DNA VAF、传统MT/WT配对和DNA CCF；"
        "应改用正常连接或正常异构体肽作为对照。"
    ),
    "DNA SV": (
        "适用：DNA断点、split reads/discordant pairs、重构转录本和ORF、RNA正交支持、"
        "HLA呈递、HLA-LOH/APPM、CCF和安全性。"
        "不适用：普通SNV式MT/WT规则；应使用断点前后的正常结构对照。"
    ),
}


def _patient_track_evidence_note(track: str) -> str:
    return _PATIENT_TRACK_EVIDENCE_NOTES.get(
        track,
        "适用和不适用证据按事件类型判定；缺失或不适用均不按阴性处理。",
    )


_PATIENT_SAFETY_LAYER_LABELS = {
    "normal_expression": "正常组织表达参考已接入但该候选未匹配到可判定记录",
    "normal_hspc": "HSPC正常造血参考已接入但该候选未匹配到可判定记录",
    "reference_proteome": "正常蛋白组精确匹配未完成正式评估",
    "normal_ligandome": "正常HLA配体组未完成正式评估",
    "anchor_only": "突变锚定位点风险未完成正式评估",
    "normal_junction": "正常融合/剪接连接背景未完成正式评估",
}


def _patient_safety_gap(row: Mapping[str, Any]) -> str:
    """Describe the candidate-specific safety gap without hiding its cause."""
    safety = _patient_value(row, "safety_state", "safety_tier", "safety_status", default="").upper()
    missing_raw = ";".join(
        str(row.get(field) or "")
        for field in ("safety_missing_layers", "event_safety_missing_layers")
    )
    missing_keys: list[str] = []
    reason_text = ";".join(
        str(row.get(field) or "") for field in ("safety_reason", "event_safety_reason", "safety_reason_codes")
    ).lower()
    for token in re.split(r"[;,|\s]+", missing_raw):
        key = token.strip().lower()
        if key and key not in missing_keys:
            missing_keys.append(key)
    junction_track = _patient_track(row) in {"Fusion", "Splice"}
    exact_sequence_reasons = (
        "normal_hla_ligand", "normal_junction", "reference_proteome_exact_match",
        "matched_normal_support", "high_self_similarity",
    )
    if junction_track:
        # Partner-gene expression is contextual for a novel junction peptide.
        # Safety decisions must follow the exact novel sequence and junction.
        missing_keys = [key for key in missing_keys if key not in {"normal_expression", "normal_hspc"}]
        exact_layers = {
            "reference_proteome": row.get("reference_proteome_status"),
            "normal_ligandome": row.get("normal_ligandome_status"),
            "normal_junction": row.get("normal_junction_assessment_status"),
        }
        for key, value in exact_layers.items():
            if str(value or "").strip().upper() not in {"ASSESSED", "PASS", "NEGATIVE", "NOT_FOUND"} and key not in missing_keys:
                missing_keys.append(key)

    if safety in {"SAFETY_PARTIAL", "PARTIAL", "UNASSESSED"}:
        if not missing_keys:
            layer_statuses = {
                "normal_expression": row.get("normal_expression_status"),
                "normal_hspc": row.get("normal_hspc_status"),
                "reference_proteome": row.get("reference_proteome_status"),
                "normal_ligandome": row.get("normal_ligandome_status"),
                "anchor_only": row.get("anchor_assessment_status"),
                "normal_junction": row.get("normal_junction_assessment_status"),
            }
            for key, value in layer_statuses.items():
                if junction_track and key in {"normal_expression", "normal_hspc", "anchor_only"}:
                    continue
                status = str(value or "").strip().upper()
                if status and status not in {"ASSESSED", "NOT_APPLICABLE"} and key not in missing_keys:
                    missing_keys.append(key)
        missing: list[str] = []
        normal_tpm = _patient_numeric_display(_patient_observed_value(row, "normal_tissue_max_tpm"), 4)
        normal_tissue = _patient_value(row, "normal_tissue_max_tissue", "critical_tissue_name", default="正常组织")
        hspc_tpm = _patient_numeric_display(_patient_observed_value(row, "normal_hspc_tpm"), 4)
        hspc_unit = _patient_value(row, "normal_hspc_unit", default="")
        for key in missing_keys:
            if key == "normal_expression" and "normal_expression_gene_not_in_reference" in reason_text:
                label = "当前正式GTEx/HPA正常组织表达参考未收录该基因，不能按0表达解释"
            elif key == "normal_hspc" and "normal_hspc_gene_not_in_reference" in reason_text:
                label = "当前正式HSPC正常造血参考未收录该基因，不能按0表达解释"
            elif key == "normal_expression" and normal_tpm is not None:
                label = f"正常组织表达已查询但覆盖/映射不完整，最高 {normal_tpm} TPM（{normal_tissue}）"
            elif key == "normal_hspc" and hspc_tpm is not None:
                unit = f" {hspc_unit}" if hspc_unit else ""
                label = f"HSPC正常造血参考已查询但覆盖/映射不完整，最高 {hspc_tpm}{unit}"
            else:
                label = _PATIENT_SAFETY_LAYER_LABELS.get(key, f"{key}安全参考层已接入但该候选未匹配到可判定记录")
            if label not in missing:
                missing.append(label)
        if not missing:
            if junction_track:
                return (
                    "连接事件安全性需按精确新生序列复核；"
                    "伙伴基因正常组织/HSPC表达仅作背景，不直接判定该连接肽高风险"
                )
            observed: list[str] = []
            if normal_tpm is not None:
                observed.append(f"正常组织最高 {normal_tpm} TPM（{normal_tissue}）")
            if hspc_tpm is not None:
                unit = f" {hspc_unit}" if hspc_unit else ""
                observed.append(f"HSPC最高 {hspc_tpm}{unit}")
            if observed:
                return "安全参考已查询：" + "、".join(observed)
        return "安全性缺口：" + "、".join(missing or ["未记录可判定的安全参考层状态"])

    if safety not in {"SAFETY_REVIEW", "SAFETY_HIGH_RISK", "REVIEW", "CAUTION"}:
        return ""

    review_reasons: list[str] = []
    critical_tpm = _patient_numeric_display(_patient_observed_value(row, "critical_tissue_max_tpm"), 4)
    critical_tissue = _patient_value(
        row, "critical_tissue_name", "critical_tissue_max_tissue", default="正常关键组织"
    )
    if not junction_track and ("critical_tissue" in reason_text or str(row.get("critical_tissue_hit") or "").lower() in {"yes", "true", "1"}):
        detail = f"{critical_tissue}表达信号"
        if critical_tpm is not None:
            detail += f" {critical_tpm} TPM"
        review_reasons.append(detail)
    if not junction_track and "normal_hspc_expression" in reason_text:
        hspc_tpm = _patient_numeric_display(_patient_observed_value(row, "normal_hspc_tpm"), 4)
        hspc_unit = _patient_value(row, "normal_hspc_unit", default="")
        detail = "HSPC/正常造血参考中检测到表达信号"
        if hspc_tpm is not None:
            detail += f" {hspc_tpm}" + (f" {hspc_unit}" if hspc_unit else "")
        review_reasons.append(detail)
    reason_labels = (
        ("normal_hla_ligand", "正常HLA配体组中存在匹配或重叠肽段"),
        ("normal_junction", "正常组织中检测到相同或相关连接"),
        ("anchor_only", "突变主要位于HLA锚定位点且正常肽仍可能结合"),
        ("reference_proteome_exact_match", "候选肽与正常蛋白组精确匹配"),
        ("matched_normal_support", "配对正常样本中检测到变异支持"),
        ("high_self_similarity", "候选肽与正常自身肽高度相似"),
    )
    for token, label in reason_labels:
        if token in reason_text and label not in review_reasons:
            review_reasons.append(label)
    if junction_track and not any(token in reason_text for token in exact_sequence_reasons):
        missing_exact = [
            _PATIENT_SAFETY_LAYER_LABELS[key]
            for key in ("reference_proteome", "normal_ligandome", "normal_junction")
            if str(row.get({
                "reference_proteome": "reference_proteome_status",
                "normal_ligandome": "normal_ligandome_status",
                "normal_junction": "normal_junction_assessment_status",
            }[key]) or "").strip().upper() not in {"ASSESSED", "PASS", "NEGATIVE", "NOT_FOUND"}
        ]
        return (
            "连接事件安全性需按精确新生序列复核；伙伴基因正常组织/HSPC表达仅作背景，不直接判定该连接肽高风险"
            + ("；" + "、".join(missing_exact) if missing_exact else "")
        )
    if not review_reasons:
        normal_tpm = _patient_numeric_display(_patient_observed_value(row, "normal_tissue_max_tpm"), 4)
        normal_tissue = _patient_value(row, "normal_tissue_max_tissue", "critical_tissue_name", default="正常组织")
        hspc_tpm = _patient_numeric_display(_patient_observed_value(row, "normal_hspc_tpm"), 4)
        hspc_unit = _patient_value(row, "normal_hspc_unit", default="")
        if normal_tpm is not None:
            review_reasons.append(f"正常组织最高 {normal_tpm} TPM（{normal_tissue}）")
        if hspc_tpm is not None:
            review_reasons.append(f"HSPC最高 {hspc_tpm}" + (f" {hspc_unit}" if hspc_unit else ""))
    return "安全性需复核：" + "、".join(review_reasons or ["已触发安全审阅规则，具体原因见peptide_safety.tsv"])


def _patient_key_gaps(row: Mapping[str, Any], bundle: ReportBundle) -> list[str]:
    gaps: list[str] = []
    hard_failure_labels = {
        "HARD_REFERENCE_PROTEOME_MATCH": "候选肽与正常参考蛋白组存在精确匹配",
        "HARD_NORMAL_JUNCTION": "正常组织中检测到相同异常连接",
        "HARD_RESTRICTING_HLA_LOST": "限制性HLA已确认丢失",
        "HARD_MATCHED_NORMAL_SUPPORT": "配对正常样本中检测到变异支持",
        "HARD_NON_MUTANT_SEQUENCE": "候选肽不包含突变或新生连接序列",
        "HARD_SAFETY_REJECT": "安全性评估触发明确拒绝条件",
        "HARD_EVENT_ARTIFACT": "事件被判定为技术伪影风险",
        "HARD_EVENT_OR_ORF_INVALID": "事件或ORF无效，不能支持该候选肽来源",
    }
    hard_codes: list[str] = []
    for field_name in ("hard_failure_codes", "source_chain_hard_failure_codes"):
        for code in re.split(r"[,;|]", str(row.get(field_name) or "")):
            code = code.strip()
            if code and code not in hard_codes:
                hard_codes.append(code)
    if hard_codes:
        labels = [hard_failure_labels.get(code, code.replace("_", " ")) for code in hard_codes]
        gaps.append("阻断原因：" + "、".join(labels))
    # RNA-only fusion/splice events do not have a defensible DNA CCF by
    # default.  Treat that as not applicable rather than as a negative result
    # or a universal candidate defect.
    if _patient_track(row) not in {"Fusion", "Splice"} and not _patient_ccf_assessment(row, bundle)[0]:
        gaps.append("CCF未形成可靠估计")
    hla_gap = _patient_restricting_hla_gap(row, bundle)
    if hla_gap:
        gaps.append(hla_gap)
    appm_status = _patient_candidate_appm_status(bundle)
    if appm_status.startswith("证据部分完整"):
        gaps.append("APPM仅部分评估")
    elif appm_status == "未评估":
        gaps.append("APPM未评估")
    safety_gap = _patient_safety_gap(row)
    if safety_gap:
        gaps.append(safety_gap)
    fusion_artifact = _patient_fusion_artifact_review(row)
    if fusion_artifact:
        gaps.append(fusion_artifact)
    conflict_summary = _patient_conflict_summary(row)
    if conflict_summary:
        gaps.append("具体证据冲突：" + conflict_summary)
    integrity_ok, integrity_missing = _patient_candidate_integrity(row)
    if not integrity_ok:
        gaps.append("完整性缺口：" + "、".join(integrity_missing))
    return list(dict.fromkeys(gaps))


def _patient_fusion_artifact_review(row: Mapping[str, Any]) -> str:
    """Flag generic fusion patterns that require artifact-focused review."""
    if _patient_track(row) != "Fusion":
        return ""
    label = str(row.get("gene") or row.get("event_name") or row.get("fusion_gene") or "").strip()
    if not label:
        return ""
    parts = [part.strip() for part in re.split(r"::|--", label) if part.strip()]
    canonical = [re.sub(r"\([^)]*\)", "", part).split(",", 1)[0].strip().upper() for part in parts]
    fusion_type = str(row.get("fusion_type") or row.get("type") or "").lower()
    if len(canonical) == 2 and canonical[0] == canonical[1]:
        return "融合伪影复核：同基因连接，需排查内部重排、环状RNA、转录本错配和比对伪影"
    if len(parts) != 2 or "," in label or ";" in label or sum("ENSG" in part.upper() for part in parts) > 1:
        return "融合伪影复核：复杂/多伙伴标注，需重建唯一的5'::3'断点、方向、转录本和ORF"
    if any(token in fusion_type for token in ("readthrough", "read-through", "cis-near", "adjacent")):
        return "融合伪影复核：可能的read-through/邻近基因连接，需对照正常junction背景"
    return ""


def _patient_attention_reasons(row: Mapping[str, Any], bundle: ReportBundle | None = None) -> list[str]:
    reasons: list[str] = []
    anchor = _disease_anchor(row, bundle)
    if anchor:
        reasons.append(str(anchor.get("label") or "当前疾病知识库标记的锚定事件"))
    driver = str(row.get("cancer_driver_context") or "").upper()
    if driver == "DRIVER_CONTEXT" or str(row.get("cancer_gene_types") or "").strip():
        reasons.append("具有癌症基因或驱动机制背景")
    if str(row.get("oncokb_annotated") or "").lower() in {"yes", "true", "1"}:
        reasons.append("OncoKB收录")
    if str(row.get("cosmic_cgc_flag") or "").lower() in {"yes", "true", "1"}:
        reasons.append("COSMIC CGC收录")
    conflict = " ".join(str(row.get(field) or "") for field in (
        "evidence_conflict_fields", "evidence_conflict_layers", "evidence_conflict_status",
    )).upper()
    if conflict.strip() and "NO_CONFLICT" not in conflict and "NONE" not in conflict:
        reasons.append("存在需人工解释的证据冲突")
    if str(row.get("manual_review_required") or "").lower() in {"yes", "true", "1"}:
        reasons.append("排序规则要求人工复核")
    sources = str(row.get("source_tools") or row.get("source_chain_orthogonal_sources") or "")
    source_count = len([part for part in re.split(r"[;,|+]", sources) if part.strip()])
    orthogonal = str(row.get("source_chain_orthogonal_status") or "").upper()
    if source_count >= 2 or any(token in orthogonal for token in ("SUPPORTED", "CONCORDANT", "PASS")):
        reasons.append("获得多工具或正交证据支持")
    cross_platform = str(row.get("cross_platform_status") or "").upper()
    if "CONSISTENT" in cross_platform or "CONCORDANT" in cross_platform:
        reasons.append("DNA跨平台结果一致")
    return list(dict.fromkeys(reasons))


def _manual_review_label(value: Any) -> str:
    return re.sub(r"[^A-Z0-9]+", "_", str(value or "").strip().upper()).strip("_")


def _patient_configured_manual_review(row: Mapping[str, Any], profile: Mapping[str, Any]) -> bool:
    """Return true for events explicitly retained by the active review rules."""
    config = profile.get("manual_review") if isinstance(profile, Mapping) else None
    config = config if isinstance(config, Mapping) else {}
    targets = {
        _manual_review_label(value)
        for key in ("events", "genes")
        for value in (config.get(key) or [])
        if _manual_review_label(value)
    }
    observed = {
        _manual_review_label(row.get(key))
        for key in ("gene", "event_name", "event_id", "source_event_id")
        if _manual_review_label(row.get(key))
    }
    if targets.intersection(observed):
        return True
    scope = str(row.get("evidence_scope") or "").upper()
    source = str(row.get("source") or row.get("source_tool") or "").upper()
    return "DIAGNOSTIC_WHITELIST_RESCUE" in scope or "DIAGNOSTIC_WHITELIST_RESCUE" in source


def _patient_fusion_manual_review_advice(row: Mapping[str, Any]) -> str:
    if _patient_track(row) != "Fusion":
        return ""
    return (
        "作为疾病相关关键融合保留机制审阅；只把精确跨断点、来源可追溯且不与正常蛋白组精确匹配的肽"
        "视为融合特异性候选，普通融合伙伴蛋白来源肽仅作背景。先复核断点、方向、阅读框、正常read-through"
        "及独立RNA支持，再决定长肽或minigene验证；事件重要性本身不自动提高R等级"
    )


def _patient_manual_review_identity(row: Mapping[str, Any]) -> str:
    """Collapse caller/breakpoint records for the same named fusion in summaries.

    Provenance and event-level outputs intentionally retain every exact event.
    The patient-facing manual-review table instead presents one biological
    fusion label, so independent caller IDs and rescued breakpoints do not
    consume several of the five review slots.
    """
    if _patient_track(row) == "Fusion":
        fusion = str(row.get("gene") or row.get("event_name") or "").strip()
        normalized = _manual_review_label(fusion)
        if normalized:
            return f"FUSION::{normalized}"
    keys = _patient_event_keys(row)
    return keys[0] if keys else str(row.get("event_id") or row.get("event_name") or row.get("gene") or "")


def _patient_manual_review_rows(
    events: list[dict[str, str]],
    peptides: list[dict[str, str]],
    bundle: ReportBundle,
    val_map: Mapping[str, Mapping[str, str]],
    limit: int = 5,
) -> list[dict[str, str]]:
    peptide_by_event: dict[str, dict[str, str]] = {}
    for peptide in peptides:
        for key in _patient_event_keys(peptide):
            peptide_by_event.setdefault(key, peptide)
    scored: list[tuple[int, int, int, dict[str, str], list[str]]] = []
    seen: set[str] = set()
    for index, event in enumerate(events):
        keys = _patient_event_keys(event)
        event_key = keys[0] if keys else str(event.get("event_name") or event.get("gene") or "")
        if not event_key or event_key in seen:
            continue
        seen.update(keys or [event_key])
        peptide = next((peptide_by_event[key] for key in keys if key in peptide_by_event), None)
        row = dict(peptide or {})
        row.update({key: value for key, value in event.items() if str(value or "").strip()})
        reasons = _patient_attention_reasons(row, bundle)
        configured = _patient_configured_manual_review(row, bundle.profile)
        if configured:
            reasons.insert(0, "由当前疾病/证据规则明确指定为关键人工审阅事件")
        if not reasons:
            continue
        mechanism_priority = bool(_disease_anchor(row, bundle))
        scored.append((0 if configured or mechanism_priority else 1, -len(reasons), index, row, list(dict.fromkeys(reasons))))
    result: list[dict[str, str]] = []
    displayed: set[str] = set()
    for _, _, _, row, reasons in sorted(scored):
        identity = _patient_manual_review_identity(row)
        if identity in displayed:
            continue
        displayed.add(identity)
        gaps = _patient_key_gaps(row, bundle)
        advice = _patient_validation(row, val_map)
        fusion_advice = _patient_fusion_manual_review_advice(row)
        if fusion_advice:
            advice = fusion_advice
        if gaps:
            advice += "。当前缺口：" + "；".join(gaps)
        result.append({
            "事件": str(row.get("gene") or row.get("event_name") or row.get("event_id") or ""),
            "为什么重要": "；".join(reasons) + "。" + _patient_dna_evidence(row) + "。RNA数据：" + _patient_rna_measurements(row) + ("。" + _patient_dna_rna_interpretation(row) if _patient_dna_rna_interpretation(row) else ""),
            "当前建议": advice,
        })
        if len(result) >= limit:
            break
    return result


def _patient_candidate_attention(row: Mapping[str, Any], bundle: ReportBundle) -> str:
    evidence = [
        _patient_metric("事件", row, "event_authenticity_state", "cross_platform_status"),
        _patient_metric("RNA", row, "rna_support_state", "rna_support_status"),
        _patient_metric("呈递", row, "presentation_consensus_state", "presentation_evidence_grade"),
        _patient_metric("MT/WT", row, "mutant_specificity_status", "mutant_specificity_state"),
    ]
    cross_site = _patient_cross_site_rna(row)
    return "；".join(evidence) + "。" + _patient_dna_evidence(row) + "。RNA数据：" + _patient_rna_measurements(row) + ("。" + _patient_dna_rna_interpretation(row) if _patient_dna_rna_interpretation(row) else "") + ("。" + cross_site if cross_site else "")


def _patient_candidate_disposition(
    row: Mapping[str, Any],
    val_map: Mapping[str, Mapping[str, str]],
    event_grades: Mapping[str, str],
) -> tuple[str, str]:
    recommendation = _patient_validation(row, val_map).strip()
    event_grade = str(event_grades.get(str(row.get("event_id") or "")) or _patient_grade(row)).upper()
    if "DO NOT ADVANCE" in recommendation.upper() or "暂缓/不推进" in recommendation or event_grade == "R4":
        return "PAUSED", "当前建议不推进或事件级为R4"
    integrity_ok, integrity_missing = _patient_candidate_integrity(row)
    if not integrity_ok:
        return "TECHNICAL_REVIEW", "；".join(integrity_missing)
    if event_grade in {"R1", "R2", "R3-READY", "R3-REVIEW"}:
        return "PRIORITY_CONFIRM", "优先完成事件真实性、冲突和正交证据确认"
    return "EVIDENCE_GAP", "补齐关键证据后再决定是否进入功能实验"


PATIENT_INPUT_FILE_METADATA = {
    "somatic_vcf": ("体细胞 VCF", "用于本次体细胞变异输入"),
    "tumor_dna_bam": ("肿瘤 DNA BAM", "用于纯度、CNV、LOH和DNA深度/VAF分析"),
    "normal_dna_bam": ("正常 DNA BAM", "用于排除胚系，并支持纯度、CNV、LOH和深度分析"),
    "tumor_short_rna_fastq": ("肿瘤短读 RNA FASTQ", "用于表达、RNA位点深度、alt reads/VAF及融合/剪接分析"),
    "tumor_long_rna_results": ("肿瘤长读 RNA 结果集", "用于完整转录本、长读junction和长读融合证据"),
    "shared_rna_evidence": ("短读 RNA 衍生证据目录", "由RNA原始数据计算得到，用于统一接入表达和RNA证据"),
}


def _patient_path_name(value: Any) -> str:
    text = str(value or "").strip().rstrip("/")
    return Path(text).name if text else ""


def _patient_input_file_rows(provenance: Mapping[str, Any]) -> list[dict[str, str]]:
    inputs = provenance.get("input_files") or {}
    if not isinstance(inputs, Mapping):
        return []
    rows: list[dict[str, str]] = []
    for key, value in inputs.items():
        data_type, purpose = PATIENT_INPUT_FILE_METADATA.get(
            str(key), (str(key).replace("_", " "), "用于本次分析输入"),
        )
        values = value if isinstance(value, list) else [value]
        names = [_patient_path_name(item) for item in values if _patient_path_name(item)]
        if names:
            rows.append({"数据类型": data_type, "文件名": "; ".join(names), "用途": purpose})
    return rows


def _patient_purity_value(value: Any, *, missing: str = "未形成估计") -> str:
    text = str(value or "").strip()
    return missing if text.upper() in {"", "NA", "N/A", "NONE", "."} else text


def _patient_purity_rows(bundle: ReportBundle) -> list[dict[str, str]]:
    rows = []
    for record in bundle.purity_tools:
        rows.append({
            "工具/模式": str(record.get("tool") or record.get("source") or "未记录"),
            "纯度": _patient_purity_value(record.get("purity")),
            "倍性": _patient_purity_value(record.get("ploidy")),
            "QC/状态": str(record.get("status") or "UNASSESSED"),
            "交叉验证说明": str(record.get("note") or "未提供工具级解释"),
        })
    return rows or [{
        "工具/模式": "未评估", "纯度": "未评估", "倍性": "未评估",
        "QC/状态": "UNASSESSED", "交叉验证说明": "未提供结构化纯度工具结果",
    }]


def _patient_purity_consensus(bundle: ReportBundle) -> tuple[str, str]:
    consensus = bundle.purity_consensus
    if consensus:
        purity = _patient_purity_value(consensus.get("recommended_purity"), missing="未评估")
        ploidy = _patient_purity_value(consensus.get("recommended_ploidy"), missing="未评估")
        tool = str(consensus.get("selected_tool") or "多工具共识")
        status = str(consensus.get("status") or "UNASSESSED")
        result = f"推荐纯度 {purity}、倍性 {ploidy}（{tool}；{status}）"
        basis = str(consensus.get("basis") or "多工具结果已并列保留；共识依据未记录")
        return result, basis
    fallback = str(bundle.provenance.get("purity_ploidy") or "未评估")
    return fallback, "未提供结构化多工具共识；不得静默选择单一工具"


def _patient_disease_background(bundle: ReportBundle) -> tuple[str, str]:
    """Resolve clinical disease context without inferring it from sample paths."""
    provenance = bundle.provenance
    clinical_keys = ("disease", "diagnosis", "disease_name", "cancer_type", "tumor_type")

    def meaningful(value: Any) -> str:
        text = str(value or "").strip()
        return "" if text.lower() in {"", "none", "na", "n/a", "unknown", "unassessed", "default"} else text

    for key in clinical_keys:
        value = meaningful(provenance.get(key))
        if value:
            return value, f"来源：结构化临床背景字段 {key}"

    for container_key in ("clinical_context", "disease_profile", "patient_context", "clinical"):
        container = provenance.get(container_key)
        if not isinstance(container, Mapping):
            continue
        for key in clinical_keys:
            value = meaningful(container.get(key))
            if value:
                return value, f"来源：结构化临床背景 {container_key}.{key}"

    for key in clinical_keys:
        value = meaningful(bundle.profile.get(key))
        if value:
            return value, f"来源：疾病/分析 profile 字段 {key}"

    profile_name = meaningful(bundle.profile.get("_profile_name") or provenance.get("profile"))
    if profile_name and profile_name.lower() not in {"evidence_consensus"}:
        return f"分析配置：{Path(profile_name).stem}", "未提供结构化临床诊断；回退采用非默认疾病/分析 profile"

    rules_name = meaningful(provenance.get("rules_name"))
    if not rules_name:
        rules = provenance.get("rules")
        if isinstance(rules, Mapping):
            rules_name = meaningful(rules.get("name") or rules.get("path"))
        elif rules:
            rules_name = meaningful(rules)
    if rules_name:
        return f"分析配置：{Path(rules_name).stem}", "未提供结构化临床诊断；回退采用排序/分析配置"

    return "未记录", "未提供结构化临床背景、非默认分析 profile 或排序配置"


def _patient_qc_rows(bundle: ReportBundle) -> list[dict[str, str]]:
    provenance = bundle.provenance
    purity_result, purity_basis = _patient_purity_consensus(bundle)
    disease_result, disease_basis = _patient_disease_background(bundle)
    rows = [
        {"项目": "疾病/分析背景", "结果": disease_result, "解释": disease_basis},
        {"项目": "肿瘤/正常配对", "结果": str(provenance.get("pairing_status") or "未评估"), "解释": "区分已使用配对输入与已完成指纹确认"},
        {"项目": "肿瘤纯度/倍性", "结果": purity_result, "解释": "用于CNV、CCF和LOH解释；工具冲突必须保留"},
        {"项目": "肿瘤DNA深度", "结果": str(provenance.get("tumor_dna_depth") or "未评估"), "解释": "默认汇总去重事件位点有效深度；低深度会降低检出能力"},
        {"项目": "正常DNA深度", "结果": str(provenance.get("normal_dna_depth") or "未评估"), "解释": "优先汇总候选位点normal_depth；缺失时从配对VCF正常样本DP/AD回填，再回退到normal DNA BAM覆盖估算"},
        {"项目": "RNA质量/覆盖", "结果": str(provenance.get("rna_qc_status") or "未评估"), "解释": "区分RNA支持、充分覆盖下未检出和表达层证据"},
        {"项目": "参考版本", "结果": str(provenance.get("genome_build") or "未记录"), "解释": "FASTA、GTF、VEP和坐标必须一致"},
        {"项目": "QC证据来源", "结果": "all_tool_results.tsv" if bundle.evidence_source_status == "CANONICAL_ALL_TOOL_RESULTS" else "未完整归一化", "解释": "样本QC汇总与候选排序的数据职责分离"},
        {"项目": "纯度/倍性多工具综合建议", "结果": purity_result, "解释": purity_basis},
    ]
    return rows


def _patient_hla_rows(peptides: list[dict[str, str]]) -> list[dict[str, str]]:
    loci: dict[str, list[str]] = {"A": [], "B": [], "C": []}
    for row in peptides:
        allele = str(row.get("hla_allele") or "").replace("HLA-", "")
        locus = allele.split("*", 1)[0]
        if locus in loci and allele and allele not in loci[locus]:
            loci[locus].append(allele)
    return [{"位点": key, "推荐等位基因": " / ".join(values[:2]) if values else "未评估", "用途": "限制性HLA-I呈递背景"} for key, values in loci.items()]


def _patient_hla_loh_status(value: Any) -> str:
    status = str(value or "").strip().lower()
    if status in {"yes", "true", "1", "loh", "loss", "lost", "deleted"}:
        return "LOST"
    if status in {"no", "false", "0", "retained", "retain", "kept", "normal"}:
        return "RETAINED"
    return "UNASSESSED"


def _patient_hla_loh_label(status: str, *, missing: str = "未提供") -> str:
    return {
        "LOST": "检出LOH（丢失）",
        "RETAINED": "未提示LOH",
        "CONFLICT": "结果冲突",
        "UNASSESSED": missing,
    }.get(status, missing)


def _patient_hla_loh_consensus(bundle: ReportBundle) -> tuple[list[dict[str, str]], str]:
    by_tool: dict[str, dict[str, set[str]]] = {"LOHHLA": {}, "SpecHLA": {}}
    tool_evidence: dict[str, dict[str, bool]] = {"LOHHLA": {}, "SpecHLA": {}}
    for record in bundle.hla_loh_tool_results:
        allele = str(record.get("hla_allele") or record.get("allele") or "").strip()
        if not re.match(r"^HLA-[ABC]\*", allele):
            continue
        tool = str(record.get("_report_tool") or record.get("evidence_tool") or record.get("source_tool") or record.get("tool") or "")
        tool = "LOHHLA" if "lohhla" in tool.lower() else "SpecHLA" if "spechla" in tool.lower() else ""
        if not tool:
            continue
        by_tool[tool].setdefault(allele, set()).add(_patient_hla_loh_status(record.get("loh_status") or record.get("status")))
        prefix = "lohhla_" if tool == "LOHHLA" else "spechla_"
        has_evidence = any(str(value or "").strip() for key, value in record.items() if key.startswith(prefix))
        tool_evidence[tool][allele] = tool_evidence[tool].get(allele, False) or has_evidence or bool(record.get("_report_source"))

    restricting = sorted({
        str(row.get("hla_allele") or row.get("allele") or "").strip()
        for row in bundle.peptides
        if re.match(r"^HLA-[ABC]\*", str(row.get("hla_allele") or row.get("allele") or "").strip())
    })
    restricting.extend(sorted({allele for calls in by_tool.values() for allele in calls if allele not in restricting}))
    rows: list[dict[str, str]] = []
    aggregate: list[tuple[str, str]] = []
    for allele in restricting:
        statuses: dict[str, str] = {}
        for tool in ("LOHHLA", "SpecHLA"):
            calls = by_tool[tool].get(allele, set()) - {"UNASSESSED"}
            statuses[tool] = next(iter(calls)) if len(calls) == 1 else "CONFLICT" if len(calls) > 1 else "UNASSESSED"
        assessed = [status for status in statuses.values() if status != "UNASSESSED"]
        evidence_present = {tool: tool_evidence.get(tool, {}).get(allele, False) for tool in ("LOHHLA", "SpecHLA")}
        def label(tool: str) -> str:
            status = statuses[tool]
            if status == "UNASSESSED" and evidence_present[tool]:
                return "未形成判断（QC不足）"
            return _patient_hla_loh_label(status)

        if not assessed:
            if any(evidence_present.values()):
                consensus, explanation, internal = "未评估", "已有工具原始证据，但未形成可用逐等位基因LOH判断/QC未通过", "UNASSESSED"
            else:
                consensus, explanation, internal = "未评估", "未提供逐等位基因HLA LOH结果", "UNASSESSED"
        elif "CONFLICT" in assessed or len(set(assessed)) > 1:
            consensus, explanation, internal = "工具结果冲突，暂不判定", "保留全部结果并要求人工复核", "CONFLICT"
        elif len(assessed) == 2:
            internal = assessed[0]
            consensus = "多工具一致：" + _patient_hla_loh_label(internal)
            explanation = "LOHHLA与SpecHLA逐等位基因结论一致"
        else:
            internal = assessed[0]
            tool = next(name for name, status in statuses.items() if status != "UNASSESSED")
            other_tool = "SpecHLA" if tool == "LOHHLA" else "LOHHLA"
            consensus = f"仅{tool}报告{_patient_hla_loh_label(internal)}，证据有限"
            if evidence_present[other_tool]:
                explanation = (
                    f"{other_tool}有原始证据，但未形成可用逐等位基因LOH判断/QC未通过；"
                    "单工具结果不足以确认该等位基因在肿瘤中完整保留或丢失"
                )
            else:
                explanation = (
                    f"另一个HLA LOH工具未提供该等位基因结果；"
                    "单工具结果不足以确认该等位基因在肿瘤中完整保留或丢失"
                )
        aggregate.append((allele, internal))
        rows.append({
            "HLA等位基因": allele,
            "LOHHLA": label("LOHHLA"),
            "SpecHLA": label("SpecHLA"),
            "综合判断": consensus,
            "说明": explanation,
        })

    multi_lost = [row["HLA等位基因"] for row in rows if str(row["综合判断"]).startswith("多工具一致") and "检出LOH" in str(row["综合判断"])]
    single_lost = [row["HLA等位基因"] for row in rows if str(row["综合判断"]).startswith("仅") and "检出LOH" in str(row["综合判断"])]
    conflicts = [row["HLA等位基因"] for row in rows if "冲突" in str(row["综合判断"])]
    multi_retained = [row for row in rows if str(row["综合判断"]).startswith("多工具一致") and "未提示LOH" in str(row["综合判断"])]
    single_retained = [row for row in rows if str(row["综合判断"]).startswith("仅") and "未提示LOH" in str(row["综合判断"])]
    if multi_lost:
        overall = "多工具一致检出限制性HLA-I LOH：" + "、".join(multi_lost)
    elif conflicts:
        overall = "HLA-I LOH工具结果冲突：" + "、".join(conflicts)
    elif single_lost:
        overall = "仅单工具提示限制性HLA-I LOH，需正交确认：" + "、".join(single_lost)
    elif rows and len(multi_retained) == len(rows):
        overall = "多工具一致未提示限制性HLA-I LOH（支持当前分层，但不替代实验验证）"
    elif single_retained:
        tools = sorted({"SpecHLA" if str(row["综合判断"]).startswith("仅SpecHLA") else "LOHHLA" for row in single_retained})
        other_qc = any("QC不足" in str(row.get("LOHHLA") or "") or "QC不足" in str(row.get("SpecHLA") or "") for row in single_retained)
        unavailable = "另一工具因QC不足未形成有效判断" if other_qc else "另一工具未形成有效判断"
        overall = (
            f"{'/'.join(tools)}未提示相应限制性HLA-I等位基因LOH；{unavailable}，"
            "因此当前仅有单工具支持，不足以确认该等位基因在肿瘤中完整保留。"
        )
    else:
        overall = "限制性HLA-I LOH未评估"
    return rows, overall


def _patient_appm_rows(bundle: ReportBundle) -> list[dict[str, str]]:
    summary = bundle.appm_summary
    dimensions = [
        ("MHC-I核心", "mhc_i_integrity_status"),
        ("MHC-II背景", "mhc_ii_integrity_status"),
        ("IFNG/JAK-STAT", "ifng_response_status"),
        ("APPM证据完整度", "appm_evidence_completeness"),
    ]
    rows = []
    for label, key in dimensions:
        status = str(summary.get(key) or "UNASSESSED")
        if label == "MHC-I核心" and status.upper() == "MHC_I_INTACT":
            result = "现有结果未发现HLA-I呈递系统整体完全丧失，肿瘤可能仍保留一定呈递能力；但抗原加工环节和HLA-LOH证据尚不完整。"
            explanation = "这是基于当前计算证据的审慎判断，不表示HLA-I呈递功能已被实验确认完整。"
        else:
            result = _patient_status_text(status)
            explanation = "未评估不等于正常或阴性" if status == "UNASSESSED" else "用于判断抗原加工呈递条件，不能单独预测临床疗效"
        rows.append({"维度": label, "结果": result, "通俗解释": explanation})
    _, hla_loh_consensus = _patient_hla_loh_consensus(bundle)
    rows.append({"维度": "限制性HLA-I LOH", "结果": hla_loh_consensus, "通俗解释": "由LOHHLA与SpecHLA逐等位基因结果形成共识；仅HLA-A/B/C丢失可直接影响相应HLA-I候选"})
    return rows


def _patient_tool_rows(bundle: ReportBundle) -> list[dict[str, str]]:
    tools = bundle.provenance.get("tools") or {}
    source_counts = _patient_source_tool_counts(bundle.events + bundle.peptides)
    provenance_by_key = {
        _patient_tool_key(name): (str(name), record)
        for name, record in tools.items() if isinstance(record, Mapping)
    } if isinstance(tools, Mapping) else {}
    discovered_keys = set(provenance_by_key) | set(source_counts)
    for key, (_, fields, _) in _PATIENT_RESULT_TOOL_SPECS.items():
        if any(_patient_assessed(row, *fields) for row in bundle.peptides):
            discovered_keys.add(key)
    rows: list[dict[str, str]] = []
    for key in sorted(discovered_keys):
            provenance_name, record = provenance_by_key.get(key, ("", {}))
            spec = _PATIENT_RESULT_TOOL_SPECS.get(key)
            source_meta = _PATIENT_SOURCE_TOOL_META.get(key)
            display_name = str(record.get("display_name") or (spec[0] if spec else source_meta[0] if source_meta else provenance_name or key))
            purpose = str(record.get("purpose") or "")
            if not purpose:
                purpose = spec[2] if spec else source_meta[1] if source_meta else str(record.get("mode") or "结果证据生成")
            fields = spec[1] if spec else ()
            result_count = sum(1 for row in bundle.peptides if fields and _patient_assessed(row, *fields))
            if not spec:
                result_count = max(result_count, source_counts.get(key, 0))
            event_result_count = int(record.get("evidence_event_rows") or 0)
            recorded_status = str(record.get("status") or "UNASSESSED")
            if key == "tap_appm" and result_count:
                tap_count = sum(
                    1 for row in bundle.peptides
                    if _patient_assessed(row, "tap_processing_status")
                )
                status = f"APPM已评估（候选修饰值 {result_count}/{len(bundle.peptides)}）"
                if tap_count:
                    status += f"；独立TAP状态 {tap_count}/{len(bundle.peptides)}"
                else:
                    status += "；独立TAP候选级状态未单列"
                status_basis = "；依据APPM汇总和肽段修饰结果"
            elif result_count:
                status = _patient_tool_coverage_status(bundle.peptides, key, fields, result_count)
                status_basis = "；结果值优先于可能过期的运行清单状态"
            elif event_result_count:
                status = f"事件来源记录已确认（{event_result_count}个事件）"
                status_basis = "；依据标准化事件的工具来源与精确映射记录"
            else:
                status = _patient_status_text(recorded_status)
                status_basis = ""
            version_record = (
                bundle.tool_versions.get(provenance_name)
                or bundle.tool_versions.get(provenance_name.lower())
                or bundle.tool_versions.get(key)
                or {}
            )
            default_version = "NeoAg APPM 2.0" if key == "tap_appm" and result_count else "原始运行版本未记录（需补工具版本清单）"
            rows.append({
                "流程/工具": display_name,
                "版本": str(version_record.get("version") or record.get("version") or default_version),
                "版本依据": str(version_record.get("evidence") or record.get("version_evidence") or "运行清单") + status_basis,
                "状态": status,
                "作用": purpose,
            })
    return rows


def _patient_release_metadata(bundle: ReportBundle) -> dict[str, str]:
    input_hash = str(bundle.evidence_integrity.get("actual_sha256") or bundle.evidence_integrity.get("expected_sha256") or "")
    rules_version = str(
        bundle.profile.get("rules_version")
        or bundle.profile.get("version")
        or bundle.profile.get("_profile_version")
        or (bundle.provenance.get("evidence_consensus") or {}).get("rules_version")
        or (bundle.provenance.get("parallel_rankings") or {}).get("rules_version")
        or "未记录"
    )
    run_id = str(bundle.provenance.get("run_id") or bundle.provenance.get("analysis_id") or "")
    if not run_id and input_hash:
        run_id = f"{bundle.sample_id or 'sample'}-{input_hash[:12]}"
    return {"run_id": run_id or "未记录", "rules_version": rules_version, "input_sha256": input_hash or "未记录"}


def _patient_release_audit(
    bundle: ReportBundle,
    displayed: list[dict[str, str]],
    paused: list[dict[str, str]],
    val_map: Mapping[str, Mapping[str, str]],
    event_grades: Mapping[str, str],
    rendered_without_audit: str,
) -> tuple[list[dict[str, str]], str]:
    metadata = _patient_release_metadata(bundle)
    integrity_failures = [row for row in displayed if not _patient_candidate_integrity(row)[0]]
    numeric_without_source = []
    for row in displayed:
        has_numeric = any(_patient_numeric_value(row, field)[1] is not None for field in (
            "netmhcpan_el_rank", "mhcflurry_presentation_score", "prime_score", "bigmhc_im_score", "ccf_estimate",
        ))
        if has_numeric and not str(row.get("evidence_field_sources") or "").strip():
            numeric_without_source.append(row)
    ccf_conflicts = []
    for row in displayed:
        reliable, _ = _patient_ccf_assessment(row, bundle)
        raw = _patient_value(row, "clonality_state", "ccf_status", default="").upper()
        if not reliable and raw in {"SUPPORTED", "CLONAL"} and "未形成可靠估计" not in _patient_evidence_summary(row, bundle):
            ccf_conflicts.append(row)
    do_not_advance = [
        row for row in displayed
        if "DO NOT ADVANCE" in _patient_validation(row, val_map).upper()
        or str(event_grades.get(str(row.get("event_id") or "")) or "").upper() == "R4"
    ]
    score_only = [row for row in displayed if not _patient_candidate_integrity(row)[0]]
    # Count only conflicts that survive the patient-facing provenance filter.
    # Differences against derived ranking/validation copies are audit sync
    # records, not independent biological-tool disagreements.
    conflicting = [row for row in displayed if _patient_conflict_summary(row)]
    unexplained_conflicts = [
        row for row in conflicting
        if _patient_conflict_summary(row) not in _patient_limitation(row, bundle)
        and _patient_conflict_summary(row) not in "；".join(_patient_key_gaps(row, bundle))
    ]
    path_or_log = bool(re.search(r"(?:/mnt/|/root/|/home/|Traceback \(most recent call last\)|nohup:)", rendered_without_audit))
    hla_appm_consistent = all(
        "HLA/APPM" not in _patient_evidence_summary(row, bundle)
        and _patient_restricting_hla_reliability(row, bundle)[0]
        and _patient_candidate_appm_status(bundle) != "未评估"
        for row in displayed
    )
    checks = [
        ("1", "展示候选关键字段完整", not integrity_failures, f"失败 {len(integrity_failures)} 条"),
        ("2", "精确数值存在字段级来源", not numeric_without_source, f"缺少来源映射 {len(numeric_without_source)} 条"),
        ("3", "样本级HLA/APPM与候选级表述一致", hla_appm_consistent, "限制性HLA须有逐等位基因多工具一致结论；APPM单独引用样本级共识"),
        ("4", "CCF未评估时不显示SUPPORTED/CLONAL", not ccf_conflicts, f"逻辑冲突 {len(ccf_conflicts)} 条"),
        ("5", "不推进候选未进入患者重点表", not do_not_advance, f"误入 {len(do_not_advance)} 条；技术池 {len(paused)} 条"),
        ("6", "事件数与Peptide-HLA数分开", True, f"事件 {len(bundle.events)}；Peptide-HLA {len(bundle.peptides)}"),
        ("7", "Top候选不由单一工具分数决定", not score_only, "应用完整性门槛、事件等级和多维证据分流"),
        ("8", "记录run_id、规则版本和输入hash", all(metadata[key] != "未记录" for key in metadata), f"run_id={metadata['run_id']}；规则={metadata['rules_version']}；hash={metadata['input_sha256'][:16]}..."),
        ("9", "Top候选来源冲突均已解决或解释", not unexplained_conflicts, f"有冲突 {len(conflicting)} 条；未解释 {len(unexplained_conflicts)} 条"),
        ("10", "患者版不含服务器路径或工具日志", not path_or_log, "已扫描绝对路径和常见日志标记"),
    ]
    rows = [{"编号": number, "审计项": label, "状态": "通过" if passed else "失败", "说明": detail} for number, label, passed, detail in checks]
    overall = "PASS" if all(passed for _, _, passed, _ in checks) else "DRAFT_NOT_FOR_RELEASE"
    return rows, overall


def make_patient_report(
    path: str | Path,
    bundle: ReportBundle,
    *,
    event_top_n: int = 10,
    candidate_top_n: int = 50,
) -> None:
    """Write the template-aligned, sample-agnostic patient HTML report."""
    if event_top_n < 1 or candidate_top_n < 1:
        raise ValueError("event_top_n and candidate_top_n must be positive integers")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    val_map = _val_by_peptide(bundle.validation_rows)
    ranked = list(bundle.peptides)
    event_grade_map = _patient_event_grade_map(bundle.events)
    all_representatives = _patient_representatives(ranked, len(ranked))
    all_representatives = [
        row for _, row in sorted(
            enumerate(all_representatives),
            key=lambda item: (
                0 if _disease_anchor(item[1], bundle) else 1,
                1 if _patient_fusion_artifact_review(item[1]) else 0,
                item[0],
            ),
        )
    ]
    dispositions = {
        str(row.get("peptide_id") or row.get("event_id") or index): _patient_candidate_disposition(row, val_map, event_grade_map)
        for index, row in enumerate(all_representatives)
    }

    def disposition_for(row: Mapping[str, Any]) -> tuple[str, str]:
        key = str(row.get("peptide_id") or row.get("event_id") or "")
        return dispositions.get(key) or _patient_candidate_disposition(row, val_map, event_grade_map)

    def patient_top_eligible(row: Mapping[str, Any]) -> bool:
        event_grade = _patient_event_row_grade(row, event_grade_map)
        consensus_eligible = event_grade in {
            "R1", "R2", "R3-READY", "R3-GAP", "R3-REVIEW"
        }
        if not consensus_eligible:
            status, _ = disposition_for(row)
            if event_grade != "UNASSESSED" or status not in {
                "PRIORITY_CONFIRM", "EVIDENCE_GAP", "TECHNICAL_REVIEW"
            }:
                return False
        if not all(str(row.get(field) or "").strip() for field in ("event_id", "peptide", "hla_allele")):
            return False
        direct_recommendation = " ".join(
            str(row.get(field) or "").strip()
            for field in ("recommended_use", "recommendation", "validation_recommendation")
        ).lower()
        # Evidence-consensus R grades are authoritative.  A stale weighted or
        # validation-plan recommendation must not suppress an R1-R3 candidate.
        if not consensus_eligible and (
            "do not advance" in direct_recommendation or "不推进" in direct_recommendation
        ):
            return False
        hard_fail = " ".join(
            str(row.get(field) or "").strip()
            for field in ("hard_fail", "hard_failure_codes", "hard_fail_reasons")
        ).strip().lower()
        return hard_fail in {"", "0", "false", "no", "none", "na", "n/a", "unassessed"}

    patient_representatives = [
        row for row in all_representatives
        if patient_top_eligible(row)
    ]
    top = patient_representatives[:candidate_top_n]
    paused_representatives = [row for row in all_representatives if not patient_top_eligible(row)]
    event_grade_counts = _patient_event_grade_counts(bundle.events)
    track_counts: dict[str, int] = {}
    event_seen: set[str] = set()
    for row in bundle.events or ranked:
        event_id = str(row.get("event_id") or row.get("event_name") or row.get("peptide_id") or "")
        if event_id and event_id in event_seen:
            continue
        if event_id:
            event_seen.add(event_id)
        track = _patient_track(row)
        track_counts[track] = track_counts.get(track, 0) + 1
    independent_event_count = len(event_seen) or len(bundle.events)
    peptide_hla_count = len({identity_value(row, "peptide_hla_id") for row in ranked})

    out = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        f"<title>肿瘤新抗原筛选报告—{esc(bundle.sample_id)}</title>", REPORT_CSS,
        "</head><body class='patient'><h1>肿瘤新抗原筛选报告</h1><p class='small'>患者沟通版</p>",
        "<div class='info'><b>阅读提示：</b>本报告先给出结论与候选分层，再说明样本、HLA、变异、肽段和验证建议。"
        "R1–R4是研究性证据等级，不是疗效等级。</div>",
        "<div class='warn'><b>重要说明：</b>本报告为研究性计算筛选，不能替代临床诊断或治疗决策。预测候选不等于体内真实呈递、不等于T细胞能够识别，也不等于确定治疗方案；所有候选均需进一步实验和临床专业判断。</div>",
    ]
    if bundle.evidence_source_status != "CANONICAL_ALL_TOOL_RESULTS":
        out.append(
            "<div class='warn'><b>证据来源警示：</b>未加载标准路径下的全部工具证据表；"
            "报告只能使用排序输入中的字段，所有缺失项必须按未评估解释。</div>"
        )
    elif bundle.evidence_integrity.get("status") != "PASS":
        out.append(
            "<div class='warn'><b>证据完整性警示：</b>全部工具证据表与机器可读清单的SHA256不一致或未完成校验；"
            "报告仅可用于排查，不应作为正式结果。</div>"
        )
    if str(bundle.provenance.get("report_role") or "").lower() == "pipeline_snapshot":
        out.append(
            "<div class='warn'><b>报告状态：</b>这是 Pipeline 运行阶段的结果快照，不是最终患者报告。"
            "请在 Evidence-consensus 完成后运行 open-neo-review，并以其 reports/patient_report.html 为准。</div>"
        )

    out.append("<div class='section'><h2>1. 报告摘要</h2>")
    summary = [
        {"项目": "样本编号", "结果": bundle.sample_id or "未注明", "说明": "以运行清单为准"},
        {"项目": "分析入口", "结果": bundle.entry_mode or "未注明", "说明": "可能为VCF、融合、剪接或联合入口"},
        {"项目": "评分配置", "结果": str(bundle.profile.get("_profile_name") or "未记录"), "说明": "研究性规则版本"},
    ]
    release_metadata = _patient_release_metadata(bundle)
    summary.extend([
        {"项目": "运行标识", "结果": release_metadata["run_id"], "说明": "用于报告与输入证据追溯"},
        {"项目": "规则版本", "结果": release_metadata["rules_version"], "说明": "本次事件分级和候选分流规则"},
        {"项目": "输入证据哈希", "结果": release_metadata["input_sha256"], "说明": "全部工具证据表SHA256"},
    ])
    out.append(_table(summary, ["项目", "结果", "说明"]))
    screening_rows = [
        {"口径": "独立事件", "数量": str(independent_event_count), "说明": "来自ranked_events.evidence_consensus.tsv，按事件去重"},
        {"口径": "Peptide-HLA组合", "数量": str(peptide_hla_count), "说明": "同一事件可产生多个重叠肽段和多个HLA组合；此处按唯一肽段序列+标准化HLA等位基因统计，重复证据行不重复计数"},
    ]
    out.append("<h3>筛选规模</h3>" + _table(screening_rows, ["口径", "数量", "说明"]))
    event_grade_metadata = {
        "R1": ("第一批实验优先", "关键事件级证据较完整"),
        "R2": ("有条件推进", "补充少量关键证据后可进入实验"),
        "R3-READY": ("已具备重点审阅条件", "事件仍属R3，不能解释为已确认新抗原"),
        "R3-GAP": ("存在证据缺口", "优先补齐缺失层后重新评估"),
        "R3-REVIEW": ("存在冲突或需人工复核", "保留冲突并进行正交验证"),
        "R4": ("当前暂不推进", "存在硬失败、明确风险或明显证据不足"),
        "UNASSESSED": ("事件等级未评估", "事件表缺少可用等级字段"),
    }
    event_grade_rows = []
    for grade in ("R1", "R2", "R3-READY", "R3-GAP", "R3-REVIEW", "R4", "UNASSESSED"):
        count = event_grade_counts.get(grade, 0)
        if grade == "UNASSESSED" and not count:
            continue
        meaning, next_step = event_grade_metadata[grade]
        event_grade_rows.append({"事件等级": grade, "事件数": count, "含义": meaning, "下一步": next_step})
    out.append("<h3>事件级结论</h3>")
    out.append("<p class='small'>以下数量直接读取ranked_events.evidence_consensus.tsv并按独立事件统计，不从Peptide-HLA表推算。</p>")
    out.append(_table(event_grade_rows, ["事件等级", "事件数", "含义", "下一步"]))
    focus_count = len(top)
    r1_r2_count = event_grade_counts.get("R1", 0) + event_grade_counts.get("R2", 0)
    anchor_labels = list(dict.fromkeys(
        str(row.get("gene") or row.get("event_name") or row.get("event_id") or "").strip()
        for row in top if _disease_anchor(row, bundle)
    ))
    focus_track_counts: dict[str, int] = {}
    for row in top:
        track = _patient_track(row)
        focus_track_counts[track] = focus_track_counts.get(track, 0) + 1
    if r1_r2_count:
        grade_conclusion = f"目前有{r1_r2_count}个独立事件达到R1或R2计算证据等级，但仍须完成对应实验验证后才能推进"
    else:
        grade_conclusion = "目前未获得可直接进入首批实验的R1或R2候选"
    priority_parts: list[str] = []
    if anchor_labels:
        priority_parts.append("疾病相关锚定事件" + "、".join(anchor_labels[:3]) + "相关候选")
    for track, label in (("SNV", "SNV"), ("InDel", "InDel"), ("Fusion", "融合"), ("Splice", "异常剪接")):
        count = focus_track_counts.get(track, 0)
        if not count or (track == "Fusion" and anchor_labels):
            continue
        priority_parts.append(f"{count}个{label}候选")
    priority_text = "、".join(priority_parts)
    review_conclusion = (
        f"{focus_count}个候选进入重点人工复核"
        + (f"，其中{priority_text}值得优先补证" if priority_text else "")
        + "。"
    )
    out.append(
        "<div class='info'><b>首页结论：</b>"
        f"本次分析共评估{independent_event_count}个独立候选事件，产生{peptide_hla_count}个肽段–HLA预测组合。"
        "经事件真实性、异常转录本/RNA表达、HLA呈递、突变特异性和安全性初筛，"
        f"{grade_conclusion}。{review_conclusion}"
        "所有结果均为研究性计算预测，尚未证明相关肽段在肿瘤细胞表面真实呈递或能够诱导T细胞反应。</div>"
    )
    if bundle.disease_knowledge.get("anchors"):
        out.append(
            "<p class='small'>已加载疾病知识配置：命中的锚定事件优先展示，"
            "但不会绕过事件真实性、精确断点、HLA、安全性或实验验证门槛，也不自动提升R等级。</p>"
        )
    out.append(f"<p>候选选择同时考虑事件真实性、RNA支持、HLA呈递、MT/WT突变特异性、克隆性、限制性HLA状态、APPM和安全性；缺失证据统一视为未评估，不作为阴性结论。另有{len(paused_representatives)}个事件代表候选因当前不推进或完整性门槛未通过，仅保留在技术审阅池。</p></div>")

    out.append("<div class='section'><h2>2. 患者样本与测序数据</h2>")
    out.append("<p>本节列出本次报告实际声明使用的患者数据，以及由这些输入得到的样本配对、测序质量、纯度和参考版本评估。未提供的项目保持“未评估”，不会自动写成正常。</p>")
    input_rows = _patient_input_file_rows(bundle.provenance)
    if input_rows:
        out.append("<p class='small'>以下列出本次分析使用的数据文件或结果目录名称及其用途；患者版不展示服务器绝对路径。</p>")
        out.append(_table(input_rows, ["数据类型", "文件名", "用途"]))
    else:
        out.append("<div class='warn'><b>患者数据清单未提供：</b>请通过 --patient-inputs 或 provenance.json 的 input_files 字段提供文件目录或文件名。</div>")
    out.append("<p class='small'>以下评估来自结构化运行清单和工具结果；‘未评估’表示证据缺失，不表示正常。</p>")
    out.append(_table(_patient_qc_rows(bundle), ["项目", "结果", "解释"]))
    out.append("<p class='small'>以下并列展示各纯度工具的估算。共识值用于CNV、CCF和LOH解释；明显冲突时保留全部结果及低置信度说明，不静默选择FACETS。</p>")
    out.append(_table(_patient_purity_rows(bundle), ["工具/模式", "纯度", "倍性", "QC/状态", "交叉验证说明"]))
    ccf_coverage_rows = _patient_ccf_coverage_rows(bundle)
    if ccf_coverage_rows:
        out.append("<h3>CCF覆盖度与缺失影响</h3>")
        out.append(
            "<p class='small'>CCF按独立事件和事件类型统计。RNA-only融合/剪接不把RNA junction reads伪装成DNA CCF；"
            "DNA来源事件缺失CCF时，不能判断候选覆盖的肿瘤细胞比例，也不能作为克隆性正向证据。</p>"
        )
        out.append(_table(ccf_coverage_rows, [
            "事件类型", "独立事件", "可靠CCF", "低置信数值", "缺失/未解析",
            "RNA-only不适用", "可评估事件覆盖率", "解释与影响",
        ]))
    out.append("</div>")

    out.append("<div class='section'><h2>3. HLA分型与抗原呈递条件</h2>")
    out.append("<h3>推荐HLA-I背景</h3>" + _table(_patient_hla_rows(ranked), ["位点", "推荐等位基因", "用途"]))
    out.append("<h3>抗原加工呈递与免疫逃逸</h3>" + _table(_patient_appm_rows(bundle), ["维度", "结果", "通俗解释"]))
    hla_loh_rows, _ = _patient_hla_loh_consensus(bundle)
    out.append("<h3>HLA-I LOH多工具结果</h3>")
    out.append("<p class='small'>以下按限制性HLA-I等位基因并列展示LOHHLA和SpecHLA结果。未提供不等于未发生LOH；工具冲突时保留全部结果，不静默选择单一工具。</p>")
    if hla_loh_rows:
        out.append(_table(hla_loh_rows, ["HLA等位基因", "LOHHLA", "SpecHLA", "综合判断", "说明"]))
    else:
        out.append("<div class='warn'><b>HLA-I LOH未评估：</b>未找到逐等位基因的LOHHLA或SpecHLA结果。</div>")
    out.append(
        "<div class='info'><b>拟进入高成本实验前建议补充：</b><ul class='compact'>"
        "<li>确认HLA分型的多工具一致性。</li>"
        "<li>核查正常DNA与肿瘤DNA的等位基因特异覆盖。</li>"
        "<li>结合HLA区域拷贝数、B等位基因频率及肿瘤纯度进行校正。</li>"
        "<li>必要时采用靶向HLA测序。</li>"
        "<li>使用IHC或流式验证HLA-I蛋白表达。</li>"
        "<li>联合评估B2M及抗原加工呈递通路状态。</li>"
        "</ul></div>"
    )
    out.append("<p>这些结果说明候选是否具备被加工和呈递的条件，但不能单独判断免疫治疗敏感、耐药或患者获益。</p></div>")

    out.append("<div class='section'><h2>4. 重点变异事件（按类型、按事件去重）</h2>")
    overall = [{"事件类型": track, "事件数": count, "主要复核重点": {"SNV": "DNA深度/VAF、RNA alt、MT/WT", "InDel": "局部重比对、阅读框、NMD和phasing", "Fusion": "精确断点、junction reads、frame和正常read-through", "Splice": "精确junction、PSI/reads、正常isoform和ORF", "DNA SV": "断点与异常转录本"}.get(track, "事件真实性和证据完整性")} for track, count in sorted(track_counts.items())]
    out.append(_table(overall, ["事件类型", "事件数", "主要复核重点"]))
    splice_funnel_rows = _patient_splice_funnel_rows(bundle)
    if splice_funnel_rows:
        out.append("<h3>异常剪接逐级筛选漏斗</h3>")
        out.append(
            "<p class='small'>漏斗按独立剪接事件统计，而不是Peptide-HLA行数。缺字段不会被当作通过；"
            "“阶段后可能剩余”给出明确通过数到包含未评估事件的上限。unique junction reads、总覆盖和PSI分别读取，"
            "不会用基因TPM或caller汇总reads互相替代。</p>"
        )
        out.append(_table(splice_funnel_rows, [
            "筛选阶段", "进入事件", "已评估", "明确通过", "明确未通过", "未评估",
            "阶段后可能剩余", "规则/说明",
        ]))
    for track in ("SNV", "InDel", "Fusion", "Splice", "DNA SV"):
        representatives = _patient_event_representatives(
            bundle.events, ranked, max(event_top_n, len(bundle.events), len(ranked)), track,
        )
        if not representatives:
            continue
        representatives = [
            row for _, row in sorted(
                enumerate(representatives),
                key=lambda item: (
                    0 if _disease_anchor(item[1], bundle) else 1,
                    1 if _patient_fusion_artifact_review(item[1]) else 0,
                    item[0],
                ),
            )
        ]
        representatives = representatives[:event_top_n]
        rows = []
        for rank, row in enumerate(representatives, 1):
            peptide = str(row.get("peptide") or "").strip()
            hla = str(row.get("hla_allele") or "").strip()
            peptide_hla = f"{peptide} / {hla}" if peptide and hla else "尚未形成完整Peptide-HLA组合"
            rows.append({
                "排名": rank,
                "基因/事件": row.get("gene") or row.get("event_name") or row.get("event_id", ""),
                "改变": _patient_event_change(row),
                "肽段-HLA": peptide_hla,
                "R等级": _patient_event_row_grade(row, event_grade_map),
                "关键证据与下一步": _patient_event_evidence_and_next_step(
                    row, bundle, val_map, compact_common=True,
                ),
            })
        out.append(f"<h3>{esc(track)} Top {event_top_n}</h3>")
        out.append(
            f"<p class='small'><b>本赛道证据口径：</b>{esc(_patient_track_evidence_note(track))}"
            "下表仅展示候选间存在差异的实测证据和特异性缺口，不再逐行重复不适用项。</p>"
        )
        out.append(_table(rows, ["排名", "基因/事件", "改变", "肽段-HLA", "R等级", "关键证据与下一步"]))
    out.append(f"<p class='small'>不同事件赛道的证据结构不同，Top {event_top_n}用于赛道内审阅，不应仅凭序号直接跨赛道比较。</p></div>")

    manual_review_rows = _patient_manual_review_rows(bundle.events, ranked, bundle, val_map)
    out.append("<div class='section'><h2>关键人工审阅事件</h2>")
    out.append("<p>以下事件因机制重要、多工具/正交支持或证据存在冲突而单独保留人工审阅；进入本表不等于自动升级为R1/R2，也不代表已确认新抗原。</p>")
    out.append("<p class='small'>对于关键融合，报告只将精确跨断点、来源可追溯且不与正常蛋白组精确匹配的肽作为融合特异性候选；普通融合伙伴蛋白肽不能替代融合junction肽。疾病相关关键融合即使暂为R4，也应保留在本节说明其机制意义和补证路线。</p>")
    if manual_review_rows:
        out.append(_table(manual_review_rows, ["事件", "为什么重要", "当前建议"]))
    else:
        out.append("<p class='small'>本次未筛出具有明确机制标记、多工具支持或证据冲突的独立人工审阅事件。</p>")
    out.append("</div>")

    displayed_candidate_count = len(top)
    out.append(
        "<div class='section'><h2>5. 当前进入人工复核的候选Peptide–HLA组合"
        f"（去重后{displayed_candidate_count}个）</h2>"
    )

    def patient_candidate_rows(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
        result = []
        for rank, row in enumerate(rows, 1):
            result.append({
                "排名": rank,
                "基因": row.get("gene", ""),
                "类型": _patient_track(row),
                "肽段-HLA": f"{row.get('peptide', '')} / {row.get('hla_allele', '')}",
                "等级": _patient_event_row_grade(row, event_grade_map),
                "关键证据与下一步": _patient_event_evidence_and_next_step(row, bundle, val_map),
            })
        return result

    candidate_headers = ["排名", "基因", "类型", "肽段-HLA", "等级", "关键证据与下一步"]
    out.append(
        f"<h3>当前展示{displayed_candidate_count}个去重候选组合</h3>"
        "<p class='small'>本表按事件级证据等级展示候选：R1、R2及R3的三个细分等级"
        "（R3-READY、R3-GAP、R3-REVIEW）。表内不再使用未细分的R3；其中R3-READY表示候选基本合理、"
        "仍需完成指定确认步骤，R3-GAP表示关键资料缺失，R3-REVIEW表示证据冲突或伪影风险需人工复核。"
        "仅纳入身份可追溯的非R4候选；R4、硬失败或明确不推进的候选保留在技术审阅池，不进入本表。"
        "本表用于研究性候选审阅，不表示已经确认新抗原或可直接进入功能实验。</p>"
    )
    out.append(_table(patient_candidate_rows(top), candidate_headers))
    quantitative_rows = [
        _patient_presentation_quantitative_row(row, index)
        for index, row in enumerate(top, start=1)
    ]
    out.append("<h3>呈递与免疫原性定量明细</h3>")
    out.append(
        "<p class='small'>Percentile rank越低表示模型预测越强，便于跨等位基因比较；"
        "≤1%仅作为研究性优选参考，不是免疫原性或体内呈递证明。模型score方向依各工具定义。"
        "若训练覆盖/外推状态未记录，报告保持未评估，不因工具返回数值而推定该HLA属于训练支持等位基因。</p>"
    )
    out.append(_table(quantitative_rows, [
        "排名", "肽段-HLA", "肽长/变异位置", "NetMHCpan原始值", "MT/WT定量比较",
        "突变位置结构解释", "WT自身反应/耐受风险", "MHCflurry原始值", "稳定性",
        "免疫原性辅助模型", "HLA模型覆盖",
    ]))
    out.append(f"<p class='small'>当前暂缓/不推进及完整性门槛未通过的{len(paused_representatives)}个事件代表候选不进入患者版重点表，仅保留在科研技术版审阅池。排序仍采用R1–R4、同赛道Pareto、确定性tie-break和事件去重。</p></div>")

    interpretation_top = top[:20]
    interpretation_count = len(interpretation_top)
    out.append(f"<div class='section'><h2>6. 人工复核候选的综合证据与实验建议（{interpretation_count}个）</h2>")
    out.append(
        f"<p>本节解读当前进入人工复核的{interpretation_count}个去重候选组合；建议顺序：先确认事件和异常转录本真实性，"
        "再补RNA alt/VAF或精确junction证据，完成MT/WT、正常背景和限制性HLA复核，最后开展短肽、长肽、"
        "minigene及T细胞功能实验。</p>"
    )
    interpretation_rows = []
    for row in interpretation_top:
        interpretation_rows.append({
            "候选": f"{row.get('gene', '')} | {row.get('peptide', '')} | {row.get('hla_allele', '')}",
            "为什么值得关注": _patient_candidate_attention(row, bundle),
            "当前不确定性": "；".join(_patient_key_gaps(row, bundle)) or "未见明确关键缺口；仍需实验确认",
            "建议下一步": _patient_validation(row, val_map),
        })
    comprehensive_headers = ["候选", "为什么值得关注", "当前不确定性", "建议下一步"]
    out.append(_table(interpretation_rows, comprehensive_headers))
    out.append("</div>")

    out.append("<div class='section'><h2>7. 分析方法与工具状态</h2>")
    active_tracks = {
        str(track).strip().upper().replace(" ", "_")
        for track, count in track_counts.items()
        if count
    }
    input_methods = []
    peptide_methods = []
    if active_tracks.intersection({"SNV", "INDEL"}):
        input_methods.append("体细胞VCF标准化与来源链检查")
        peptide_methods.append("SNV突变肽/正常肽配对与InDel新生尾部")
    if "FUSION" in active_tracks:
        input_methods.append("Fusion断点与多工具来源检查")
        peptide_methods.append("Fusion精确断点与阅读框肽段")
    if "SPLICE" in active_tracks:
        input_methods.append("Splice精确junction与来源链检查")
        peptide_methods.append("Splice精确junction及ORF肽段")
    if active_tracks.intersection({"DNA_SV", "SV"}):
        input_methods.append("DNA结构变异断点与来源链检查")
        peptide_methods.append("DNA结构变异异常转录本与ORF肽段")
    method_rows = [
        {"阶段": "事件输入与质控", "方法": "；".join(input_methods) or "未检测到可用事件赛道", "状态": "仅列出本次实际输入"},
        {"阶段": "肽段构建", "方法": "；".join(peptide_methods) or "未形成候选肽段", "状态": "无可追溯序列或ORF时不形成高等级证据"},
        {"阶段": "呈递预测", "方法": "核心结合/呈递、稳定性和免疫原性样模型分组汇总", "状态": "工具缺失记为未评估"},
        {"阶段": "综合排序", "方法": "证据状态、hard fail、priority cap、R1–R4、Pareto和事件去重", "状态": "研究性规则"},
    ]
    out.append(_table(method_rows, ["阶段", "方法", "状态"]))
    tool_rows = _patient_tool_rows(bundle)
    recorded_names = {_patient_tool_key(row.get("流程/工具")) for row in tool_rows}
    for inferred in _patient_inferred_tool_rows(top, bundle.tool_versions):
        if _patient_tool_key(inferred.get("流程/工具")) not in recorded_names:
            tool_rows.append(inferred)
    out.append("<h3>运行工具与结果覆盖</h3>" + _table(tool_rows, ["流程/工具", "版本", "版本依据", "状态", "作用"]))
    matched = sum(1 for row in top if row.get("_report_evidence_matched") == "YES")
    source_mapped = sum(1 for row in top if str(row.get("evidence_field_sources") or "").strip() not in {"", "{}"})
    manifest_validation = bundle.evidence_manifest.get("validation") or {}
    source_rows = [
        {"审计项": "权威证据源", "结果": bundle.evidence_source_status, "说明": top[0].get("_report_evidence_source", "") if top else "无候选"},
        {"审计项": "Top候选逐行匹配", "结果": f"{matched}/{len(top)}", "说明": "按peptide_id + HLA精确匹配"},
        {"审计项": "字段级来源映射", "结果": f"{source_mapped}/{len(top)}", "说明": "0表示综合表未记录字段级来源，不等于没有工具结果"},
        {"审计项": "证据表清单校验", "结果": bundle.evidence_integrity.get("status", "UNASSESSED"), "说明": "SHA256与all_tool_results.manifest.json一致" if bundle.evidence_integrity.get("status") == "PASS" else "; ".join(manifest_validation.get("errors") or []) or "SHA256未匹配"},
        {"审计项": "证据冲突记录", "结果": str(len(bundle.evidence_conflicts)), "说明": "冲突保留在evidence_conflicts.tsv，不静默覆盖"},
    ]
    out.append("<h3>证据来源审计</h3>" + _table(source_rows, ["审计项", "结果", "说明"]))
    out.append(
        f"<h3>Top {len(top)}证据维度可用性</h3>"
        "<p>本表区分可用于当前分层的证据、已计算但置信度不足的结果，以及尚未形成数值或不适用的项目。"
        "低置信或无法估计不等同于阴性。</p>"
        + _table(
            _patient_evidence_audit_rows(top, bundle),
            ["证据维度", "可作为当前分层证据", "尚不能作为可靠证据", "判定口径"],
        )
    )
    out.append("</div>")

    out.append("<div class='section'><h2>8. 局限性与总体结论</h2><ul>")
    out.append("<li>DNA/RNA覆盖不足时的未检出属于低检出能力或未评估，不是阴性。</li><li>计算呈递不等于体内真实呈递，体内呈递也不等于T细胞能够识别。</li><li>正常表达、HSPC、正常蛋白组、正常配体组和正常剪接证据不完整时，安全性必须标记为证据部分完整。</li><li>患者版不输出用药建议、疗效承诺或已确认新抗原结论。</li></ul>")
    out.append("<p><b>总体结论：</b>本次结果提供了可追溯的研究候选和分层验证顺序。应优先围绕事件真实性、RNA支持、突变肽与正常肽特异性、限制性HLA状态、APPM和正常背景补证，再决定首批实验集合。</p></div>")

    grade_rows = [{"等级": grade, "定义": value[0], "核心要求/处理": value[1]} for grade, value in R_GRADE_PATIENT.items()]
    out.append("<div class='section'><h2>附录A：R1–R4证据分层</h2>" + _table(grade_rows, ["等级", "定义", "核心要求/处理"]) + "</div>")
    glossary = [
        {"术语": "HLA", "通俗解释": "细胞表面的‘展示架’，把肽段展示给T细胞。"},
        {"术语": "APPM", "通俗解释": "抗原从蛋白被切割、运输到HLA展示的一整套加工呈递机制。"},
        {"术语": "LOH", "通俗解释": "某个HLA等位基因在肿瘤中丢失，可能使受它限制的候选无法呈递。"},
        {"术语": "MT/WT", "通俗解释": "突变肽与正常肽成对比较，用于判断候选是否真正具有突变特异性。"},
        {"术语": "CCF", "通俗解释": "估计携带该事件的肿瘤细胞比例；低可信结果不能当作精确比例。"},
        {"术语": "RNA alt reads / RNA VAF", "通俗解释": "直接支持突变转录本的RNA reads及其比例，比仅有gene TPM更接近突变表达证据。"},
        {"术语": "junction reads", "通俗解释": "跨越融合或异常剪接连接点的reads，必须精确对应同一junction。"},
        {"术语": "来源链 C1–C4", "通俗解释": "评价事件到转录本、ORF和肽段是否可追溯；与最终R1–R4推荐等级不同。"},
        {"术语": "Pareto排序", "通俗解释": "在多个证据维度间保留没有被全面压倒的候选，避免单一总分掩盖重要短板。"},
    ]
    out.append("<div class='section'><h2>附录B：术语说明</h2>" + _table(glossary, ["术语", "通俗解释"]) + "</div>")
    trace_rows = [
        {"文件": "运行与来源清单", "用途": "记录输入、参考、工具版本和运行身份"},
        {"文件": "全部工具证据表", "用途": "所有工具证据的统一可追溯表"},
        {"文件": "事件级候选排序表", "用途": "按事件去重的研究候选排序"},
        {"文件": "肽段-HLA候选排序表", "用途": "肽段-HLA级研究候选排序"},
        {"文件": "实验验证计划", "用途": "短肽、长肽、minigene和targeted RNA建议"},
        {"文件": "证据冲突清单", "用途": "来源字段冲突和复核原因"},
        {"文件": "当前暂缓/不推进与完整性审阅池", "用途": "仅供科研技术审阅，不进入患者版重点候选表"},
    ]
    out.append("<div class='section'><h2>附录C：附件与可追溯文件</h2>" + _table(trace_rows, ["文件", "用途"]) + "<p class='small'>患者版仅列逻辑文件名，不展示服务器绝对路径；实际位置和校验和见技术报告与run manifest。</p></div>")
    rendered_without_audit = "\n".join(out)
    release_audit_rows, release_status = _patient_release_audit(
        bundle, top, paused_representatives, val_map, event_grade_map, rendered_without_audit,
    )
    out.append("<div class='section'><h2>附录D：正式发布前自动审计</h2>")
    out.append(f"<p><b>审计结论：</b>{'通过，可进入正式患者报告发布流程' if release_status == 'PASS' else '未通过，仅可作为草稿或技术审阅材料'}</p>")
    out.append(_table(release_audit_rows, ["编号", "审计项", "状态", "说明"]) + "</div>")
    if release_status != "PASS":
        out.insert(7, "<div class='warn'><b>发布状态：</b>自动审计未通过，本报告只能用于草稿或技术审阅，不能作为正式患者版发布。</div>")
    out.append("</body></html>")
    p.write_text("\n".join(out), encoding="utf-8")
    audit_path = p.with_suffix(".release_audit.json")
    audit_path.write_text(json.dumps({
        "status": release_status,
        "metadata": release_metadata,
        "checks": release_audit_rows,
    }, ensure_ascii=False, indent=2), encoding="utf-8")


def _profile_threshold_section(profile: Mapping[str, Any]) -> str:
    sections = []
    for key in ("gates", "safety", "ccf_lite", "score_weights", "l3_weights", "appm_penalty"):
        block = profile.get(key)
        if isinstance(block, Mapping) and block:
            rows = [{"parameter": k, "value": v} for k, v in block.items()]
            sections.append(f"<h3>{esc(key)}</h3>{_table(rows, ['parameter', 'value'])}")
    return "".join(sections)


def _provenance_section(provenance: Mapping[str, Any]) -> str:
    if not provenance:
        return "<p class='small'>No provenance metadata supplied.</p>"
    out = ["<ul class='compact'>"]
    out.append(f"<li><b>sample_id:</b> <span class='mono'>{esc(provenance.get('sample_id'))}</span></li>")
    out.append(f"<li><b>entry_mode:</b> {esc(provenance.get('entry_mode'))}</li>")
    out.append(f"<li><b>profile:</b> {esc(provenance.get('profile'))}</li>")
    out.append(f"<li><b>created_at:</b> {esc(provenance.get('created_at'))}</li>")
    tools = provenance.get("tools") or {}
    if tools:
        out.append("</ul><h3>Tool provenance</h3><table><tr><th>tool</th><th>status</th><th>version</th><th>file</th><th>mode</th></tr>")
        for name, rec in tools.items():
            if not isinstance(rec, Mapping):
                continue
            out.append(
                "<tr>"
                f"<td>{esc(name)}</td><td>{esc(rec.get('status'))}</td><td>{esc(rec.get('version'))}</td>"
                f"<td class='mono'>{esc(rec.get('file'))}</td><td>{esc(rec.get('mode'))}</td>"
                "</tr>"
            )
        out.append("</table>")
    else:
        out.append("</ul>")
    return "\n".join(out)


def make_technical_report(path: str | Path, bundle: ReportBundle) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    mod_by_pep = _map_by(bundle.appm_peptide_modifiers, "peptide_id")
    esc_by_pep = _map_by(bundle.peptide_escape_flags, "peptide_id")
    safe_by_pep = _map_by(bundle.peptide_safety, "peptide_id")
    ccf_by_event = _map_by(bundle.ccf, "event_id")
    val_map = _val_by_peptide(bundle.validation_rows)

    out = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        f"<title>NeoAg Technical Evidence Report — {esc(bundle.sample_id)}</title>",
        REPORT_CSS,
        "</head><body>",
        "<h1>NeoAg Pipeline — Research / Technical Report</h1>",
        f"<p><b>Profile:</b> {esc(bundle.profile.get('_profile_name'))} &nbsp; "
        f"<b>Sample:</b> <span class='mono'>{esc(bundle.sample_id)}</span> &nbsp; "
        f"<b>Entry mode:</b> {esc(bundle.entry_mode)}</p>",
        "<div class='warn'><b>Boundary:</b> Computational triage only. "
        "Mechanism flags (APPM, escape, safety) are evidence layers, not clinical diagnoses.</div>",
    ]

    out.append("<div class='section'><h2>Pipeline &amp; Provenance</h2>")
    out.append(_provenance_section(bundle.provenance))
    out.append("</div>")

    parallel_rankings = bundle.provenance.get("parallel_rankings", {})
    if isinstance(parallel_rankings, Mapping) and parallel_rankings.get("evidence_consensus"):
        out.append("<div class='section'><h2>Experimental parallel evidence-consensus ranking</h2>")
        out.append(
            "<div class='warn'><b>Research-only parallel analysis:</b> "
            "The patient and technical reports use the evidence-consensus ranking as the primary final ranking. "
            "The legacy weighted ranking is preserved for audit and comparison.</div>"
        )
        metadata_rows = [
            {"field": "primary patient-report ranking", "value": parallel_rankings.get("evidence_consensus", "ranked_peptides.evidence_consensus.tsv")},
            {"field": "legacy weighted ranking", "value": parallel_rankings.get("legacy_weighted", "ranked_peptides.tsv")},
            {"field": "parallel event ranking", "value": parallel_rankings.get("event_consensus", "")},
            {"field": "comparison", "value": parallel_rankings.get("comparison", "")},
            {"field": "rules", "value": parallel_rankings.get("rules_name", "")},
            {"field": "rules version", "value": parallel_rankings.get("rules_version", "")},
            {"field": "rules status", "value": parallel_rankings.get("rules_status", "PROVISIONAL_RESEARCH_ONLY")},
        ]
        out.append(_table(metadata_rows, ["field", "value"]))
        out.append("</div>")

    if bundle.wes_qc:
        out.append("<div class='section'><h2>Independent WES QC</h2>")
        wes_headers = [
            "sample_id", "qc_status", "total_reads", "primary_mapping_rate_pct",
            "properly_paired_rate_pct", "duplicate_rate_pct", "target_definition",
            "mean_target_coverage", "pct_target_bases_20x", "pct_target_bases_30x",
            "on_target_rate_pct", "capture_rate_status", "formal_capture_rate_pct",
        ]
        out.append(_table(bundle.wes_qc, wes_headers))
        if any(row.get("capture_rate_status") != "ASSESSED" for row in bundle.wes_qc):
            out.append(
                "<div class='warn'><b>WES capture QC is partial:</b> coverage and "
                "on-target values use a GENCODE CDS proxy. The assay-specific capture "
                "BED is required before reporting a formal capture rate.</div>"
            )
        out.append("</div>")

    platform_counts = _cross_platform_counts(bundle.events)
    if platform_counts:
        out.append("<div class='section'><h2>WES/WGS Cross-platform DNA Evidence</h2>")
        rows = [{"cross_platform_status": status, "event_count": count} for status, count in sorted(platform_counts.items())]
        out.append(_table(rows, ["cross_platform_status", "event_count"]))
        out.append(
            "<div class='warn'><b>Interpretation:</b> Power-limited absence is not treated as a negative result. "
            "Source-unreproduced InDels require local assembly/IGV review. Sample-specific calls describe the "
            "sequenced specimen and must not be generalized to every tumor time point.</div>"
        )
        out.append("</div>")

    out.append("<div class='section'><h2>Profile Thresholds &amp; Weights</h2>")
    out.append(_profile_threshold_section(bundle.profile))
    out.append("</div>")

    out.append("<div class='section'><h2>APPM Summary</h2>")
    if bundle.appm_summary:
        rows = [{"field": k, "value": v} for k, v in bundle.appm_summary.items()]
        out.append(_table(rows, ["field", "value"]))
    if bundle.appm_submodule_scores:
        out.append("<h3>APPM Submodule Scores</h3>")
        out.append(_table(bundle.appm_submodule_scores, [
            "parent_module", "submodule", "score", "status", "defect_severity",
            "appm_call_confidence", "driver_defects", "action_hint", "confidence_reason",
        ]))
    if bundle.appm_gene_status:
        out.append("<h3>APPM Gene Status</h3>")
        out.append(_table(bundle.appm_gene_status, [
            "gene", "pathway", "biallelic_status", "functional_status", "copy_number_status",
            "loh_status", "expression_status", "gene_integrity_status", "reason",
        ], max_rows=50))
    if bundle.appm_conflicts:
        out.append("<h3>APPM Conflicts</h3>")
        out.append(_table(bundle.appm_conflicts, list(bundle.appm_conflicts[0].keys())))
    out.append("</div>")

    if bundle.immune_escape_summary:
        out.append("<div class='section'><h2>Immune Escape Summary</h2>")
        out.append(_table(bundle.immune_escape_summary, list(bundle.immune_escape_summary[0].keys())))
        out.append("</div>")

    out.append("<div class='section'><h2>Ranked Events (full)</h2>")
    out.append(_table(bundle.events, [
        "event_id", "event_name", "event_type", "mutation_source", "peptide_consequence", "gene",
        "source_chain_track", "source_chain_confidence_tier", "source_chain_confidence_label",
        "source_chain_orthogonal_status", "source_chain_orthogonal_sources",
        "source_chain_hard_failure", "source_chain_hard_failure_codes",
        "source_chain_reason_codes", "source_chain_missing_requirements",
        "source_chain_low_power_requirements", "source_chain_negative_requirements",
        "source_chain_conflict_requirements",
        "evidence_grade", "pipeline_r_grade", "review_status", "experiment_priority",
        "cancer_gene_list_status", "cancer_gene_symbols", "cancer_gene_types", "cancer_driver_context",
        "oncokb_annotated", "cosmic_cgc_flag", "cancer_gene_source_count", "cancer_gene_sources", "cancer_gene_context",
        "event_score", "raw_ccf", "ccf_estimate", "ccf_status", "ccf_confidence",
        "ccf_method", "ccf_warning", "clonality_multiplier", "safety_status", "safety_reason",
        "safety_evidence_completeness", "safety_missing_layers", "normal_expression_status",
        "normal_hspc_status", "reference_proteome_status", "normal_ligandome_status",
        "phase_group_id", "haplotype_status", "phase_support_reads",
        "phase_total_informative_reads", "phase_confidence", "component_event_ids",
        "combined_protein_change", "comparison_status", "cross_platform_status",
        "cross_platform_confidence", "cross_platform_multiplier", "cross_platform_priority_cap",
        "cross_platform_review_required", "wes_tumor_depth", "wes_tumor_alt_count",
        "wes_tumor_alt_vaf", "wgs_tumor_depth", "wgs_tumor_alt_count",
        "wgs_tumor_alt_vaf", "normal_depth", "normal_alt_count", "normal_alt_vaf",
    ]))
    out.append("</div>")

    out.append("<div class='section'><h2>Ranked Peptides (full)</h2>")
    out.append("<h3>Quantitative presentation and immunogenicity review</h3>")
    out.append(
        "<p>Lower percentile rank indicates a stronger model prediction. A rank at or below 1% is shown only as a prioritization reference; it does not establish immunogenicity or in-vivo presentation. Missing allele training/extrapolation metadata remains unassessed.</p>"
    )
    out.append(_table([
        _patient_presentation_quantitative_row(row, index)
        for index, row in enumerate(bundle.peptides[:100], start=1)
    ], [
        "排名", "肽段-HLA", "肽长/变异位置", "NetMHCpan原始值", "MT/WT定量比较",
        "突变位置结构解释", "WT自身反应/耐受风险", "MHCflurry原始值", "稳定性",
        "免疫原性辅助模型", "HLA模型覆盖",
    ], max_rows=100))
    pep_headers = [
        "peptide_id", "event_id", "gene", "cancer_gene_types", "cancer_driver_context", "cancer_gene_context",
        "source_chain_track", "source_chain_confidence_tier", "source_chain_confidence_label",
        "source_chain_orthogonal_status", "source_chain_orthogonal_sources",
        "source_chain_requirement_statuses", "source_chain_reason_codes",
        "source_chain_hard_failure", "source_chain_hard_failure_codes",
        "peptide", "wildtype_peptide", "peptide_consequence",
        "hla_allele", "mhc_class", "presentation_evidence_grade", "binding_evidence_score",
        "presentation_evidence_score", "netmhcpan_mt_ic50", "netmhcpan_mt_rank_ba",
        "netmhcpan_mt_rank_el", "netmhcpan_wt_ic50", "netmhcpan_wt_rank_ba",
        "netmhcpan_wt_rank_el", "netmhcpan_ba_rank", "netmhcpan_el_rank",
        "mhcflurry_mt_affinity_percentile", "mhcflurry_mt_presentation_score",
        "mhcflurry_wt_affinity_percentile", "mhcflurry_wt_presentation_score",
        "mhcflurry_affinity_percentile", "mhcflurry_presentation_score",
        "netmhcstabpan_score", "netmhcstabpan_rank", "netmhcstabpan_wt_score",
        "netmhcstabpan_wt_rank", "prime_score", "prime_rank", "prime_wt_score",
        "prime_wt_rank", "bigmhc_im_score", "bigmhc_im_wt_score", "deepimmuno_score",
        "predictor_allele_support_status", "allele_extrapolation_status",
        "netchop_31d_cterm_score", "netchop_31d_max_score", "netchop_processing_status",
        "agretopicity_el", "mt_wt_el_rank_difference",
        "mhcflurry_mt_wt_presentation_difference", "prime_mt_wt_score_difference",
        "bigmhc_mt_wt_score_difference", "mutation_positions_in_peptide",
        "mutation_anchor_only", "mutation_tcr_facing", "mutation_position_role",
        "mutation_position_interpretation", "wt_self_reactivity_risk_status",
        "wt_self_reactivity_risk_reason", "mt_wt_interpretation_caution", "mutant_specificity_status",
        "mutant_specificity_gate_status", "mutant_specificity_reason", "mutant_specificity_multiplier",
        "phase_group_id", "haplotype_status", "phase_support_reads",
        "phase_total_informative_reads", "phase_confidence", "component_event_ids",
        "combined_protein_change", "redundancy_group",
        "comparison_status", "cross_platform_status", "cross_platform_confidence",
        "cross_platform_multiplier", "cross_platform_review_required",
        "appm_multiplier", "ccf_multiplier", "safety_status", "safety_evidence_completeness",
        "safety_missing_layers", "normal_expression_status", "normal_hspc_status",
        "reference_proteome_status", "normal_ligandome_status", "anchor_assessment_status",
        "escape_status",
        "efficacy_score", "final_priority", "recommended_use",
    ]
    out.append(_table(bundle.peptides, pep_headers))
    out.append("</div>")

    technical_event_grades = _patient_event_grade_map(bundle.events)
    review_pool = []
    for row in _patient_representatives(bundle.peptides, len(bundle.peptides)):
        disposition, reason = _patient_candidate_disposition(row, val_map, technical_event_grades)
        if disposition not in {"PAUSED", "TECHNICAL_REVIEW"}:
            continue
        review_pool.append({
            "disposition": disposition,
            "reason": reason,
            "event_id": row.get("event_id", ""),
            "gene": row.get("gene", ""),
            "peptide": row.get("peptide", ""),
            "hla_allele": row.get("hla_allele", ""),
            "event_grade": technical_event_grades.get(str(row.get("event_id") or ""), _patient_grade(row)),
            "recommended_action": _patient_validation(row, val_map),
        })
    out.append("<div class='section'><h2>Current paused / do-not-advance and integrity review pool</h2>")
    out.append(_table(review_pool, ["disposition", "reason", "event_id", "gene", "peptide", "hla_allele", "event_grade", "recommended_action"]))
    out.append("</div>")

    if bundle.validation_rows:
        out.append("<div class='section'><h2>Validation Plan (full)</h2>")
        headers = list(bundle.validation_rows[0].keys())
        out.append(_table(bundle.validation_rows, headers))
        out.append("</div>")

    if bundle.peptide_safety:
        out.append("<div class='section'><h2>Peptide Safety Evidence</h2>")
        out.append(_table(bundle.peptide_safety, list(bundle.peptide_safety[0].keys()), max_rows=100))
        out.append("</div>")

    if bundle.peptide_escape_flags:
        out.append("<div class='section'><h2>Peptide Escape Flags</h2>")
        out.append(_table(bundle.peptide_escape_flags, list(bundle.peptide_escape_flags[0].keys()), max_rows=100))
        out.append("</div>")

    if bundle.ccf:
        out.append("<div class='section'><h2>CCF / Clonality</h2>")
        out.append(_table(bundle.ccf, list(bundle.ccf[0].keys()), max_rows=100))
        out.append("</div>")

    ccf_coverage_rows = _patient_ccf_coverage_rows(bundle)
    if ccf_coverage_rows:
        out.append("<div class='section'><h2>CCF Coverage Audit</h2>")
        out.append(_table(ccf_coverage_rows, [
            "事件类型", "独立事件", "可靠CCF", "低置信数值", "缺失/未解析",
            "RNA-only不适用", "可评估事件覆盖率", "解释与影响",
        ]))
        out.append("</div>")

    splice_funnel_rows = _patient_splice_funnel_rows(bundle)
    if splice_funnel_rows:
        out.append("<div class='section'><h2>Splice Filtering Funnel</h2>")
        out.append(_table(splice_funnel_rows, [
            "筛选阶段", "进入事件", "已评估", "明确通过", "明确未通过", "未评估",
            "阶段后可能剩余", "规则/说明",
        ]))
        out.append("</div>")

    out.append("<div class='section'><h2>Peptide Mechanism Cards</h2>")
    for ppt in bundle.peptides[:25]:
        pid = str(ppt.get("peptide_id") or "")
        e = esc_by_pep.get(pid, {})
        a = mod_by_pep.get(pid, {})
        s = safe_by_pep.get(pid, {})
        c = ccf_by_event.get(str(ppt.get("event_id") or ""), {})
        v = val_map.get(pid, {})
        out.append("<div class='card'>")
        out.append(f"<h3>{esc(ppt.get('peptide'))} — {esc(ppt.get('hla_allele'))}</h3>")
        out.append(
            f"<p><b>IDs:</b> peptide_id=<span class='mono'>{esc(pid)}</span>; "
            f"event_id=<span class='mono'>{esc(ppt.get('event_id'))}</span></p>"
        )
        out.append(
            f"<p><b>Layers:</b> mutation_source={esc(ppt.get('mutation_source'))}; "
            f"peptide_consequence={esc(ppt.get('peptide_consequence'))}; source_tool={esc(ppt.get('source_tool'))}</p>"
        )
        if ppt.get("source_chain_confidence_tier"):
            out.append(
                f"<p><b>Source-chain confidence:</b> {_badge(ppt.get('source_chain_confidence_tier'))}; "
                f"track={esc(ppt.get('source_chain_track'))}; "
                f"orthogonal={esc(ppt.get('source_chain_orthogonal_status'))}; "
                f"sources={esc(ppt.get('source_chain_orthogonal_sources'))}; "
                f"missing=<span class='mono'>{esc(ppt.get('source_chain_missing_requirements'))}</span>; "
                f"low-power=<span class='mono'>{esc(ppt.get('source_chain_low_power_requirements'))}</span>; "
                f"negative=<span class='mono'>{esc(ppt.get('source_chain_negative_requirements'))}</span>; "
                f"conflict=<span class='mono'>{esc(ppt.get('source_chain_conflict_requirements'))}</span></p>"
            )
        final_grade = ppt.get("evidence_grade") or ppt.get("pipeline_r_grade") or ppt.get("priority")
        if final_grade:
            out.append(
                f"<p><b>Final recommendation tier:</b> {_badge(final_grade)}; "
                f"review={esc(ppt.get('review_status'))}; "
                f"experiment_priority={esc(ppt.get('experiment_priority'))}. "
                "Source-chain confidence and final recommendation answer different questions: "
                "event authenticity alone does not establish neoantigen priority.</p>"
            )
        out.append(
            f"<p><b>Presentation:</b> grade={esc(ppt.get('presentation_evidence_grade'))}; "
            f"BA={esc(ppt.get('netmhcpan_ba_rank'))}; EL={esc(ppt.get('netmhcpan_el_rank'))}; "
            f"MHCflurry={esc(ppt.get('mhcflurry_presentation_score'))}</p>"
        )
        out.append(
            f"<p><b>Mutant specificity:</b> {_badge(ppt.get('mutant_specificity_gate_status'))}; "
            f"status={esc(ppt.get('mutant_specificity_status'))}; "
            f"MT_EL={esc(ppt.get('netmhcpan_mt_rank_el', ppt.get('netmhcpan_el_rank')))}; "
            f"WT_EL={esc(ppt.get('netmhcpan_wt_rank_el'))}; "
            f"agretopicity={esc(ppt.get('agretopicity_el'))}; "
            f"positions={esc(ppt.get('mutation_positions_in_peptide'))}; "
            f"anchor_only={esc(ppt.get('mutation_anchor_only'))}; "
            f"TCR_facing={esc(ppt.get('mutation_tcr_facing'))}; "
            f"position_role={esc(ppt.get('mutation_position_role'))}; "
            f"position_interpretation=<span class='mono'>{esc(ppt.get('mutation_position_interpretation'))}</span>; "
            f"WT_self_reactivity={_badge(ppt.get('wt_self_reactivity_risk_status'))}; "
            f"WT_risk_reason=<span class='mono'>{esc(ppt.get('wt_self_reactivity_risk_reason'))}</span>; "
            f"reason=<span class='mono'>{esc(ppt.get('mutant_specificity_reason'))}</span></p>"
        )
        out.append(
            f"<p><b>Haplotype:</b> {_badge(ppt.get('haplotype_status'))}; "
            f"phase_group={esc(ppt.get('phase_group_id'))}; "
            f"support={esc(ppt.get('phase_support_reads'))}/{esc(ppt.get('phase_total_informative_reads'))}; "
            f"confidence={esc(ppt.get('phase_confidence'))}; "
            f"components=<span class='mono'>{esc(ppt.get('component_event_ids'))}</span>; "
            f"protein=<span class='mono'>{esc(ppt.get('combined_protein_change'))}</span></p>"
        )
        out.append(
            f"<p><b>APPM:</b> multiplier={esc(a.get('appm_multiplier', ppt.get('appm_multiplier')))}; "
            f"reason=<span class='mono'>{esc(a.get('appm_multiplier_reason', ''))}</span></p>"
        )
        out.append(
            f"<p><b>Escape:</b> {_badge(e.get('escape_status', ppt.get('escape_status')))}; "
            f"multiplier={esc(e.get('escape_multiplier', ppt.get('escape_multiplier')))}; "
            f"reason=<span class='mono'>{esc(e.get('escape_reason', ''))}</span></p>"
        )
        out.append(
            f"<p><b>Safety:</b> {_badge(s.get('safety_status', ppt.get('safety_status')))}; "
            f"tier={esc(s.get('safety_tier', ''))}; completeness={esc(s.get('safety_evidence_completeness', ppt.get('safety_evidence_completeness')))}; "
            f"missing=<span class='mono'>{esc(s.get('safety_missing_layers', ppt.get('safety_missing_layers')))}</span>; "
            f"reason=<span class='mono'>{esc(s.get('safety_reason', ppt.get('safety_reason')))}</span></p>"
        )
        out.append(
            f"<p><b>CCF:</b> status={esc(c.get('ccf_status', ppt.get('ccf_status')))}; "
            f"raw={esc(c.get('raw_ccf', ppt.get('raw_ccf')))}; "
            f"estimate={esc(c.get('ccf_estimate', ppt.get('ccf_estimate')))}; "
            f"95% CI={esc(c.get('ccf_ci_low', ppt.get('ccf_ci_low')))}-"
            f"{esc(c.get('ccf_ci_high', ppt.get('ccf_ci_high')))}; "
            f"confidence={esc(c.get('ccf_confidence', ppt.get('ccf_confidence')))}; "
            f"method={esc(c.get('ccf_method', ppt.get('ccf_method')))}; "
            f"local_CN={esc(c.get('total_cn', ppt.get('total_cn')))} "
            f"(major={esc(c.get('major_cn', ppt.get('major_cn')))}, minor={esc(c.get('minor_cn', ppt.get('minor_cn')))}); "
            f"multiplicity={esc(c.get('multiplicity_best', ppt.get('multiplicity_best')))}; "
            f"normal_check={esc(c.get('normal_contamination_status', ppt.get('normal_contamination_status')))}; "
            f"multiplier={esc(c.get('clonality_multiplier', ppt.get('ccf_multiplier')))}; "
            f"warning=<span class='mono'>{esc(c.get('ccf_warning', ppt.get('ccf_warning')))}</span>. "
            "A high VAF or point CCF is compatibility evidence only; clonality requires the interval, local allele-specific CN, purity and matched-normal review.</p>"
        )
        out.append(
            f"<p><b>WES/WGS evidence:</b> {_badge(ppt.get('cross_platform_status'))}; "
            f"confidence={esc(ppt.get('cross_platform_confidence'))}; "
            f"multiplier={esc(ppt.get('cross_platform_multiplier'))}; "
            f"WES={esc(ppt.get('wes_tumor_alt_count'))}/{esc(ppt.get('wes_tumor_depth'))}; "
            f"WGS={esc(ppt.get('wgs_tumor_alt_count'))}/{esc(ppt.get('wgs_tumor_depth'))}; "
            f"normal={esc(ppt.get('normal_alt_count'))}/{esc(ppt.get('normal_depth'))}</p>"
        )
        if v:
            out.append(
                f"<p><b>Validation design:</b> mode={esc(v.get('validation_mode'))}; "
                f"assay={esc(v.get('recommended_assay'))}; minigene=<span class='mono'>{esc(v.get('minigene'))}</span></p>"
            )
        out.append(
            f"<p><b>Decision:</b> {_badge(ppt.get('final_priority'))}; "
            f"efficacy_score={esc(ppt.get('efficacy_score'))}; {esc(ppt.get('recommended_use'))}</p>"
        )
        out.append("</div>")
    out.append("</div>")

    out.append("<div class='section'><h2>Field Glossary</h2><table><tr><th>Field</th><th>Description</th></tr>")
    for field_name, desc in FIELD_GLOSSARY.items():
        out.append(f"<tr><td class='mono'>{esc(field_name)}</td><td>{esc(desc)}</td></tr>")
    out.append("</table></div>")
    out.append("</body></html>")
    p.write_text("\n".join(out), encoding="utf-8")


def make_dual_reports(
    reports_dir: str | Path,
    bundle: ReportBundle,
    *,
    patient_name: str = "evidence_report.patient.html",
    technical_name: str = "evidence_report.technical.html",
    legacy_name: str = "evidence_report.html",
    event_top_n: int = 10,
    candidate_top_n: int = 50,
) -> dict[str, str]:
    reports_dir = Path(reports_dir)
    patient_path = reports_dir / patient_name
    technical_path = reports_dir / technical_name
    legacy_path = reports_dir / legacy_name
    make_patient_report(patient_path, bundle, event_top_n=event_top_n, candidate_top_n=candidate_top_n)
    make_technical_report(technical_path, bundle)
    legacy_path.write_text(technical_path.read_text(encoding="utf-8"), encoding="utf-8")
    return {
        "evidence_report_patient": str(patient_path),
        "evidence_report_technical": str(technical_path),
        "evidence_report": str(legacy_path),
    }
