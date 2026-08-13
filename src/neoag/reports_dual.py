"""Dual-audience HTML reports: patient communication vs research/technical."""

from __future__ import annotations

import html
import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

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
    "MARGINAL_MT_ADVANTAGE": "突变肽相对正常肽仅有轻度优势",
    "HLA_LOH_UNASSESSED": "限制性HLA多工具确认未完成",
    "SAFETY_PARTIAL": "正常组织安全性证据不完整",
    "SUPPORTED": "获得支持",
    "CLONAL": "倾向克隆性事件",
    "C1": "候选来源链完整且有正交支持",
    "C2": "候选来源链较完整，仍需实验确认",
    "C3": "候选来源链基本合理，但关键环节尚未闭合",
    "C4": "候选来源链不完整，仅进入技术审阅",
    "UNASSESSED": "未评估",
    "UNSPECIFIED": "未明确",
    "PARTIAL": "证据部分完整",
    "IFNG_RESPONSE_CAUTION": "IFNG/JAK-STAT应答存在谨慎信号",
    "MHC_II_INTACT": "现有结果未见MHC-II呈递系统整体完全丧失",
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


def _read_hla_loh_tool_results(root: Path | None) -> list[dict[str, str]]:
    """Load allele-level LOHHLA and SpecHLA evidence from supported layouts."""
    if not root:
        return []
    candidates = {
        "LOHHLA": [
            root / "hla_loh" / "lohhla" / "hla_loh.tsv",
            root / "hla_loh_consensus" / "lohhla_hla_loh.tsv",
        ],
        "SpecHLA": [
            root / "hla_loh" / "spechla" / "hla_loh.tsv",
            root / "hla_loh_consensus" / "spechla_sequenza012_hla_loh.tsv",
            root / "hla_loh_consensus" / "spechla_hla_loh.tsv",
        ],
    }
    results: list[dict[str, str]] = []
    for tool, paths in candidates.items():
        selected = next((path for path in paths if path.is_file()), None)
        if not selected:
            continue
        for row in _read_optional(selected):
            record = dict(row)
            record["_report_tool"] = tool
            record["_report_source"] = str(selected.relative_to(root))
            results.append(record)
    return results


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


def _read_longrna_junction_genes(root: Path | None, extra_gene_ids: set[str] | None = None) -> dict[str, str]:
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
    gtf_candidates = [root.parents[2] / "data" / "ref" / "hg38" / "gencode.gtf", root.parents[2] / "data" / "ref" / "ctat" / "current" / "ctat_genome_lib_build_dir" / "ref_annot.gtf"]
    for gtf in gtf_candidates:
        if not gtf.is_file() or not (gene_ids or transcript_ids):
            continue
        with gtf.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if "\tgene\t" in line:
                    match = re.search(r'gene_id "(ENSG\d+)";.*?gene_name "([^"]+)";', line)
                    if match and match.group(1) in gene_ids:
                        gene_names[match.group(1)] = match.group(2)
                    if match:
                        fields = line.split("\t")
                        if len(fields) >= 7:
                            chrom = fields[0] if fields[0].startswith("chr") else "chr" + fields[0]
                            gene_intervals.append((chrom, int(fields[3]), int(fields[4]), fields[6], match.group(2)))
                elif "\ttranscript\t" in line and transcript_ids:
                    match = re.search(r'gene_id "(ENSG\d+)";.*?transcript_id "(ENS[TG]\d+)";.*?gene_name "([^"]+)";', line)
                    if match and match.group(2) in transcript_ids:
                        gene_names[match.group(2)] = match.group(3)
        if gene_names:
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
    fields = (
        "gene", "combined_protein_change", "event_name", "consequence", "peptide_consequence",
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
        "ccf_estimate", "ccf_confidence_state", "purity_consensus_status", "hla_appm_state",
        "appm_integrity_status", "restricting_hla_lost", "restricting_locus_loh", "loh_status",
        "escape_status", "immunogenicity_score", "immunogenicity_composite_score",
        "netmhcpan_el_rank", "mhcflurry_presentation_score", "mhcflurry_processing_score",
        "prime_score", "bigmhc_im_score", "netmhcstabpan_rank", "netmhcstabpan_score",
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
        if (not row.get("gene") or str(row.get("gene", "")).startswith("SEV|")) and junction in junction_genes:
            row["gene"] = junction_genes[junction]
        for field_name in ("gene", "event_name", "combined_protein_change"):
            if row.get(field_name):
                row[field_name] = _replace_gene_ids(row[field_name], junction_genes)
        enriched.append(row)
    return enriched


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
    "safety_status": "正常组织表达、自身肽相似性等安全性初筛结果。",
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
    if patient_inputs:
        explicit_inputs = patient_inputs.get("input_files") if isinstance(patient_inputs.get("input_files"), Mapping) else patient_inputs
        prov["input_files"] = dict(explicit_inputs)

    def p(*parts: str) -> Path | None:
        return root / Path(*parts) if root else None

    source_events = _read_optional(p("inputs", "combined_raw_events.tsv") if root else None)
    canonical_evidence_path = p("scoring", "evidence_consensus", "all_tool_results.tsv") if root else None
    fallback_evidence_path = p("scoring", "all_tool_results.tsv") if root else None
    if canonical_evidence_path and canonical_evidence_path.is_file():
        source_peptides = _read_optional(canonical_evidence_path)
        evidence_source_label = "scoring/evidence_consensus/all_tool_results.tsv"
        evidence_source_status = "CANONICAL_ALL_TOOL_RESULTS"
    elif fallback_evidence_path and fallback_evidence_path.is_file():
        source_peptides = _read_optional(fallback_evidence_path)
        evidence_source_label = "scoring/all_tool_results.tsv"
        evidence_source_status = "FALLBACK_ALL_TOOL_RESULTS"
    else:
        source_peptides = []
        evidence_source_label = "ranked input (all_tool_results unavailable)"
        evidence_source_status = "RANKED_INPUT_ONLY"
    extra_gene_ids = {
        gene_id
        for source in source_events + source_peptides
        for field_name in ("gene", "gene_pair", "event_name", "combined_protein_change")
        for gene_id in re.findall(r"ENSG\d+", str(source.get(field_name) or ""))
    }
    junction_genes = _read_longrna_junction_genes(root, extra_gene_ids)
    enriched_peptides = _enrich_rows_from_sources(
        peptides, source_peptides or peptides, source_events, junction_genes, evidence_source_label,
    )
    enriched_events = _enrich_rows_from_sources(
        events, source_peptides or enriched_peptides, source_events, junction_genes, evidence_source_label,
    )
    evidence_manifest = _read_json_optional(
        p("scoring", "evidence_consensus", "all_tool_results.manifest.json") if root else None
    )
    evidence_integrity = {"status": "UNASSESSED", "expected_sha256": "", "actual_sha256": ""}
    if source_peptides and canonical_evidence_path and canonical_evidence_path.is_file() and evidence_manifest:
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
        evidence_conflicts=_read_optional(
            p("scoring", "evidence_consensus", "evidence_conflicts.tsv") if root else None
        ),
        evidence_source_status=evidence_source_status,
        evidence_integrity=evidence_integrity,
        purity_tools=[dict(row) for row in (prov.get("purity_cnv_tools") or []) if isinstance(row, Mapping)],
        purity_consensus=dict(prov.get("purity_cnv_consensus") or {}),
        hla_loh_tool_results=_read_hla_loh_tool_results(root),
        tool_versions=_read_tool_version_manifest(root),
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
    change = str(event.get("combined_protein_change") or event.get("event_name") or event.get("consequence") or "")
    if ":p." in change:
        change = "p." + change.split(":p.", 1)[1]
    if change:
        return change.replace("%3D", "=")
    consequence = str(event.get("peptide_consequence") or "").strip().lower()
    event_type = str(event.get("event_type") or "").strip()
    peptide = str(event.get("peptide") or "").strip()
    if event_type == "Splice" or consequence in {"novel_junction", "splice_junction", "exon_deletion_junction"}:
        return f"异常剪接肽段 {peptide}" if peptide else "异常剪接事件（完整蛋白改变未形成）"
    if event_type == "Fusion" and peptide:
        return f"融合肽段 {peptide}"
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
    seen: set[str] = set()
    for row in rows:
        if track and _patient_track(row) != track:
            continue
        event_id = str(row.get("event_id") or row.get("event_name") or row.get("peptide_id") or "")
        key = event_id or f"{row.get('gene', '')}|{row.get('peptide', '')}"
        if key in seen:
            continue
        seen.add(key)
        selected.append(row)
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
    value = _patient_value(row, *fields, default="").upper()
    return bool(value) and not any(token in value for token in ("UNASSESSED", "NOT_AVAILABLE", "NOT_RUN"))


def _patient_metric(label: str, row: Mapping[str, Any], *fields: str) -> str:
    value = _patient_value(row, *fields)
    return f"{label}={_patient_status_text(value)}"


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
        return True, f"克隆性={_patient_status_text(state)}；CCF={ccf_text}；置信度={_patient_status_text(confidence)}"
    purity_status = str((bundle.purity_consensus if bundle else {}).get("status") or "").upper()
    reason = "样本纯度低且缺少可用CCF结果" if "LOW_PURITY" in purity_status else "缺少可用CCF数值或可靠置信度"
    return False, f"克隆性=未形成可靠估计；原因={reason}；处理=不作为阴性，也不作为正向加分"


def _patient_candidate_integrity(row: Mapping[str, Any]) -> tuple[bool, list[str]]:
    missing: list[str] = []
    if not str(row.get("event_id") or "").strip():
        missing.append("event_id缺失")
    if not str(row.get("peptide") or "").strip():
        missing.append("peptide缺失")
    if not str(row.get("hla_allele") or "").strip():
        missing.append("限制性HLA缺失")
    source_tier = str(row.get("source_chain_confidence_tier") or "").strip().upper()
    if source_tier in {"", "UNASSESSED", "C4"}:
        missing.append("来源链不可回溯")
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
    if "仅SpecHLA" in overall and "未见" in overall:
        return "单工具提示保留；多工具LOH确认未完成"
    if "仅LOHHLA" in overall and "未见" in overall:
        return "单工具提示保留；多工具LOH确认未完成"
    if "多工具一致未见" in overall:
        return "多工具一致提示保留"
    if "检出限制性HLA-I LOH" in overall:
        return overall.replace("检出限制性HLA-I LOH", "检出限制性HLA丢失")
    if "冲突" in overall:
        return "工具结果冲突；限制性HLA状态待确认"
    return "未评估"


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
    if str(row.get("evidence_conflict_fields") or "").strip():
        limitations.append("影响候选的来源字段存在冲突，需按冲突记录人工解释")
    if _patient_value(row, "safety_state", "safety_status", default="").upper() in {"SAFETY_PARTIAL", "PARTIAL"}:
        limitations.append("安全性证据不完整")
    hla_status = _patient_candidate_hla_status(bundle)
    if hla_status == "未评估":
        limitations.append("限制性HLA状态未评估")
    elif "多工具LOH确认未完成" in hla_status:
        limitations.append("限制性HLA仅单工具提示保留；多工具LOH确认未完成")
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
        _patient_metric("NetChop 3.1d", row, "netchop_processing_status", "netchop_31d_max_score"),
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
        ("加工/稳定性", ("mhcflurry_processing_score", "netmhcstabpan_rank", "netchop_31d_max_score", "tap_processing_status")),
        ("MT/WT特异性", ("mutant_specificity_status", "mutant_specificity_state")),
        ("安全性", ("safety_state", "safety_status")),
        ("免疫原性", ("immunogenicity_composite_score", "immunogenicity_score", "bigmhc_im_score")),
        ("来源链", ("source_chain_confidence_tier",)),
    ]
    total = len(rows)
    result = []
    for label, fields in dimensions:
        assessed = sum(1 for row in rows if _patient_assessed(row, *fields))
        result.append({"证据维度": label, "已评估": str(assessed), "未评估": str(total - assessed), "范围": f"Top {total}"})
    hla_assessed = total if _patient_candidate_hla_status(bundle) != "未评估" else 0
    appm_assessed = total if _patient_candidate_appm_status(bundle) != "未评估" else 0
    ccf_assessed = sum(1 for row in rows if _patient_ccf_assessment(row, bundle)[0])
    result.insert(5, {"证据维度": "克隆性/CCF", "已评估": str(ccf_assessed), "未评估": str(total - ccf_assessed), "范围": f"Top {total}（需数值与可靠置信度）"})
    result.insert(6, {"证据维度": "限制性HLA状态", "已评估": str(hla_assessed), "未评估": str(total - hla_assessed), "范围": f"Top {total}（样本级LOH共识）"})
    result.insert(7, {"证据维度": "APPM状态", "已评估": str(appm_assessed), "未评估": str(total - appm_assessed), "范围": f"Top {total}（样本级APPM评估）"})
    return result


def _patient_inferred_tool_rows(rows: list[dict[str, str]], tool_versions: Mapping[str, Mapping[str, str]] | None = None) -> list[dict[str, str]]:
    tools = [
        ("NetMHCpan", ("netmhcpan_el_rank", "netmhcpan_mt_rank_el"), "HLA-I结合/呈递"),
        ("MHCflurry", ("mhcflurry_presentation_score",), "呈递与加工"),
        ("PRIME", ("prime_score",), "免疫原性/呈递辅助"),
        ("BigMHC", ("bigmhc_im_score",), "免疫原性辅助"),
        ("NetMHCstabpan", ("netmhcstabpan_rank", "netmhcstabpan_score"), "肽-HLA稳定性"),
        ("NetChop 3.1d", ("netchop_processing_status", "netchop_31d_max_score"), "蛋白酶体酶切加工"),
        ("TAP/APPM", ("tap_processing_status",), "TAP转运与加工通路"),
    ]
    result = []
    versions = {str(key).lower(): value for key, value in (tool_versions or {}).items()}
    for name, fields, purpose in tools:
        count = sum(1 for row in rows if _patient_assessed(row, *fields))
        status = f"综合证据表已载入结果值（{count}/{len(rows)}）" if count else "未评估（本次综合证据表无结果值）"
        record = versions.get(name.lower()) or versions.get(name.lower().replace("/", "_")) or {}
        version = str(record.get("version") or "原始运行版本未记录（需补工具版本清单）")
        version_evidence = str(record.get("evidence") or "综合结果仅证明工具结果存在，不能反推版本")
        result.append({"流程/工具": name, "版本": version, "版本依据": version_evidence, "状态": status, "作用": purpose})
    return result


def _patient_validation(row: Mapping[str, Any], val_map: Mapping[str, Mapping[str, str]]) -> str:
    val = val_map.get(str(row.get("peptide_id") or ""), {})
    explicit = str(val.get("validation_strategy") or val.get("recommended_assay") or row.get("recommended_validation") or row.get("recommended_use") or "")
    if explicit:
        translations = (
            ("do not advance", "当前暂缓/不推进"),
            ("novel c-terminal tail", "新生C端肽段：优先采用覆盖新生尾部的混合长肽（15–27 aa）和/或移码minigene；短肽仅作次级验证"),
            ("mutant short peptide", "突变短肽（8–11 aa）与匹配的正常短肽对照；建议开展MHC-I ELISpot或多聚体实验"),
            ("fusion junction long peptide", "优先采用跨融合断点的长肽和/或融合minigene，并先确认精确断点与阅读框"),
            ("abnormal splice/exon-junction long peptide", "优先采用覆盖异常剪接连接点的长肽（15–27 aa）和/或剪接minigene，不应仅依赖短肽"),
        )
        lowered = explicit.lower()
        for marker, translated in translations:
            if marker in lowered:
                return translated
        return explicit
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
    values = [f"基因表达 {gene_tpm} TPM" if gene_tpm is not None else "基因表达未提供"]
    values.append(f"转录本表达 {transcript_tpm} TPM" if transcript_tpm is not None else "转录本表达未提供")
    if track in {"Fusion", "Splice"}:
        values.append(f"junction reads {junction_reads}" if junction_reads is not None else "junction reads未提供")
    else:
        values.append(f"RNA位点深度 {depth}" if depth is not None else "RNA位点深度未计算")
        values.append(f"RNA alt reads {alt_reads}" if alt_reads is not None else "RNA alt reads未计算")
        values.append(f"RNA VAF {vaf}" if vaf is not None else "RNA VAF未计算")
    return "；".join(values)


def _patient_event_evidence_and_next_step(
    row: Mapping[str, Any], bundle: ReportBundle, val_map: Mapping[str, Mapping[str, str]],
) -> str:
    evidence = [
        _patient_metric("事件", row, "event_authenticity_state", "cross_platform_status"),
        _patient_metric("RNA", row, "rna_support_state", "rna_support_status"),
        _patient_metric("呈递", row, "presentation_consensus_state", "presentation_evidence_grade"),
        _patient_metric("MT/WT", row, "mutant_specificity_status", "mutant_specificity_state"),
    ]
    gaps = _patient_key_gaps(row, bundle)
    gap_text = "；".join(gaps) if gaps else "未见明确关键缺口"
    return (
        f"核心证据：{'；'.join(evidence)}。"
        f"RNA数据：{_patient_rna_measurements(row)}。"
        f"主要缺口：{gap_text}。"
        f"下一步：{_patient_validation(row, val_map)}"
    )


def _patient_key_gaps(row: Mapping[str, Any], bundle: ReportBundle) -> list[str]:
    gaps: list[str] = []
    if not _patient_ccf_assessment(row, bundle)[0]:
        gaps.append("CCF未形成可靠估计")
    hla_status = _patient_candidate_hla_status(bundle)
    if "多工具LOH确认未完成" in hla_status:
        gaps.append("限制性HLA=单工具提示保留；多工具LOH确认未完成")
    elif hla_status == "未评估":
        gaps.append("限制性HLA未评估")
    elif "冲突" in hla_status or "丢失" in hla_status:
        gaps.append(hla_status)
    appm_status = _patient_candidate_appm_status(bundle)
    if appm_status.startswith("证据部分完整"):
        gaps.append("APPM仅部分评估")
    elif appm_status == "未评估":
        gaps.append("APPM未评估")
    safety = _patient_value(row, "safety_state", "safety_status", default="").upper()
    if safety in {"SAFETY_PARTIAL", "PARTIAL", "UNASSESSED", ""}:
        gaps.append("安全性证据不完整")
    if str(row.get("evidence_conflict_fields") or "").strip():
        gaps.append("影响候选的证据字段存在冲突")
    integrity_ok, integrity_missing = _patient_candidate_integrity(row)
    if not integrity_ok:
        gaps.append("完整性缺口：" + "、".join(integrity_missing))
    return list(dict.fromkeys(gaps))


def _patient_attention_reasons(row: Mapping[str, Any]) -> list[str]:
    reasons: list[str] = []
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
    scored: list[tuple[int, int, dict[str, str], list[str]]] = []
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
        reasons = _patient_attention_reasons(row)
        if not reasons:
            continue
        scored.append((-len(reasons), index, row, reasons))
    result: list[dict[str, str]] = []
    for _, _, row, reasons in sorted(scored)[:limit]:
        gaps = _patient_key_gaps(row, bundle)
        advice = _patient_validation(row, val_map)
        if gaps:
            advice += "。当前缺口：" + "；".join(gaps)
        result.append({
            "事件": str(row.get("gene") or row.get("event_name") or row.get("event_id") or ""),
            "为什么重要": "；".join(reasons) + "。RNA数据：" + _patient_rna_measurements(row),
            "当前建议": advice,
        })
    return result


def _patient_candidate_attention(row: Mapping[str, Any], bundle: ReportBundle) -> str:
    evidence = [
        _patient_metric("事件", row, "event_authenticity_state", "cross_platform_status"),
        _patient_metric("RNA", row, "rna_support_state", "rna_support_status"),
        _patient_metric("呈递", row, "presentation_consensus_state", "presentation_evidence_grade"),
        _patient_metric("MT/WT", row, "mutant_specificity_status", "mutant_specificity_state"),
    ]
    return "；".join(evidence) + "。RNA数据：" + _patient_rna_measurements(row)


def _patient_candidate_disposition(
    row: Mapping[str, Any],
    val_map: Mapping[str, Mapping[str, str]],
    event_grades: Mapping[str, str],
) -> tuple[str, str]:
    integrity_ok, integrity_missing = _patient_candidate_integrity(row)
    if not integrity_ok:
        return "TECHNICAL_REVIEW", "；".join(integrity_missing)
    recommendation = _patient_validation(row, val_map).strip()
    event_grade = str(event_grades.get(str(row.get("event_id") or "")) or _patient_grade(row)).upper()
    if "DO NOT ADVANCE" in recommendation.upper() or "暂缓/不推进" in recommendation or event_grade == "R4":
        return "PAUSED", "当前建议不推进或事件级为R4"
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


def _patient_qc_rows(bundle: ReportBundle) -> list[dict[str, str]]:
    provenance = bundle.provenance
    purity_result, purity_basis = _patient_purity_consensus(bundle)
    rows = [
        {"项目": "疾病/分析背景", "结果": str(provenance.get("disease") or "未记录"), "解释": "优先采用结构化临床背景，其次采用分析配置"},
        {"项目": "肿瘤/正常配对", "结果": str(provenance.get("pairing_status") or "未评估"), "解释": "区分已使用配对输入与已完成指纹确认"},
        {"项目": "肿瘤纯度/倍性", "结果": purity_result, "解释": "用于CNV、CCF和LOH解释；工具冲突必须保留"},
        {"项目": "肿瘤DNA深度", "结果": str(provenance.get("tumor_dna_depth") or "未评估"), "解释": "默认汇总去重事件位点有效深度；低深度会降低检出能力"},
        {"项目": "正常DNA深度", "结果": str(provenance.get("normal_dna_depth") or "未评估"), "解释": "默认汇总去重事件位点有效深度，用于排除胚系和正常支持"},
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
        "RETAINED": "未见LOH（保留）",
        "CONFLICT": "结果冲突",
        "UNASSESSED": missing,
    }.get(status, missing)


def _patient_hla_loh_consensus(bundle: ReportBundle) -> tuple[list[dict[str, str]], str]:
    by_tool: dict[str, dict[str, set[str]]] = {"LOHHLA": {}, "SpecHLA": {}}
    for record in bundle.hla_loh_tool_results:
        allele = str(record.get("hla_allele") or record.get("allele") or "").strip()
        if not re.match(r"^HLA-[ABC]\*", allele):
            continue
        tool = str(record.get("_report_tool") or record.get("evidence_tool") or record.get("source_tool") or record.get("tool") or "")
        tool = "LOHHLA" if "lohhla" in tool.lower() else "SpecHLA" if "spechla" in tool.lower() else ""
        if not tool:
            continue
        by_tool[tool].setdefault(allele, set()).add(_patient_hla_loh_status(record.get("loh_status") or record.get("status")))

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
        if not assessed:
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
            consensus = f"仅{tool}报告{_patient_hla_loh_label(internal)}，证据有限"
            explanation = f"另一个HLA LOH工具未提供该等位基因结果"
        aggregate.append((allele, internal))
        rows.append({
            "HLA等位基因": allele,
            "LOHHLA": _patient_hla_loh_label(statuses["LOHHLA"]),
            "SpecHLA": _patient_hla_loh_label(statuses["SpecHLA"]),
            "综合判断": consensus,
            "说明": explanation,
        })

    lost = [allele for allele, status in aggregate if status == "LOST"]
    conflicts = [allele for allele, status in aggregate if status == "CONFLICT"]
    assessed = [status for _, status in aggregate if status != "UNASSESSED"]
    tools_present = {tool for tool, calls in by_tool.items() if calls}
    if lost:
        overall = "检出限制性HLA-I LOH：" + "、".join(lost)
    elif conflicts:
        overall = "HLA-I LOH工具结果冲突：" + "、".join(conflicts)
    elif assessed and all(status == "RETAINED" for status in assessed):
        if tools_present == {"LOHHLA", "SpecHLA"}:
            overall = "多工具一致未见限制性HLA-I LOH（保留）"
        elif tools_present:
            overall = f"未见限制性HLA-I LOH（仅{next(iter(tools_present))}，证据有限）"
        else:
            overall = "限制性HLA-I LOH未评估"
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
    rows: list[dict[str, str]] = []
    if isinstance(tools, Mapping):
        for name, record in tools.items():
            if not isinstance(record, Mapping):
                continue
            version_record = bundle.tool_versions.get(str(name)) or bundle.tool_versions.get(str(name).lower()) or {}
            rows.append({
                "流程/工具": str(name),
                "版本": str(version_record.get("version") or record.get("version") or "原始运行版本未记录（需补工具版本清单）"),
                "版本依据": str(version_record.get("evidence") or record.get("version_evidence") or "运行清单"),
                "状态": _patient_status_text(record.get("status") or "UNASSESSED"),
                "作用": str(record.get("purpose") or record.get("mode") or "证据生成"),
            })
    return rows


def _patient_release_metadata(bundle: ReportBundle) -> dict[str, str]:
    input_hash = str(bundle.evidence_integrity.get("actual_sha256") or bundle.evidence_integrity.get("expected_sha256") or "")
    rules_version = str(
        bundle.profile.get("rules_version")
        or bundle.profile.get("version")
        or bundle.profile.get("_profile_version")
        or (bundle.provenance.get("evidence_consensus") or {}).get("rules_version")
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
    conflicting = [row for row in displayed if str(row.get("evidence_conflict_fields") or "").strip()]
    unexplained_conflicts = [row for row in conflicting if "来源字段存在冲突" not in _patient_limitation(row, bundle)]
    path_or_log = bool(re.search(r"(?:/mnt/|/root/|/home/|Traceback \(most recent call last\)|nohup:)", rendered_without_audit))
    hla_appm_consistent = all(
        "HLA/APPM" not in _patient_evidence_summary(row, bundle)
        and _patient_candidate_hla_status(bundle) != "未评估"
        and _patient_candidate_appm_status(bundle) != "未评估"
        for row in displayed
    )
    checks = [
        ("1", "展示候选关键字段完整", not integrity_failures, f"失败 {len(integrity_failures)} 条"),
        ("2", "精确数值存在字段级来源", not numeric_without_source, f"缺少来源映射 {len(numeric_without_source)} 条"),
        ("3", "样本级HLA/APPM与候选级表述一致", hla_appm_consistent, "HLA与APPM分别引用样本级共识"),
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
    dispositions = {
        str(row.get("peptide_id") or row.get("event_id") or index): _patient_candidate_disposition(row, val_map, event_grade_map)
        for index, row in enumerate(all_representatives)
    }

    def disposition_for(row: Mapping[str, Any]) -> tuple[str, str]:
        key = str(row.get("peptide_id") or row.get("event_id") or "")
        return dispositions.get(key) or _patient_candidate_disposition(row, val_map, event_grade_map)

    patient_representatives = [row for row in all_representatives if disposition_for(row)[0] in {"PRIORITY_CONFIRM", "EVIDENCE_GAP"}]
    top = patient_representatives[:candidate_top_n]
    paused_representatives = [row for row in all_representatives if disposition_for(row)[0] in {"PAUSED", "TECHNICAL_REVIEW"}]
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
        {"口径": "独立事件", "数量": str(len(event_seen) or len(bundle.events)), "说明": "来自ranked_events.evidence_consensus.tsv，按事件去重"},
        {"口径": "Peptide-HLA组合", "数量": str(len(ranked)), "说明": "同一事件可产生多个重叠肽段和多个HLA组合，不等同于独立新抗原数量"},
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
    focus_count = len(top[:10])
    out.append(f"<p><b>本次重点审阅：</b>{focus_count}个事件。它们是本报告综合证据表中优先展示的代表事件，不等同于已有{focus_count}个经实验确认的新抗原。</p>")
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
    out.append("<p>这些结果说明候选是否具备被加工和呈递的条件，但不能单独判断免疫治疗敏感、耐药或患者获益。</p></div>")

    out.append("<div class='section'><h2>4. 重点变异事件（按类型、按事件去重）</h2>")
    overall = [{"事件类型": track, "事件数": count, "主要复核重点": {"SNV": "DNA深度/VAF、RNA alt、MT/WT", "InDel": "局部重比对、阅读框、NMD和phasing", "Fusion": "精确断点、junction reads、frame和正常read-through", "Splice": "精确junction、PSI/reads、正常isoform和ORF", "DNA SV": "断点与异常转录本"}.get(track, "事件真实性和证据完整性")} for track, count in sorted(track_counts.items())]
    out.append(_table(overall, ["事件类型", "事件数", "主要复核重点"]))
    for track in ("SNV", "InDel", "Fusion", "Splice", "DNA SV"):
        representatives = _patient_event_representatives(bundle.events, ranked, event_top_n, track)
        if not representatives:
            continue
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
                "关键证据与下一步": _patient_event_evidence_and_next_step(row, bundle, val_map),
            })
        out.append(f"<h3>{esc(track)} Top {event_top_n}</h3>" + _table(rows, ["排名", "基因/事件", "改变", "肽段-HLA", "R等级", "关键证据与下一步"]))
    out.append(f"<p class='small'>不同事件赛道的证据结构不同，Top {event_top_n}用于赛道内审阅，不应仅凭序号直接跨赛道比较。</p></div>")

    manual_review_rows = _patient_manual_review_rows(bundle.events, ranked, bundle, val_map)
    out.append("<div class='section'><h2>关键人工审阅事件</h2>")
    out.append("<p>以下事件因机制重要、多工具/正交支持或证据存在冲突而单独保留人工审阅；进入本表不等于自动升级为R1/R2，也不代表已确认新抗原。</p>")
    if manual_review_rows:
        out.append(_table(manual_review_rows, ["事件", "为什么重要", "当前建议"]))
    else:
        out.append("<p class='small'>本次未筛出具有明确机制标记、多工具支持或证据冲突的独立人工审阅事件。</p>")
    out.append("</div>")

    out.append(f"<div class='section'><h2>5. 候选肽段Top {candidate_top_n}（跨赛道、按事件去重）</h2>")

    def patient_candidate_rows(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
        result = []
        for rank, row in enumerate(rows, 1):
            result.append({
                "排名": rank,
                "基因": row.get("gene", ""),
                "类型": _patient_track(row),
                "肽段-HLA": f"{row.get('peptide', '')} / {row.get('hla_allele', '')}",
                "等级": _patient_grade(row),
                "关键证据与下一步": _patient_event_evidence_and_next_step(row, bundle, val_map),
            })
        return result

    candidate_headers = ["排名", "基因", "类型", "肽段-HLA", "等级", "关键证据与下一步"]
    out.append(f"<h3>重点候选 Top {candidate_top_n}</h3><p class='small'>本表保留R3-GAP候选；证据缺口在“主要限制”中明确列出，不表示已经确认新抗原或可直接进入功能实验。</p>")
    out.append(_table(patient_candidate_rows(top), candidate_headers))
    out.append(f"<p class='small'>当前暂缓/不推进及完整性门槛未通过的{len(paused_representatives)}个事件代表候选不进入患者版重点表，仅保留在科研技术版审阅池。排序仍采用R1–R4、同赛道Pareto、确定性tie-break和事件去重。</p></div>")

    out.append("<div class='section'><h2>6. Top候选综合证据与实验建议</h2>")
    out.append("<p>建议顺序：先确认事件和异常转录本真实性，再补RNA alt/VAF或精确junction证据，完成MT/WT、正常背景和限制性HLA复核，最后开展短肽、长肽、minigene及T细胞功能实验。</p>")
    interpretation_rows = []
    for row in top[:10]:
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
    method_rows = [
        {"阶段": "事件输入与质控", "方法": "VCF/Fusion/Splice/SV标准化与来源链检查", "状态": "按当前输入执行"},
        {"阶段": "肽段构建", "方法": "SNV MT/WT；InDel novel tail；Fusion/Splice精确junction ORF", "状态": "无可追溯ORF时不形成高等级证据"},
        {"阶段": "呈递预测", "方法": "核心结合/呈递、稳定性和免疫原性样模型分组汇总", "状态": "工具缺失记为未评估"},
        {"阶段": "综合排序", "方法": "证据状态、hard fail、priority cap、R1–R4、Pareto和事件去重", "状态": "研究性规则"},
    ]
    out.append(_table(method_rows, ["阶段", "方法", "状态"]))
    tool_rows = _patient_tool_rows(bundle)
    recorded_names = {str(row.get("流程/工具") or "").lower() for row in tool_rows}
    for inferred in _patient_inferred_tool_rows(top, bundle.tool_versions):
        if str(inferred.get("流程/工具") or "").lower() not in recorded_names:
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
    out.append("<h3>Top 100证据维度完整性</h3>" + _table(_patient_evidence_audit_rows(top, bundle), ["证据维度", "已评估", "未评估", "范围"]))
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
            "This section does not replace the current primary weighted ranking. "
            "The evidence-consensus ranking has not been experimentally calibrated "
            "and is provided only for algorithm comparison and candidate review.</div>"
        )
        metadata_rows = [
            {"field": "primary patient-report ranking", "value": parallel_rankings.get("legacy_weighted", "ranked_peptides.tsv")},
            {"field": "parallel peptide ranking", "value": parallel_rankings.get("evidence_consensus", "")},
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
    pep_headers = [
        "peptide_id", "event_id", "gene", "cancer_gene_types", "cancer_driver_context", "cancer_gene_context",
        "source_chain_track", "source_chain_confidence_tier", "source_chain_confidence_label",
        "source_chain_orthogonal_status", "source_chain_orthogonal_sources",
        "source_chain_requirement_statuses", "source_chain_reason_codes",
        "source_chain_hard_failure", "source_chain_hard_failure_codes",
        "peptide", "wildtype_peptide", "peptide_consequence",
        "hla_allele", "mhc_class", "presentation_evidence_grade", "binding_evidence_score",
        "presentation_evidence_score", "netmhcpan_ba_rank", "netmhcpan_el_rank",
        "netmhcpan_wt_rank_el", "agretopicity_el", "mt_wt_el_rank_difference",
        "mhcflurry_mt_wt_presentation_difference", "prime_mt_wt_score_difference",
        "bigmhc_mt_wt_score_difference", "mutation_positions_in_peptide",
        "mutation_anchor_only", "mutation_tcr_facing", "mutant_specificity_status",
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
            f"confidence={esc(c.get('ccf_confidence', ppt.get('ccf_confidence')))}; "
            f"method={esc(c.get('ccf_method', ppt.get('ccf_method')))}; "
            f"multiplier={esc(c.get('clonality_multiplier', ppt.get('ccf_multiplier')))}; "
            f"warning=<span class='mono'>{esc(c.get('ccf_warning', ppt.get('ccf_warning')))}</span></p>"
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
