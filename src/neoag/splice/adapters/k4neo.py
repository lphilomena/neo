"""k4neo normal-background plugin with exact cts_id provenance.

k4neo output is imported only through a project-generated query map. Database
non-detection is kept distinct from locus-coverage-aware matched-normal
non-detection and cannot independently establish the strongest normal-safety
grade.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

from neoag.splice.identifiers import stable_id

from .base import as_float_text, as_int, get, read_delimited, row_hash, source_record_id


def _query_map(path: str | Path) -> dict[str, dict[str, str]]:
    mapping: dict[str, dict[str, str]] = {}
    for row in read_delimited(path):
        qid = get(row, "query_id", "cts_id", "query_name")
        if qid:
            mapping[qid] = row
    return mapping


def _normal_row(
    *,
    query: dict[str, str], sample_id: str, source_name: str, source_type: str,
    tissue: str, stage: str, study_id: str, count: int, total: int,
    prevalence: float | None, critical_tissues: set[str], source_file: str,
    source_record_id_value: str, uniqueness_rate: str = "",
) -> dict[str, str]:
    detected = (count > 0) or (prevalence is not None and prevalence > 0)
    if detected:
        assessment = "NORMAL_DETECTED"
        reason = "k4neo detected the exact query sequence in healthy samples."
        detection = "DETECTED"
    elif total > 0 or prevalence == 0:
        assessment = "NOT_DETECTED_KMER_SCREEN"
        reason = "No k4neo healthy-sample hit was reported for the exact query; this is a sequence-index screen, not locus coverage."
        detection = "NOT_DETECTED"
    else:
        assessment = "UNASSESSED"
        reason = "k4neo output did not provide a resolvable healthy-sample denominator or sample rate."
        detection = "UNASSESSED"
    critical = tissue.casefold() in critical_tissues if tissue else False
    return {
        "normal_background_id": stable_id("NBG", source_name, query.get("query_id"), tissue, stage, study_id, source_type),
        "splice_event_id": query.get("splice_event_id", ""), "junction_id": query.get("junction_id", ""),
        "origin_peptide_id": query.get("origin_peptide_id", ""), "query_id": query.get("query_id", ""),
        "sample_id": sample_id, "normal_source": source_name, "normal_source_type": source_type,
        "normal_tissue": tissue, "critical_tissue": "true" if critical else "false",
        "developmental_stage": stage, "study_id": study_id, "detection_status": detection,
        "coverage_status": "SEQUENCE_INDEX_QUERIED" if assessment != "UNASSESSED" else "UNASSESSED",
        "junction_reads": "", "sample_count": str(count) if total or count else "",
        "total_samples": str(total) if total else "",
        "sample_prevalence": "" if prevalence is None else f"{prevalence:.12g}",
        "kmer_prevalence": "" if prevalence is None else f"{prevalence:.12g}",
        "uniqueness_rate": uniqueness_rate, "assessment_status": assessment,
        "assessment_reason": reason, "source_file": source_file,
        "source_record_id": source_record_id_value, "evidence_conflict_status": "NONE",
    }


def parse_k4neo(
    *,
    sample_id: str,
    query_map: str | Path,
    healthy_sample_rate: Iterable[str | Path] = (),
    annotated: Iterable[str | Path] = (),
    uniqueness: Iterable[str | Path] = (),
    source_tool_version: str = "UNASSESSED",
    critical_tissues: Iterable[str] = (),
    strict: bool = False,
) -> dict[str, list[dict[str, str]]]:
    mapping = _query_map(query_map)
    critical = {x.casefold() for x in critical_tissues if x}
    result: dict[str, list[dict[str, str]]] = {"normal_background": [], "tool_evidence": [], "conflicts": []}

    def unresolved(path: Path, row_no: int, row: dict[str, str], qid: str) -> None:
        rid = source_record_id("k4neo", path, row_no, row)
        result["conflicts"].append({
            "conflict_id": stable_id("CFL", "k4neo", rid, qid),
            "entity_type": "K4NEO_RESULT", "entity_id": rid, "sample_id": sample_id,
            "conflict_type": "K4NEO_QUERY_ID_UNRESOLVED", "field_name": "cts_id",
            "observed_values": qid, "source_tools": "k4neo", "source_record_ids": rid,
            "severity": "ERROR" if strict else "WARNING", "resolution_status": "UNRESOLVED",
            "resolution_reason": "cts_id must match a query_id in the exact project k4neo query map.",
        })

    for value in healthy_sample_rate:
        p = Path(value)
        for row_no, row in enumerate(read_delimited(p), start=2):
            qid = get(row, "cts_id", "query_id")
            query = mapping.get(qid)
            if not query:
                unresolved(p, row_no, row, qid)
                continue
            rid = source_record_id("k4neo", p, row_no, row)
            rate_text = as_float_text(get(row, "sample_rate", "sample_prevalence"))
            rate = float(rate_text) if rate_text else None
            nrow = _normal_row(
                query=query, sample_id=sample_id, source_name="k4neo",
                source_type="K4NEO_HEALTHY_SAMPLE_RATE", tissue=get(row, "tissue"),
                stage=get(row, "developmental_stage"), study_id=get(row, "study_id"),
                count=0, total=0, prevalence=rate, critical_tissues=critical,
                source_file=str(p), source_record_id_value=rid,
            )
            result["normal_background"].append(nrow)
            result["tool_evidence"].append({
                "evidence_id": stable_id("EVD", "k4neo", qid, rid), "entity_type": "SEQUENCE_QUERY",
                "entity_id": qid, "sample_id": sample_id, "evidence_group": "NORMAL_BACKGROUND_KMER",
                "evidence_type": nrow["assessment_status"], "source_tool": "k4neo",
                "source_tool_version": source_tool_version, "source_file": str(p),
                "source_row_number": str(row_no), "source_record_id": rid,
                "provided_value": rate_text, "verified_value": nrow["assessment_status"],
                "resolution_status": "MAPPED_EXACT_QUERY_ID", "resolution_reason": nrow["assessment_reason"],
                "raw_payload_sha256": row_hash(row),
            })

    for value in annotated:
        p = Path(value)
        for row_no, row in enumerate(read_delimited(p), start=2):
            qid = get(row, "cts_id", "query_id")
            query = mapping.get(qid)
            if not query:
                unresolved(p, row_no, row, qid)
                continue
            disease = get(row, "disease").casefold()
            if disease and disease not in {"healthy", "normal", "non-tumor", "nontumor"}:
                # Tumour hits are informative for prevalence but are not normal-safety evidence.
                continue
            rid = source_record_id("k4neo", p, row_no, row)
            count, total = as_int(get(row, "count"), 0), as_int(get(row, "total"), 0)
            prevalence = (count / total) if total > 0 else None
            nrow = _normal_row(
                query=query, sample_id=sample_id, source_name="k4neo",
                source_type="K4NEO_ANNOTATED_HEALTHY_STUDY", tissue=get(row, "tissue"),
                stage=get(row, "developmental_stage"), study_id=get(row, "study_id"),
                count=count, total=total, prevalence=prevalence, critical_tissues=critical,
                source_file=str(p), source_record_id_value=rid,
            )
            result["normal_background"].append(nrow)
            result["tool_evidence"].append({
                "evidence_id": stable_id("EVD", "k4neo", qid, rid), "entity_type": "SEQUENCE_QUERY",
                "entity_id": qid, "sample_id": sample_id, "evidence_group": "NORMAL_BACKGROUND_KMER",
                "evidence_type": nrow["assessment_status"], "source_tool": "k4neo",
                "source_tool_version": source_tool_version, "source_file": str(p),
                "source_row_number": str(row_no), "source_record_id": rid,
                "provided_value": f"{count}/{total}", "verified_value": nrow["assessment_status"],
                "resolution_status": "MAPPED_EXACT_QUERY_ID", "resolution_reason": nrow["assessment_reason"],
                "raw_payload_sha256": row_hash(row),
            })

    for value in uniqueness:
        p = Path(value)
        for row_no, row in enumerate(read_delimited(p), start=2):
            qid = get(row, "cts_id", "query_id", "query_cts_id")
            query = mapping.get(qid)
            if not query:
                unresolved(p, row_no, row, qid)
                continue
            rid = source_record_id("k4neo", p, row_no, row)
            unique = as_float_text(get(row, "uniqueness_rate", "unique_rate", "uniq_rate", "fraction_unique"))
            result["normal_background"].append({
                "normal_background_id": stable_id("NBG", "k4neo-uniq", qid, rid),
                "splice_event_id": query.get("splice_event_id", ""), "junction_id": query.get("junction_id", ""),
                "origin_peptide_id": query.get("origin_peptide_id", ""), "query_id": qid,
                "sample_id": sample_id, "normal_source": "k4neo-uniq",
                "normal_source_type": "K4NEO_SEQUENCE_UNIQUENESS", "normal_tissue": "",
                "critical_tissue": "false", "developmental_stage": "", "study_id": "",
                "detection_status": "NOT_APPLICABLE", "coverage_status": "SEQUENCE_UNIQUENESS_ASSESSED",
                "junction_reads": "", "sample_count": "", "total_samples": "",
                "sample_prevalence": "", "kmer_prevalence": "", "uniqueness_rate": unique,
                "assessment_status": "UNIQUENESS_ASSESSED", "assessment_reason": "k4neo uniqueness result for exact query_id.",
                "source_file": str(p), "source_record_id": rid, "evidence_conflict_status": "NONE",
            })
    return result
