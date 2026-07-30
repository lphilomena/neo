"""EasyQuant exact-query adapter.

Only query names present in the project-generated query map are accepted.  The
adapter never resolves an EasyQuant row by gene name, sequence similarity or
row order.
"""
from __future__ import annotations

from pathlib import Path

from neoag.splice.identifiers import stable_id

from .base import as_float_text, as_int, get, read_delimited, row_hash, source_record_id


def parse_easyquant(
    path: str | Path,
    *,
    sample_id: str,
    query_map: str | Path,
    source_tool_version: str = "UNASSESSED",
    strict: bool = False,
) -> dict[str, list[dict[str, str]]]:
    p = Path(path)
    map_path = Path(query_map)
    mapping = {row.get("query_id") or row.get("query_name"): row for row in read_delimited(map_path)}
    result: dict[str, list[dict[str, str]]] = {
        "targeted_quantification": [], "causal_links": [], "tool_evidence": [], "conflicts": [],
    }
    for row_no, row in enumerate(read_delimited(p), start=2):
        rid = source_record_id("EasyQuant", p, row_no, row)
        name = get(row, "name", "query_id", "sequence_name")
        query = mapping.get(name)
        if not query:
            result["conflicts"].append({
                "conflict_id": stable_id("CFL", "EasyQuant", rid, name),
                "entity_type": "TARGETED_QUANTIFICATION", "entity_id": rid, "sample_id": sample_id,
                "conflict_type": "EASYQUANT_QUERY_ID_UNRESOLVED", "field_name": "name",
                "observed_values": name, "source_tools": "EasyQuant", "source_record_ids": rid,
                "severity": "ERROR" if strict else "WARNING", "resolution_status": "UNRESOLVED",
                "resolution_reason": "The EasyQuant name must match a query_id in the exact project query map.",
            })
            continue
        junction_reads = as_int(get(row, "junc", "junction_reads"), 0)
        spanning = as_int(get(row, "span", "spanning_pairs"), 0)
        if junction_reads > 0:
            support = "TARGETED_REQUANT_SUPPORTED"
        elif spanning > 0:
            support = "TARGETED_SPANNING_ONLY"
        else:
            support = "TARGETED_REQUANT_NEGATIVE"
        qid = query.get("query_id", name)
        tqid = stable_id("TQ", qid, rid)
        result["targeted_quantification"].append({
            "targeted_quant_id": tqid, "query_id": qid, "sample_id": sample_id,
            "query_name": name, "splice_event_id": query.get("splice_event_id", ""),
            "junction_id": query.get("junction_id", ""), "variant_id": query.get("variant_id", ""),
            "position_1based": get(row, "pos", "position", default=query.get("position_1based", "")),
            "junction_reads": str(junction_reads), "spanning_pairs": str(spanning),
            "max_anchor": str(as_int(get(row, "anch", "max_anchor"), 0)),
            "left_reads": str(as_int(get(row, "a", "left_reads"), 0)),
            "right_reads": str(as_int(get(row, "b", "right_reads"), 0)),
            "interval": get(row, "interval"), "within_interval": str(as_int(get(row, "within_interval"), 0)),
            "coverage_percent": as_float_text(get(row, "coverage_perc", "coverage_percent")),
            "coverage_mean": as_float_text(get(row, "coverage_mean")),
            "coverage_median": as_float_text(get(row, "coverage_median")),
            "support_status": support, "source_tool": "EasyQuant",
            "source_tool_version": source_tool_version, "source_file": str(p),
            "source_record_id": rid, "mapping_status": "MAPPED_EXACT_QUERY_ID",
            "evidence_conflict_status": "NONE",
        })
        entity_id = query.get("junction_id") or query.get("splice_event_id") or qid
        result["tool_evidence"].append({
            "evidence_id": stable_id("EVD", "EasyQuant", qid, rid),
            "entity_type": "JUNCTION" if query.get("junction_id") else "SEQUENCE_QUERY",
            "entity_id": entity_id, "sample_id": sample_id, "evidence_group": "TARGETED_REQUANT",
            "evidence_type": support, "source_tool": "EasyQuant", "source_tool_version": source_tool_version,
            "source_file": str(p), "source_row_number": str(row_no), "source_record_id": rid,
            "provided_value": f"junc={junction_reads};span={spanning}",
            "verified_value": str(junction_reads), "resolution_status": "MAPPED_EXACT_QUERY_ID",
            "resolution_reason": f"query_id={qid}", "raw_payload_sha256": row_hash(row),
        })
        variant_id = query.get("variant_id", "")
        junction_id = query.get("junction_id", "")
        event_id = query.get("splice_event_id", "")
        if variant_id and junction_id and event_id:
            link = stable_id("DCL", variant_id, junction_id, event_id)
            result["causal_links"].append({
                "causal_link_id": link, "variant_id": variant_id, "junction_id": junction_id,
                "splice_event_id": event_id, "sample_id": sample_id, "gene": "", "gene_id": "",
                "transcript_id": query.get("transcript_hypothesis_id", ""),
                "causal_status": support if junction_reads > 0 else "DNA_RNA_CIS_SUPPORTED",
                "prediction_status": "UNASSESSED", "rna_junction_status": "EXACT_RNA_SUPPORTED" if junction_reads > 0 else "UNASSESSED",
                "targeted_requant_status": support, "pvacsplice_status": "UNASSESSED",
                "junction_reads": "", "easyquant_junction_reads": str(junction_reads),
                "easyquant_spanning_pairs": str(spanning), "spliceai_score": "", "pangolin_score": "",
                "mmsplice_score": "", "ci_spliceai_score": "", "source_tools": "EasyQuant",
                "source_tool_versions": source_tool_version, "source_files": str(p), "source_record_ids": rid,
                "link_resolution_status": "RESOLVED_EXACT_QUERY_TO_VARIANT_AND_JUNCTION",
                "resolution_reason": f"EasyQuant query {qid} was generated from the exact causal context.",
                "evidence_conflict_status": "NONE",
            })
    return result
