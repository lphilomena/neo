"""Evidence-group-aware consensus for splice-derived peptide origins."""
from __future__ import annotations

from collections import defaultdict
from typing import Iterable

from neoag.splice.identifiers import stable_id


def _split(value: str) -> set[str]:
    return {x for x in str(value or "").split(";") if x}


def _grade_num(value: str) -> int:
    try:
        return int(value[1:])
    except Exception:
        return 0


def build_consensus(tables: dict[str, list[dict[str, str]]], *, sample_id: str) -> list[dict[str, str]]:
    evidence_by_entity: dict[str, set[str]] = defaultdict(set)
    for row in tables.get("tool_evidence", []):
        evidence_by_entity[row.get("entity_id", "")].add(row.get("evidence_group", ""))
    links_by_event: dict[str, set[str]] = defaultdict(set)
    for row in tables.get("event_junction_links", []):
        links_by_event[row.get("splice_event_id", "")].add(row.get("junction_id", ""))
    exact_rna_junctions = {
        row.get("entity_id", "") for row in tables.get("tool_evidence", [])
        if row.get("evidence_group") == "RNA_JUNCTION"
        and row.get("resolution_status", "").startswith("RESOLVED")
        and str(row.get("verified_value", "")).strip() not in {"", "0", "0.0"}
    }
    event_groups: dict[str, set[str]] = defaultdict(set)
    for event_id, jids in links_by_event.items():
        for jid in jids:
            event_groups[event_id].update(evidence_by_entity.get(jid, set()))
        event_groups[event_id].update(evidence_by_entity.get(event_id, set()))
    orfs = {row.get("orf_id", ""): row for row in tables.get("orfs", [])}
    origins = tables.get("peptide_origins", [])
    presentations_by_origin: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in tables.get("presentation", []):
        presentations_by_origin[row.get("origin_peptide_id", "")].append(row)
    normal_by_entity: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in tables.get("normal_background", []):
        for key in (row.get("origin_peptide_id", ""), row.get("splice_event_id", ""), row.get("junction_id", "")):
            if key:
                normal_by_entity[key].append(row)

    protein_generators: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for orf in tables.get("orfs", []):
        key = (orf.get("splice_event_id", ""), orf.get("protein_sequence_sha256", ""), orf.get("frame_status", ""))
        if key[1]:
            protein_generators[key].add(orf.get("source_generator", ""))

    rows: list[dict[str, str]] = []
    for origin in origins:
        event_id = origin.get("splice_event_id", "")
        por = origin.get("origin_peptide_id", "")
        oid = origin.get("orf_id", "")
        orf = orfs.get(oid, {})
        groups = set(event_groups.get(event_id, set()))
        groups.update(evidence_by_entity.get(oid, set()))
        groups.update(evidence_by_entity.get(por, set()))
        event_jids = links_by_event.get(event_id, set()) | _split(origin.get("junction_ids", ""))
        has_exact_rna = bool(event_jids & exact_rna_junctions)
        has_event_model = bool(groups & {"SPLICE_GRAPH", "INTRON_RETENTION", "DNA_CAUSAL", "LONG_READ"})
        if "LONG_READ" in groups or "DNA_CAUSAL" in groups:
            event_grade = "E3"
        elif has_exact_rna and has_event_model:
            event_grade = "E2"
        elif has_exact_rna:
            event_grade = "E1"
        else:
            event_grade = "E0"

        generator_key = (event_id, orf.get("protein_sequence_sha256", ""), orf.get("frame_status", ""))
        generators = {x for x in protein_generators.get(generator_key, set()) if x}
        if "PROTEIN_VALIDATION" in groups or "LONG_READ" in groups:
            orf_grade = "O3"
        elif len(generators) >= 2:
            orf_grade = "O2"
        elif orf and orf.get("orf_validity_status") not in {"", "INVALID", "UNRESOLVED"}:
            orf_grade = "O1"
        else:
            orf_grade = "O0"

        normal_rows: list[dict[str, str]] = []
        for key in [por, event_id, *event_jids]:
            normal_rows.extend(normal_by_entity.get(key, []))
        detected_critical = any(
            row.get("assessment_status") == "NORMAL_DETECTED" and row.get("critical_tissue", "").lower() == "true"
            for row in normal_rows
        )
        detected_any = any(row.get("assessment_status") == "NORMAL_DETECTED" for row in normal_rows)
        adequate_negative_sources = {
            (row.get("normal_source_type", ""), row.get("normal_source", ""))
            for row in normal_rows if row.get("assessment_status") == "NOT_DETECTED_ADEQUATE_COVERAGE"
        }
        if detected_critical or detected_any:
            normal_grade = "N0"
            normal_status = "NORMAL_DETECTED"
        elif len(adequate_negative_sources) >= 2:
            normal_grade = "N3"
            normal_status = "MULTI_SOURCE_NOT_DETECTED_ADEQUATE_COVERAGE"
        elif len(adequate_negative_sources) == 1:
            normal_grade = "N2"
            normal_status = "NOT_DETECTED_ADEQUATE_COVERAGE"
        else:
            normal_grade = "N1"
            normal_status = "NORMAL_BACKGROUND_INCOMPLETE"

        presentations = presentations_by_origin.get(por, [])
        presentation_grade = "P1" if presentations else "P0"
        hard: list[str] = []
        caps: list[str] = []
        if not orf or orf_grade == "O0":
            hard.append("HARD_ORF_INVALID")
        crosses_status = str(origin.get("crosses_junction", "")).lower()
        novel_status = str(origin.get("contains_novel_aa", "")).lower()
        if crosses_status == "false" and novel_status == "false":
            hard.append("HARD_NO_NOVEL_AMINO_ACID")
        elif crosses_status != "true" and novel_status != "true":
            caps.append("CAP_PEPTIDE_NOVELTY_UNRESOLVED_R3")
        if normal_grade == "N0":
            hard.append("HARD_NORMAL_BACKGROUND_DETECTED")
        if not event_id:
            hard.append("HARD_PEPTIDE_ORIGIN_UNRESOLVED")
        if len(generators) < 2:
            caps.append("CAP_SINGLE_PEPTIDE_GENERATOR_R2")
        if normal_grade == "N1":
            caps.append("CAP_NORMAL_BACKGROUND_INCOMPLETE_R3")
        if event_grade in {"E0"}:
            caps.append("CAP_SPLICE_EVENT_RNA_UNCONFIRMED_R3")
        if presentation_grade == "P0":
            caps.append("CAP_PRESENTATION_UNASSESSED_R3")

        if hard:
            tier = "R4"
        elif _grade_num(event_grade) >= 2 and _grade_num(orf_grade) >= 2 and _grade_num(normal_grade) >= 2 and presentation_grade == "P1":
            tier = "R1"
        elif _grade_num(event_grade) >= 1 and _grade_num(orf_grade) >= 1 and presentation_grade == "P1":
            tier = "R2"
        else:
            tier = "R3"
        if any(code in caps for code in {
            "CAP_NORMAL_BACKGROUND_INCOMPLETE_R3", "CAP_SPLICE_EVENT_RNA_UNCONFIRMED_R3",
            "CAP_PRESENTATION_UNASSESSED_R3", "CAP_PEPTIDE_NOVELTY_UNRESOLVED_R3",
        }):
            if tier in {"R1", "R2"}:
                tier = "R3"
        elif "CAP_SINGLE_PEPTIDE_GENERATOR_R2" in caps and tier == "R1":
            tier = "R2"
        cap = max((code.rsplit("R", 1)[-1] for code in caps if "_R" in code), default="", key=lambda x: int(x))
        rows.append({
            "consensus_id": stable_id("CON", event_id, por), "splice_event_id": event_id,
            "origin_peptide_id": por, "peptide_id": origin.get("peptide_id", ""), "sample_id": sample_id,
            "event_evidence_grade": event_grade, "orf_evidence_grade": orf_grade,
            "normal_safety_grade": normal_grade, "presentation_grade": presentation_grade,
            "independent_evidence_groups": ";".join(sorted(x for x in groups if x)),
            "independent_translation_generators": ";".join(sorted(generators)),
            "event_consensus_status": "RNA_EVENT_SUPPORTED" if has_exact_rna else "RNA_EVENT_UNCONFIRMED",
            "orf_consensus_status": "MULTI_GENERATOR_CONSENSUS" if len(generators) >= 2 else "SINGLE_GENERATOR",
            "normal_background_status": normal_status, "final_evidence_tier": tier,
            "priority_cap": f"R{cap}" if cap else "", "consensus_reason": (
                f"event={event_grade}; orf={orf_grade}; normal={normal_grade}; presentation={presentation_grade}; "
                f"groups={','.join(sorted(groups)) or 'NONE'}"
            ),
            "hard_fail_codes": ";".join(hard), "cap_codes": ";".join(caps),
        })
    return rows
