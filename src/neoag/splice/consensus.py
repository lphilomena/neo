"""Evidence-chain-aware consensus for splice-derived peptide origins."""
from __future__ import annotations

from collections import defaultdict

from neoag.splice.identifiers import stable_id

_RNA_GENERATORS = {"ImmunoPepper", "moPepGen"}
_CAUSAL_PRIORITY = {
    "UNASSESSED": 0, "DNA_PREDICTION_ONLY": 1, "DNA_RNA_CIS_SUPPORTED": 2,
    "TARGETED_SPANNING_ONLY": 3, "PVACSPLICE_SUPPORTED": 4,
    "TARGETED_REQUANT_SUPPORTED": 5,
}


def _split(value: str) -> set[str]:
    return {x for x in str(value or "").split(";") if x}


def _grade_num(value: str) -> int:
    try:
        return int(value[1:])
    except Exception:
        return 0


def _resolved_evidence(row: dict[str, str]) -> bool:
    status = str(row.get("resolution_status", "")).upper()
    if not status or any(token in status for token in ("UNRESOLVED", "UNVERIFIED", "REJECTED", "AMBIGUOUS")):
        return False
    return str(row.get("verified_value", "")).strip() not in {"", "0", "0.0"}


def build_consensus(tables: dict[str, list[dict[str, str]]], *, sample_id: str) -> list[dict[str, str]]:
    evidence_by_entity: dict[str, set[str]] = defaultdict(set)
    for row in tables.get("tool_evidence", []):
        if _resolved_evidence(row):
            evidence_by_entity[row.get("entity_id", "")].add(row.get("evidence_group", ""))
    links_by_event: dict[str, set[str]] = defaultdict(set)
    for row in tables.get("event_junction_links", []):
        links_by_event[row.get("splice_event_id", "")].add(row.get("junction_id", ""))
    exact_rna_junctions = {
        row.get("entity_id", "") for row in tables.get("tool_evidence", [])
        if row.get("evidence_group") == "RNA_JUNCTION"
        and row.get("resolution_status", "") == "RESOLVED_EXACT"
        and _resolved_evidence(row)
    }
    rna_sources_by_junction: dict[str, set[str]] = defaultdict(set)
    for row in tables.get("tool_evidence", []):
        if row.get("entity_id", "") in exact_rna_junctions and row.get("evidence_group") == "RNA_JUNCTION" and _resolved_evidence(row):
            rna_sources_by_junction[row.get("entity_id", "")].add(
                row.get("source_assay_id", "") or "RNA_ASSAY_UNRESOLVED"
            )
    event_groups: dict[str, set[str]] = defaultdict(set)
    for event_id, jids in links_by_event.items():
        for jid in jids:
            event_groups[event_id].update(evidence_by_entity.get(jid, set()))
        event_groups[event_id].update(evidence_by_entity.get(event_id, set()))
    for transcript in tables.get("transcripts", []):
        event_groups[transcript.get("splice_event_id", "")].update(
            evidence_by_entity.get(transcript.get("transcript_hypothesis_id", ""), set())
        )
    orfs = {row.get("orf_id", ""): row for row in tables.get("orfs", [])}
    origins = tables.get("peptide_origins", [])
    presentations_by_origin: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in tables.get("presentation", []):
        presentations_by_origin[row.get("origin_peptide_id", "")].append(row)

    normal_by_entity: dict[str, list[dict[str, str]]] = defaultdict(list)
    query_to_entities: dict[str, set[str]] = defaultdict(set)
    for q in tables.get("sequence_queries", []):
        for key in (q.get("origin_peptide_id", ""), q.get("splice_event_id", ""), q.get("junction_id", "")):
            if q.get("query_id") and key:
                query_to_entities[q["query_id"]].add(key)
    for row in tables.get("normal_background", []):
        keys = [row.get("origin_peptide_id", ""), row.get("splice_event_id", ""), row.get("junction_id", "")]
        keys.extend(query_to_entities.get(row.get("query_id", ""), set()))
        for key in keys:
            if key:
                normal_by_entity[key].append(row)

    # Exact full/partial ORF consensus requires the same protein sequence and
    # frame from independent generators, and excludes peptide-only products.
    protein_generators: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for orf in tables.get("orfs", []):
        validity = orf.get("orf_validity_status", "")
        if validity.endswith("_ONLY"):
            continue
        key = (orf.get("splice_event_id", ""), orf.get("protein_sequence_sha256", ""), orf.get("frame_status", ""))
        if key[1]:
            protein_generators[key].add(orf.get("source_generator", ""))

    def _resolved_peptide_origin(origin: dict[str, str]) -> bool:
        status = origin.get("origin_status", "").upper()
        conflict = origin.get("evidence_conflict_status", "").upper()
        if not status or any(token in status for token in ("UNRESOLVED", "AMBIGUOUS", "INVALID", "REJECTED")):
            return False
        return conflict in {"", "NONE", "RESOLVED"}

    peptide_generators: dict[tuple[str, str], set[str]] = defaultdict(set)
    for origin in origins:
        generator = origin.get("source_generator", "")
        if generator in _RNA_GENERATORS and _resolved_peptide_origin(origin):
            peptide_generators[(origin.get("splice_event_id", ""), origin.get("peptide_sequence", ""))].add(generator)

    causal_by_event: dict[str, list[dict[str, str]]] = defaultdict(list)
    causal_by_junction: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in tables.get("causal_links", []):
        causal_by_event[row.get("splice_event_id", "")].append(row)
        if row.get("junction_id"):
            causal_by_junction[row.get("junction_id", "")].append(row)

    chain_by_origin_type: dict[tuple[str, str], dict[str, str]] = {}
    for row in tables.get("evidence_chains", []):
        chain_by_origin_type[(row.get("origin_peptide_id", ""), row.get("chain_type", ""))] = row

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
        independent_rna_sources = {
            source for jid in event_jids for source in rna_sources_by_junction.get(jid, set()) if source
        }
        causal_rows = list(causal_by_event.get(event_id, []))
        for jid in event_jids:
            causal_rows.extend(causal_by_junction.get(jid, []))
        seen_causal: set[str] = set()
        causal_rows = [r for r in causal_rows if not (r.get("causal_link_id", "") in seen_causal or seen_causal.add(r.get("causal_link_id", "")))]
        causal_status = max(
            (row.get("causal_status", "UNASSESSED") for row in causal_rows),
            key=lambda x: _CAUSAL_PRIORITY.get(x, 0), default="UNASSESSED",
        )
        has_causal = _CAUSAL_PRIORITY.get(causal_status, 0) >= 2
        has_event_model = bool(groups & {"SPLICE_GRAPH", "INTRON_RETENTION", "LONG_READ"})
        if "LONG_READ" in groups or has_causal:
            event_grade = "E3"
        elif has_exact_rna and has_event_model:
            event_grade = "E2"
        elif has_exact_rna:
            event_grade = "E1"
        else:
            event_grade = "E0"

        generator_key = (event_id, orf.get("protein_sequence_sha256", ""), orf.get("frame_status", ""))
        exact_orf_generators = {x for x in protein_generators.get(generator_key, set()) if x}
        exact_peptide_generators = peptide_generators.get((event_id, origin.get("peptide_sequence", "")), set())
        if groups & {"PROTEIN_VALIDATION", "LIGANDOME", "LONG_READ"}:
            orf_grade = "O3"
            translation_level = "DIRECT_OR_LONGREAD_VALIDATION"
        elif len(exact_orf_generators) >= 2:
            orf_grade = "O2"
            translation_level = "EXACT_ORF"
        elif orf and orf.get("orf_validity_status") not in {"", "INVALID", "UNRESOLVED"}:
            orf_grade = "O1"
            translation_level = "EXACT_PEPTIDE" if len(exact_peptide_generators) >= 2 else "SINGLE_GENERATOR"
        else:
            orf_grade = "O0"
            translation_level = "UNRESOLVED"

        normal_rows: list[dict[str, str]] = []
        for key in [por, event_id, *event_jids]:
            normal_rows.extend(normal_by_entity.get(key, []))
        detected_statuses = {
            "DETECTED_MATCHED_NORMAL", "DETECTED_CRITICAL_TISSUE",
            "DETECTED_BROAD_NORMAL", "LOW_LEVEL_NONCRITICAL_NORMAL",
            "NORMAL_DETECTED",
        }
        seen: set[str] = set()
        normal_rows = [
            row for row in normal_rows
            if not (
                row.get("normal_background_id", "") in seen
                or seen.add(row.get("normal_background_id", ""))
            )
        ]
        detected_critical = any(
            row.get("assessment_status") == "DETECTED_CRITICAL_TISSUE"
            or (
                row.get("assessment_status") == "NORMAL_DETECTED"
                and row.get("critical_tissue", "").lower() == "true"
            )
            for row in normal_rows
        )
        detected_any = any(row.get("assessment_status") in detected_statuses for row in normal_rows)
        adequate_negative_sources = {
            (row.get("normal_source_type", ""), row.get("normal_source", ""))
            for row in normal_rows if row.get("assessment_status") == "NOT_DETECTED_ADEQUATE_COVERAGE"
        }
        kmer_negative_sources = {
            (row.get("normal_source_type", ""), row.get("normal_source", ""))
            for row in normal_rows if row.get("assessment_status") == "NOT_DETECTED_KMER_SCREEN"
        }
        if detected_critical or detected_any:
            normal_grade = "N0"
            normal_status = "NORMAL_DETECTED_CRITICAL" if detected_critical else "NORMAL_DETECTED"
        elif adequate_negative_sources and kmer_negative_sources:
            normal_grade = "N3"
            normal_status = "LOCUS_AND_KMER_NOT_DETECTED"
        elif len(adequate_negative_sources) >= 2:
            normal_grade = "N3"
            normal_status = "MULTI_SOURCE_NOT_DETECTED_ADEQUATE_COVERAGE"
        elif len(adequate_negative_sources) == 1:
            normal_grade = "N2"
            normal_status = "NOT_DETECTED_ADEQUATE_COVERAGE"
        elif kmer_negative_sources:
            normal_grade = "N1"
            normal_status = "K4NEO_NEGATIVE_ONLY"
        else:
            normal_grade = "N1"
            normal_status = "NORMAL_BACKGROUND_INCOMPLETE"

        presentations = presentations_by_origin.get(por, [])
        presentation_grade = "P1" if presentations else "P0"
        rna_chain = chain_by_origin_type.get((por, "RNA_DRIVEN"), {})
        dna_chain = chain_by_origin_type.get((por, "DNA_CAUSAL"), {})
        normal_chain = chain_by_origin_type.get((por, "NORMAL_BACKGROUND"), {})

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
        if len(exact_orf_generators) < 2:
            if len(exact_peptide_generators) >= 2:
                caps.append("CAP_DUAL_GENERATOR_PEPTIDE_ONLY_R2")
            else:
                caps.append("CAP_SINGLE_PEPTIDE_GENERATOR_R2")
        if normal_grade == "N1":
            caps.append("CAP_K4NEO_ONLY_R3" if kmer_negative_sources else "CAP_NORMAL_BACKGROUND_INCOMPLETE_R3")
        if event_grade == "E0":
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
        r3_caps = {
            "CAP_NORMAL_BACKGROUND_INCOMPLETE_R3", "CAP_K4NEO_ONLY_R3",
            "CAP_SPLICE_EVENT_RNA_UNCONFIRMED_R3", "CAP_PRESENTATION_UNASSESSED_R3",
            "CAP_PEPTIDE_NOVELTY_UNRESOLVED_R3",
        }
        if any(code in caps for code in r3_caps) and tier in {"R1", "R2"}:
            tier = "R3"
        elif any(code in caps for code in {"CAP_SINGLE_PEPTIDE_GENERATOR_R2", "CAP_DUAL_GENERATOR_PEPTIDE_ONLY_R2"}) and tier == "R1":
            tier = "R2"
        cap = max((code.rsplit("R", 1)[-1] for code in caps if "_R" in code), default="", key=lambda x: int(x))
        all_generators = exact_orf_generators | exact_peptide_generators
        rows.append({
            "consensus_id": stable_id("CON", event_id, por), "splice_event_id": event_id,
            "origin_peptide_id": por, "peptide_id": origin.get("peptide_id", ""), "sample_id": sample_id,
            "event_evidence_grade": event_grade, "orf_evidence_grade": orf_grade,
            "normal_safety_grade": normal_grade, "presentation_grade": presentation_grade,
            "rna_driven_chain_status": rna_chain.get("chain_status", "RNA_DRIVEN_UNASSESSED"),
            "dna_causal_chain_status": dna_chain.get("chain_status", causal_status),
            "normal_background_chain_status": normal_chain.get("chain_status", normal_status),
            "translation_consensus_level": translation_level,
            "independent_evidence_groups": ";".join(sorted(x for x in groups if x)),
            "independent_rna_sources": ";".join(sorted(independent_rna_sources)),
            "independent_translation_generators": ";".join(sorted(exact_orf_generators)),
            "independent_peptide_generators": ";".join(sorted(exact_peptide_generators)),
            "event_consensus_status": "DNA_CAUSAL_AND_RNA_SUPPORTED" if has_causal and has_exact_rna else ("RNA_EVENT_SUPPORTED" if has_exact_rna else "RNA_EVENT_UNCONFIRMED"),
            "orf_consensus_status": "MULTI_GENERATOR_EXACT_ORF" if len(exact_orf_generators) >= 2 else ("MULTI_GENERATOR_EXACT_PEPTIDE" if len(exact_peptide_generators) >= 2 else "SINGLE_GENERATOR"),
            "normal_background_status": normal_status, "final_evidence_tier": tier,
            "priority_cap": f"R{cap}" if cap else "", "consensus_reason": (
                f"event={event_grade}; orf={orf_grade}; translation={translation_level}; normal={normal_grade}; "
                f"presentation={presentation_grade}; causal={causal_status}; generators={','.join(sorted(all_generators)) or 'NONE'}"
            ),
            "hard_fail_codes": ";".join(hard), "cap_codes": ";".join(caps),
        })
    return rows


def consensus_reason_conflicts(consensus_rows: list[dict[str, str]], *, sample_id: str) -> list[dict[str, str]]:
    conflicts: list[dict[str, str]] = []
    for row in consensus_rows:
        entity_id = row.get("origin_peptide_id", "") or row.get("splice_event_id", "")
        for field_name, conflict_type, severity in (
            ("hard_fail_codes", "CONSENSUS_HARD_FAIL", "ERROR"),
            ("cap_codes", "CONSENSUS_PRIORITY_CAP", "WARNING"),
        ):
            for code in _split(row.get(field_name, "")):
                conflict = {
                    "entity_type": "PEPTIDE_ORIGIN", "entity_id": entity_id, "sample_id": sample_id,
                    "conflict_type": conflict_type, "field_name": field_name,
                    "observed_values": code, "source_tools": "NeoAgEvidenceConsensus",
                    "source_record_ids": row.get("consensus_id", ""), "severity": severity,
                    "resolution_status": "RULE_APPLIED",
                    "resolution_reason": f"{code} contributed to final tier {row.get('final_evidence_tier', '')}.",
                }
                conflict["conflict_id"] = stable_id("CFL", conflict)
                conflicts.append(conflict)
    return conflicts
