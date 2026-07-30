"""Independent evidence-chain construction for NeoAg v0.5.1."""
from __future__ import annotations

from collections import defaultdict
from typing import Iterable

from neoag.splice.identifiers import stable_id

_RNA_GENERATORS = {"ImmunoPepper", "moPepGen"}
_CAUSAL_PRIORITY = {
    "UNASSESSED": 0,
    "DNA_PREDICTION_ONLY": 1,
    "DNA_RNA_CIS_SUPPORTED": 2,
    "PVACSPLICE_SUPPORTED": 3,
    "TARGETED_SPANNING_ONLY": 3,
    "TARGETED_REQUANT_SUPPORTED": 4,
}


def _tokens(value: str) -> set[str]:
    return {x for x in str(value or "").split(";") if x}


def _chain_row(
    *,
    event_id: str,
    origin: dict[str, str],
    sample_id: str,
    chain_type: str,
    status: str,
    strength: str,
    groups: Iterable[str] = (),
    tools: Iterable[str] = (),
    entities: Iterable[str] = (),
    evidence_ids: Iterable[str] = (),
    limits: Iterable[str] = (),
    conflict: str = "NONE",
    reason: str = "",
) -> dict[str, str]:
    origin_id = origin.get("origin_peptide_id", "")
    return {
        "evidence_chain_id": stable_id("ECH", event_id, origin_id, chain_type),
        "splice_event_id": event_id, "origin_peptide_id": origin_id,
        "peptide_id": origin.get("peptide_id", ""), "sample_id": sample_id,
        "chain_type": chain_type, "chain_status": status, "chain_strength": strength,
        "independent_source_groups": ";".join(sorted(set(x for x in groups if x))),
        "source_tools": ";".join(sorted(set(x for x in tools if x))),
        "supporting_entity_ids": ";".join(sorted(set(x for x in entities if x))),
        "supporting_evidence_ids": ";".join(sorted(set(x for x in evidence_ids if x))),
        "limiting_reasons": ";".join(sorted(set(x for x in limits if x))),
        "conflict_status": conflict, "chain_reason": reason,
    }


def build_evidence_chains(tables: dict[str, list[dict[str, str]]], *, sample_id: str) -> list[dict[str, str]]:
    origins = tables.get("peptide_origins", [])
    event_links: dict[str, set[str]] = defaultdict(set)
    for row in tables.get("event_junction_links", []):
        event_links[row.get("splice_event_id", "")].add(row.get("junction_id", ""))

    evidence_by_entity: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in tables.get("tool_evidence", []):
        evidence_by_entity[row.get("entity_id", "")].append(row)

    origin_by_event_peptide: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for origin in origins:
        origin_by_event_peptide[(origin.get("splice_event_id", ""), origin.get("peptide_sequence", ""))].append(origin)

    orfs = {row.get("orf_id", ""): row for row in tables.get("orfs", [])}
    orf_generators: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for orf in tables.get("orfs", []):
        key = (orf.get("splice_event_id", ""), orf.get("protein_sequence_sha256", ""), orf.get("frame_status", ""))
        if key[1]:
            orf_generators[key].add(orf.get("source_generator", ""))

    causal_by_event: dict[str, list[dict[str, str]]] = defaultdict(list)
    causal_by_junction: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in tables.get("causal_links", []):
        causal_by_event[row.get("splice_event_id", "")].append(row)
        if row.get("junction_id"):
            causal_by_junction[row.get("junction_id", "")].append(row)

    normal_by_entity: dict[str, list[dict[str, str]]] = defaultdict(list)
    query_to_entities: dict[str, set[str]] = defaultdict(set)
    for q in tables.get("sequence_queries", []):
        qid = q.get("query_id", "")
        for key in (q.get("origin_peptide_id", ""), q.get("splice_event_id", ""), q.get("junction_id", "")):
            if qid and key:
                query_to_entities[qid].add(key)
    for row in tables.get("normal_background", []):
        keys = [row.get("origin_peptide_id", ""), row.get("splice_event_id", ""), row.get("junction_id", "")]
        keys.extend(query_to_entities.get(row.get("query_id", ""), set()))
        for key in keys:
            if key:
                normal_by_entity[key].append(row)

    rows: list[dict[str, str]] = []
    for origin in origins:
        event_id = origin.get("splice_event_id", "")
        por = origin.get("origin_peptide_id", "")
        jids = event_links.get(event_id, set()) | _tokens(origin.get("junction_ids", ""))

        # RNA-driven chain -------------------------------------------------
        same_peptide = origin_by_event_peptide.get((event_id, origin.get("peptide_sequence", "")), [])
        def _resolved_peptide_origin(row: dict[str, str]) -> bool:
            """Accept exact/partial translated peptide products, but never unresolved origins.

            ImmunoPepper deliberately labels short-read local translations as
            ``PARTIAL_TRANSLATED_CANDIDATE`` rather than claiming a full-length
            ORF.  They are nevertheless valid for *peptide-level* cross-generator
            agreement when the event and peptide sequence match exactly.
            """
            status = row.get("origin_status", "").upper()
            conflict = row.get("evidence_conflict_status", "").upper()
            if not status or any(token in status for token in ("UNRESOLVED", "AMBIGUOUS", "INVALID", "REJECTED")):
                return False
            return conflict in {"", "NONE", "RESOLVED"}

        peptide_generators = {
            row.get("source_generator", "") for row in same_peptide
            if row.get("source_generator", "") in _RNA_GENERATORS
            and _resolved_peptide_origin(row)
        }
        current_orf = orfs.get(origin.get("orf_id", ""), {})
        exact_orf_generators = {
            x for x in orf_generators.get((event_id, current_orf.get("protein_sequence_sha256", ""), current_orf.get("frame_status", "")), set())
            if x in _RNA_GENERATORS
        }
        exact_rna_evidence = []
        for jid in jids:
            exact_rna_evidence.extend([
                ev for ev in evidence_by_entity.get(jid, [])
                if ev.get("evidence_group") == "RNA_JUNCTION"
                and ev.get("resolution_status", "").startswith("RESOLVED")
                and str(ev.get("verified_value", "")) not in {"", "0", "0.0"}
            ])
        if len(exact_orf_generators) >= 2 and exact_rna_evidence:
            rna_status, rna_strength = "DUAL_GENERATOR_EXACT_ORF", "STRONG"
            limits: list[str] = []
        elif len(peptide_generators) >= 2 and exact_rna_evidence:
            rna_status, rna_strength = "DUAL_GENERATOR_EXACT_PEPTIDE", "MODERATE_STRONG"
            limits = ["ORF_LEVEL_CONSENSUS_NOT_ESTABLISHED"]
        elif origin.get("source_generator") in _RNA_GENERATORS and exact_rna_evidence:
            rna_status, rna_strength = "SINGLE_GENERATOR_RNA_SUPPORTED", "MODERATE"
            limits = ["SECOND_RNA_DRIVEN_GENERATOR_MISSING"]
        elif exact_rna_evidence:
            rna_status, rna_strength = "RNA_EVENT_ONLY", "WEAK"
            limits = ["RNA_DRIVEN_TRANSLATION_MISSING"]
        else:
            rna_status, rna_strength = "RNA_DRIVEN_UNCONFIRMED", "UNASSESSED"
            limits = ["EXACT_RNA_JUNCTION_SUPPORT_MISSING"]
        rows.append(_chain_row(
            event_id=event_id, origin=origin, sample_id=sample_id, chain_type="RNA_DRIVEN",
            status=rna_status, strength=rna_strength, groups=["RNA_JUNCTION", "RNA_DRIVEN_TRANSLATION"],
            tools=peptide_generators | exact_orf_generators | {ev.get("source_tool", "") for ev in exact_rna_evidence},
            entities=[event_id, por, origin.get("orf_id", ""), *jids],
            evidence_ids=[ev.get("evidence_id", "") for ev in exact_rna_evidence], limits=limits,
            reason=(
                f"exact_rna_junctions={len(exact_rna_evidence)}; peptide_generators={','.join(sorted(peptide_generators)) or 'NONE'}; "
                f"exact_orf_generators={','.join(sorted(exact_orf_generators)) or 'NONE'}"
            ),
        ))

        # DNA-causal chain -------------------------------------------------
        causal_rows = list(causal_by_event.get(event_id, []))
        for jid in jids:
            causal_rows.extend(causal_by_junction.get(jid, []))
        seen_causal: set[str] = set()
        causal_rows = [r for r in causal_rows if not (r.get("causal_link_id", "") in seen_causal or seen_causal.add(r.get("causal_link_id", "")))]
        best = max((row.get("causal_status", "UNASSESSED") for row in causal_rows), key=lambda x: _CAUSAL_PRIORITY.get(x, 0), default="UNASSESSED")
        if best == "TARGETED_REQUANT_SUPPORTED":
            strength = "STRONG"
        elif best in {"PVACSPLICE_SUPPORTED", "TARGETED_SPANNING_ONLY"}:
            strength = "MODERATE_STRONG"
        elif best == "DNA_RNA_CIS_SUPPORTED":
            strength = "MODERATE"
        elif best == "DNA_PREDICTION_ONLY":
            strength = "WEAK"
        else:
            strength = "UNASSESSED"
        rows.append(_chain_row(
            event_id=event_id, origin=origin, sample_id=sample_id, chain_type="DNA_CAUSAL",
            status=best, strength=strength, groups=["DNA_CAUSAL", "TARGETED_REQUANT"],
            tools=[tool for row in causal_rows for tool in _tokens(row.get("source_tools", ""))],
            entities=[row.get("causal_link_id", "") for row in causal_rows],
            evidence_ids=[ev.get("evidence_id", "") for row in causal_rows for ev in evidence_by_entity.get(row.get("causal_link_id", ""), [])],
            limits=[] if _CAUSAL_PRIORITY.get(best, 0) >= 2 else ["DNA_RNA_CAUSAL_LINK_INCOMPLETE"],
            reason=f"strongest_exact_causal_status={best}; causal_links={len(causal_rows)}",
        ))

        # Normal-background chain -----------------------------------------
        normal_rows: list[dict[str, str]] = []
        for key in [por, event_id, *jids]:
            normal_rows.extend(normal_by_entity.get(key, []))
        # Remove repeated row objects reached through multiple entities.
        seen = set()
        normal_rows = [r for r in normal_rows if not (r.get("normal_background_id", "") in seen or seen.add(r.get("normal_background_id", "")))]
        detected_critical = [r for r in normal_rows if r.get("assessment_status") == "NORMAL_DETECTED" and r.get("critical_tissue", "").lower() == "true"]
        detected = [r for r in normal_rows if r.get("assessment_status") == "NORMAL_DETECTED"]
        coverage_negative = {
            (r.get("normal_source_type", ""), r.get("normal_source", ""))
            for r in normal_rows if r.get("assessment_status") == "NOT_DETECTED_ADEQUATE_COVERAGE"
        }
        kmer_negative = {
            (r.get("normal_source_type", ""), r.get("normal_source", ""))
            for r in normal_rows if r.get("assessment_status") == "NOT_DETECTED_KMER_SCREEN"
        }
        if detected_critical:
            nstatus, nstrength, nlimits = "NORMAL_DETECTED_CRITICAL", "HARD_NEGATIVE", []
        elif detected:
            nstatus, nstrength, nlimits = "NORMAL_DETECTED", "NEGATIVE", []
        elif len(coverage_negative) >= 1 and len(kmer_negative) >= 1:
            nstatus, nstrength, nlimits = "LOCUS_AND_KMER_NEGATIVE", "STRONG", []
        elif len(coverage_negative) >= 2:
            nstatus, nstrength, nlimits = "MULTI_SOURCE_LOCUS_NEGATIVE", "STRONG", []
        elif len(coverage_negative) == 1:
            nstatus, nstrength, nlimits = "SINGLE_SOURCE_LOCUS_NEGATIVE", "MODERATE", ["SECOND_NORMAL_SOURCE_MISSING"]
        elif kmer_negative:
            nstatus, nstrength, nlimits = "K4NEO_NEGATIVE_ONLY", "WEAK_MODERATE", ["LOCUS_COVERAGE_NORMAL_EVIDENCE_MISSING"]
        else:
            nstatus, nstrength, nlimits = "NORMAL_BACKGROUND_UNASSESSED", "UNASSESSED", ["NORMAL_BACKGROUND_INCOMPLETE"]
        rows.append(_chain_row(
            event_id=event_id, origin=origin, sample_id=sample_id, chain_type="NORMAL_BACKGROUND",
            status=nstatus, strength=nstrength, groups=["NORMAL_BACKGROUND_LOCUS", "NORMAL_BACKGROUND_KMER"],
            tools=[r.get("normal_source", "") for r in normal_rows],
            entities=[r.get("normal_background_id", "") for r in normal_rows],
            evidence_ids=[], limits=nlimits,
            conflict="NORMAL_DETECTED" if detected else "NONE",
            reason=(
                f"detected={len(detected)}; coverage_negative_sources={len(coverage_negative)}; "
                f"kmer_negative_sources={len(kmer_negative)}"
            ),
        ))
    return rows
