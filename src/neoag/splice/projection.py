"""Backward-compatible projection from formal splice entities to NeoAg catalogs."""
from __future__ import annotations

from collections import defaultdict

from neoag.model_layers import enrich_event_layers, enrich_peptide_layers
from neoag.schemas import EVENT_FIELDS as LEGACY_EVENT_FIELDS, PEPTIDE_FIELDS as LEGACY_PEPTIDE_FIELDS, RNA_JUNCTION_EVIDENCE_FIELDS
from neoag.utils import safe_id


def _split(value: str) -> list[str]:
    return [x for x in str(value or "").split(";") if x]


def project_legacy(tables: dict[str, list[dict[str, str]]], *, sample_id: str, disease_profile: str = "default") -> dict[str, list[dict[str, str]]]:
    junctions = {row.get("junction_id", ""): row for row in tables.get("junctions", [])}
    consensus_by_origin = {row.get("origin_peptide_id", ""): row for row in tables.get("consensus", [])}
    event_consensus: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in tables.get("consensus", []):
        event_consensus[row.get("splice_event_id", "")].append(row)
    presentations: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in tables.get("presentation", []):
        presentations[row.get("origin_peptide_id", "")].append(row)

    raw_events: list[dict[str, str]] = []
    for event in tables.get("events", []):
        event_id = event.get("splice_event_id", "")
        jids = _split(event.get("junction_ids", ""))
        jrows = [junctions[j] for j in jids if j in junctions]
        exact_counts = [int(float(row.get("total_split_reads") or 0)) for row in jrows]
        verified = min(exact_counts) if jids and len(exact_counts) == len(jids) else 0
        positive = sum(1 for x in exact_counts if x > 0)
        if verified > 0:
            support_status = "SUPPORTED_ALL_EVENT_JUNCTIONS_EXACT"
        elif positive:
            support_status = "PARTIAL_EVENT_JUNCTION_SUPPORT"
        else:
            support_status = "NO_EXACT_RNA_JUNCTION_SUPPORT"
        cands = event_consensus.get(event_id, [])
        best = sorted(cands, key=lambda x: x.get("final_evidence_tier", "R9"))[0] if cands else {}
        canonical = jids[0] if len(jids) == 1 else ""
        j0 = junctions.get(canonical, {})
        row = {
            "event_id": event_id, "splice_event_id": event_id, "sample_id": sample_id,
            "disease_profile": disease_profile, "event_type": "Splice", "splice_event_type": event.get("event_type", ""),
            "mutation_source": "RNA_splice", "peptide_consequence": "splice_derived",
            "evidence_scope": "FORMAL_SPLICE_PROVENANCE_V050", "priority_cap": best.get("priority_cap", ""),
            "gene": event.get("gene", ""), "event_name": event.get("event_type", ""),
            "genome_build": event.get("genome_build", ""), "canonical_junction_id": canonical,
            "junction_ids": event.get("junction_ids", ""), "source_junction_id": "",
            "junction_chrom": j0.get("chrom", ""), "junction_start": j0.get("intron_start_1based", ""),
            "junction_end": j0.get("intron_end_1based", ""), "junction_strand": j0.get("strand", ""),
            "junction_donor": j0.get("donor_1based", ""), "junction_acceptor": j0.get("acceptor_1based", ""),
            "junction_coordinate_system": "intron_1based_closed", "junction_resolution_status": event.get("event_resolution_status", ""),
            "junction_resolution_reason": "Formal Splice Provenance Layer event projection.",
            "junction_match_status": "EXACT_EVENT_JUNCTION_SET", "junction_match_method": "canonical_event_junction_links",
            "junction_support_status": support_status, "junction_support_conflict": "",
            "junction_support_reason": "Conservative minimum across all exact event junctions; partial sets are not promoted.",
            "provided_rna_junction_reads": str(verified), "rna_junction_reads": str(verified),
            "rna_junction_source": ";".join(sorted({x for jr in jrows for x in _split(jr.get("source_tools", "")) if x in {"RegTools", "STAR", "STAR-SJ"}})),
            "normal_junction_status": best.get("normal_background_status", "UNASSESSED"),
            "rna_frame_status": "UNASSESSED", "event_confidence": event.get("event_confidence", ""),
            "event_expression": event.get("event_expression", ""), "rna_support_status": support_status,
            "rna_evidence_completeness": "COMPLETE" if verified > 0 else "PARTIAL",
            "tumor_specificity": "0.5", "source_file": event.get("source_files", ""),
            "source_record_id": event.get("source_record_ids", ""), "source_tools": event.get("source_tools", ""),
            "source_records": event.get("source_record_ids", ""), "provenance_record_count": event.get("provenance_record_count", "1"),
            "evidence_conflict_status": event.get("evidence_conflict_status", "NONE"),
            "splice_event_evidence_grade": best.get("event_evidence_grade", ""),
            "orf_evidence_grade": best.get("orf_evidence_grade", ""),
            "normal_safety_grade": best.get("normal_safety_grade", ""),
            "splice_consensus_tier": best.get("final_evidence_tier", ""),
            "source": "SpliceProvenanceLayer_v0.5.1",
        }
        raw_events.append(enrich_event_layers(row))

    raw_peptides: list[dict[str, str]] = []
    for origin in tables.get("peptide_origins", []):
        por = origin.get("origin_peptide_id", "")
        preds = presentations.get(por, []) or [{}]
        consensus = consensus_by_origin.get(por, {})
        jids = _split(origin.get("junction_ids", ""))
        jrows = [junctions[j] for j in jids if j in junctions]
        counts = [int(float(x.get("total_split_reads") or 0)) for x in jrows]
        verified = min(counts) if jids and len(counts) == len(jids) else 0
        for pred_idx, pred in enumerate(preds, start=1):
            peptide = pred.get("epitope_sequence") or origin.get("peptide_sequence", "")
            hla = pred.get("hla_allele", "")
            peptide_row = {
                "peptide_id": safe_id(f"{por}|{hla}|{peptide}|{pred_idx}"),
                "event_id": origin.get("splice_event_id", ""), "splice_event_id": origin.get("splice_event_id", ""),
                "transcript_hypothesis_id": origin.get("transcript_hypothesis_id", ""), "orf_id": origin.get("orf_id", ""),
                "origin_peptide_id": por, "sample_id": sample_id, "event_type": "Splice",
                "splice_event_type": "", "mutation_source": "RNA_splice", "peptide_consequence": "splice_derived",
                "evidence_scope": "FORMAL_SPLICE_PROVENANCE_V050", "priority_cap": consensus.get("priority_cap", ""),
                "gene": origin.get("gene", ""), "peptide": peptide, "wildtype_peptide": origin.get("wildtype_peptide", ""),
                "genome_build": junctions.get(jids[0], {}).get("genome_build", "") if jids else "",
                "canonical_junction_id": jids[0] if len(jids) == 1 else "", "junction_ids": origin.get("junction_ids", ""),
                "junction_coordinate_system": "intron_1based_closed", "junction_resolution_status": "RESOLVED" if jids else "UNRESOLVED",
                "junction_match_status": "EXACT_ORIGIN_LINK", "junction_match_method": "peptide_origin_to_event_to_junction",
                "junction_support_status": "SUPPORTED_ALL_EVENT_JUNCTIONS_EXACT" if verified > 0 else "NO_OR_PARTIAL_EXACT_SUPPORT",
                "provided_rna_junction_reads": str(verified), "rna_junction_reads": str(verified),
                "rna_junction_source": ";".join(sorted({x for jr in jrows for x in _split(jr.get("source_tools", "")) if x in {"RegTools", "STAR", "STAR-SJ"}})),
                "rna_frame_status": "UNASSESSED", "crosses_junction": origin.get("crosses_junction", ""),
                "contains_novel_aa": origin.get("contains_novel_aa", ""), "hla_allele": hla,
                "mhc_class": pred.get("mhc_class", ""), "source_tool": origin.get("source_generator", ""),
                "source_file": origin.get("source_file", ""), "source_record_id": origin.get("source_record_id", ""),
                "source_tools": origin.get("source_generator", ""), "source_records": origin.get("source_record_id", ""),
                "provenance_record_count": "1", "evidence_conflict_status": origin.get("evidence_conflict_status", "NONE"),
                "generation_status": origin.get("origin_status", ""),
                "binding_rank": pred.get("best_percentile", ""), "presentation_score": pred.get("presentation_score", ""),
                "immunogenicity_score": pred.get("immunogenicity_score", ""),
                "normal_junction_assessment_status": consensus.get("normal_background_status", "UNASSESSED"),
                "splice_event_evidence_grade": consensus.get("event_evidence_grade", ""),
                "orf_evidence_grade": consensus.get("orf_evidence_grade", ""),
                "normal_safety_grade": consensus.get("normal_safety_grade", ""),
                "independent_translation_generators": consensus.get("independent_translation_generators", ""),
                "splice_consensus_tier": consensus.get("final_evidence_tier", ""),
                "final_priority": consensus.get("final_evidence_tier", ""),
                "recommended_use": "EXPERIMENTAL_VALIDATION" if consensus.get("final_evidence_tier") in {"R1", "R2"} else "REVIEW",
            }
            raw_peptides.append(enrich_peptide_layers(peptide_row))

    rna_evidence: list[dict[str, str]] = []
    for link in tables.get("event_junction_links", []):
        junction = junctions.get(link.get("junction_id", ""), {})
        reads = junction.get("total_split_reads", "0")
        rna_evidence.append({
            "evidence_id": link.get("event_junction_link_id", ""), "event_id": link.get("splice_event_id", ""),
            "peptide_id": "", "sample_id": sample_id, "gene": "", "gene_pair": "",
            "junction_reads": reads, "junction_source": junction.get("source_tools", ""),
            "mutation_source": "RNA_splice", "peptide_consequence": "splice_derived",
            "rna_frame_status": "UNASSESSED",
            "rna_support_status": "SUPPORTED_EXACT_JUNCTION" if str(reads) not in {"", "0", "0.0"} else "NO_EXACT_READ_SUPPORT",
            "rna_evidence_completeness": "COMPLETE" if str(reads) not in {"", "0", "0.0"} else "PARTIAL",
            "rna_evidence_score": "1" if str(reads) not in {"", "0", "0.0"} else "0",
        })
    return {"raw_events": raw_events, "raw_peptides": raw_peptides, "rna_junction_evidence": rna_evidence}
