"""pVACbind output adapter with exact FASTA-index provenance mapping."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any

from neoag.splice.identifiers import link_id, peptide_id, peptide_origin_id, sequence_sha256, stable_id

from .base import as_float_text, as_int, clean, get, infer_mhc_class, read_delimited, row_hash, source_record_id


def _boundary_crossing(offsets: str, start: int, end: int) -> str:
    """Assess whether a protein interval crosses an event-local AA boundary.

    Boundary tokens are emitted by the ImmunoPepper adapter. ``N-M`` denotes
    a boundary between residues; ``N`` denotes a codon that spans the junction.
    """
    if not offsets or start <= 0 or end <= 0:
        return "UNASSESSED"
    assessed = False
    for token in str(offsets).split(";"):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            left_text, right_text = token.split("-", 1)
            try:
                left, right = int(left_text), int(right_text)
            except ValueError:
                continue
            assessed = True
            if start <= left and end >= right:
                return "true"
        else:
            try:
                residue = int(token)
            except ValueError:
                continue
            assessed = True
            if start <= residue <= end:
                return "true"
    return "false" if assessed else "UNASSESSED"


def _inherit_parent_state(parents: list[dict[str, str]], field: str) -> str:
    values = {str(row.get(field, "")).strip() for row in parents if str(row.get(field, "")).strip()}
    if len(values) == 1:
        return next(iter(values))
    for preferred in ("true", "false", "UNASSESSED"):
        if preferred in values and all(value in {preferred, "UNASSESSED"} for value in values):
            return preferred
    return "UNASSESSED"


def _read_map(path: str | Path | None, bundle: dict[str, list[dict[str, str]]] | None) -> dict[str, list[dict[str, str]]]:
    rows: list[dict[str, str]] = []
    if path:
        rows.extend(read_delimited(path))
    elif bundle:
        rows.extend(bundle.get("pvacbind_fasta_map", []))
    index: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = get(row, "index", "Index", "fasta_id", "id")
        if key:
            index[key].append(row)
    return index


def parse_pvacbind(
    path: str | Path,
    *,
    sample_id: str,
    fasta_map: str | Path | None = None,
    source_tool_version: str = "UNASSESSED",
    entity_bundle: dict[str, list[dict[str, str]]] | None = None,
    result_scope: str = "ALL_EPITOPES",
) -> dict[str, list[dict[str, str]]]:
    p = Path(path)
    result = {
        "peptide_origins": [], "peptide_origin_links": [], "presentation": [],
        "tool_evidence": [], "conflicts": [],
    }
    map_index = _read_map(fasta_map, entity_bundle)
    events = {row["splice_event_id"]: row for row in (entity_bundle or {}).get("events", [])}
    orfs = {row["orf_id"]: row for row in (entity_bundle or {}).get("orfs", [])}
    transcripts = {row["transcript_hypothesis_id"]: row for row in (entity_bundle or {}).get("transcripts", [])}
    existing_origins: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    origins_by_orf: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in (entity_bundle or {}).get("peptide_origins", []):
        existing_origins[(row.get("orf_id", ""), row.get("peptide_sequence", ""))].append(row)
        origins_by_orf[row.get("orf_id", "")].append(row)

    for row_no, row in enumerate(read_delimited(p), start=2):
        record_id = source_record_id("pVACbind", p, row_no, row)
        index_value = get(row, "Index", "index", "Mutation", "sequence_name")
        mappings = map_index.get(index_value, [])
        epitope = get(row, "Epitope Seq", "epitope_seq", "MT Epitope Seq", "peptide", "sequence").upper()
        hla = get(row, "HLA Allele", "hla_allele", "Allele")
        mhc_class = infer_mhc_class(hla)
        if not epitope or not hla:
            result["conflicts"].append({
                "entity_type": "PRESENTATION", "entity_id": record_id, "sample_id": sample_id,
                "conflict_type": "PVACBIND_REQUIRED_FIELD_MISSING", "field_name": "HLA Allele/Epitope Seq",
                "observed_values": f"HLA={hla or 'MISSING'};epitope={epitope or 'MISSING'}",
                "source_tools": "pVACbind", "source_record_ids": record_id, "severity": "ERROR",
                "resolution_status": "UNRESOLVED",
                "resolution_reason": "A presentation record requires both an HLA allele and a non-empty epitope sequence.",
            })
            continue
        if len(mappings) != 1:
            result["conflicts"].append({
                "entity_type": "PRESENTATION", "entity_id": record_id, "sample_id": sample_id,
                "conflict_type": "PVACBIND_INDEX_AMBIGUOUS" if mappings else "PVACBIND_INDEX_UNRESOLVED",
                "field_name": "Index", "observed_values": index_value, "source_tools": "pVACbind",
                "source_record_ids": record_id, "severity": "ERROR", "resolution_status": "UNRESOLVED",
                "resolution_reason": f"FASTA map matches={len(mappings)}; no sequence-only fallback was used.",
            })
            result["tool_evidence"].append({
                "entity_type": "PRESENTATION", "entity_id": record_id, "sample_id": sample_id,
                "evidence_group": "PRESENTATION", "evidence_type": "PVACBIND_UNMAPPED_RESULT",
                "source_tool": "pVACbind", "source_tool_version": source_tool_version,
                "source_file": str(p), "source_row_number": str(row_no), "source_record_id": record_id,
                "provided_value": f"{index_value}|{hla}|{epitope}", "verified_value": "",
                "resolution_status": "UNRESOLVED", "resolution_reason": f"FASTA map matches={len(mappings)}",
                "raw_payload_sha256": row_hash(row),
            })
            continue
        mapping = mappings[0]
        oid = get(mapping, "orf_id")
        sth = get(mapping, "transcript_hypothesis_id")
        event_id = get(mapping, "splice_event_id")
        gene = get(mapping, "gene")
        orf = orfs.get(oid, {})
        transcript = transcripts.get(sth, {})
        event = events.get(event_id, {})
        protein = orf.get("protein_sequence", "")
        map_sequence_hash = get(mapping, "sequence_sha256")
        computed_sequence_hash = sequence_sha256(protein)
        provenance_errors: list[str] = []
        mapped_sample = get(mapping, "sample_id")
        if mapped_sample and mapped_sample != sample_id:
            provenance_errors.append(f"FASTA map sample mismatch: {mapped_sample} != {sample_id}")
        if not oid or not orf:
            provenance_errors.append(f"ORF missing: {oid or 'EMPTY'}")
        if not sth or not transcript:
            provenance_errors.append(f"transcript hypothesis missing: {sth or 'EMPTY'}")
        if not event_id or not event:
            provenance_errors.append(f"splice event missing: {event_id or 'EMPTY'}")
        if orf and orf.get("transcript_hypothesis_id", "") != sth:
            provenance_errors.append("ORF-to-transcript foreign key mismatch")
        if orf and orf.get("splice_event_id", "") != event_id:
            provenance_errors.append("ORF-to-event foreign key mismatch")
        if transcript and transcript.get("splice_event_id", "") != event_id:
            provenance_errors.append("transcript-to-event foreign key mismatch")
        if not protein:
            provenance_errors.append("ORF protein sequence missing")
        if not map_sequence_hash:
            provenance_errors.append("FASTA map sequence_sha256 missing")
        elif computed_sequence_hash and map_sequence_hash != computed_sequence_hash:
            provenance_errors.append("FASTA map sequence_sha256 does not match ORF sequence")
        orf_sequence_hash = orf.get("protein_sequence_sha256", "")
        if orf and not orf_sequence_hash:
            provenance_errors.append("stored ORF protein_sequence_sha256 missing")
        elif orf_sequence_hash and computed_sequence_hash and orf_sequence_hash != computed_sequence_hash:
            provenance_errors.append("stored ORF protein_sequence_sha256 does not match ORF sequence")
        if provenance_errors:
            result["conflicts"].append({
                "entity_type": "PRESENTATION", "entity_id": record_id, "sample_id": sample_id,
                "conflict_type": "PVACBIND_PROVENANCE_CHAIN_INVALID", "field_name": "Index/ORF/transcript/event/sequence_sha256",
                "observed_values": f"Index={index_value};orf={oid};transcript={sth};event={event_id};map_sha256={map_sequence_hash}",
                "source_tools": "pVACbind", "source_record_ids": record_id, "severity": "ERROR",
                "resolution_status": "UNRESOLVED", "resolution_reason": "; ".join(provenance_errors),
            })
            result["tool_evidence"].append({
                "entity_type": "PRESENTATION", "entity_id": record_id, "sample_id": sample_id,
                "evidence_group": "PRESENTATION", "evidence_type": "PVACBIND_PROVENANCE_REJECTED",
                "source_tool": "pVACbind", "source_tool_version": source_tool_version,
                "source_file": str(p), "source_row_number": str(row_no), "source_record_id": record_id,
                "provided_value": f"{index_value}|{hla}|{epitope}", "verified_value": "",
                "resolution_status": "UNRESOLVED", "resolution_reason": "; ".join(provenance_errors),
                "raw_payload_sha256": row_hash(row),
            })
            continue
        position = as_int(get(row, "Sub-peptide Position", "sub_peptide_position", "position"), 0)
        occurrences = [i for i in range(len(protein)) if epitope and protein.startswith(epitope, i)] if protein else []
        stated_matches = bool(position and protein and protein[position - 1: position - 1 + len(epitope)] == epitope)
        if position and protein and not stated_matches:
            if len(occurrences) == 1:
                corrected = occurrences[0] + 1
                result["conflicts"].append({
                    "entity_type": "PRESENTATION", "entity_id": record_id, "sample_id": sample_id,
                    "conflict_type": "PVACBIND_POSITION_CORRECTED", "field_name": "Sub-peptide Position",
                    "observed_values": f"stated={position};exact_sequence_position={corrected}",
                    "source_tools": "pVACbind", "source_record_ids": record_id, "severity": "WARNING",
                    "resolution_status": "RESOLVED_EXACT_SEQUENCE",
                    "resolution_reason": "The reported position did not match the exact mapped ORF; a unique exact epitope occurrence was used.",
                })
                position = corrected
            else:
                result["conflicts"].append({
                    "entity_type": "PRESENTATION", "entity_id": record_id, "sample_id": sample_id,
                    "conflict_type": "PVACBIND_EPITOPE_ORF_MISMATCH", "field_name": "Epitope Seq",
                    "observed_values": f"Index={index_value};position={position};epitope={epitope};occurrences={len(occurrences)}",
                    "source_tools": "pVACbind", "source_record_ids": record_id, "severity": "ERROR",
                    "resolution_status": "UNRESOLVED",
                    "resolution_reason": "The epitope could not be placed uniquely and exactly in the FASTA-mapped ORF.",
                })
                continue
        elif not position and protein and epitope:
            if len(occurrences) == 1:
                position = occurrences[0] + 1
            elif len(occurrences) != 1:
                result["conflicts"].append({
                    "entity_type": "PRESENTATION", "entity_id": record_id, "sample_id": sample_id,
                    "conflict_type": "PVACBIND_EPITOPE_POSITION_UNRESOLVED", "field_name": "Sub-peptide Position",
                    "observed_values": f"Index={index_value};epitope={epitope};occurrences={len(occurrences)}",
                    "source_tools": "pVACbind", "source_record_ids": record_id, "severity": "ERROR",
                    "resolution_status": "UNRESOLVED",
                    "resolution_reason": "No reported position and no unique exact epitope occurrence in the mapped ORF.",
                })
                continue
        origin_matches = existing_origins.get((oid, epitope), [])
        if len(origin_matches) == 1:
            origin = origin_matches[0]
            por = origin["origin_peptide_id"]
            pid = origin["peptide_id"]
        elif len(origin_matches) > 1:
            result["conflicts"].append({
                "entity_type": "PRESENTATION", "entity_id": record_id, "sample_id": sample_id,
                "conflict_type": "PVACBIND_PEPTIDE_ORIGIN_AMBIGUOUS", "field_name": "origin_peptide_id",
                "observed_values": ";".join(sorted(origin.get("origin_peptide_id", "") for origin in origin_matches)),
                "source_tools": "pVACbind", "source_record_ids": record_id, "severity": "ERROR",
                "resolution_status": "UNRESOLVED",
                "resolution_reason": "Multiple existing peptide origins share this ORF and epitope; presentation evidence was withheld.",
            })
            continue
        else:
            pid = peptide_id(epitope)
            por = peptide_origin_id(
                orf_id_value=oid, splice_event_id_value=event_id, peptide_sequence=epitope,
                protein_start=position or "", protein_end=(position + len(epitope) - 1) if position else "",
            )
            parent_origins = origins_by_orf.get(oid, [])
            crossing_states = {
                _boundary_crossing(parent.get("junction_offset_in_peptide", ""), position, position + len(epitope) - 1)
                for parent in parent_origins
                if parent.get("junction_offset_in_peptide")
            }
            crosses = "true" if "true" in crossing_states else ("false" if crossing_states == {"false"} else "UNASSESSED")
            result["peptide_origins"].append({
                "origin_peptide_id": por, "peptide_id": pid, "orf_id": oid,
                "transcript_hypothesis_id": sth, "splice_event_id": event_id, "sample_id": sample_id,
                "gene": gene, "peptide_sequence": epitope, "peptide_length": str(len(epitope)),
                "protein_start": str(position or ""),
                "protein_end": str(position + len(epitope) - 1) if position else "",
                "crosses_junction": crosses, "junction_ids": transcript.get("junction_chain", ""),
                "junction_offset_in_peptide": ";".join(sorted({p.get("junction_offset_in_peptide", "") for p in parent_origins if p.get("junction_offset_in_peptide")})),
                "contains_novel_aa": _inherit_parent_state(parent_origins, "contains_novel_aa"),
                "novel_aa_positions": "", "wildtype_counterpart_status": _inherit_parent_state(parent_origins, "wildtype_counterpart_status") or "WT_UNRESOLVED",
                "wildtype_peptide": "", "reference_proteome_match": _inherit_parent_state(parent_origins, "reference_proteome_match"),
                "generator_group": "PRESENTATION_DERIVED", "source_generator": "pVACbind",
                "source_generator_version": source_tool_version, "source_file": str(p),
                "source_record_id": record_id, "origin_status": "EPITOPE_FROM_MAPPED_ORF",
                "evidence_conflict_status": "ORIGIN_AMBIGUOUS" if len(origin_matches) > 1 else "NONE",
            })
            result["peptide_origin_links"].append({
                "peptide_origin_link_id": link_id("POL", pid, por, oid, sth, event_id),
                "peptide_id": pid, "origin_peptide_id": por, "orf_id": oid,
                "transcript_hypothesis_id": sth, "splice_event_id": event_id,
                "sample_id": sample_id, "link_status": "RESOLVED",
            })
        presentation_id = stable_id("PRE", por, hla, epitope, index_value, result_scope)
        presentation = {
            "presentation_id": presentation_id, "origin_peptide_id": por, "peptide_id": pid,
            "orf_id": oid, "transcript_hypothesis_id": sth, "splice_event_id": event_id,
            "sample_id": sample_id, "index": index_value, "hla_allele": hla, "mhc_class": mhc_class,
            "epitope_sequence": epitope, "epitope_length": str(len(epitope)),
            "sub_peptide_position": str(position or ""),
            "median_ic50": as_float_text(get(row, "Median IC50 Score", "median_ic50", "Median MT IC50 Score")),
            "best_ic50": as_float_text(get(row, "Best IC50 Score", "best_ic50", "Best MT IC50 Score")),
            "median_percentile": as_float_text(get(row, "Median Percentile", "median_percentile", "Median MT Percentile")),
            "best_percentile": as_float_text(get(row, "Best Percentile", "best_percentile", "Best MT Percentile")),
            "median_binding_percentile": as_float_text(get(row, "Median IC50 Percentile", "Median Binding Percentile", "median_binding_percentile")),
            "best_binding_percentile": as_float_text(get(row, "Best IC50 Percentile", "Best Binding Percentile", "best_binding_percentile")),
            "median_presentation_percentile": as_float_text(get(row, "Median Presentation Percentile", "median_presentation_percentile")),
            "best_presentation_percentile": as_float_text(get(row, "Best Presentation Percentile", "best_presentation_percentile")),
            "median_immunogenicity_percentile": as_float_text(get(row, "Median Immunogenicity Percentile", "median_immunogenicity_percentile")),
            "best_immunogenicity_percentile": as_float_text(get(row, "Best Immunogenicity Percentile", "best_immunogenicity_percentile")),
            "presentation_score": as_float_text(get(row, "MHCflurryEL Presentation Score", "presentation_score", "Best MHCflurryEL Presentation Score")),
            "immunogenicity_score": as_float_text(get(row, "Immunogenicity Score", "immunogenicity_score", "BigMHC_IM Score")),
            "prediction_methods": get(row, "Best IC50 Score Method", "Best Percentile Method", "prediction_methods"),
            "aggregate_tier": get(row, "Tier", "tier", "Evaluation"), "result_scope": result_scope,
            "source_tool": "pVACbind", "source_tool_version": source_tool_version,
            "source_file": str(p), "source_record_id": record_id,
            "mapping_status": "MAPPED_EXACT_FASTA_INDEX", "evidence_conflict_status": "NONE",
        }
        result["presentation"].append(presentation)
        result["tool_evidence"].append({
            "entity_type": "PEPTIDE_ORIGIN", "entity_id": por, "sample_id": sample_id,
            "evidence_group": "PRESENTATION", "evidence_type": "PVACBIND_PREDICTION",
            "source_tool": "pVACbind", "source_tool_version": source_tool_version,
            "source_file": str(p), "source_row_number": str(row_no), "source_record_id": record_id,
            "provided_value": f"{hla}|{epitope}", "verified_value": presentation_id,
            "resolution_status": "MAPPED_EXACT_FASTA_INDEX",
            "resolution_reason": f"Index={index_value}; HLA={hla}; epitope={epitope}",
            "raw_payload_sha256": row_hash(row),
        })
    return result
