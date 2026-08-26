"""splice2neo adapter for exact mutation-to-junction causal provenance.

The adapter consumes a headered TSV exported from a splice2neo R workflow. It
recognises the documented ``junc_id``, transcript/context sequence and peptide
context columns while accepting common explicit aliases. No gene-only or
nearest-locus linkage is performed.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from neoag.splice.coordinates import normalize_genome_build
from neoag.splice.identifiers import (
    link_id, orf_id, peptide_id, peptide_origin_id, sequence_sha256,
    splice_event_id, stable_id, transcript_hypothesis_id,
)

from ..sequence_queries import clean_dna, make_sequence_query
from ..variants import ExactEventIndex, canonical_variant_id, parse_junction_token, parse_variant_token, variant_from_row, variant_type
from .base import as_float_text, as_int, clean, get, join_tokens, read_delimited, row_hash, source_record_id, truth_text

_AA = set("ACDEFGHIKLMNPQRSTVWY")


def _aa(value: Any) -> str:
    seq = "".join(clean(value).split()).upper().replace("*", "")
    return seq if seq and set(seq) <= _AA else ""


def _score(row: Mapping[str, Any], *aliases: str) -> str:
    return as_float_text(get(row, *aliases))


def parse_splice2neo(
    path: str | Path,
    *,
    sample_id: str,
    genome_build: str = "GRCh38",
    source_tool_version: str = "UNASSESSED",
    entity_bundle: Mapping[str, list[Mapping[str, str]]] | None = None,
    strict: bool = False,
) -> dict[str, list[dict[str, str]]]:
    p = Path(path)
    build = normalize_genome_build(genome_build)
    result: dict[str, list[dict[str, str]]] = {
        "junctions": [], "events": [], "event_junction_links": [], "variants": [],
        "causal_links": [], "transcripts": [], "orfs": [], "peptide_origins": [],
        "peptide_origin_links": [], "sequence_queries": [], "tool_evidence": [], "conflicts": [],
    }
    event_index = ExactEventIndex.from_tables(entity_bundle or {})
    for row_no, row in enumerate(read_delimited(p), start=2):
        rid = source_record_id("splice2neo", p, row_no, row)
        junc_token = get(row, "junc_id", "junction_id", "junction", "canonical_junction_id")
        junction = parse_junction_token(junc_token, genome_build=build)
        if junction is None:
            result["conflicts"].append({
                "conflict_id": stable_id("CFL", "splice2neo", rid, "junction"),
                "entity_type": "SPLICE2NEO_ROW", "entity_id": rid, "sample_id": sample_id,
                "conflict_type": "SPLICE2NEO_JUNCTION_UNRESOLVED", "field_name": "junc_id",
                "observed_values": junc_token, "source_tools": "splice2neo", "source_record_ids": rid,
                "severity": "ERROR" if strict else "WARNING", "resolution_status": "UNRESOLVED",
                "resolution_reason": "A canonical or chr:start-end:strand junction is required.",
            })
            if strict:
                continue
            else:
                continue
        jid = junction.junction_id
        gene = get(row, "gene", "gene_name", "gene_symbol", "GENE_SYMBOL")
        gene_id = get(row, "gene_id", "ensembl_gene_id")
        event_type = get(row, "event_type", "splice_event_type", "type", default="MUTATION_ASSOCIATED_SPLICE").upper()
        resolved_event, _ = event_index.resolve([jid], gene=gene)
        event_id = resolved_event or splice_event_id(
            genome_build=build, event_type=event_type, strand=junction.strand,
            junction_ids=[jid], gene=gene or gene_id,
        )
        result["junctions"].append({
            "junction_id": jid, "sample_id": sample_id, "genome_build": build,
            "chrom": junction.chrom, "intron_start_1based": str(junction.intron_start_1based),
            "intron_end_1based": str(junction.intron_end_1based), "strand": junction.strand,
            "donor_1based": str(junction.donor_1based), "acceptor_1based": str(junction.acceptor_1based),
            "splice_motif": get(row, "splice_motif"), "annotation_status": "SPLICE2NEO_ASSOCIATED",
            "unique_split_reads": str(as_int(get(row, "junction_reads", "junc_reads", "read_count", "score"), 0)),
            "multi_split_reads": "0", "total_split_reads": str(as_int(get(row, "junction_reads", "junc_reads", "read_count", "score"), 0)),
            "max_overhang": get(row, "max_overhang", "anchor"),
            "source_coordinate_systems": "intron_1based_closed", "source_tools": "splice2neo",
            "source_tool_versions": source_tool_version, "source_files": str(p), "source_record_ids": rid,
            "provenance_record_count": "1", "junction_resolution_status": "RESOLVED_EXACT",
            "evidence_conflict_status": "NONE",
        })
        result["events"].append({
            "splice_event_id": event_id, "sample_id": sample_id, "genome_build": build,
            "event_type": event_type, "gene": gene, "gene_id": gene_id, "strand": junction.strand,
            "junction_ids": jid, "reference_junction_ids": get(row, "reference_junction_ids"),
            "alternative_junction_ids": jid, "affected_exons": get(row, "affected_exons", "exons"),
            "annotation_status": "SPLICE2NEO_MUTATION_ASSOCIATED", "cryptic_exon_status": get(row, "cryptic_exon_status", default="UNASSESSED"),
            "psi": _score(row, "psi"), "delta_psi": _score(row, "delta_psi", "dpsi"),
            "qvalue": _score(row, "qvalue", "FDR"), "outlier_score": _score(row, "outlier_score"),
            "event_expression": get(row, "event_expression", "junction_reads", "junc_reads"),
            "event_confidence": get(row, "event_confidence", default="MUTATION_ASSOCIATED"),
            "reference_path_status": get(row, "reference_path_status", default="UNASSESSED"),
            "cohort_analysis_status": get(row, "cohort_analysis_status", default="NOT_APPLICABLE"),
            "source_tools": "splice2neo", "source_tool_versions": source_tool_version,
            "source_files": str(p), "source_record_ids": rid, "provenance_record_count": "1",
            "event_resolution_status": "RESOLVED_EXACT_JUNCTION",
            "evidence_conflict_status": "NONE",
        })
        result["event_junction_links"].append({
            "event_junction_link_id": link_id("EJL", event_id, jid, "splice2neo"),
            "splice_event_id": event_id, "junction_id": jid, "sample_id": sample_id,
            "path_id": get(row, "tx_mod_id", "cts_id"), "path_role": "ALTERNATIVE",
            "edge_index": "1", "junction_role": "MUTATION_ASSOCIATED",
            "source_tool": "splice2neo", "source_record_id": rid, "link_status": "RESOLVED_EXACT",
        })

        parsed_variant = variant_from_row(row, genome_build=build)
        if not parsed_variant:
            token = get(row, "variant", "variant_key", "mut_id", "mutation_id")
            parsed_variant = parse_variant_token(token, genome_build=build) if token else None
        variant_id = ""
        if parsed_variant:
            vb, vc, vp, vr, va = parsed_variant
            variant_id = canonical_variant_id(vb, vc, vp, vr, va)
            scores = {
                "spliceai_score": _score(row, "spliceai_score", "SpliceAI", "spliceai"),
                "pangolin_score": _score(row, "pangolin_score", "Pangolin", "pangolin"),
                "mmsplice_score": _score(row, "mmsplice_score", "MMSplice", "mmsplice"),
                "ci_spliceai_score": _score(row, "ci_spliceai_score", "CI-SpliceAI", "ci_spliceai"),
            }
            tx_id = get(row, "tx_id", "transcript_id", "Transcript")
            result["variants"].append({
                "variant_id": variant_id, "sample_id": sample_id, "genome_build": vb,
                "chrom": vc, "pos_1based": str(vp), "ref": vr, "alt": va,
                "variant_type": variant_type(vr, va), "gene": gene, "gene_id": gene_id,
                "transcript_ids": tx_id, "hgvsc": get(row, "hgvsc", "HGVSc"),
                "hgvsp": get(row, "hgvsp", "HGVSp"), **scores,
                "source_tools": "splice2neo", "source_tool_versions": source_tool_version,
                "source_files": str(p), "source_record_ids": rid,
                "variant_resolution_status": "RESOLVED_EXACT", "evidence_conflict_status": "NONE",
            })
            reads = as_int(get(row, "junction_reads", "junc_reads", "read_count", "score"), 0)
            prediction_present = any(value for value in scores.values())
            causal_status = "DNA_RNA_CIS_SUPPORTED" if reads > 0 else "DNA_PREDICTION_ONLY"
            result["causal_links"].append({
                "causal_link_id": stable_id("DCL", variant_id, jid, event_id),
                "variant_id": variant_id, "junction_id": jid, "splice_event_id": event_id,
                "sample_id": sample_id, "gene": gene, "gene_id": gene_id, "transcript_id": tx_id,
                "causal_status": causal_status,
                "prediction_status": "PREDICTED_SPLICE_EFFECT" if prediction_present else "UNASSESSED",
                "rna_junction_status": "EXACT_RNA_SUPPORTED" if reads > 0 else "UNASSESSED",
                "targeted_requant_status": "UNASSESSED", "pvacsplice_status": "UNASSESSED",
                "junction_reads": str(reads), "easyquant_junction_reads": "", "easyquant_spanning_pairs": "",
                **scores, "source_tools": "splice2neo", "source_tool_versions": source_tool_version,
                "source_files": str(p), "source_record_ids": rid,
                "link_resolution_status": "RESOLVED_EXACT_VARIANT_AND_JUNCTION",
                "resolution_reason": "Variant and junction were supplied in the same splice2neo record.",
                "evidence_conflict_status": "NONE",
            })
            result["tool_evidence"].append({
                "evidence_id": stable_id("EVD", "splice2neo", variant_id, jid, rid),
                "entity_type": "CAUSAL_LINK", "entity_id": stable_id("DCL", variant_id, jid, event_id),
                "sample_id": sample_id, "evidence_group": "DNA_CAUSAL",
                "evidence_type": causal_status, "source_tool": "splice2neo",
                "source_tool_version": source_tool_version, "source_file": str(p),
                "source_row_number": str(row_no), "source_record_id": rid,
                "provided_value": f"{variant_id}|{jid}", "verified_value": event_id,
                "resolution_status": "RESOLVED_EXACT_VARIANT_AND_JUNCTION",
                "resolution_reason": "No gene-only or nearest-locus link was used.",
                "raw_payload_sha256": row_hash(row),
            })

        tx_id = get(row, "tx_id", "transcript_id", "Transcript")
        tx_mod_id = get(row, "tx_mod_id", "modified_transcript_id", default=tx_id)
        context = clean_dna(get(row, "cts_seq", "context_sequence", "modified_transcript_sequence"))
        position_raw = get(row, "cts_junc_pos", "context_junction_position", "junc_pos_context")
        # splice2neo can emit comma-separated positions for complex contexts. For
        # targeted evidence, only a unique position can be imported automatically.
        positions = [x.strip() for x in position_raw.split(",") if x.strip().isdigit()]
        query_position = int(positions[0]) if len(positions) == 1 else 0
        sth = ""
        if tx_id or context:
            sth = transcript_hypothesis_id(
                splice_event_id_value=event_id, reference_transcript_id=tx_id,
                junction_chain=[jid], path_role="MUTATION_ASSOCIATED_CONTEXT",
                sequence_sha256=sequence_sha256(context), source_generator="splice2neo",
            )
            result["transcripts"].append({
                "transcript_hypothesis_id": sth, "splice_event_id": event_id, "sample_id": sample_id,
                "gene": gene, "gene_id": gene_id, "reference_transcript_id": tx_id,
                "mane_status": get(row, "mane_status", "MANE Select", default="UNASSESSED"),
                "path_id": tx_mod_id, "path_role": "MUTATION_ASSOCIATED_CONTEXT",
                "exon_chain": get(row, "exon_chain"), "junction_chain": jid,
                "cds_start": get(row, "cds_start"), "cds_stop": get(row, "cds_stop"),
                "cds_phase_before_event": get(row, "cds_phase_before_event"),
                "cds_phase_after_event": get(row, "cds_phase_after_event"),
                "frame_status": get(row, "frame_status", "cds_description", default="UNASSESSED"),
                "translation_start_source": "SPLICE2NEO_CONTEXT",
                "transcript_expression_tpm": get(row, "transcript_expression_tpm", "tpm"),
                "full_length_status": "CONTEXT_SEQUENCE" if context else "STRUCTURE_ONLY",
                "long_read_support": "UNASSESSED", "nucleotide_sequence_sha256": sequence_sha256(context),
                "source_generator": "splice2neo", "source_generator_version": source_tool_version,
                "source_file": str(p), "source_record_id": rid,
                "hypothesis_status": "RESOLVED_CONTEXT" if context else "PARTIAL",
                "evidence_conflict_status": "NONE",
            })
        if context:
            query = make_sequence_query(
                sample_id=sample_id, query_type="TARGETED_JUNCTION_AND_NORMAL_SCREEN",
                nucleotide_sequence=context, position_1based=query_position,
                query_length=get(row, "query_length", default="62"), splice_event_id=event_id,
                junction_id=jid, variant_id=variant_id, transcript_hypothesis_id=sth,
                sequence_scope="SPLICE2NEO_CONTEXT", source_generator="splice2neo",
                source_file=str(p), source_record_id=rid, query_name=get(row, "cts_id", default=""),
            )
            if not query_position:
                query["query_status"] = "K4NEO_READY_EASYQUANT_AMBIGUOUS_POSITION"
            result["sequence_queries"].append(query)

        protein_context = _aa(get(row, "peptide_context", "protein_context", "translated_context"))
        explicit_peptide = _aa(get(row, "neo_peptide", "peptide", "epitope"))
        junc_in_orf = truth_text(get(row, "junc_in_orf", "junction_in_orf"))
        if protein_context and sth:
            frame = get(row, "frame_status", "cds_description", default="FRAME_RESOLVED" if junc_in_orf == "true" else "UNASSESSED")
            oid = orf_id(
                transcript_hypothesis_id_value=sth, protein_sequence_sha256=sequence_sha256(protein_context),
                frame_status=frame,
            )
            result["orfs"].append({
                "orf_id": oid, "transcript_hypothesis_id": sth, "splice_event_id": event_id,
                "sample_id": sample_id, "gene": gene, "protein_sequence": protein_context,
                "protein_sequence_sha256": sequence_sha256(protein_context),
                "protein_length": str(len(protein_context)), "orf_start": "", "orf_stop": "",
                "frame_status": frame, "frameshift_status": get(row, "frameshift_status", default="UNASSESSED"),
                "novel_aa_start": get(row, "novel_aa_start"), "novel_aa_end": get(row, "novel_aa_end"),
                "premature_stop_status": get(row, "premature_stop_status", default="UNASSESSED"),
                "nmd_risk": get(row, "nmd_risk", default="UNASSESSED"), "nmd_reason": get(row, "nmd_reason"),
                "orf_validity_status": "VALID_TRANSLATED_CONTEXT" if junc_in_orf == "true" else "PARTIAL_TRANSLATED_CONTEXT",
                "source_generator": "splice2neo", "source_generator_version": source_tool_version,
                "source_file": str(p), "source_record_id": rid, "evidence_conflict_status": "NONE",
            })
            result["tool_evidence"].append({
                "evidence_id": stable_id("EVD", "splice2neo", oid, rid), "entity_type": "ORF",
                "entity_id": oid, "sample_id": sample_id, "evidence_group": "DNA_CAUSAL_TRANSLATION",
                "evidence_type": "SPLICE2NEO_TRANSLATED_CONTEXT", "source_tool": "splice2neo",
                "source_tool_version": source_tool_version, "source_file": str(p),
                "source_row_number": str(row_no), "source_record_id": rid,
                "provided_value": protein_context, "verified_value": sequence_sha256(protein_context),
                "resolution_status": "RESOLVED_EVENT_TRANSCRIPT_ORF",
                "resolution_reason": f"junction_in_orf={junc_in_orf}", "raw_payload_sha256": row_hash(row),
            })
            if explicit_peptide and explicit_peptide in protein_context:
                start = protein_context.find(explicit_peptide) + 1
                pid = peptide_id(explicit_peptide)
                por = peptide_origin_id(
                    orf_id_value=oid, splice_event_id_value=event_id, peptide_sequence=explicit_peptide,
                    protein_start=start, protein_end=start + len(explicit_peptide) - 1,
                )
                result["peptide_origins"].append({
                    "origin_peptide_id": por, "peptide_id": pid, "orf_id": oid,
                    "transcript_hypothesis_id": sth, "splice_event_id": event_id, "sample_id": sample_id,
                    "gene": gene, "peptide_sequence": explicit_peptide, "peptide_length": str(len(explicit_peptide)),
                    "protein_start": str(start), "protein_end": str(start + len(explicit_peptide) - 1),
                    "crosses_junction": get(row, "crosses_junction", default="UNASSESSED"),
                    "junction_ids": jid, "required_junction_ids": jid,
                    "junction_offset_in_peptide": get(row, "junction_offset_in_peptide"),
                    "contains_novel_aa": get(row, "contains_novel_aa", default="true"),
                    "novel_aa_positions": get(row, "novel_aa_positions"),
                    "wildtype_counterpart_status": get(row, "wildtype_counterpart_status", default="UNRESOLVED"),
                    "wildtype_peptide": get(row, "wildtype_peptide"), "reference_proteome_match": get(row, "reference_proteome_match", default="UNASSESSED"),
                    "generator_group": "DNA_CAUSAL", "source_generator": "splice2neo",
                    "source_generator_version": source_tool_version, "source_file": str(p),
                    "source_record_id": rid, "origin_status": "RESOLVED_EXACT_SEQUENCE",
                    "evidence_conflict_status": "NONE",
                })
                result["peptide_origin_links"].append({
                    "peptide_origin_link_id": link_id("POL", pid, por), "peptide_id": pid,
                    "origin_peptide_id": por, "orf_id": oid, "transcript_hypothesis_id": sth,
                    "splice_event_id": event_id, "sample_id": sample_id, "link_status": "RESOLVED_EXACT",
                })
    return result
