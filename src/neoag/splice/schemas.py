"""Canonical v0.5.0 Splice Provenance Layer table schemas."""
from __future__ import annotations

SPLICE_PROVENANCE_SCHEMA_VERSION = "0.5.0"

JUNCTION_FIELDS = [
    "junction_id", "sample_id", "genome_build", "chrom",
    "intron_start_1based", "intron_end_1based", "strand",
    "donor_1based", "acceptor_1based", "splice_motif", "annotation_status",
    "unique_split_reads", "multi_split_reads", "total_split_reads", "max_overhang",
    "source_coordinate_systems", "source_tools", "source_tool_versions",
    "source_files", "source_record_ids", "provenance_record_count",
    "junction_resolution_status", "evidence_conflict_status",
]

EVENT_FIELDS = [
    "splice_event_id", "sample_id", "genome_build", "event_type", "gene", "gene_id", "strand",
    "junction_ids", "reference_junction_ids", "alternative_junction_ids", "affected_exons",
    "annotation_status", "cryptic_exon_status", "psi", "delta_psi", "qvalue", "outlier_score",
    "event_expression", "event_confidence", "reference_path_status", "cohort_analysis_status",
    "source_tools", "source_tool_versions", "source_files", "source_record_ids",
    "provenance_record_count", "event_resolution_status", "evidence_conflict_status",
]

EVENT_JUNCTION_LINK_FIELDS = [
    "event_junction_link_id", "splice_event_id", "junction_id", "sample_id",
    "path_id", "path_role", "edge_index", "junction_role", "source_tool",
    "source_record_id", "link_status",
]

TRANSCRIPT_HYPOTHESIS_FIELDS = [
    "transcript_hypothesis_id", "splice_event_id", "sample_id", "gene", "gene_id",
    "reference_transcript_id", "mane_status", "path_id", "path_role", "exon_chain",
    "junction_chain", "cds_start", "cds_stop", "cds_phase_before_event", "cds_phase_after_event",
    "frame_status", "translation_start_source", "transcript_expression_tpm",
    "full_length_status", "long_read_support", "nucleotide_sequence_sha256",
    "source_generator", "source_generator_version", "source_file", "source_record_id",
    "hypothesis_status", "evidence_conflict_status",
]

ORF_FIELDS = [
    "orf_id", "transcript_hypothesis_id", "splice_event_id", "sample_id", "gene",
    "protein_sequence", "protein_sequence_sha256", "protein_length", "orf_start", "orf_stop",
    "frame_status", "frameshift_status", "novel_aa_start", "novel_aa_end",
    "premature_stop_status", "nmd_risk", "nmd_reason", "orf_validity_status",
    "source_generator", "source_generator_version", "source_file", "source_record_id",
    "evidence_conflict_status",
]

PEPTIDE_ORIGIN_FIELDS = [
    "origin_peptide_id", "peptide_id", "orf_id", "transcript_hypothesis_id", "splice_event_id",
    "sample_id", "gene", "peptide_sequence", "peptide_length", "protein_start", "protein_end",
    "crosses_junction", "junction_ids", "junction_offset_in_peptide", "contains_novel_aa",
    "novel_aa_positions", "wildtype_counterpart_status", "wildtype_peptide",
    "reference_proteome_match", "generator_group", "source_generator", "source_generator_version",
    "source_file", "source_record_id", "origin_status", "evidence_conflict_status",
]

PEPTIDE_ORIGIN_LINK_FIELDS = [
    "peptide_origin_link_id", "peptide_id", "origin_peptide_id", "orf_id",
    "transcript_hypothesis_id", "splice_event_id", "sample_id", "link_status",
]

PRESENTATION_FIELDS = [
    "presentation_id", "origin_peptide_id", "peptide_id", "orf_id", "transcript_hypothesis_id",
    "splice_event_id", "sample_id", "index", "hla_allele", "mhc_class", "epitope_sequence",
    "epitope_length", "sub_peptide_position", "median_ic50", "best_ic50",
    "median_percentile", "best_percentile", "median_binding_percentile", "best_binding_percentile",
    "median_presentation_percentile", "best_presentation_percentile",
    "median_immunogenicity_percentile", "best_immunogenicity_percentile",
    "presentation_score", "immunogenicity_score", "prediction_methods",
    "aggregate_tier", "result_scope", "source_tool", "source_tool_version", "source_file",
    "source_record_id", "mapping_status", "evidence_conflict_status",
]

NORMAL_BACKGROUND_FIELDS = [
    "normal_background_id", "splice_event_id", "junction_id", "origin_peptide_id", "sample_id",
    "normal_source", "normal_source_type", "normal_tissue", "critical_tissue",
    "detection_status", "coverage_status", "junction_reads", "sample_prevalence",
    "kmer_prevalence", "assessment_status", "assessment_reason", "source_file",
    "source_record_id", "evidence_conflict_status",
]

TOOL_EVIDENCE_FIELDS = [
    "evidence_id", "entity_type", "entity_id", "sample_id", "evidence_group", "evidence_type",
    "source_tool", "source_tool_version", "source_assay_id", "source_file", "source_row_number", "source_record_id",
    "provided_value", "verified_value", "resolution_status", "resolution_reason", "raw_payload_sha256",
]

CONSENSUS_FIELDS = [
    "consensus_id", "splice_event_id", "origin_peptide_id", "peptide_id", "sample_id",
    "event_evidence_grade", "orf_evidence_grade", "normal_safety_grade", "presentation_grade",
    "independent_evidence_groups", "independent_rna_sources", "independent_translation_generators",
    "event_consensus_status", "orf_consensus_status", "normal_background_status",
    "final_evidence_tier", "priority_cap", "consensus_reason", "hard_fail_codes", "cap_codes",
]

CONFLICT_FIELDS = [
    "conflict_id", "entity_type", "entity_id", "sample_id", "conflict_type", "field_name",
    "observed_values", "source_tools", "source_record_ids", "severity", "resolution_status",
    "resolution_reason",
]

QC_FIELDS = ["metric", "value", "status", "detail"]

PVACBIND_FASTA_MAP_FIELDS = [
    "index", "orf_id", "transcript_hypothesis_id", "splice_event_id", "sample_id", "gene",
    "sequence_type", "sequence_sha256", "source_generator", "source_record_id",
]

OUTPUT_FILENAMES = {
    "junctions": "splice_junctions.tsv",
    "events": "splice_events.tsv",
    "event_junction_links": "splice_event_junction_links.tsv",
    "transcripts": "splice_transcript_hypotheses.tsv",
    "orfs": "splice_orfs.tsv",
    "peptide_origins": "splice_peptide_origins.tsv",
    "peptide_origin_links": "splice_peptide_origin_links.tsv",
    "presentation": "splice_pvacbind_predictions.tsv",
    "normal_background": "splice_normal_background.tsv",
    "tool_evidence": "splice_tool_evidence.long.tsv",
    "consensus": "splice_consensus.tsv",
    "conflicts": "splice_conflicts.tsv",
    "qc": "splice_qc.tsv",
    "pvacbind_fasta": "splice_pvacbind_input.fasta",
    "pvacbind_fasta_map": "splice_pvacbind_fasta_map.tsv",
    "raw_events": "raw_events.tsv",
    "raw_peptides": "raw_peptides.tsv",
    "rna_junction_evidence": "rna_junction_evidence.tsv",
    "manifest": "provenance_manifest.json",
}

TABLE_FIELDS = {
    "junctions": JUNCTION_FIELDS,
    "events": EVENT_FIELDS,
    "event_junction_links": EVENT_JUNCTION_LINK_FIELDS,
    "transcripts": TRANSCRIPT_HYPOTHESIS_FIELDS,
    "orfs": ORF_FIELDS,
    "peptide_origins": PEPTIDE_ORIGIN_FIELDS,
    "peptide_origin_links": PEPTIDE_ORIGIN_LINK_FIELDS,
    "presentation": PRESENTATION_FIELDS,
    "normal_background": NORMAL_BACKGROUND_FIELDS,
    "tool_evidence": TOOL_EVIDENCE_FIELDS,
    "consensus": CONSENSUS_FIELDS,
    "conflicts": CONFLICT_FIELDS,
    "qc": QC_FIELDS,
    "pvacbind_fasta_map": PVACBIND_FASTA_MAP_FIELDS,
}
