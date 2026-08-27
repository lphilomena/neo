"""Canonical v0.5.1 Splice Provenance Layer table schemas.

v0.5.1 adds three independently auditable evidence chains:

* RNA-driven translation (ImmunoPepper + moPepGen)
* DNA-causal splicing (splice2neo + EasyQuant + pVACsplice)
* normal-background sequence screening (coverage-aware panels + k4neo)
"""
from __future__ import annotations

SPLICE_PROVENANCE_SCHEMA_VERSION = "0.5.3-splicemutr-normal-p0"

JUNCTION_FIELDS = [
    "junction_id", "sample_id", "genome_build", "chrom",
    "intron_start_1based", "intron_end_1based", "strand",
    "donor_1based", "acceptor_1based", "splice_motif", "annotation_status",
    "unique_split_reads", "multi_split_reads", "total_split_reads", "max_overhang",
    "source_coordinate_systems", "source_tools", "source_tool_versions",
    "source_files", "source_record_ids", "provenance_record_count",
    "junction_resolution_status", "evidence_conflict_status",
]

JUNCTION_READ_QC_FIELDS = [
    "junction_read_qc_id", "junction_id", "sample_id", "source_tool",
    "source_tool_version", "source_assay_id", "source_file", "source_record_id",
    "resolution_status", "unique_split_reads", "total_split_reads",
    "unique_fragment_starts", "max_overhang", "mapping_quality",
    "multimapping_fraction", "tumor_psi", "caller_filter_status", "qc_policy",
    "min_unique_split_reads", "min_unique_fragment_starts", "min_overhang",
    "min_mapping_quality", "max_multimapping_fraction", "min_tumor_psi",
    "qc_status", "failed_checks", "missing_checks", "qc_reason",
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
    "crosses_junction", "junction_ids", "required_junction_ids",
    "junction_offset_in_peptide", "contains_novel_aa", "structural_novelty_status",
    "tumor_specificity_status", "cohort_analysis_status",
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

VARIANT_FIELDS = [
    "variant_id", "sample_id", "genome_build", "chrom", "pos_1based", "ref", "alt",
    "variant_type", "gene", "gene_id", "transcript_ids", "hgvsc", "hgvsp",
    "spliceai_score", "pangolin_score", "mmsplice_score", "ci_spliceai_score",
    "source_tools", "source_tool_versions", "source_files", "source_record_ids",
    "variant_resolution_status", "evidence_conflict_status",
]

CAUSAL_LINK_FIELDS = [
    "causal_link_id", "variant_id", "junction_id", "splice_event_id", "sample_id",
    "gene", "gene_id", "transcript_id", "causal_status", "prediction_status",
    "rna_junction_status", "targeted_requant_status", "pvacsplice_status",
    "junction_reads", "easyquant_junction_reads", "easyquant_spanning_pairs",
    "spliceai_score", "pangolin_score", "mmsplice_score", "ci_spliceai_score",
    "source_tools", "source_tool_versions", "source_files", "source_record_ids",
    "link_resolution_status", "resolution_reason", "evidence_conflict_status",
]

SEQUENCE_QUERY_FIELDS = [
    "query_id", "sample_id", "query_name", "query_type", "splice_event_id", "junction_id",
    "variant_id", "transcript_hypothesis_id", "orf_id", "origin_peptide_id",
    "nucleotide_sequence", "sequence_sha256", "position_1based", "query_length",
    "sequence_scope", "source_generator", "source_file", "source_record_id",
    "query_status", "evidence_conflict_status",
]

TARGETED_QUANT_FIELDS = [
    "targeted_quant_id", "query_id", "sample_id", "query_name", "splice_event_id",
    "junction_id", "variant_id", "position_1based", "junction_reads", "spanning_pairs",
    "max_anchor", "left_reads", "right_reads", "interval", "within_interval",
    "coverage_percent", "coverage_mean", "coverage_median", "support_status",
    "source_tool", "source_tool_version", "source_file", "source_record_id",
    "mapping_status", "evidence_conflict_status",
]

PVACSPLICE_PREDICTION_FIELDS = [
    "pvacsplice_prediction_id", "causal_link_id", "variant_id", "junction_id",
    "splice_event_id", "origin_peptide_id", "peptide_id", "sample_id", "chrom",
    "variant_start_1based", "variant_stop_1based", "ref", "alt", "junction_score",
    "junction_anchor", "transcript_id", "gene", "gene_id", "hla_allele", "mhc_class",
    "epitope_sequence", "epitope_length", "protein_position", "best_ic50",
    "best_percentile", "best_binding_percentile", "best_presentation_percentile",
    "presentation_score", "immunogenicity_score", "prediction_methods", "aggregate_tier",
    "source_tool", "source_tool_version", "source_file", "source_record_id",
    "mapping_status", "evidence_conflict_status",
]

NORMAL_BACKGROUND_FIELDS = [
    "normal_background_id", "splice_event_id", "junction_id", "origin_peptide_id", "query_id",
    "sample_id", "normal_source", "normal_source_type", "normal_tissue", "critical_tissue",
    "developmental_stage", "study_id", "detection_status", "coverage_status",
    "junction_reads", "normal_total_junction_reads", "sample_count", "total_samples",
    "sample_prevalence", "normal_tissue_count", "normal_tissues",
    "source_dataset", "reference_release",
    "kmer_prevalence", "uniqueness_rate", "assessment_status", "assessment_reason",
    "source_file", "source_record_id", "evidence_conflict_status",
]

TOOL_EVIDENCE_FIELDS = [
    "evidence_id", "entity_type", "entity_id", "sample_id", "evidence_group", "evidence_type",
    "source_tool", "source_tool_version", "source_assay_id", "source_file", "source_row_number", "source_record_id",
    "provided_value", "verified_value", "resolution_status", "resolution_reason", "raw_payload_sha256",
]

EVIDENCE_CHAIN_FIELDS = [
    "evidence_chain_id", "splice_event_id", "origin_peptide_id", "peptide_id", "sample_id",
    "chain_type", "chain_status", "chain_strength", "independent_source_groups",
    "source_tools", "supporting_entity_ids", "supporting_evidence_ids",
    "limiting_reasons", "conflict_status", "chain_reason",
]

CONSENSUS_FIELDS = [
    "consensus_id", "splice_event_id", "origin_peptide_id", "peptide_id", "sample_id",
    "event_evidence_grade", "orf_evidence_grade", "normal_safety_grade", "presentation_grade",
    "independent_evidence_groups", "independent_rna_sources", "independent_translation_generators",
    "rna_driven_chain_status", "dna_causal_chain_status", "normal_background_chain_status",
    "translation_consensus_level", "independent_peptide_generators",
    "required_junction_ids", "required_junction_qc_status",
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
    "junction_read_qc": "splice_junction_read_qc.tsv",
    "events": "splice_events.tsv",
    "event_junction_links": "splice_event_junction_links.tsv",
    "transcripts": "splice_transcript_hypotheses.tsv",
    "orfs": "splice_orfs.tsv",
    "peptide_origins": "splice_peptide_origins.tsv",
    "peptide_origin_links": "splice_peptide_origin_links.tsv",
    "variants": "splice_variants.tsv",
    "causal_links": "splice_causal_links.tsv",
    "sequence_queries": "splice_sequence_queries.tsv",
    "targeted_quantification": "splice_targeted_quantification.tsv",
    "pvacsplice_predictions": "splice_pvacsplice_predictions.tsv",
    "presentation": "splice_pvacbind_predictions.tsv",
    "normal_background": "splice_normal_background.tsv",
    "tool_evidence": "splice_tool_evidence.long.tsv",
    "evidence_chains": "splice_evidence_chains.tsv",
    "consensus": "splice_consensus.tsv",
    "conflicts": "splice_conflicts.tsv",
    "qc": "splice_qc.tsv",
    "pvacbind_fasta": "splice_pvacbind_input.fasta",
    "pvacbind_fasta_map": "splice_pvacbind_fasta_map.tsv",
    "easyquant_input": "splice_easyquant_input.tsv",
    "easyquant_query_map": "splice_easyquant_query_map.tsv",
    "k4neo_input": "splice_k4neo_input.tsv",
    "k4neo_query_map": "splice_k4neo_query_map.tsv",
    "raw_events": "raw_events.tsv",
    "raw_peptides": "raw_peptides.tsv",
    "rna_junction_evidence": "rna_junction_evidence.tsv",
    "manifest": "provenance_manifest.json",
}

TABLE_FIELDS = {
    "junctions": JUNCTION_FIELDS,
    "junction_read_qc": JUNCTION_READ_QC_FIELDS,
    "events": EVENT_FIELDS,
    "event_junction_links": EVENT_JUNCTION_LINK_FIELDS,
    "transcripts": TRANSCRIPT_HYPOTHESIS_FIELDS,
    "orfs": ORF_FIELDS,
    "peptide_origins": PEPTIDE_ORIGIN_FIELDS,
    "peptide_origin_links": PEPTIDE_ORIGIN_LINK_FIELDS,
    "variants": VARIANT_FIELDS,
    "causal_links": CAUSAL_LINK_FIELDS,
    "sequence_queries": SEQUENCE_QUERY_FIELDS,
    "targeted_quantification": TARGETED_QUANT_FIELDS,
    "pvacsplice_predictions": PVACSPLICE_PREDICTION_FIELDS,
    "presentation": PRESENTATION_FIELDS,
    "normal_background": NORMAL_BACKGROUND_FIELDS,
    "tool_evidence": TOOL_EVIDENCE_FIELDS,
    "evidence_chains": EVIDENCE_CHAIN_FIELDS,
    "consensus": CONSENSUS_FIELDS,
    "conflicts": CONFLICT_FIELDS,
    "qc": QC_FIELDS,
    "pvacbind_fasta_map": PVACBIND_FASTA_MAP_FIELDS,
    "easyquant_query_map": SEQUENCE_QUERY_FIELDS,
    "k4neo_query_map": SEQUENCE_QUERY_FIELDS,
}
