# Layer 1: mutation_source | Layer 2: peptide_consequence | Layer 3: scored in ranked_peptides (l3_*)
EVENT_FIELDS = [
    "event_id","splice_event_id","sample_id","disease_profile","event_type","mutation_source","peptide_consequence",
    "evidence_scope","priority_cap","wes_confidence_tier",
    "gene","event_name","splice_event_type","junction_ids",
    "genome_build","canonical_junction_id","source_junction_id",
    "junction_chrom","junction_start","junction_end","junction_strand","junction_donor","junction_acceptor",
    "junction_coordinate_system","junction_resolution_status","junction_resolution_reason",
    "junction_match_status","junction_match_method","junction_support_status","junction_support_conflict","junction_support_reason",
    "provided_rna_junction_reads",
    "cancer_gene_list_status","cancer_gene_symbols","cancer_gene_types",
    "cancer_driver_context","oncokb_annotated","cosmic_cgc_flag",
    "cancer_gene_source_count","cancer_gene_sources","cancer_gene_match_basis","cancer_gene_context",
    "chrom","pos","ref","alt","transcript_id","consequence",
    "normal_junction_status",
    "rna_junction_reads","rna_junction_source","rna_frame_status",
    "event_confidence","event_expression","gene_expression_tpm","transcript_expression_tpm",
    "expression_evidence_status","rna_support_status","rna_evidence_completeness","rna_evidence_score",
    "driver_relevance",
    "tumor_vaf","tumor_depth","tumor_alt_count",
    "rna_vaf","rna_alt_reads","rna_depth","rna_vaf_source",
    "clonality","persistence","tumor_specificity",
    "raw_ccf","ccf_estimate","ccf_status","ccf_confidence","ccf_warning","ccf_method",
    "ccf_resolution","ccf_resolution_reason","clonality_multiplier",
    "phase_group_id","haplotype_status","phase_support_reads","phase_total_informative_reads",
    "phase_confidence","component_event_ids","combined_protein_change","redundancy_group",
    "cross_platform_variant_key","comparison_status","cross_platform_status",
    "cross_platform_confidence","cross_platform_multiplier","cross_platform_priority_cap",
    "cross_platform_review_required","source_vcf_tumor_ad","source_vcf_tumor_af",
    "wes_tumor_depth","wes_tumor_alt_count","wes_tumor_alt_vaf",
    "wgs_tumor_depth","wgs_tumor_alt_count","wgs_tumor_alt_vaf",
    "normal_depth","normal_alt_count","normal_alt_vaf",
    "other_zero_alt_probability_at_source_pileup_vaf",
    "normal_tissue_max_tpm","normal_tissue_max_tissue","critical_tissue_max_tpm",
    "critical_tissue_name","normal_hspc_tpm","normal_hspc_unit","critical_tissue_hit",
    "normal_expression_status","normal_hspc_status","reference_proteome_status",
    "normal_ligandome_status","anchor_assessment_status","safety_evidence_completeness",
    "safety_missing_layers",
    "safety_status","safety_reason","appm_mhc_i_integrity","appm_mhc_ii_integrity",
    "source_file","source_row_number","source_record_id","source_tools","source_records",
    "provenance_record_count","evidence_conflict_status",
    "splice_event_evidence_grade","orf_evidence_grade","normal_safety_grade","splice_consensus_tier",
    "event_score","source"
]

PEPTIDE_FIELDS = [
    "peptide_id","event_id","splice_event_id","transcript_hypothesis_id","orf_id","origin_peptide_id","sample_id","event_type","mutation_source","peptide_consequence",
    "evidence_scope","priority_cap","wes_confidence_tier",
    "gene","peptide","wildtype_peptide","splice_event_type","junction_ids",
    "genome_build","canonical_junction_id","source_junction_id",
    "junction_chrom","junction_start","junction_end","junction_strand","junction_donor","junction_acceptor",
    "junction_coordinate_system","junction_resolution_status","junction_resolution_reason",
    "junction_match_status","junction_match_method","junction_support_status","junction_support_conflict","junction_support_reason",
    "provided_rna_junction_reads",
    "cancer_gene_list_status","cancer_gene_symbols","cancer_gene_types",
    "cancer_driver_context","oncokb_annotated","cosmic_cgc_flag",
    "cancer_gene_source_count","cancer_gene_sources","cancer_gene_match_basis","cancer_gene_context",
    "crosses_junction","contains_novel_aa",
    "rna_junction_reads","rna_junction_source","rna_frame_status",
    "gene_expression_tpm","transcript_expression_tpm","expression_evidence_status",
    "rna_support_status","rna_evidence_completeness","rna_evidence_score",
    "rna_vaf","rna_alt_reads","rna_depth","rna_vaf_source",
    "hla_allele","mhc_class","source_tool",
    "source_file","source_row_number","source_record_id","source_tools","source_records",
    "provenance_record_count","evidence_conflict_status","generation_status",
    "binding_rank","el_rank","presentation_score","immunogenicity_score",
    "wildtype_binding_rank","self_similarity_score","normal_hla_ligand_overlap",
    "netmhcpan_mt_ic50","netmhcpan_mt_rank_ba","netmhcpan_mt_rank_el",
    "netmhcpan_wt_ic50","netmhcpan_wt_rank_ba","netmhcpan_wt_rank_el",
    "netmhcpan_ba_rank","netmhcpan_el_rank",
    "netmhcstabpan_score","netmhcstabpan_rank",
    "netchop_31d_max_score","netchop_31d_mean_score","netchop_31d_cleavage_sites","netchop_processing_status",
    "mhcflurry_affinity_percentile","mhcflurry_processing_score","mhcflurry_presentation_score",
    "mhcflurry_wt_affinity_percentile","mhcflurry_wt_processing_score","mhcflurry_wt_presentation_score",
    "binding_evidence_score","presentation_evidence_score","presentation_evidence_grade",
    "iedb_immunogenicity_score","immunogenicity_resolved",
    "prime_score","prime_rank","bigmhc_im_score","deepimmuno_score",
    "prime_wt_score","prime_wt_rank","bigmhc_im_wt_score",
    "immunogenicity_composite_score","immunogenicity_source",
    "presentation_gate_status","presentation_gate_reason","presentation_gate_multiplier",
    "appm_multiplier","appm_multiplier_reason","appm_integrity_status",
    "appm_evidence_completeness","appm_review_required","appm_action",
    "raw_ccf","ccf_estimate","ccf_status","ccf_confidence","ccf_warning","ccf_method",
    "ccf_resolution","ccf_resolution_reason","ccf_multiplier",
    "safety_tier","safety_status","safety_reason","safety_multiplier","review_required",
    "reference_proteome_exact_match","normal_ligand_tissue","mutation_anchor_only",
    "normal_tissue_max_tpm","normal_tissue_max_tissue","critical_tissue_max_tpm",
    "critical_tissue_name","normal_hspc_tpm","normal_hspc_unit",
    "normal_expression_status","normal_hspc_status","reference_proteome_status",
    "normal_ligandome_status","anchor_assessment_status","normal_junction_assessment_status",
    "safety_evidence_completeness","safety_missing_layers","safety_priority_cap",
    "mutation_positions_in_peptide","mutation_tcr_facing",
    "phase_group_id","haplotype_status","phase_support_reads","phase_total_informative_reads",
    "phase_confidence","component_event_ids","combined_protein_change","redundancy_group",
    "cross_platform_variant_key","comparison_status","cross_platform_status",
    "cross_platform_confidence","cross_platform_multiplier","cross_platform_priority_cap",
    "cross_platform_review_required","source_vcf_tumor_ad","source_vcf_tumor_af",
    "wes_tumor_depth","wes_tumor_alt_count","wes_tumor_alt_vaf",
    "wgs_tumor_depth","wgs_tumor_alt_count","wgs_tumor_alt_vaf",
    "normal_depth","normal_alt_count","normal_alt_vaf",
    "other_zero_alt_probability_at_source_pileup_vaf",
    "agretopicity_el","mt_wt_el_rank_difference","mhcflurry_mt_wt_presentation_difference",
    "prime_mt_wt_score_difference","bigmhc_mt_wt_score_difference",
    "mutant_specificity_status","mutant_specificity_gate_status","mutant_specificity_reason",
    "mutant_specificity_multiplier","mutant_specificity_priority_cap",
    "escape_status","escape_flag","escape_reason","resistance_risk","escape_action","escape_multiplier","restricting_hla_lost",
    "l3_event_confidence_score","l3_expression_score","l3_clonality_score","l3_tumor_specificity_score",
    "l3_hla_binding_score","l3_hla_presentation_score","l3_rna_support_score","l3_rna_junction_support_score",
    "l3_normal_tissue_safety_score","l3_apm_integrity_score","l3_immunogenicity_score",
    "immunology_composite_score",
    "splice_event_evidence_grade","orf_evidence_grade","normal_safety_grade",
    "independent_translation_generators","splice_consensus_tier",
    "efficacy_score","final_priority","recommended_use"
]

PRESENTATION_FIELDS = [
    "peptide_id","event_id","sample_id","peptide","hla_allele","mhc_class",
    "netmhcpan_ba_rank","netmhcpan_el_rank",
    "netmhcpan_mt_rank_ba","netmhcpan_mt_rank_el","netmhcpan_wt_rank_ba","netmhcpan_wt_rank_el",
    "netmhcstabpan_score","netmhcstabpan_rank",
    "netchop_31d_max_score","netchop_31d_mean_score","netchop_31d_cleavage_sites","netchop_processing_status",
    "mhcflurry_affinity_percentile","mhcflurry_processing_score","mhcflurry_presentation_score",
    "mhcflurry_wt_affinity_percentile","mhcflurry_wt_processing_score","mhcflurry_wt_presentation_score",
    "iedb_immunogenicity_score",
    "prime_score","prime_rank","bigmhc_im_score","deepimmuno_score",
    "prime_wt_score","prime_wt_rank","bigmhc_im_wt_score",
    "immunogenicity_composite_score","immunogenicity_source",
    "binding_evidence_score","presentation_evidence_score",
    "evidence_completeness","presentation_evidence_grade"
]

# Standard Project B intermediate layer (multi-entry A–F → unified event ranking)
EXPRESSION_EVIDENCE_FIELDS = [
    "event_id", "sample_id", "gene", "transcript_id", "event_expression",
    "gene_expression_tpm", "transcript_expression_tpm", "expression_tpm",
    "expression_evidence_status", "expression_source", "transcript_expression_source",
    "mutation_source", "peptide_consequence",
]

RNA_JUNCTION_EVIDENCE_FIELDS = [
    "evidence_id", "event_id", "peptide_id", "sample_id", "gene", "gene_pair",
    "junction_reads", "junction_source", "mutation_source", "peptide_consequence",
    "rna_alt_reads", "rna_ref_reads", "rna_depth", "rna_vaf", "rna_vaf_source",
    "rna_frame_status", "rna_support_status", "rna_evidence_completeness", "rna_evidence_score",
    "targeted_validation_status", "targeted_validation_source",
    "targeted_validation_method",
]

FUSION_EVIDENCE_FIELDS = [
    "evidence_id", "event_id", "sample_id", "bpid", "ftid", "fusion_gene",
    "breakpoint1", "breakpoint2", "fusion_type", "frame_status", "bp1_frame", "bp2_frame",
    "exon_boundary", "neo_peptide_sequence", "fusion_protein_sequence",
    "rna_junction_reads", "rna_spanning_reads", "anchor_size",
    "caller_support_frac", "caller_prob", "caller_pass", "tools_detected",
    "filter_status", "filter_reason", "source_file",
]

SAFETY_EVIDENCE_FIELDS = [
    "evidence_id", "level", "event_id", "peptide_id", "sample_id", "gene", "peptide",
    "safety_status", "safety_reason", "normal_tissue_max_tpm", "normal_hspc_tpm", "normal_hspc_unit",
    "critical_tissue_hit", "normal_hla_ligand_overlap",
]

# Multi-entry modes (see docs/INPUT_ARCHITECTURE.md)
INPUT_MODES = {
    "snv_indel": "A — annotated VCF / pVACseq + HLA + expression",
    "fusion": "B — EasyFuse / STAR-Fusion / Arriba / AGFusion + HLA + RNA support",
    "splice_junction": "C — annotated VCF + RegTools/junction + RNA + HLA",
    "sv": "D — SV/BND VCF + GTF/reference + HLA + RNA junction",
    "peptide_only": "E — peptide table/FASTA + HLA + optional evidence",
    "e2e": "F — WES/WGS/WTS BAM/FASTQ end-to-end (optional)",
    "intermediates": "Pre-built raw_events + raw_peptides passthrough",
    "pvac": "Legacy alias for snv_indel (+ optional fusion/splice pVAC outputs)",
}

STANDARD_INTERMEDIATE_PATHS = {
    "raw_events": "parsed/raw_events.tsv",
    "raw_peptides": "parsed/raw_peptides.tsv",
    "splice_junctions": "parsed/splice_junctions.tsv",
    "splice_tool_evidence": "parsed/splice_tool_evidence.long.tsv",
    "splice_peptide_provenance": "parsed/splice_peptide_provenance.long.tsv",
    "splice_conflicts": "parsed/splice_conflicts.tsv",
    "splice_consensus": "parsed/splice_consensus.tsv",
    "splice_consensus_provenance": "parsed/splice_consensus_provenance.tsv",
    "splice_consensus_conflicts": "parsed/splice_consensus_conflicts.tsv",
    "junction_aliases": "parsed/junction_aliases.tsv",
    "splice_layer_junctions": "parsed/splice/splice_junctions.tsv",
    "splice_layer_events": "parsed/splice/splice_events.tsv",
    "splice_layer_event_junction_links": "parsed/splice/splice_event_junction_links.tsv",
    "splice_layer_transcripts": "parsed/splice/splice_transcript_hypotheses.tsv",
    "splice_layer_orfs": "parsed/splice/splice_orfs.tsv",
    "splice_layer_peptide_origins": "parsed/splice/splice_peptide_origins.tsv",
    "splice_layer_consensus": "parsed/splice/splice_consensus.tsv",
    "presentation_evidence": "presentation/presentation_evidence.tsv",
    "expression_evidence": "parsed/expression_evidence.tsv",
    "rna_junction_evidence": "parsed/rna_junction_evidence.tsv",
    "fusion_evidence": "parsed/fusion_evidence.tsv",
    "ccf_2": "clonality/ccf_2.tsv",
    "ccf_lite": "clonality/ccf_lite.tsv",
    "safety_evidence": "safety/safety_evidence.tsv",
}

# v0.4.4 strict splice-junction registry/provenance schemas.
from .splice.registry import (
    JUNCTION_ENTITY_FIELDS as SPLICE_JUNCTION_FIELDS,
    SPLICE_CONFLICT_FIELDS,
    SPLICE_PEPTIDE_PROVENANCE_FIELDS,
    SPLICE_TOOL_EVIDENCE_FIELDS,
)



PEPTIDE_SAFETY_FIELDS = [
    "peptide_id","event_id","sample_id","event_type","mutation_source","peptide_consequence",
    "gene","peptide","hla_allele","mhc_class",
    "matched_normal_status","normal_alt_reads","normal_vaf","tumor_only_flag",
    "reference_proteome_exact_match","reference_match_gene","reference_match_protein","reference_match_position",
    "normal_hla_ligand_exact_match","normal_ligand_tissue","normal_ligand_hla","normal_ligand_source_protein",
    "normal_tissue_max_tpm","normal_tissue_max_tissue","critical_tissue_max_tpm",
    "normal_hspc_tpm","normal_hspc_unit","critical_tissue_hit",
    "normal_junction_seen","normal_junction_source","normal_junction_max_reads","normal_junction_tissue",
    "wildtype_peptide","mt_binding_rank","wt_binding_rank","mt_wt_fold_change",
    "mutation_position_in_peptide","mutation_anchor_only","anchor_risk_status",
    "closest_self_peptide","closest_self_gene","closest_self_similarity","closest_self_hla_binding_rank","closest_self_normal_expression_tpm",
    "safety_tier","safety_status","safety_reason","safety_multiplier","review_required",
    "normal_expression_status","normal_hspc_status","reference_proteome_status",
    "normal_ligandome_status","anchor_assessment_status","normal_junction_assessment_status",
    "safety_evidence_completeness","safety_missing_layers","safety_priority_cap"
]

EVENT_SAFETY_FIELDS = [
    "event_id","sample_id","gene","event_type","mutation_source",
    "normal_expression_status","normal_junction_status","matched_normal_status",
    "event_safety_status","event_safety_reason","normal_hspc_status",
    "reference_proteome_status","normal_ligandome_status","anchor_assessment_status",
    "normal_tissue_max_tpm","normal_tissue_max_tissue","critical_tissue_max_tpm",
    "critical_tissue_name","normal_hspc_tpm","normal_hspc_unit",
    "safety_evidence_completeness","safety_missing_layers"
]

IMMUNE_ESCAPE_SUMMARY_FIELDS = [
    "sample_id","mhc_i_escape_status","mhc_ii_escape_status","ifng_response_status",
    "cytotoxic_killing_resistance_status","hla_loh_status","lost_hla_alleles",
    "lost_hla_i_alleles","lost_hla_ii_alleles","unclassified_lost_hla_alleles",
    "b2m_biallelic_loss","jak1_biallelic_loss","jak2_biallelic_loss",
    "tap_defect","nlrc5_defect","ciita_defect",
    "overall_immune_escape_risk","mechanism_summary","evidence_completeness","interpretation"
]

PEPTIDE_ESCAPE_FIELDS = [
    "peptide_id","event_id","sample_id","peptide","hla_allele","mhc_class",
    "restricting_hla_lost","lost_hla_alleles","b2m_status","hla_class_i_global_status",
    "jak_stat_status","tap_processing_status","nlrc5_status","ciita_status",
    "escape_status","escape_reason","escape_multiplier","priority_cap"
]

# Evidence provenance fields used by standard evidence TSV writers.
EVIDENCE_PROVENANCE_FIELDS = [
    "evidence_source",
    "evidence_tool",
    "evidence_tool_version",
    "evidence_mode",
    "evidence_file",
    "evidence_status",
]


DIAGNOSTIC_FUSION_RESCUE_FIELDS = [
    "rescue_id",
    "sample_id",
    "fusion_gene",
    "fusion_gene_raw",
    "fusion_gene_normalized",
    "gene5",
    "gene3",
    "breakpoint1",
    "breakpoint2",
    "ftid",
    "fusion_type",
    "frame_status",
    "neo_peptide_sequence",
    "neo_peptide_sequence_bp",
    "fusion_protein_sequence",
    "rna_junction_reads",
    "rna_spanning_reads",
    "anchor_size",
    "star_detected",
    "fusioncatcher_detected",
    "arriba_detected",
    "tools_detected",
    "tool_count",
    "prediction_class",
    "prediction_prob",
    "easyfuse_pass_status",
    "diagnostic_whitelist_status",
    "diagnostic_relevance",
    "rescue_reason",
    "peptide_generation_status",
    "source_file",
    "notes",
]

TOOL_PROVENANCE_TOOLS = (
    "pvacseq",
    "pvacfuse",
    "netmhcpan",
    "mhcflurry",
    "netmhcstabpan",
    "prime",
    "bigmhc_im",
    "deepimmuno",
    "iedb",
    "vep",
    "lohhla",
    "spechla",
    "facets",
    "appm_lite",
    "ccf_lite",
)

APPM_LITE_FIELDS = [
    "sample_id",
    "pathway",
    "gene",
    "mutation_status",
    "mutation_consequence",
    "expression_tpm",
    "expression_status",
    "copy_number_status",
    "loh_status",
    "risk_flag",
    "risk_reason",
] + EVIDENCE_PROVENANCE_FIELDS

APPM_SUMMARY_FIELDS = [
    "sample_id",
    "mhc_i_integrity_score",
    "mhc_ii_integrity_score",
    "hla_i_loh_flag",
    "hla_loh_alleles",
    "hla_i_loh_alleles",
    "hla_ii_loh_flag",
    "hla_ii_loh_alleles",
    "hla_loh_unclassified_alleles",
    "b2m_risk",
    "tap_risk",
    "nlrc5_risk",
    "ciita_risk",
    "expression_assessment_status",
    "appm_overall_status",
] + EVIDENCE_PROVENANCE_FIELDS

NETMHCPAN_EVIDENCE_FIELDS = [
    "sample_id", "peptide", "hla_allele", "peptide_hla_key",
    "netmhcpan_ba_score", "netmhcpan_ba_rank", "netmhcpan_el_score", "netmhcpan_el_rank",
    "source_file",
] + EVIDENCE_PROVENANCE_FIELDS

MHCFLURRY_EVIDENCE_FIELDS = [
    "sample_id", "peptide", "hla_allele", "peptide_hla_key",
    "mhcflurry_affinity", "mhcflurry_affinity_percentile",
    "mhcflurry_processing_score", "mhcflurry_presentation_score", "source_file",
] + EVIDENCE_PROVENANCE_FIELDS

NETMHCSTABPAN_EVIDENCE_FIELDS = [
    "sample_id", "peptide", "hla_allele", "peptide_hla_key",
    "netmhcstabpan_score", "netmhcstabpan_rank", "source_file",
] + EVIDENCE_PROVENANCE_FIELDS

PRIME_EVIDENCE_FIELDS = [
    "sample_id", "peptide", "hla_allele", "prime_score", "prime_rank", "source_file",
] + EVIDENCE_PROVENANCE_FIELDS

BIGMHC_IM_EVIDENCE_FIELDS = [
    "sample_id", "peptide", "hla_allele", "bigmhc_im_score", "source_file",
] + EVIDENCE_PROVENANCE_FIELDS

DEEPIMMUNO_EVIDENCE_FIELDS = [
    "sample_id", "peptide", "hla_allele", "deepimmuno_score", "source_file",
] + EVIDENCE_PROVENANCE_FIELDS

IEDB_IMMUNOGENICITY_FIELDS = [
    "sample_id", "peptide", "hla_allele", "iedb_immunogenicity_score", "source_file",
] + EVIDENCE_PROVENANCE_FIELDS

PURITY_EVIDENCE_FIELDS = ["sample_id", "purity"] + EVIDENCE_PROVENANCE_FIELDS

CNV_SEGMENT_FIELDS = ["chrom", "start", "end", "total_cn"] + EVIDENCE_PROVENANCE_FIELDS

HLA_LOH_EVIDENCE_FIELDS = [
    "hla_allele",
    "loh_status",
    "call_rule",
    "call_qc",
    "lohhla_pval",
    "lohhla_unpaired_pval",
    "lohhla_pval_unique",
    "lohhla_unpaired_pval_unique",
    "lohhla_copy_number_with_baf",
    "lohhla_cn_lower",
    "lohhla_cn_upper",
    "lohhla_mismatch_sites",
    "lohhla_prop_supportive_sites",
    "lohhla_loss_allele_raw",
    "lohhla_kept_allele_raw",
    "spechla_loh_raw",
    "spechla_copyratio",
    "spechla_allele_frequency",
    "spechla_purity",
    "spechla_ploidy",
    "spechla_het_num",
    "spechla_loss_hla_raw",
    "spechla_kept_hla_raw",
] + EVIDENCE_PROVENANCE_FIELDS

CCF_LITE_FIELDS = [
    "event_id", "gene", "chrom", "pos", "tumor_vaf", "tumor_depth", "tumor_alt_count",
    "purity", "total_copy_number", "mutation_multiplicity_assumption", "ccf_estimate",
    "ccf_status", "clonality_multiplier", "ccf_confidence", "ccf_warning",
] + EVIDENCE_PROVENANCE_FIELDS

IMMUNE_ESCAPE_EVENT_FIELDS = [
    "event_id", "sample_id", "gene", "pathway", "mechanism", "gene_status",
    "loss_status", "loss_mechanism", "risk_level", "evidence",
    "resistance_risk", "peptide_action",
] + EVIDENCE_PROVENANCE_FIELDS

PEPTIDE_ESCAPE_FLAG_FIELDS = [
    "peptide_id", "event_id", "sample_id", "peptide", "hla_allele", "mhc_class",
    "restricting_hla_lost", "global_mhc_i_escape", "ifng_escape", "mhc_ii_escape",
    "escape_flag", "escape_status", "escape_risk", "resistance_risk",
    "escape_reason", "escape_action", "escape_multiplier", "priority_cap",
] + EVIDENCE_PROVENANCE_FIELDS

# Candidate source-chain confidence fields (SNV/InDel/Fusion/Splice C1-C4).
# These fields describe whether the peptide can be traced to a valid upstream
# event/transcript/ORF. They are distinct from peptide-HLA recommendation R1-R4.
SOURCE_CHAIN_FIELDS = [
    "source_chain_track",
    "source_chain_confidence_tier",
    "source_chain_confidence_label",
    "source_chain_confidence_grade",
    "source_chain_rule_version",
    "source_chain_integration_mode",
    "source_chain_orthogonal_status",
    "source_chain_orthogonal_sources",
    "source_chain_hard_failure",
    "source_chain_hard_failure_codes",
    "source_chain_reason_codes",
    "source_chain_confidence_reason_codes",
    "event_authenticity_status",
    "orthogonal_confirmation_status",
    "transcript_orf_status",
    "novel_sequence_status",
    "phasing_status",
    "normal_background_status",
    "source_chain_supported_requirements",
    "source_chain_missing_requirements",
    "source_chain_low_power_requirements",
    "source_chain_negative_requirements",
    "source_chain_conflict_requirements",
    "source_chain_not_applicable_requirements",
    "source_chain_requirement_count",
    "source_chain_supported_count",
    "source_chain_unassessed_count",
    "source_chain_low_power_count",
    "source_chain_negative_count",
    "source_chain_conflict_count",
    "source_chain_not_applicable_count",
    "source_chain_requirement_statuses",
    "source_chain_requirement_details",
]

SOURCE_CHAIN_REQUIREMENT_FIELDS = [
    "sample_id",
    "event_id",
    "peptide_id",
    "gene",
    "source_chain_track",
    "source_chain_confidence_tier",
    "requirement_name",
    "requirement_label",
    "requirement_applicability",
    "requirement_status",
    "requirement_value",
    "requirement_core",
    "fatal_if_negative",
    "fatal_if_conflict",
    "reason_code",
    "reason",
    "source_fields",
    "requirement_conflict",
    "rule_version",
]
