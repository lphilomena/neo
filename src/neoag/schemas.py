# Layer 1: mutation_source | Layer 2: peptide_consequence | Layer 3: scored in ranked_peptides (l3_*)
EVENT_FIELDS = [
    "event_id","sample_id","disease_profile","event_type","mutation_source","peptide_consequence",
    "evidence_scope","priority_cap","wes_confidence_tier",
    "gene","event_name",
    "cancer_gene_list_status","cancer_gene_symbols","cancer_gene_types",
    "cancer_driver_context","oncokb_annotated","cosmic_cgc_flag",
    "cancer_gene_source_count","cancer_gene_sources","cancer_gene_match_basis","cancer_gene_context",
    "chrom","pos","ref","alt","transcript_id","consequence",
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
    "event_score","source"
]

PEPTIDE_FIELDS = [
    "peptide_id","event_id","sample_id","event_type","mutation_source","peptide_consequence",
    "evidence_scope","priority_cap","wes_confidence_tier",
    "gene","peptide","wildtype_peptide",
    "cancer_gene_list_status","cancer_gene_symbols","cancer_gene_types",
    "cancer_driver_context","oncokb_annotated","cosmic_cgc_flag",
    "cancer_gene_source_count","cancer_gene_sources","cancer_gene_match_basis","cancer_gene_context",
    "crosses_junction","contains_novel_aa","rna_junction_reads","rna_junction_source","rna_frame_status",
    "gene_expression_tpm","transcript_expression_tpm","expression_evidence_status",
    "rna_support_status","rna_evidence_completeness","rna_evidence_score",
    "rna_vaf","rna_alt_reads","rna_depth","rna_vaf_source",
    "hla_allele","mhc_class","source_tool",
    "binding_rank","el_rank","presentation_score","immunogenicity_score",
    "wildtype_binding_rank","self_similarity_score","normal_hla_ligand_overlap",
    "netmhcpan_mt_ic50","netmhcpan_mt_rank_ba","netmhcpan_mt_rank_el",
    "netmhcpan_wt_ic50","netmhcpan_wt_rank_ba","netmhcpan_wt_rank_el",
    "netmhcpan_ba_rank","netmhcpan_el_rank",
    "netmhcstabpan_score","netmhcstabpan_rank",
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
    "efficacy_score","final_priority","recommended_use"
]

PRESENTATION_FIELDS = [
    "peptide_id","event_id","sample_id","peptide","hla_allele","mhc_class",
    "netmhcpan_ba_rank","netmhcpan_el_rank",
    "netmhcpan_mt_rank_ba","netmhcpan_mt_rank_el","netmhcpan_wt_rank_ba","netmhcpan_wt_rank_el",
    "netmhcstabpan_score","netmhcstabpan_rank",
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
    "presentation_evidence": "presentation/presentation_evidence.tsv",
    "expression_evidence": "parsed/expression_evidence.tsv",
    "rna_junction_evidence": "parsed/rna_junction_evidence.tsv",
    "fusion_evidence": "parsed/fusion_evidence.tsv",
    "ccf_2": "clonality/ccf_2.tsv",
    "ccf_lite": "clonality/ccf_lite.tsv",
    "safety_evidence": "safety/safety_evidence.tsv",
}


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

HLA_LOH_EVIDENCE_FIELDS = ["hla_allele", "loh_status"] + EVIDENCE_PROVENANCE_FIELDS

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
