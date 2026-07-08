SV_EVENT_FULL_FIELDS = [
    "sv_event_id", "event_id", "sample_id", "evidence_scope", "svtype",
    "chrom1", "pos1", "strand1", "chrom2", "pos2", "strand2", "bnd_alt",
    "cipos", "ciend", "svlen", "inserted_sequence",
    "callers", "caller_count", "record_ids",
    "tumor_sr", "tumor_pe", "tumor_alt_support", "tumor_local_depth",
    "normal_sr", "normal_pe", "normal_alt_support", "normal_local_depth",
    "sniffles_support", "sniffles_coverage", "sniffles_precise", "sniffles_rnames_count",
    "sv_vaf_like", "pon_overlap", "population_sv_overlap", "blacklist_overlap",
    "breakpoint_precision_bp",
    "breakend1_capture_status", "breakend2_capture_status",
    "breakend1_capture_distance_bp", "breakend2_capture_distance_bp",
    "capture_interpretability", "capture_filter_status", "capture_filter_reason",
    "event_confidence_tier", "event_confidence_score", "wes_confidence_tier", "priority_cap",
    "gene1", "gene2", "transcript1", "transcript2", "exon1", "exon2",
    "cds_phase1", "cds_phase2", "effect_class", "fusion_in_frame", "frameshift",
    "protein_sequence_id", "junction_aa_position", "rna_junction_reads",
    "rna_support_status", "final_sv_confidence", "reconstruction_status", "reconstruction_reason",
    "filter_status", "filter_reason",
]

SV_PROTEIN_FIELDS = [
    "protein_sequence_id", "event_id", "sample_id", "gene", "transcript_id",
    "protein_type", "protein_sequence", "wt_protein_sequence", "wt_prefix_aa",
    "novel_aa", "junction_aa_position", "novel_start_aa", "frameshift_start_aa",
    "stop_gain_position", "in_frame", "reconstruction_method", "reconstruction_confidence",
    "reconstruction_reason",
]

SV_EVENT_TO_PEPTIDE_FIELDS = [
    "peptide_id", "event_id", "sample_id", "protein_sequence_id", "peptide", "wildtype_peptide",
    "hla_allele", "mhc_class", "peptide_length", "peptide_start_aa", "peptide_end_aa",
    "crosses_junction", "contains_novel_aa", "novel_aa_count", "wildtype_match_status",
    "reference_proteome_match", "normal_hla_ligand_overlap",
]

SV_VALIDATION_DESIGN_FIELDS = [
    "event_id", "sample_id", "protein_sequence_id", "gene", "protein_type",
    "design_type", "long_peptide", "long_peptide_length", "minigene_aa",
    "junction_aa_position", "novel_start_aa", "design_reason",
]
