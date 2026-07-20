/*
 * EVAL_IMMUNOGENICITY — Shared immunogenicity evaluation sub-workflow.
 *
 * Composes the core scoring chain used by both SNV/Indel and SV pipelines:
 *   presentation → APPM → CCF → safety → immune escape → scoring → report
 *
 * Included by:
 *   - workflows/neoag_v03_rc.nf  (SNV/Indel scoring chain)
 *   - workflows/sv_score_v03.nf   (SV scoring chain)
 */

include { BUILD_PRESENTATION } from '../../../modules/build_presentation/main.nf'
include { APPM_2 } from '../../../modules/appm_2/main.nf'
include { CCF_2 } from '../../../modules/ccf_2/main.nf'
include { PEPTIDE_SAFETY } from '../../../modules/peptide_safety/main.nf'
include { IMMUNE_ESCAPE } from '../../../modules/immune_escape/main.nf'
include { SCORE_V041 } from '../../../modules/score_v041/main.nf'
include { VALIDATION_PLAN_V03 } from '../../../modules/validation_plan/main.nf'
include { REPORT_V041 } from '../../../modules/report_v041/main.nf'
include { WORKFLOW_PROVENANCE } from '../../../modules/workflow_provenance/main.nf'

workflow EVAL_IMMUNOGENICITY {
    take:
        sample_id
        profile_name
        raw_events
        raw_peptides
        netmhcpan_evidence
        mhcflurry_evidence
        vep_appm_file
        expression_file
        hla_loh_file
        purity_file
        cnv_file
        normal_expression_file
        normal_hla_ligands_file

    main:
        BUILD_PRESENTATION(
            raw_peptides,
            netmhcpan_evidence,
            mhcflurry_evidence,
            profile_name
        )

        APPM_2(
            sample_id,
            profile_name,
            vep_appm_file,
            expression_file,
            hla_loh_file,
            cnv_file,
            raw_peptides,
            purity_file
        )

        CCF_2(
            raw_events,
            profile_name,
            purity_file,
            cnv_file
        )

        PEPTIDE_SAFETY(
            raw_events,
            raw_peptides,
            profile_name,
            normal_expression_file,
            normal_hla_ligands_file
        )

        IMMUNE_ESCAPE(
            sample_id,
            raw_peptides,
            profile_name,
            vep_appm_file,
            cnv_file,
            expression_file,
            hla_loh_file,
            APPM_2.out.appm_gene_status,
            APPM_2.out.appm_pathway_status,
            CCF_2.out.ccf_2
        )

        SCORE_V041(
            raw_events,
            raw_peptides,
            BUILD_PRESENTATION.out.presentation_evidence,
            APPM_2.out.appm_summary,
            CCF_2.out.ccf_2,
            normal_expression_file,
            normal_hla_ligands_file,
            PEPTIDE_SAFETY.out.peptide_safety,
            IMMUNE_ESCAPE.out.peptide_escape_flags,
            APPM_2.out.appm_peptide_modifiers,
            profile_name
        )

        VALIDATION_PLAN_V03(SCORE_V041.out.ranked_peptides)

        REPORT_V041(
            SCORE_V041.out.ranked_events,
            SCORE_V041.out.ranked_peptides,
            APPM_2.out.appm_summary,
            VALIDATION_PLAN_V03.out.validation_plan,
            APPM_2.out.appm_gene_status,
            APPM_2.out.appm_module_scores,
            APPM_2.out.appm_submodule_scores,
            APPM_2.out.appm_conflicts,
            APPM_2.out.appm_peptide_modifiers,
            IMMUNE_ESCAPE.out.immune_escape_summary,
            IMMUNE_ESCAPE.out.peptide_escape_flags,
            PEPTIDE_SAFETY.out.peptide_safety,
            CCF_2.out.ccf_2,
            profile_name
        )

        version_files = APPM_2.out.versions.mix(CCF_2.out.versions, PEPTIDE_SAFETY.out.versions, IMMUNE_ESCAPE.out.versions, SCORE_V041.out.versions, REPORT_V041.out.versions).collect()
        WORKFLOW_PROVENANCE(version_files)

    emit:
        ranked_events = SCORE_V041.out.ranked_events
        ranked_peptides = SCORE_V041.out.ranked_peptides
        report = REPORT_V041.out.report
        report_v041 = REPORT_V041.out.report
        workflow_provenance = WORKFLOW_PROVENANCE.out.workflow_provenance
}
