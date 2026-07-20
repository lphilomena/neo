nextflow.enable.dsl=2

include { RUN_BINDING_PREDICTORS } from '../modules/run_binding_predictors/main.nf'
include { PARSE_NETMHCPAN } from '../modules/parse_netmhcpan/main.nf'
include { PARSE_MHCFLURRY } from '../modules/parse_mhcflurry/main.nf'
include { EVAL_IMMUNOGENICITY } from './subworkflows/immunogenicity/main.nf'

/*
 * Score SV Phase 1 raw tables with NetMHCpan + MHCflurry + score_v03.
 *
 * Required params:
 *   --sample_id SVMINI
 *   --profile_name sv_wgs_phase1
 *   --raw_events path/to/parsed/raw_events.tsv
 *   --raw_peptides path/to/parsed/raw_peptides.tsv
 *   --outdir results/SVMINI_sv_scored
 */

workflow NEOAG_SV_SCORE_V03 {
    take:
        sample_id
        profile_name
        raw_events
        raw_peptides
        vep_appm_file
        expression_file
        hla_loh_file
        purity_file
        cnv_file
        normal_expression_file
        normal_hla_ligands_file
        binding_stub

    main:
        RUN_BINDING_PREDICTORS(sample_id, raw_peptides, binding_stub)
        PARSE_NETMHCPAN(sample_id, RUN_BINDING_PREDICTORS.out.netmhcpan)
        PARSE_MHCFLURRY(sample_id, RUN_BINDING_PREDICTORS.out.mhcflurry)

        EVAL_IMMUNOGENICITY(
            sample_id,
            profile_name,
            raw_events,
            raw_peptides,
            PARSE_NETMHCPAN.out.netmhcpan_evidence,
            PARSE_MHCFLURRY.out.mhcflurry_evidence,
            vep_appm_file,
            expression_file,
            hla_loh_file,
            purity_file,
            cnv_file,
            normal_expression_file,
            normal_hla_ligands_file
        )

    emit:
        ranked_events = EVAL_IMMUNOGENICITY.out.ranked_events
        ranked_peptides = EVAL_IMMUNOGENICITY.out.ranked_peptides
        report = EVAL_IMMUNOGENICITY.out.report
        report_v041 = EVAL_IMMUNOGENICITY.out.report_v041
        workflow_provenance = EVAL_IMMUNOGENICITY.out.workflow_provenance
}
