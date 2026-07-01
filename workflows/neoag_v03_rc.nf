include { PARSE_PVAC } from '../modules/parse_pvac/main.nf'
include { PARSE_NETMHCPAN } from '../modules/parse_netmhcpan/main.nf'
include { PARSE_MHCFLURRY } from '../modules/parse_mhcflurry/main.nf'
include { BUILD_PRESENTATION } from '../modules/build_presentation/main.nf'
include { APPM_2 } from '../modules/appm_2/main.nf'
include { CCF_2 } from '../modules/ccf_2/main.nf'
include { PEPTIDE_SAFETY } from '../modules/peptide_safety/main.nf'
include { IMMUNE_ESCAPE } from '../modules/immune_escape/main.nf'
include { SCORE_V041 } from '../modules/score_v041/main.nf'
include { VALIDATION_PLAN_V03 } from '../modules/validation_plan/main.nf'
include { REPORT_V041 } from '../modules/report_v041/main.nf'
include { WORKFLOW_PROVENANCE } from '../modules/workflow_provenance/main.nf'

workflow NEOAG_V03_RC {
  take:
    sample_id
    profile_name
    pvac_files
    netmhcpan_file
    mhcflurry_file
    vep_appm_file
    expression_file
    hla_loh_file
    purity_file
    cnv_file
    normal_expression_file
    normal_hla_ligands_file

  main:
    PARSE_PVAC(sample_id, profile_name, pvac_files)
    PARSE_NETMHCPAN(sample_id, netmhcpan_file)
    PARSE_MHCFLURRY(sample_id, mhcflurry_file)

    BUILD_PRESENTATION(
      PARSE_PVAC.out.raw_peptides,
      PARSE_NETMHCPAN.out.netmhcpan_evidence,
      PARSE_MHCFLURRY.out.mhcflurry_evidence,
      profile_name
    )

    APPM_2(
      sample_id,
      profile_name,
      vep_appm_file,
      expression_file,
      hla_loh_file,
      cnv_file,
      PARSE_PVAC.out.raw_peptides,
      purity_file
    )

    CCF_2(
      PARSE_PVAC.out.raw_events,
      profile_name,
      purity_file,
      cnv_file
    )

    PEPTIDE_SAFETY(
      PARSE_PVAC.out.raw_events,
      PARSE_PVAC.out.raw_peptides,
      profile_name,
      normal_expression_file,
      normal_hla_ligands_file
    )

    IMMUNE_ESCAPE(
      sample_id,
      PARSE_PVAC.out.raw_peptides,
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
      PARSE_PVAC.out.raw_events,
      PARSE_PVAC.out.raw_peptides,
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
