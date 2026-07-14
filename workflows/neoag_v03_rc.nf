/*
 * NEOAG_V03_RC — Core scoring chain sub-workflow.
 *
 * Included by main.nf, main_full.nf, and phase1 workflows.
 * Uses v2 modules (APPM_2, CCF_2) which supersede the legacy lite versions
 * (appm_lite, ccf_lite). The lite modules are deprecated and kept only for
 * backward compatibility; new workflows should use the _2 variants exclusively.
 *
 * Supports two entry modes (matching run_v03 in pipeline_v03.py):
 *   a) From pvac files  → PARSE_PVAC extracts raw_events + raw_peptides
 *   b) From raw files   → raw_events / raw_peptides used directly
 *                          (skip PARSE_PVAC, matching run-full's
 *                           build_raw_intermediates path)
 */

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
    raw_events_file             // optional: pre-built from upstream build_raw_intermediates
    raw_peptides_file           // optional: pre-built from upstream build_raw_intermediates

  main:
    // --- Entry point: pvac parsing or pre-built raw files -------------------
    // When raw_events_file + raw_peptides_file are provided (non-empty channel),
    // use them directly — matching run_v03's raw_events/raw_peptides path.
    // Otherwise parse from pvac files.
    PARSE_PVAC(sample_id, profile_name, pvac_files)

    // If caller provided raw files, prefer them over PARSE_PVAC output.
    // Use mix() + first() to select the provided file when available.
    raw_events   = raw_events_file.mix(PARSE_PVAC.out.raw_events).first()
    raw_peptides = raw_peptides_file.mix(PARSE_PVAC.out.raw_peptides).first()

    PARSE_NETMHCPAN(sample_id, netmhcpan_file)
    PARSE_MHCFLURRY(sample_id, mhcflurry_file)

    BUILD_PRESENTATION(
      raw_peptides,
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
      CCF_2.out.ccf_lite
    )

    SCORE_V041(
      raw_events,
      raw_peptides,
      BUILD_PRESENTATION.out.presentation_evidence,
      APPM_2.out.appm_summary,
      CCF_2.out.ccf_lite,
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
      CCF_2.out.ccf_lite,
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
