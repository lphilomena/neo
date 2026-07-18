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
include { EVAL_IMMUNOGENICITY } from './subworkflows/immunogenicity/main.nf'

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
