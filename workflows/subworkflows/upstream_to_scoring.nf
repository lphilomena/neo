/*
 * UPSTREAM_TO_SCORING — Shared tail: upstream output → scoring chain.
 *
 * Takes the RUN_UPSTREAM output directory and wires every downstream channel
 * into NEOAG_V03_RC.  This is the common suffix shared by all top-level
 * workflows (main_fromVCF, main_fromVCF_nohla, main_all_nohla, main_all).
 *
 * Included by:
 *   - workflows/main_fromVCF.nf
 *   - workflows/main_fromVCF_nohla.nf
 *   - workflows/main_all_nohla.nf
 *   - workflows/main_all.nf
 */

include { NEOAG_V03_RC } from '../../workflows/neoag_v03_rc.nf'

workflow UPSTREAM_TO_SCORING {
  take:
    sample_id
    profile_name
    upstream_dir                // path channel from RUN_UPSTREAM.out.upstream_dir
    normal_expression           // val: path string or ''
    normal_hla_ligands          // val: path string or ''

  main:
    // --- Derive input channels from upstream output --------------------------
    // pvac files are optional: when the TOML doesn't enable pvac tools,
    // run-upstream produces raw_events/raw_peptides directly (via
    // build_raw_intermediates).  Filter out non-existent pvac files so
    // PARSE_PVAC gets an empty channel and Nextflow skips it.
    pvac_ch = upstream_dir
      .map { d -> file("${d}/tools/pvacseq_aggregated.tsv") }
      .mix(upstream_dir.map { d -> file("${d}/tools/pvacfuse_aggregated.tsv") })
      .filter { it.exists() }
      .collect()

    netmhcpan_ch    = upstream_dir.map { d -> file("${d}/tools/netmhcpan.xls") }
    mhcflurry_ch    = upstream_dir.map { d -> file("${d}/tools/mhcflurry.csv") }
    vep_ch          = upstream_dir.map { d -> file("${d}/tools/vep_appm.tsv") }
    hla_loh_ch      = upstream_dir.map { d -> file("${d}/tools/hla_loh.tsv") }
    purity_ch       = upstream_dir.map { d -> file("${d}/tools/purity.tsv") }
    expression_ch   = upstream_dir.map { d -> file("${d}/tools/expression.tsv") }
    cnv_ch          = upstream_dir.map { d -> file("${d}/tools/cnv.tsv") }

    raw_events_ch   = upstream_dir.map { d -> file("${d}/parsed/raw_events.tsv") }
    raw_peptides_ch = upstream_dir.map { d -> file("${d}/parsed/raw_peptides.tsv") }

    // Normal reference files: use user-specified paths, or empty channel.
    normal_expression_ch  = normal_expression
      ? Channel.value(file(normal_expression))
      : Channel.empty()
    normal_hla_ligands_ch = normal_hla_ligands
      ? Channel.value(file(normal_hla_ligands))
      : Channel.empty()

    // --- Scoring chain -------------------------------------------------------
    NEOAG_V03_RC(
      sample_id,
      profile_name,
      pvac_ch,
      netmhcpan_ch,
      mhcflurry_ch,
      vep_ch,
      expression_ch,
      hla_loh_ch,
      purity_ch,
      cnv_ch,
      normal_expression_ch,
      normal_hla_ligands_ch,
      raw_events_ch,
      raw_peptides_ch,
    )

  emit:
    ranked_events      = NEOAG_V03_RC.out.ranked_events
    ranked_peptides    = NEOAG_V03_RC.out.ranked_peptides
    report             = NEOAG_V03_RC.out.report
    report_v041        = NEOAG_V03_RC.out.report_v041
    workflow_provenance = NEOAG_V03_RC.out.workflow_provenance
}
