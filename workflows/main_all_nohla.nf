/*
 * NeoAg — Complete End-to-End Pipeline (Manual HLA + Variant Calling + Scoring)
 *
 * Based on main_all.nf but uses manually-configured HLA alleles instead of
 * OptiType de-novo typing.  Suitable when HLA typing is already known or
 * performed externally before running this pipeline.
 *
 * Pipeline stages:
 *   1. GATK_MUTECT2              — somatic variant calling (auto-indexes BAMs)
 *   2. GATK_FILTER_MUTECT_CALLS  — filter raw variants
 *   3. SNV_WRITE_RUN_CONFIG      — generate TOML config with manual HLA + VCF paths
 *   4. RUN_UPSTREAM              — pVACtools / NetMHCpan / MHCflurry / VEP
 *   5. NEOAG_V03_RC              — scoring chain → ranked report
 *
 * Required parameters:
 *   --normal_bam          Normal sample BAM (for variant calling)
 *   --tumor_bam           Tumor BAM (for variant calling)
 *   --hla_alleles         Comma-separated HLA alleles,
 *                         e.g. "HLA-A*02:01,HLA-B*07:02,HLA-C*07:02"
 *   --reference_fasta     Reference genome FASTA
 *                         (default from NEOAG_REFERENCE_FASTA env var)
 *   --sample_id           Sample identifier
 *
 * Optional parameters:
 *   --tumor_sample_name   Tumor sample name in BAM (default: TUMOR)
 *   --normal_sample_name  Normal sample name in BAM (default: NORMAL)
 *   --profile_name        Scoring profile (default: default)
 *   --outdir              Output directory (default: results/all_nohla)
 *   --upstream_stub       Run upstream in stub mode (default: false)
 *   --normal_expression   Normal tissue expression TSV
 *   --normal_hla_ligands  Normal HLA ligands TSV
 *
 * Usage:
 *   bin/neoag-nextflow run workflows/main_all_nohla.nf \
 *     --normal_bam /path/to/normal.bam \
 *     --tumor_bam /path/to/tumor.bam \
 *     --hla_alleles "HLA-A*02:01,HLA-B*07:02,HLA-C*07:02" \
 *     --sample_id SAMPLE001 \
 *     -c conf/main_full.config
 *
 * Usage (with env var for reference):
 *   NEOAG_REFERENCE_FASTA=/path/to/hg38.fa \
 *   bin/neoag-nextflow run workflows/main_all_nohla.nf \
 *     --normal_bam /path/to/normal.bam \
 *     --tumor_bam /path/to/tumor.bam \
 *     --hla_alleles "HLA-A*02:01,HLA-B*07:02,HLA-C*07:02" \
 *     --sample_id SAMPLE001 \
 *     -c conf/main_full.config
 */

nextflow.enable.dsl=2

include { GATK_MUTECT2 } from '../modules/gatk_mutect2/main.nf'
include { GATK_FILTER_MUTECT_CALLS } from '../modules/gatk_filter_mutect_calls/main.nf'
include { SNV_WRITE_RUN_CONFIG } from '../modules/snv_write_run_config/main.nf'
include { RUN_UPSTREAM } from '../modules/run_upstream/main.nf'
include { NEOAG_V03_RC } from './neoag_v03_rc.nf'

// --- Defaults -----------------------------------------------------------------
params.sample_id = params.sample_id ?: 'SAMPLE001'
params.profile_name = params.profile_name ?: 'default'
params.outdir = params.outdir ?: 'results/all_nohla'
params.normal_bam = params.normal_bam ?: ''
params.tumor_bam = params.tumor_bam ?: ''
params.hla_alleles = params.hla_alleles ?: ''
params.reference_fasta = params.reference_fasta ?: ''
params.tumor_sample_name = params.tumor_sample_name ?: 'TUMOR'
params.normal_sample_name = params.normal_sample_name ?: 'NORMAL'
params.strict_mode = params.strict_mode ?: false
params.upstream_stub = params.upstream_stub ?: false
params.normal_expression = params.normal_expression ?: ''
params.normal_hla_ligands = params.normal_hla_ligands ?: ''

// --- Workflow -----------------------------------------------------------------
workflow {

  // --- Input validation ----------------------------------------------------
  if (!params.normal_bam)   { error "Missing --normal_bam" }
  if (!params.tumor_bam)    { error "Missing --tumor_bam" }
  if (!params.hla_alleles)  { error "Missing --hla_alleles. Provide comma-separated HLA alleles, e.g. --hla_alleles \"HLA-A*02:01,HLA-B*07:02,HLA-C*07:02\"" }

  // --- Reference resolution ------------------------------------------------
  ref_fasta = params.reference_fasta
    ? file(params.reference_fasta)
    : file(System.getenv('NEOAG_REFERENCE_FASTA') ?: '/dev/null')
  if (!ref_fasta.exists()) {
    error "Missing reference FASTA. Set --reference_fasta or NEOAG_REFERENCE_FASTA."
  }
  ref_fai  = file("${ref_fasta}.fai")
  def fasta_str = ref_fasta.toString()
  // Derive .dict path: .chr.fa → .dict, .chr.fasta → .dict, .fa → .dict
  if (fasta_str.endsWith('.chr.fa'))       { fasta_str = fasta_str[0..-8] + '.dict' }
  else if (fasta_str.endsWith('.chr.fasta')) { fasta_str = fasta_str[0..-11] + '.dict' }
  else if (fasta_str.endsWith('.fa'))        { fasta_str = fasta_str[0..-4] + '.dict' }
  else if (fasta_str.endsWith('.fasta'))     { fasta_str = fasta_str[0..-7] + '.dict' }
  ref_dict = file(fasta_str)

  // --- Resolve sample metadata ---------------------------------------------
  sample_id         = params.sample_id
  profile_name      = params.profile_name
  tumor_sample_name  = params.tumor_sample_name
  normal_sample_name = params.normal_sample_name

  normal_bam = file(params.normal_bam)
  tumor_bam  = file(params.tumor_bam)

  // --- Parse HLA alleles from manual input ---------------------------------
  // Same format as OptiType output (list of strings), but sourced from --hla_alleles.
  hla_alleles_list = Channel.value(
    params.hla_alleles.split(',')*.trim()
  )

  // === Stage 1 — Somatic variant calling (GATK Mutect2, auto-indexes BAMs) ==
  GATK_MUTECT2(
    sample_id,
    tumor_bam,
    normal_bam,
    ref_fasta,
    ref_fai,
    ref_dict,
    '',            // intervals_bed — empty = WGS mode (no capture intervals)
    tumor_sample_name,
    normal_sample_name,
  )

  // === Stage 2 — Filter variants ============================================
  GATK_FILTER_MUTECT_CALLS(
    sample_id,
    GATK_MUTECT2.out.raw_vcf,
    ref_fasta,
    ref_fai,
    ref_dict,
    '',   // germline_resource — none for now
    '',   // panel_of_normals — none for now
  )

  // === Stage 3 — Generate run config (TOML with manual HLA) ================
  // Uses manually-provided HLA alleles instead of OptiType output.
  SNV_WRITE_RUN_CONFIG(
    sample_id,
    profile_name,
    GATK_FILTER_MUTECT_CALLS.out.filtered_vcf,
    hla_alleles_list,
    tumor_sample_name,
    normal_sample_name,
    params.upstream_stub,
  )

  // === Stage 4 — Run upstream tools ========================================
  RUN_UPSTREAM(sample_id, SNV_WRITE_RUN_CONFIG.out.run_config)
  upstream_dir = RUN_UPSTREAM.out.upstream_dir

  // === Stage 5 — Derive channels from upstream output ======================
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

  normal_expression_ch  = params.normal_expression
    ? Channel.value(file(params.normal_expression))
    : Channel.empty()
  normal_hla_ligands_ch = params.normal_hla_ligands
    ? Channel.value(file(params.normal_hla_ligands))
    : Channel.empty()

  // === Stage 6 — NeoAg scoring chain =======================================
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
}
