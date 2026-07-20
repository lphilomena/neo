/*
 * NeoAg — Complete End-to-End Pipeline + Quality Control (QC)
 *
 * Extends main_all.nf with two parallel QC subworkflows:
 *
 *   1. OPTITYPE                       — de-novo HLA typing from normal blood BAM
 *   2. GATK_MUTECT2                   — somatic variant calling
 *   3. GATK_FILTER_MUTECT_CALLS       — filter raw variants
 *   4. SNV_WRITE_RUN_CONFIG           — generate TOML config with HLA + VCF paths
 *   5. RUN_UPSTREAM                   — pVACtools / NetMHCpan / MHCflurry / VEP
 *   6. UPSTREAM_TO_SCORING            — scoring chain (shared with main_fromVCF.nf)
 *
 *   [PARALLEL] LOH_CHECK              — LOHHLA + SpecHLA HLA LOH detection
 *   [PARALLEL] PURITY_CHECK           — FACETS + PURPLE tumor purity estimation
 *
 * The QC subworkflows start once HLA typing (Stage 1) completes and run in
 * parallel with Stages 2–6.  Results are written to:
 *     <outdir>/qc/quality_control.tsv
 *
 * Required parameters:
 *   --normal_bam          Normal blood/sample BAM (for HLA typing + variant calling)
 *   --tumor_bam           Tumor BAM (for variant calling)
 *   --reference_fasta     Reference genome FASTA (default from NEOAG_REFERENCE_FASTA)
 *   --dbsnp_vcf           dbSNP/common SNP VCF (for FACETS snp-pileup, optional)
 *   --sample_id           Sample identifier
 *
 * Optional parameters:
 *   --tumor_sample_name   Tumor sample name in BAM (default: TUMOR)
 *   --normal_sample_name  Normal sample name in BAM (default: NORMAL)
 *   --profile_name        Scoring profile (default: default)
 *   --outdir              Output directory
 *   --skip_qc             Skip QC subworkflows (default: false)
 *
 * Usage:
 *   bin/neoag-nextflow run workflows/main_all_qc.nf \
 *     --normal_bam /path/to/normal_blood.bam \
 *     --tumor_bam /path/to/tumor.bam \
 *     --dbsnp_vcf /path/to/dbsnp_chr.vcf.gz \
 *     --sample_id sunbinbin \
 *     --outdir /home/na/project/working/result_sunbinbin \
 *     -c conf/main_full.config
 *
 * Production (Docker):
 *   NEOAG_RUNNER_MODE=docker \
 *   bin/neoag-nextflow run workflows/main_all_qc.nf \
 *     --normal_bam /path/to/normal_blood.bam \
 *     --tumor_bam /path/to/tumor.bam \
 *     --dbsnp_vcf /path/to/dbsnp_chr.vcf.gz \
 *     --sample_id SAMPLE001 \
 *     -c conf/main_full.config -profile docker
 */

nextflow.enable.dsl=2

include { OPTITYPE } from '../modules/optitype/main.nf'
include { GATK_MUTECT2 } from '../modules/gatk_mutect2/main.nf'
include { GATK_FILTER_MUTECT_CALLS } from '../modules/gatk_filter_mutect_calls/main.nf'
include { SNV_WRITE_RUN_CONFIG } from '../modules/snv_write_run_config/main.nf'
include { RUN_UPSTREAM } from '../modules/run_upstream/main.nf'
include { UPSTREAM_TO_SCORING } from './subworkflows/upstream_to_scoring.nf'
include { LOH_CHECK } from './subworkflows/loh_check.nf'
include { PURITY_CHECK } from './subworkflows/purity_check.nf'

// --- Defaults -----------------------------------------------------------------
params.sample_id = params.sample_id ?: 'SAMPLE001'
params.profile_name = params.profile_name ?: 'default'
params.outdir = params.outdir ?: 'results/all_qc'
params.normal_bam = params.normal_bam ?: ''
params.tumor_bam = params.tumor_bam ?: ''
params.reference_fasta = params.reference_fasta ?: ''
params.dbsnp_vcf = params.dbsnp_vcf ?: ''
params.tumor_sample_name = params.tumor_sample_name ?: 'TUMOR'
params.normal_sample_name = params.normal_sample_name ?: 'NORMAL'
params.strict_mode = params.strict_mode ?: false
params.upstream_stub = params.upstream_stub ?: false
params.normal_expression = params.normal_expression ?: ''
params.normal_hla_ligands = params.normal_hla_ligands ?: ''
params.skip_qc = params.skip_qc ?: false

// --- Workflow -----------------------------------------------------------------
workflow {

  // --- Input validation ----------------------------------------------------
  if (!params.normal_bam) { error "Missing --normal_bam" }
  if (!params.tumor_bam)  { error "Missing --tumor_bam" }

  // --- Reference resolution ------------------------------------------------
  ref_fasta = params.reference_fasta
    ? file(params.reference_fasta)
    : file(System.getenv('NEOAG_REFERENCE_FASTA') ?: '/dev/null')
  if (!ref_fasta.exists()) {
    error "Missing reference FASTA. Set --reference_fasta or NEOAG_REFERENCE_FASTA."
  }
  ref_fai  = file("${ref_fasta}.fai")
  def fasta_str = ref_fasta.toString()
  if (fasta_str.endsWith('.chr.fa'))       { fasta_str = fasta_str[0..-8] + '.dict' }
  else if (fasta_str.endsWith('.chr.fasta')) { fasta_str = fasta_str[0..-11] + '.dict' }
  else if (fasta_str.endsWith('.fa'))        { fasta_str = fasta_str[0..-4] + '.dict' }
  else if (fasta_str.endsWith('.fasta'))     { fasta_str = fasta_str[0..-7] + '.dict' }
  ref_dict = file(fasta_str)

  // --- dbSNP VCF resolution (optional — for FACETS snp-pileup) -------------
  dbsnp_vcf = params.dbsnp_vcf
    ? file(params.dbsnp_vcf)
    : file(System.getenv('NEOAG_DBSNP_VCF') ?: '/dev/null')

  // --- Resolve sample metadata ---------------------------------------------
  sample_id         = params.sample_id
  profile_name      = params.profile_name
  tumor_sample_name  = params.tumor_sample_name
  normal_sample_name = params.normal_sample_name

  normal_bam = file(params.normal_bam)
  tumor_bam  = file(params.tumor_bam)

  // === Stage 1 — HLA typing (OptiType on normal blood BAM) ==================
  OPTITYPE(
    sample_id,
    normal_bam,
    'dna',
  )

  // === Stage 2 — Somatic variant calling (GATK Mutect2) =====================
  GATK_MUTECT2(
    sample_id,
    tumor_bam,
    normal_bam,
    ref_fasta,
    ref_fai,
    ref_dict,
    '',            // intervals_bed — empty = WGS mode
    tumor_sample_name,
    normal_sample_name,
  )

  // === Stage 3 — Filter variants ============================================
  GATK_FILTER_MUTECT_CALLS(
    sample_id,
    GATK_MUTECT2.out.raw_vcf,
    ref_fasta,
    ref_fai,
    ref_dict,
    '',   // germline_resource — none for now
    '',   // panel_of_normals — none for now
  )

  // === Stage 4 — Generate run config (TOML with OptiType HLA) ==============
  hla_alleles_list = OPTITYPE.out.hla_alleles.map { it.text.trim().tokenize('\n') }

  SNV_WRITE_RUN_CONFIG(
    sample_id,
    profile_name,
    GATK_FILTER_MUTECT_CALLS.out.filtered_vcf,
    hla_alleles_list,
    tumor_sample_name,
    normal_sample_name,
    params.upstream_stub,
  )

  // === Stage 5 — Run upstream tools ========================================
  RUN_UPSTREAM(sample_id, SNV_WRITE_RUN_CONFIG.out.run_config)

  // === Stage 6 — Scoring chain (shared sub-workflow) =======================
  UPSTREAM_TO_SCORING(
    sample_id,
    profile_name,
    RUN_UPSTREAM.out.upstream_dir,
    params.normal_expression,
    params.normal_hla_ligands,
  )

  // ==========================================================================
  //  QC Subworkflows — run in parallel with the main pipeline after Stage 1
  // ==========================================================================

  if (!params.skip_qc) {
    // --- LOH_CHECK: LOHHLA + SpecHLA HLA LOH detection ---------------------
    LOH_CHECK(
      sample_id,
      tumor_bam,
      normal_bam,
      hla_alleles_list,
      ref_fasta,
    )

    // --- PURITY_CHECK: FACETS + PURPLE tumor purity estimation -------------
    PURITY_CHECK(
      sample_id,
      tumor_bam,
      normal_bam,
      ref_fasta,
      dbsnp_vcf,
    )
  }
}
