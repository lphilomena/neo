/*
 * NeoAg — Complete End-to-End Pipeline with HLA Typing
 *
 * Extends main_fromVCF.nf with automatic HLA typing via OptiType.  Unlike
 * main_fromVCF.nf which requires pre-configured inputs.hla_alleles in the TOML,
 * this workflow computes HLA alleles directly from sequencing data.
 *
 * Pipeline stages:
 *   1. OPTITYPE          — de-novo 4-digit HLA typing from BAM/FASTQ
 *   2. MERGE_HLA_CONFIG  — inject computed HLA alleles into the run TOML
 *   3. RUN_UPSTREAM      — pVACtools / NetMHCpan / MHCflurry / VEP
 *   4. UPSTREAM_TO_SCORING — scoring chain (shared with main_fromVCF.nf)
 *
 * Required parameters:
 *   --run_config     TOML config for upstream tools (hla_alleles is auto-filled)
 *   --input_bam      Tumor BAM for HLA typing + variant calling
 *   --sample_id      Sample identifier
 *
 * Optional parameters:
 *   --input_fq1 / --input_fq2   Paired FASTQ (alternative to BAM)
 *   --seq_type       'dna' (default) or 'rna'
 *   --reference_fasta           Reference genome FASTA
 *   --profile_name              Scoring profile (default from TOML)
 *   --outdir                    Output directory
 *
 * Usage:
 *   bin/neoag-nextflow run workflows/main_fromVCF_nohla.nf \
 *     --run_config conf/run.private.toml \
 *     --input_bam tumor.bam \
 *     --sample_id SAMPLE001 \
 *     -c conf/main_full.config
 *
 * With pre-computed HLA (skip OptiType):
 *   bin/neoag-nextflow run workflows/main_fromVCF_nohla.nf \
 *     --run_config conf/run.private.toml \
 *     --skip_hla_typing \
 *     -c conf/main_full.config
 *
 * Production (Docker):
 *   NEOAG_RUNNER_MODE=docker \
 *   bin/neoag-nextflow run workflows/main_fromVCF_nohla.nf \
 *     --run_config conf/run.private.toml \
 *     --input_bam tumor.bam --sample_id SAMPLE001 \
 *     -c conf/main_full.config -profile docker
 */

nextflow.enable.dsl=2

include { OPTITYPE } from '../modules/optitype/main.nf'
include { MERGE_HLA_CONFIG } from '../modules/merge_hla_config/main.nf'
include { RUN_UPSTREAM } from '../modules/run_upstream/main.nf'
include { UPSTREAM_TO_SCORING } from './subworkflows/upstream_to_scoring.nf'

// --- Defaults -----------------------------------------------------------------
params.run_config = params.run_config ?: "${launchDir}/conf/run.stub.toml"
params.sample_id = params.sample_id ?: ''
params.profile_name = params.profile_name ?: ''
params.outdir = params.outdir ?: 'results/all'
params.input_bam = params.input_bam ?: ''
params.input_fq1 = params.input_fq1 ?: ''
params.input_fq2 = params.input_fq2 ?: ''
params.seq_type = params.seq_type ?: 'dna'
params.skip_hla_typing = params.skip_hla_typing ?: false
params.strict_mode = params.strict_mode ?: false
params.entry_mode = params.entry_mode ?: 'snv_indel'
params.normal_expression = params.normal_expression ?: ''
params.normal_hla_ligands = params.normal_hla_ligands ?: ''

// --- Workflow -----------------------------------------------------------------
workflow {

  // --- Resolve config and sample metadata ----------------------------------
  if (!params.run_config) {
    error "Missing required parameter: --run_config"
  }
  run_config = file(params.run_config)

  def toml_text = run_config.text
  def _id_match = (toml_text =~ 'id\\s*=\\s*"([^"]+)"')
  def _profile_match = (toml_text =~ 'profile\\s*=\\s*"([^"]+)"')

  sample_id   = params.sample_id    ?: (_id_match.find()      ? _id_match.group(1)      : 'SAMPLE001')
  profile_name = params.profile_name ?: (_profile_match.find() ? _profile_match.group(1) : 'default')
  entry_mode  = params.entry_mode   ?: 'snv_indel'

  // === Stage 1 — HLA typing (OptiType) =====================================
  if (!params.skip_hla_typing) {
    if (!params.input_bam && !params.input_fq1) {
      error "HLA typing requires --input_bam or --input_fq1. " +
            "Use --skip_hla_typing if HLA alleles are pre-configured in TOML."
    }
    OPTITYPE(
      sample_id,
      params.input_bam ? file(params.input_bam) : [],
      params.input_fq1 ? file(params.input_fq1) : [],
      params.input_fq2 ? file(params.input_fq2) : [],
      params.seq_type,
    )
  }

  // === Stage 2 — Merge HLA alleles into run config ========================
  if (!params.skip_hla_typing) {
    MERGE_HLA_CONFIG(sample_id, run_config, OPTITYPE.out.hla_alleles)
    merged_config = MERGE_HLA_CONFIG.out.merged_config
  } else {
    merged_config = Channel.value(run_config)
  }

  // === Stage 3 — Run upstream tools ========================================
  RUN_UPSTREAM(sample_id, merged_config)

  // === Stage 4 — Scoring chain (shared sub-workflow) =======================
  UPSTREAM_TO_SCORING(
    sample_id,
    profile_name,
    RUN_UPSTREAM.out.upstream_dir,
    params.normal_expression,
    params.normal_hla_ligands,
  )
}
