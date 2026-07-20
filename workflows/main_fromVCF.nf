/*
 * NeoAg — End-to-End Full Workflow (Nextflow equivalent of `neoag-v03 run-full`)
 *
 * Mirrors the data flow and program settings of `neoag-v03 run-full` (used by
 * tests/run_demo/run_demo.sh), orchestrating the pipeline as discrete Nextflow
 * processes for resumability, caching, and HPC/cloud portability.
 *
 *   Phase 1. RUN_UPSTREAM  — pVACtools / NetMHCpan / MHCflurry / VEP
 *   Phase 2. UPSTREAM_TO_SCORING — scoring chain (parse → evidence → scoring → report)
 *
 * Metadata (sample_id, profile_name) is auto-extracted from the TOML
 * [sample] section, exactly as run-full does.  CLI overrides:
 *   --sample_id     Override TOML sample.id
 *   --profile_name  Override TOML sample.profile
 *
 * Usage (equivalent to tests/run_demo/run_demo.sh):
 *   bin/neoag-nextflow run workflows/main_fromVCF.nf \
 *     --run_config tests/run_demo/run.sliding.private.toml \
 *     --outdir results/SAMPLE001_full \
 *     -c conf/main_full.config
 *
 * Usage (production with Docker):
 *   NEOAG_RUNNER_MODE=docker \
 *   bin/neoag-nextflow run workflows/main_fromVCF.nf \
 *     --run_config conf/run.private.toml \
 *     -c conf/main_full.config -profile docker
 */

nextflow.enable.dsl=2

include { RUN_UPSTREAM } from '../modules/run_upstream/main.nf'
include { UPSTREAM_TO_SCORING } from './subworkflows/upstream_to_scoring.nf'

// --- Defaults (overridable via CLI or config) -------------------------------
params.run_config = params.run_config ?: "${launchDir}/conf/run.stub.toml"
params.sample_id = params.sample_id ?: ''
params.profile_name = params.profile_name ?: ''
params.outdir = params.outdir ?: 'results/full'
params.entry_mode = params.entry_mode ?: 'snv_indel'
params.reference_fasta = params.reference_fasta ?: ''
params.normal_expression = params.normal_expression ?: ''
params.normal_hla_ligands = params.normal_hla_ligands ?: ''
params.immunogenicity_stub = params.immunogenicity_stub ?: false
params.strict_mode = params.strict_mode ?: false

// --- Workflow ---------------------------------------------------------------
workflow {

  // --- Resolve run config ---------------------------------------------------
  if (!params.run_config) {
    error "Missing required parameter: --run_config. Provide a TOML run configuration file."
  }
  run_config = file(params.run_config)

  // --- Extract metadata from TOML [sample] section --------------------------
  // run-full reads sample.id and sample.profile from the TOML (cli.py:804-805).
  // We replicate that here so the Nextflow invocation matches: just --run_config.
  // CLI params --sample_id / --profile_name take precedence when set explicitly.
  def toml_text = run_config.text

  def _id_match = (toml_text =~ 'id\\s*=\\s*"([^"]+)"')
  def _profile_match = (toml_text =~ 'profile\\s*=\\s*"([^"]+)"')
  def _entry_match = (toml_text =~ 'entry_mode\\s*=\\s*"([^"]+)"')

  sample_id   = params.sample_id    ?: (_id_match.find()      ? _id_match.group(1)      : 'SAMPLE001')
  profile_name = params.profile_name ?: (_profile_match.find() ? _profile_match.group(1) : 'default')
  entry_mode  = params.entry_mode   ?: (_entry_match.find()   ? _entry_match.group(1)   : 'snv_indel')

  // --- Reference resolution ------------------------------------------------
  ref_fasta = params.reference_fasta
    ? file(params.reference_fasta)
    : file(System.getenv('NEOAG_REFERENCE_FASTA') ?: '/dev/null')
  if (!ref_fasta.exists()) {
    error "Missing reference FASTA. Set --reference_fasta or NEOAG_REFERENCE_FASTA."
  }

  // --- Strict mode gate ----------------------------------------------------
  if (params.strict_mode && params.run_upstream_stub) {
    error "Strict mode enabled; --run_upstream_stub is not allowed in production."
  }

  // === Phase 1 — Run upstream tools =========================================
  RUN_UPSTREAM(sample_id, run_config)

  // === Phase 2 — Scoring chain (shared sub-workflow) =========================
  UPSTREAM_TO_SCORING(
    sample_id,
    profile_name,
    RUN_UPSTREAM.out.upstream_dir,
    params.normal_expression,
    params.normal_hla_ligands,
  )
}
