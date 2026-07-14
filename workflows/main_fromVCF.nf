/*
 * NeoAg — End-to-End Full Workflow (Nextflow equivalent of `neoag-v03 run-full`)
 *
 * Mirrors the data flow and program settings of `neoag-v03 run-full` (used by
 * tests/run_demo/run_demo.sh), orchestrating the pipeline as discrete Nextflow
 * processes for resumability, caching, and HPC/cloud portability.
 *
 *   Phase 1. RUN_UPSTREAM  — pVACtools / NetMHCpan / MHCflurry / VEP
 *   Phase 2. NEOAG_V03_RC  — scoring chain (parse → evidence → scoring → report)
 *
 * Metadata (sample_id, profile_name) is auto-extracted from the TOML
 * [sample] section, exactly as run-full does.  CLi overrides:
 *   --sample_id     Override TOML sample.id
 *   --profile_name  Override TOML sample.profile
 *
 * Usage (equivalent to tests/run_demo/run_demo.sh):
 *   bin/neoag-nextflow run workflows/main_full.nf \
 *     --run_config tests/run_demo/run.sliding.private.toml \
 *     --outdir results/SAMPLE001_full \
 *     -c conf/main_full.config
 *
 * Usage (production with SLURM + Docker):
 *   NEOAG_RUNNER_MODE=docker \
 *   bin/neoag-nextflow run workflows/main_full.nf \
 *     --run_config conf/run.private.toml \
 *     -c conf/main_full.config -profile slurm
 */

nextflow.enable.dsl=2

include { RUN_UPSTREAM } from '../modules/run_upstream/main.nf'
include { NEOAG_V03_RC } from './neoag_v03_rc.nf'

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
  ref_fai  = file("${ref_fasta}.fai")
  def fasta_str = ref_fasta.toString()
  // Derive .dict path: .chr.fa → .dict, .chr.fasta → .dict, .fa → .dict
  if (fasta_str.endsWith('.chr.fa'))       { fasta_str = fasta_str[0..-8] + '.dict' }
  else if (fasta_str.endsWith('.chr.fasta')) { fasta_str = fasta_str[0..-11] + '.dict' }
  else if (fasta_str.endsWith('.fa'))        { fasta_str = fasta_str[0..-4] + '.dict' }
  else if (fasta_str.endsWith('.fasta'))     { fasta_str = fasta_str[0..-7] + '.dict' }
  ref_dict = file(fasta_str)

  // --- Strict mode gate ----------------------------------------------------
  if (params.strict_mode && params.run_upstream_stub) {
    error "Strict mode enabled; --run_upstream_stub is not allowed in production."
  }

  // === Phase 1 — Run upstream tools =========================================
  // Equivalent to: upstream = run_upstream(config, outdir / "upstream")
  // in cmd_run_full (cli.py:803)
  RUN_UPSTREAM(sample_id, run_config)
  upstream_dir = RUN_UPSTREAM.out.upstream_dir

  // === Phase 2 — Derive input channels from upstream output ==================
  // run-full reads these keys from the upstream results dict (cli.py:807-849).
  // We construct Nextflow channels that resolve after RUN_UPSTREAM completes.
  //
  // pvac files are optional: when the TOML doesn't enable pvac tools,
  // run-upstream produces raw_events/raw_peptides directly (via
  // build_raw_intermediates).  Filter out non-existent pvac files so
  // PARSE_PVAC gets an empty channel and Nextflow skips it, matching
  // run_v03's fallback path in pipeline_v03.py:74-80.
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

  // raw_events / raw_peptides from upstream build_raw_intermediates
  // (produced when pvac tools are not enabled — matching run-full fallback)
  raw_events_ch   = upstream_dir.map { d -> file("${d}/parsed/raw_events.tsv") }
  raw_peptides_ch = upstream_dir.map { d -> file("${d}/parsed/raw_peptides.tsv") }

  // Use upstream-derived files when available, fall back to params
  // (matching run-full's fallback logic in cli.py:844-845)
  normal_expression_ch  = Channel.value(file(params.normal_expression))
  normal_hla_ligands_ch = Channel.value(file(params.normal_hla_ligands))

  // === Phase 3 — NeoAg scoring chain ========================================
  // Equivalent to: run_v03(...) in cmd_run_full (cli.py:854)
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
