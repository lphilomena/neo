nextflow.enable.dsl=2

include { RUN_UPSTREAM } from '../modules/run_upstream/main.nf'
include { NEOAG } from './neoag_rc.nf'

params.profile_name = params.profile_name ?: 'leukemia'
params.sample_id = params.sample_id ?: 'DEMO'
params.outdir = params.outdir ?: 'results/full'
params.run_config = params.run_config ?: "${projectDir}/conf/run.stub.toml"
params.run_upstream = true

params.expression = params.expression ?: "${projectDir}/data/fixtures/gene_expression.tsv"
params.cnv = params.cnv ?: "${projectDir}/data/fixtures/cnv_segments.tsv"
params.normal_expression = params.normal_expression ?: "${projectDir}/resources/normal_expression.example.tsv"
params.normal_hla_ligands = params.normal_hla_ligands ?: "${projectDir}/resources/normal_hla_ligands.example.tsv"

workflow {
  def up = "${params.outdir}/upstream"

  if (params.run_upstream) {
    RUN_UPSTREAM(params.sample_id, file(params.run_config))
  }

  pvac_ch = Channel.fromPath([
    "${up}/tools/pvacseq_aggregated.tsv",
    "${up}/tools/pvacfuse_aggregated.tsv",
  ])

  NEOAG(
    params.sample_id,
    params.profile_name,
    pvac_ch.collect(),
    Channel.fromPath("${up}/tools/netmhcpan.xls"),
    Channel.fromPath("${up}/tools/mhcflurry.csv"),
    Channel.fromPath("${up}/tools/vep_appm.tsv"),
    Channel.fromPath(params.expression),
    Channel.fromPath("${up}/tools/hla_loh.tsv"),
    Channel.fromPath("${up}/tools/purity.tsv"),
    Channel.fromPath(params.cnv),
    Channel.fromPath(params.normal_expression),
    Channel.fromPath(params.normal_hla_ligands),
  )
}
