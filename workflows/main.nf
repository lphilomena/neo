nextflow.enable.dsl=2

include { NEOAG_V03_RC } from './neoag_v03_rc.nf'

params.profile_name = params.profile_name ?: 'default'
params.sample_id = params.sample_id ?: 'SAMPLE001'
params.outdir = params.outdir ?: 'results/neoag_v03_rc'
params.strict_mode = params.strict_mode ?: false

params.pvac_files = params.pvac_files ?: ''
params.netmhcpan = params.netmhcpan ?: "${projectDir}/../assets/empty_netmhcpan.tsv"
params.mhcflurry = params.mhcflurry ?: "${projectDir}/../assets/empty_mhcflurry.tsv"
params.vep_appm = params.vep_appm ?: "${projectDir}/../assets/empty_vep.tsv"
params.expression = params.expression ?: "${projectDir}/../assets/empty_expression.tsv"
params.hla_loh = params.hla_loh ?: "${projectDir}/../assets/empty_hla_loh.tsv"
params.purity = params.purity ?: "${projectDir}/../assets/empty_purity.tsv"
params.cnv = params.cnv ?: "${projectDir}/../assets/empty_cnv.tsv"
params.normal_expression = params.normal_expression ?: "${projectDir}/../assets/empty_normal_expression.tsv"
params.normal_hla_ligands = params.normal_hla_ligands ?: "${projectDir}/../assets/empty_normal_hla_ligands.tsv"

workflow {
  if (!params.pvac_files) {
    error "Missing --pvac_files. Provide comma-separated pVACtools-like TSV paths."
  }

  pvac_ch = Channel.fromPath(params.pvac_files.split(',').collect{ it.trim() })
  netmhcpan_ch = Channel.fromPath(params.netmhcpan)
  mhcflurry_ch = Channel.fromPath(params.mhcflurry)
  vep_ch = Channel.fromPath(params.vep_appm)
  expression_ch = Channel.fromPath(params.expression)
  hla_loh_ch = Channel.fromPath(params.hla_loh)
  purity_ch = Channel.fromPath(params.purity)
  cnv_ch = Channel.fromPath(params.cnv)
  normal_expression_ch = Channel.fromPath(params.normal_expression)
  normal_hla_ligands_ch = Channel.fromPath(params.normal_hla_ligands)

  NEOAG_V03_RC(
    params.sample_id,
    params.profile_name,
    pvac_ch.collect(),
    netmhcpan_ch,
    mhcflurry_ch,
    vep_ch,
    expression_ch,
    hla_loh_ch,
    purity_ch,
    cnv_ch,
    normal_expression_ch,
    normal_hla_ligands_ch
  )
}
