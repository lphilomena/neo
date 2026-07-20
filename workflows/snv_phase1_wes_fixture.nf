nextflow.enable.dsl=2

include { SNV_WRITE_RUN_CONFIG } from '../modules/snv_write_run_config/main.nf'
include { RUN_UPSTREAM } from '../modules/run_upstream/main.nf'
include { NEOAG_V03_RC } from './neoag_v03_rc.nf'

/*
 * WES SNV Phase 1 fixture: pre-built somatic VCF → pVAC stub → score_v03.
 *
 *   nextflow run workflows/snv_phase1_wes_fixture.nf -c conf/snv_wes_demo.config
 */

params.sample_id = params.sample_id ?: 'SNVMINI'
params.profile_name = params.profile_name ?: 'default'
params.outdir = params.outdir ?: 'results/SNVMINI_snv_wes_e2e'
params.upstream_stub = params.upstream_stub != false

params.somatic_vcf = params.somatic_vcf ?: "${projectDir}/../data/fixtures_snv/mini_somatic.vcf"
params.tumor_sample_name = params.tumor_sample_name ?: 'SNVMINI_TUMOR'
params.normal_sample_name = params.normal_sample_name ?: 'SNVMINI_NORMAL'
params.hla_alleles = params.hla_alleles ?: ['HLA-A*02:01', 'HLA-B*07:02']
params.expression_tsv = params.expression_tsv ?: "${projectDir}/../data/fixtures/gene_expression.tsv"
params.cnv_tsv = params.cnv_tsv ?: "${projectDir}/../assets/empty_cnv.tsv"
params.normal_expression = params.normal_expression ?: "${projectDir}/../resources/normal_expression.example.tsv"
params.normal_hla_ligands = params.normal_hla_ligands ?: "${projectDir}/../resources/normal_hla_ligands.example.tsv"

workflow {
    filtered_ch = Channel.value(file(params.somatic_vcf))

    SNV_WRITE_RUN_CONFIG(
        params.sample_id,
        params.profile_name,
        filtered_ch,
        params.hla_alleles,
        params.tumor_sample_name,
        params.normal_sample_name,
        params.upstream_stub,
    )

    RUN_UPSTREAM(params.sample_id, SNV_WRITE_RUN_CONFIG.out.run_config)

    upstream_done = RUN_UPSTREAM.out.upstream_dir
    pvac_ch = upstream_done.map { d -> file("${d}/tools/pvacseq_aggregated.tsv") }.collect()
    net_ch = upstream_done.map { d -> file("${d}/tools/netmhcpan.xls") }
    mhc_ch = upstream_done.map { d -> file("${d}/tools/mhcflurry.csv") }
    vep_ch = upstream_done.map { d -> file("${d}/tools/vep_appm.tsv") }

    NEOAG_V03_RC(
        params.sample_id,
        params.profile_name,
        pvac_ch,
        net_ch,
        mhc_ch,
        vep_ch,
        file(params.expression_tsv),
        file("${projectDir}/../assets/empty_hla_loh.tsv"),
        file("${projectDir}/../assets/empty_purity.tsv"),
        file(params.cnv_tsv),
        file(params.normal_expression),
        file(params.normal_hla_ligands),
        Channel.empty(),
        Channel.empty(),
    )
}
