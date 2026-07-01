nextflow.enable.dsl=2

include { GATK_MUTECT2 } from '../modules/gatk_mutect2/main.nf'
include { GATK_FILTER_MUTECT_CALLS } from '../modules/gatk_filter_mutect_calls/main.nf'
include { SNV_WRITE_RUN_CONFIG } from '../modules/snv_write_run_config/main.nf'
include { RUN_UPSTREAM } from '../modules/run_upstream/main.nf'
include { NEOAG_V03_RC } from './neoag_v03_rc.nf'

/*
 * WES SNV Phase 1: tumor-normal BAM → Mutect2 → FilterMutectCalls → pVAC → score_v03.
 *
 *   nextflow run workflows/snv_phase1_wes.nf \
 *     --sample_id P001 \
 *     --tumor_bam tumor.bam --normal_bam normal.bam \
 *     --reference_fasta GRCh38.fa --intervals_bed capture.bed \
 *     --tumor_sample_name TUMOR --normal_sample_name NORMAL \
 *     -c conf/snv_wes_demo.config
 */

params.sample_id = params.sample_id ?: 'P001'
params.profile_name = params.profile_name ?: 'default'
params.outdir = params.outdir ?: 'results/P001_snv_wes'
params.upstream_stub = params.upstream_stub == true
params.threads = params.threads ?: 8

params.reference_fasta = params.reference_fasta ?: "${projectDir}/../data/fixtures_snv/mini_ref.fa"
params.intervals_bed = params.intervals_bed ?: "${projectDir}/../data/fixtures_snv/wes_capture.bed"
params.tumor_sample_name = params.tumor_sample_name ?: 'TUMOR'
params.normal_sample_name = params.normal_sample_name ?: 'NORMAL'
params.hla_alleles = params.hla_alleles ?: ['HLA-A*02:01', 'HLA-B*07:02']
params.gnomad_vcf = params.gnomad_vcf ?: ''
params.panel_of_normals = params.panel_of_normals ?: ''
params.expression_tsv = params.expression_tsv ?: "${projectDir}/../data/fixtures/gene_expression.tsv"
params.cnv_tsv = params.cnv_tsv ?: "${projectDir}/../assets/empty_cnv.tsv"
params.normal_expression = params.normal_expression ?: "${projectDir}/../resources/normal_expression.example.tsv"
params.reference_dict = params.reference_dict ?: "${params.reference_fasta}".replaceAll(/\\.fa\$/, '.dict').replaceAll(/\\.fasta\$/, '.dict')

workflow {
    sample_ch = Channel.value(params.sample_id)
    tumor_bam_ch = Channel.value(file(params.tumor_bam))
    tumor_bai_ch = Channel.value(file(params.tumor_bai ?: params.tumor_bam + '.bai'))
    normal_bam_ch = Channel.value(file(params.normal_bam))
    normal_bai_ch = Channel.value(file(params.normal_bai ?: params.normal_bam + '.bai'))
    ref_ch = Channel.value(file(params.reference_fasta))
    ref_fai_ch = Channel.value(file("${params.reference_fasta}.fai"))
    ref_dict_ch = Channel.value(file(params.reference_dict))
    bed_ch = Channel.value(file(params.intervals_bed))

    GATK_MUTECT2(
        sample_ch,
        tumor_bam_ch,
        tumor_bai_ch,
        normal_bam_ch,
        normal_bai_ch,
        ref_ch,
        ref_fai_ch,
        ref_dict_ch,
        bed_ch,
        params.tumor_sample_name,
        params.normal_sample_name,
    )

    GATK_FILTER_MUTECT_CALLS(
        sample_ch,
        GATK_MUTECT2.out.raw_vcf,
        ref_ch,
        ref_fai_ch,
        ref_dict_ch,
        params.gnomad_vcf,
        params.panel_of_normals,
    )

    SNV_WRITE_RUN_CONFIG(
        params.sample_id,
        params.profile_name,
        GATK_FILTER_MUTECT_CALLS.out.filtered_vcf,
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
    )
}
