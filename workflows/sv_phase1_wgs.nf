nextflow.enable.dsl=2

include { SV_MANTA } from '../modules/sv_manta/main.nf'
include { SV_SVABA } from '../modules/sv_svaba/main.nf'
include { SV_GRIDSS } from '../modules/sv_gridss/main.nf'
include { SV_NORMALIZE_MERGE } from '../modules/sv_normalize_merge/main.nf'
include { SV_BUILD_RAW } from '../modules/sv_build_raw/main.nf'
include { NEOAG_SV_SCORE_V03 } from './sv_score_v03.nf'

/*
 * WGS SV Phase 1 workflow (+ optional v0.3 scoring).
 *
 * Required params:
 *   --sample_id P001
 *   --tumor_bam tumor.bam --tumor_bai tumor.bam.bai
 *   --normal_bam normal.bam --normal_bai normal.bam.bai
 *   --reference_fasta GRCh38.fa
 *   --gencode_gtf gencode.gtf
 *   --hla 'HLA-A*02:01,HLA-B*07:02,...'
 *   --outdir results/P001_sv_phase1
 *
 * Scoring (default on): NetMHCpan + MHCflurry + score_v03
 *   --run_scoring false   # adapter-only mode
 *   --binding_stub true   # fixture binding for CI
 */

params.run_scoring = params.run_scoring != false
params.binding_stub = params.binding_stub == true
params.profile_name = params.profile_name ?: 'sv_wgs_phase1'
params.vep_appm = params.vep_appm ?: "${projectDir}/../assets/empty_vep.tsv"
params.hla_loh = params.hla_loh ?: "${projectDir}/../assets/empty_hla_loh.tsv"
params.purity = params.purity ?: "${projectDir}/../assets/empty_purity.tsv"
params.cnv = params.cnv ?: "${projectDir}/../assets/empty_cnv.tsv"
params.normal_expression = params.normal_expression ?: "${projectDir}/../resources/normal_expression.example.tsv"
params.normal_hla_ligands = params.normal_hla_ligands ?: "${projectDir}/../resources/normal_hla_ligands.example.tsv"

workflow {
    sample_ch = Channel.value(params.sample_id ?: 'SAMPLE001')
    tumor_bam_ch = Channel.value(file(params.tumor_bam))
    tumor_bai_ch = Channel.value(file(params.tumor_bai ?: params.tumor_bam + '.bai'))
    normal_bam_ch = Channel.value(file(params.normal_bam))
    normal_bai_ch = Channel.value(file(params.normal_bai ?: params.normal_bam + '.bai'))
    ref_ch = Channel.value(file(params.reference_fasta))
    gtf_ch = Channel.value(file(params.gencode_gtf))

    SV_MANTA(sample_ch, tumor_bam_ch, tumor_bai_ch, normal_bam_ch, normal_bai_ch, ref_ch)
    SV_SVABA(sample_ch, tumor_bam_ch, tumor_bai_ch, normal_bam_ch, normal_bai_ch, ref_ch)
    SV_GRIDSS(sample_ch, tumor_bam_ch, tumor_bai_ch, normal_bam_ch, normal_bai_ch, ref_ch)

    SV_NORMALIZE_MERGE(sample_ch, SV_MANTA.out.somatic_vcf, SV_SVABA.out.somatic_vcf, SV_GRIDSS.out.vcf)
    SV_BUILD_RAW(sample_ch, SV_NORMALIZE_MERGE.out.sv_list, ref_ch, gtf_ch)

    if (params.run_scoring) {
        NEOAG_SV_SCORE_V03(
            params.sample_id,
            params.profile_name,
            SV_BUILD_RAW.out.raw_events,
            SV_BUILD_RAW.out.raw_peptides,
            file(params.vep_appm),
            params.expression_tsv ? file(params.expression_tsv) : file("${projectDir}/../assets/empty_expression.tsv"),
            file(params.hla_loh),
            file(params.purity),
            file(params.cnv),
            file(params.normal_expression),
            file(params.normal_hla_ligands),
            params.binding_stub,
        )
    }
}
