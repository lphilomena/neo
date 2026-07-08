nextflow.enable.dsl=2

include { SV_BUILD_RAW } from '../modules/sv_build_raw/main.nf'
include { NEOAG_SV_SCORE_V03 } from './sv_score_v03.nf'

/*
 * WES SV Phase 1.5 end-to-end on pre-built SV VCFs (no BAM callers).
 *
 * Uses exome-capture-limited confidence tiers (WES_Tier1/2) and the
 * sv_wes_phase1_5 scoring profile. RNA junction evidence is strongly
 * recommended for WES samples.
 *
 * Smoke test with project fixtures:
 *
 *   nextflow run workflows/sv_phase1_5_wes.nf -c conf/sv_wes_demo.config
 */

params.sample_id = params.sample_id ?: 'SVMINI'
params.profile_name = params.profile_name ?: 'sv_wes_phase1_5'
params.outdir = params.outdir ?: 'results/SVMINI_sv_wes_e2e'
params.strict_mode = params.strict_mode ?: false
params.binding_stub = params.binding_stub != false
params.wes_mode = params.wes_mode != false

params.reference_fasta = params.reference_fasta ?: "${projectDir}/../data/fixtures_sv/mini_ref.fa"
params.gencode_gtf = params.gencode_gtf ?: "${projectDir}/../data/fixtures_sv/mini.gtf"
params.hla = params.hla ?: "${projectDir}/../data/fixtures_sv/hla.txt"
params.sv_vcf = params.sv_vcf ?: "${projectDir}/../data/fixtures_sv/mini_sv.vcf"
params.sv_callers = params.sv_callers ?: ['GRIDSS2']
params.expression_tsv = params.expression_tsv ?: "${projectDir}/../data/fixtures_sv/expression.tsv"
params.rna_junction_tsv = params.rna_junction_tsv ?: "${projectDir}/../data/fixtures_sv/rna_junctions.tsv"
params.normal_expression = params.normal_expression ?: "${projectDir}/../data/fixtures_sv/normal_expression.tsv"
params.normal_hla_ligands = params.normal_hla_ligands ?: "${projectDir}/../data/fixtures_sv/normal_hla_ligands.tsv"
params.vep_appm = params.vep_appm ?: "${projectDir}/../assets/empty_vep.tsv"
params.hla_loh = params.hla_loh ?: "${projectDir}/../assets/empty_hla_loh.tsv"
params.purity = params.purity ?: "${projectDir}/../assets/empty_purity.tsv"
params.cnv = params.cnv ?: "${projectDir}/../assets/empty_cnv.tsv"

workflow {
    if (params.strict_mode && params.binding_stub) {
        error "Strict production mode forbids binding_stub; provide real binding predictor outputs/tools or set --strict_mode false."
    }
    sample_ch = Channel.value(params.sample_id)
    ref_ch = Channel.value(file(params.reference_fasta))
    gtf_ch = Channel.value(file(params.gencode_gtf))

    sv_list = Channel.fromPath(params.sv_vcf).collect()

    SV_BUILD_RAW(sample_ch, sv_list, ref_ch, gtf_ch)

    NEOAG_SV_SCORE_V03(
        params.sample_id,
        params.profile_name,
        SV_BUILD_RAW.out.raw_events,
        SV_BUILD_RAW.out.raw_peptides,
        file(params.vep_appm),
        file(params.expression_tsv),
        file(params.hla_loh),
        file(params.purity),
        file(params.cnv),
        file(params.normal_expression),
        file(params.normal_hla_ligands),
        params.binding_stub,
    )
}
