process SV_NORMALIZE_MERGE {
    tag "${sample_id}"
    publishDir "${params.outdir}/sv_merged", mode: 'copy'

    input:
    val sample_id
    path manta_vcf
    path svaba_vcf
    path gridss_vcf

    output:
    path "${sample_id}.sv_inputs.list", emit: sv_list

    script:
    """
    # The Python SV adapter accepts multiple VCFs directly and performs its own
    # first-pass breakpoint clustering. This process records the exact VCF list
    # used by the downstream SV_BUILD_RAW process.
    printf "%s\n" ${manta_vcf} ${svaba_vcf} ${gridss_vcf} > ${sample_id}.sv_inputs.list
    """
}
