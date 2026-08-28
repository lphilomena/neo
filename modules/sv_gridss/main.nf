process SV_GRIDSS {
    tag "${sample_id}"
    cpus params.threads ?: 16
    publishDir "${params.outdir}/sv_calls/gridss", mode: 'copy'

    input:
    val sample_id
    path tumor_bam
    path tumor_bai
    path normal_bam
    path normal_bai
    path reference_fasta

    output:
    path "${sample_id}.gridss.vcf.gz", emit: vcf
    path "${sample_id}.gridss.assembly.bam", optional: true, emit: assembly_bam

    script:
    """
    gridss \
      --reference ${reference_fasta} \
      --output ${sample_id}.gridss.vcf.gz \
      --assembly ${sample_id}.gridss.assembly.bam \
      ${tumor_bam} ${normal_bam}
    """
}
