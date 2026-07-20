process SV_SVABA {
    tag "${sample_id}"
    cpus params.threads ?: 16
    container 'quay.io/biocontainers/svaba:1.2.0--h69ac913_1'
    publishDir "${params.outdir}/sv_calls/svaba", mode: 'copy'

    input:
    val sample_id
    path tumor_bam
    path tumor_bai
    path normal_bam
    path normal_bai
    path reference_fasta

    output:
    path "${sample_id}.svaba.somatic.sv.vcf", emit: somatic_vcf
    path "${sample_id}.bps.txt.gz", optional: true, emit: breakpoints

    script:
    """
    svaba run \
      -t ${tumor_bam} \
      -n ${normal_bam} \
      -G ${reference_fasta} \
      -a ${sample_id} \
      -p ${task.cpus}
    """
}
