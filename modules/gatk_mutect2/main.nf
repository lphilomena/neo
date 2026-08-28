process GATK_MUTECT2 {
    tag "${sample_id}"
    cpus params.threads ?: 8
    publishDir "${params.outdir}/calling/mutect2", mode: 'copy'

    input:
    val sample_id
    path tumor_bam
    path tumor_bai
    path normal_bam
    path normal_bai
    path reference_fasta
    path reference_fai
    path reference_dict
    path intervals_bed
    val tumor_sample_name
    val normal_sample_name

    output:
    path "${sample_id}.mutect2.raw.vcf.gz", emit: raw_vcf

    script:
    """
    gatk Mutect2 \\
      -R ${reference_fasta} \\
      -I ${tumor_bam} \\
      -I ${normal_bam} \\
      -normal ${normal_sample_name} \\
      -L ${intervals_bed} \\
      -O ${sample_id}.mutect2.raw.vcf.gz
    """
}
