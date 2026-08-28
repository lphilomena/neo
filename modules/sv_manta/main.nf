process SV_MANTA {
    tag "${sample_id}"
    cpus params.threads ?: 16
    publishDir "${params.outdir}/sv_calls/manta", mode: 'copy'

    input:
    val sample_id
    path tumor_bam
    path tumor_bai
    path normal_bam
    path normal_bai
    path reference_fasta

    output:
    path "manta/results/variants/somaticSV.vcf.gz", emit: somatic_vcf

    script:
    """
    configManta.py \
      --tumorBam ${tumor_bam} \
      --normalBam ${normal_bam} \
      --referenceFasta ${reference_fasta} \
      --runDir manta
    manta/runWorkflow.py -m local -j ${task.cpus}
    """
}
