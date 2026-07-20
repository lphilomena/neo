process SAMTOOLS_INDEX {
  tag "${bam.getName()}"
  label 'medium'
  container 'quay.io/biocontainers/samtools:1.23.1--ha83d96e_0'
  publishDir "${params.outdir}/calling/index", mode: 'copy'

  input:
    path bam

  output:
    path "${bam}.bai", emit: bai
    path "${bam}.csi", emit: csi, optional: true

  script:
  """
  samtools index '${bam}'
  """
}
