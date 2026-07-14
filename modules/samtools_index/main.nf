process SAMTOOLS_INDEX {
  tag "${bam.getName()}"
  label 'medium'
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
