process GATK_FILTER_MUTECT_CALLS {
    tag "${sample_id}"
    publishDir "${params.outdir}/calling/filter", mode: 'copy'

    input:
    val sample_id
    path raw_vcf
    path reference_fasta
    path reference_fai
    path reference_dict
    val germline_resource
    val panel_of_normals

    output:
    path "${sample_id}.mutect2.filtered.vcf.gz", emit: filtered_vcf

    script:
    def germline = germline_resource ? "--germline-resource ${germline_resource}" : ""
    def pon = panel_of_normals ? "--panel-of-normals ${panel_of_normals}" : ""
    """
    gatk FilterMutectCalls \\
      -R ${reference_fasta} \\
      -V ${raw_vcf} \\
      -O ${sample_id}.mutect2.filtered.vcf.gz \\
      ${germline} ${pon}
    """
}
