process GATK_FILTER_MUTECT_CALLS {
    tag "${sample_id}"
    container 'broadinstitute/gatk:4.6.2.0'
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
    # GATK/htsjdk looks for <ref>.dict by *replacing* the .fa/.fasta extension,
    # not by appending .dict.  E.g. assembly.chr.fa → assembly.chr.dict
    #
    # Generate the dict directly from the FASTA using GATK (always correct).
    ref_dict_gatk="\$(echo '${reference_fasta}' | sed -e 's/\\.fasta\$/.dict/' -e 's/\\.fa\$/.dict/')"
    if [ ! -f "\$ref_dict_gatk" ]; then
      echo "Generating sequence dictionary: \$ref_dict_gatk"
      gatk CreateSequenceDictionary -R '${reference_fasta}' -O "\$ref_dict_gatk"
    fi

    gatk FilterMutectCalls \\
      -R ${reference_fasta} \\
      -V ${raw_vcf} \\
      -O ${sample_id}.mutect2.filtered.vcf.gz \\
      ${germline} ${pon}
    """
}
