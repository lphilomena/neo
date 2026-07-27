process GATK_FILTER_MUTECT_CALLS {
    tag "${sample_id}"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir}/calling/filter", mode: 'copy'

    input:
    val sample_id
    path raw_vcf
    path raw_vcf_stats
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
    # When running inside the GATK Docker container, nextflow.config env.PATH
    # replaces the container's default PATH, shadowing /gatk and Python.
    # Restore the container's essential paths first.
    if [ -d "/gatk" ]; then
      export PATH="/gatk:/opt/miniconda/envs/gatk/bin:\$PATH"
    fi

    # Resolve GATK binary: /gatk/gatk inside the Docker container, gatk on conda PATH.
    if [ -x "/gatk/gatk" ]; then
      _GATK="/gatk/gatk"
    else
      _GATK="gatk"
    fi

    # GATK/htsjdk looks for <ref>.dict by *replacing* the .fa/.fasta extension,
    # not by appending .dict.  E.g. assembly.chr.fa → assembly.chr.dict
    #
    # Generate the dict directly from the FASTA using GATK (always correct).
    # Ensure input VCF is indexed (required by GATK FilterMutectCalls)
    if [ ! -f "${raw_vcf}.tbi" ] && [ ! -f "${raw_vcf}.idx" ]; then
      echo "Indexing input VCF: ${raw_vcf}"
      \$_GATK IndexFeatureFile -I ${raw_vcf}
    fi

    ref_dict_gatk="\$(echo '${reference_fasta}' | sed -e 's/\\.fasta\$/.dict/' -e 's/\\.fa\$/.dict/')"
    if [ ! -f "\$ref_dict_gatk" ]; then
      echo "Generating sequence dictionary: \$ref_dict_gatk"
      \$_GATK CreateSequenceDictionary -R '${reference_fasta}' -O "\$ref_dict_gatk"
    fi

    \$_GATK FilterMutectCalls \\
      -R ${reference_fasta} \\
      -V ${raw_vcf} \\
      --stats ${raw_vcf_stats} \\
      -O ${sample_id}.mutect2.filtered.vcf.gz \\
      ${germline} ${pon}
    """
}
