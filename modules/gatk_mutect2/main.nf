process GATK_MUTECT2 {
    tag "${sample_id}"
    label 'large'
    container 'broadinstitute/gatk:4.6.2.0'
    time '168.h'
    publishDir "${params.outdir}/calling/mutect2", mode: 'copy'

    input:
    val sample_id
    path tumor_bam
    path normal_bam
    path reference_fasta
    path reference_fai
    path reference_dict
    val intervals_bed           // path to BED file or '' for WGS (no intervals)
    val tumor_sample_name
    val normal_sample_name

    output:
    path "${sample_id}.mutect2.raw.vcf.gz",       emit: raw_vcf
    path "${sample_id}.mutect2.raw.vcf.gz.tbi",   emit: raw_vcf_idx
    path "${sample_id}.mutect2.raw.vcf.gz.stats", emit: raw_vcf_stats

    script:
    def intervals = (intervals_bed && intervals_bed != '') ? "-L ${intervals_bed}" : ""
    """
    # When running inside the GATK Docker container, nextflow.config env.PATH
    # replaces the container's default PATH, shadowing /gatk and Python.
    # Restore the container's essential paths first.
    if [ -d "/gatk" ]; then
      export PATH="/gatk:/opt/miniconda/envs/gatk/bin:\$PATH"
    fi

    # Index BAMs if needed (GATK requires .bai indices).
    # Also check next to the original file (symlink target) on NFS — avoids
    # re-indexing multi-hundred-GB BAMs on every run.
    for bam_label in '${tumor_bam}' '${normal_bam}'; do
      bam_bai="\${bam_label}.bai"
      bam_csi="\${bam_label}.csi"
      if [ ! -f "\$bam_bai" ] && [ ! -f "\$bam_csi" ]; then
        # Check NFS location (symlink target)
        target_bai="\$(readlink -f "\$bam_label" 2>/dev/null || readlink "\$bam_label").bai"
        if [ -f "\$target_bai" ]; then
          echo "Linking BAI: \$target_bai -> \$bam_bai"
          ln -sf "\$target_bai" "\$bam_bai"
        else
          echo "Indexing BAM: \$bam_label"
          samtools index "\$bam_label"
        fi
      fi
    done

    # Resolve GATK binary: /gatk/gatk inside the Docker container, gatk on conda PATH.
    if [ -x "/gatk/gatk" ]; then
      _GATK="/gatk/gatk"
    else
      _GATK="gatk"
    fi

    # GATK/htsjdk looks for <ref>.dict by *replacing* the .fa/.fasta extension,
    # NOT by appending .dict.  E.g. assembly.chr.fa → assembly.chr.dict
    # (htsjdk: ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence)
    #
    # Generate the dict directly from the FASTA using GATK's bundled Picard tool.
    # This always produces a dict with correct contig names matching the FASTA,
    # unlike linking a pre-existing assembly-level dict.
    ref_dict_gatk="\$(echo '${reference_fasta}' | sed -e 's/\\.fasta\$/.dict/' -e 's/\\.fa\$/.dict/')"
    if [ ! -f "\$ref_dict_gatk" ]; then
      echo "Generating sequence dictionary: \$ref_dict_gatk"
      \$_GATK CreateSequenceDictionary -R '${reference_fasta}' -O "\$ref_dict_gatk"
    fi

    \$_GATK Mutect2 \\
      -R ${reference_fasta} \\
      -I ${tumor_bam} \\
      -I ${normal_bam} \\
      -normal ${normal_sample_name} \\
      ${intervals} \\
      -O ${sample_id}.mutect2.raw.vcf.gz

    # Index the output VCF (required by downstream FilterMutectCalls / GATK tools)
    echo "Indexing output VCF..."
    \$_GATK IndexFeatureFile -I ${sample_id}.mutect2.raw.vcf.gz
    """
}
