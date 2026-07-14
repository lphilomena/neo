bin/neoag-nextflow run workflows/main_all_nohla.nf \
    --normal_bam /mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data/dsrct_data/dsrct/sunbinbin/wgs/sunbinbin_blood.align.bam \
    --tumor_bam /mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data/dsrct_data/dsrct/sunbinbin/wgs/sunbinbin_tumor.align.bam \
    --hla_alleles "HLA-A*01:01,HLA-B*27:05,HLA-C*02:02" \
    --sample_id sunbinbin \
    --normal_sample_name sunbinbin_blood \
    --tumor_sample_name sunbinbin_tumor \
    --outdir /home/na/project/working/result_sunbinbin \
    -c conf/main_full.config \
    -resume
