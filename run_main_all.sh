NEOAG_RUNNER_MODE=docker \
bin/neoag-nextflow run workflows/main_all.nf \
 --normal_bam /mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data/dsrct_data/dsrct/sunbinbin/wgs/sunbinbin_blood.align.bam \
 --tumor_bam /mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data/dsrct_data/dsrct/sunbinbin/wgs/sunbinbin_tumor.align.bam \
 --sample_id sunbinbin \
 --tumor_sample_name TUMOR \
 --normal_sample_name NORMAL \
 --outdir /home/na/project/working/result_sunbinbin \
 -c conf/main_full.config \
 -profile docker \
 "$@"

