NEOAG_RUNNER_MODE=docker \
bin/neoag-nextflow run workflows/main_fromVCF_nohla.nf \
    --run_config conf/run.mycase.toml \
    --input_bam /mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data/dsrct_data/dsrct/sunbinbin/wgs/sunbinbin_blood.align.bam \
    --sample_id sunbinbin \
    --outdir /home/wb/working/neoantigen \
    -c conf/main_full.config \
    -profile docker
