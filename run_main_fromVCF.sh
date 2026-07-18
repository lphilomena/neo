bin/neoag-nextflow run workflows/main_fromVCF.nf \
    --run_config conf/run.mycase.toml \
    --sample_id sunbinbin \
    --outdir /home/wb/working/neoantigen \
    -c conf/main_full.config
