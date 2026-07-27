#!/usr/bin/env bash
export NEOAG_RUNNER_MODE=docker
export NEOAG_PROFILE=docker
exec bash /mnt/disk_c/data_transfer/users/samba_wb/indev/neo/scripts/run_pipeline.sh \
  --workflow main_all_nohla \
  --normal_bam /home/wb/wbdata/dataset/bio/neoantigen/tutorial_9183/hcc1143_N_subset50K_fixedMT_v3.bam \
  --tumor_bam /home/wb/wbdata/dataset/bio/neoantigen/tutorial_9183/hcc1143_T_subset50K_fixedMT_v3.bam \
  --hla_alleles "HLA-A*02:01,HLA-B*07:02,HLA-C*07:02" \
  --sample_id test \
  --reference_fasta /mnt/zjl-bgi-zzb/peixunban/gl/data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa \
  --tumor_sample_name HCC1143_tumor \
  --normal_sample_name HCC1143_normal \
  --outdir results/test \
  --resume
