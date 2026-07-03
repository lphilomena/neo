#!/usr/bin/env bash
set -euo pipefail

PROJECT=/home/na/project/neoantigen/neoag_event_pipeline_v03_rc
CONDA_SH=/home/na/miniforge3/etc/profile.d/conda.sh
ENV=neoag-sequenza
REF=/mnt/zjl-bgi-zzb/peixunban/gl/data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa
NORMAL=/mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data/data/blood_0427/retransfer_6927_7362/ML150006927_L01_470/ML150006927_L01_470.markdup.bam
OUTROOT=${PROJECT}/results/chenxiaoliang_sequenza
GC=${OUTROOT}/reference/Homo_sapiens.GRCh38.dna.primary_assembly.chr.gc50.wig.gz
SAMTOOLS=/home/na/miniforge3/envs/${ENV}/bin/samtools
TABIX=/home/na/miniforge3/envs/${ENV}/bin/tabix

source ${CONDA_SH}
mkdir -p ${OUTROOT}/reference ${OUTROOT}/logs

run_env() {
  mamba run -n ${ENV} "$@"
}

if [[ ! -s ${GC} ]]; then
  echo "[$(date)] Generating GC wiggle: ${GC}"
  run_env sequenza-utils gc_wiggle -f ${REF} -w 50 -o ${GC}
else
  echo "[$(date)] Reusing GC wiggle: ${GC}"
fi

run_one() {
  local label=$1
  local tumor=$2
  local outdir=${OUTROOT}/${label}
  local raw=${outdir}/${label}.seqz.gz
  local bin=${outdir}/${label}.small.seqz.gz
  mkdir -p ${outdir}/sequenza_fit
  echo "[$(date)] Start ${label}"
  echo -e "sample_id\ttumor_bam\tnormal_bam\treference\tgc_wiggle\tqlimit\tmin_depth_N\thom\thet\tgc_window\tbin_window" > ${outdir}/run_parameters.tsv
  echo -e "${label}\t${tumor}\t${NORMAL}\t${REF}\t${GC}\t20\t20\t0.9\t0.25\t50\t50" >> ${outdir}/run_parameters.tsv

  if [[ ! -s ${raw} ]]; then
    echo "[$(date)] bam2seqz ${label}"
    run_env sequenza-utils bam2seqz \
      -n ${NORMAL} \
      -t ${tumor} \
      -gc ${GC} \
      -F ${REF} \
      -S ${SAMTOOLS} \
      -T ${TABIX} \
      -q 20 \
      -N 20 \
      --hom 0.9 \
      --het 0.25 \
      -o ${raw}
  else
    echo "[$(date)] Reusing raw seqz ${raw}"
  fi

  if [[ ! -s ${bin} ]]; then
    echo "[$(date)] seqz_binning ${label}"
    run_env sequenza-utils seqz_binning \
      -s ${raw} \
      -w 50 \
      -T ${TABIX} \
      -o ${bin}
  else
    echo "[$(date)] Reusing binned seqz ${bin}"
  fi

  echo "[$(date)] R fit ${label}"
  run_env Rscript ${PROJECT}/scripts/run_sequenza_fit.R ${bin} ${outdir}/sequenza_fit ${label}
  echo "[$(date)] Done ${label}"
}

run_one ML150006946_L01_137_3yearsAgo /mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data/data/tumor_3yearsAgo/ML150006946_L01_137.align.bam
run_one M1ML150017383_currentYear /mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data/data/liver_0520_WGS_shortReads/seq_liver_26052/M1ML150017383_L01_438.align.bam

echo "[$(date)] All Sequenza runs finished"
