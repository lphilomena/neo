#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export NEOAG_CONDA_BASE="${NEOAG_CONDA_BASE:-/home/na/miniforge3}"
export NEOAG_FUSION_ENV="${NEOAG_FUSION_ENV:-neoag-fusion-r36}"
export NEOAG_GATK_ENV="${NEOAG_GATK_ENV:-neoag-gatk}"
export PATH="${NEOAG_CONDA_BASE}/bin:${PATH}"
DATA_ROOT="${CHENXIAOLIANG_DATA_ROOT:-/mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data}"
NORMAL_BAM="${DATA_ROOT}/data/blood_0427/retransfer_6927_7362/ML150006927_L01_470/ML150006927_L01_470.markdup.bam"
COMMON_SNP_VCF="${ROOT}/data/facets/reference/1000G_phase1.common_af0.05.hg38.biallelic.vcf.gz"

run_one() {
  local patient_id="$1"
  local tumor_bam="$2"
  local outdir="$3"
  echo "[$(date -Is)] start ${patient_id} common_snp"
  PATIENT_ID="${patient_id}" \
  TUMOR_BAM="${tumor_bam}" \
  NORMAL_BAM="${NORMAL_BAM}" \
  FACETS_MODE=common_snp \
  COMMON_SNP_VCF="${COMMON_SNP_VCF}" \
  OUTDIR="${outdir}" \
  FACETS_TARGET_ROWS="${FACETS_TARGET_ROWS:-1000000}" \
  FACETS_CVAL_PRE="${FACETS_CVAL_PRE:-25}" \
  FACETS_CVAL_PROC="${FACETS_CVAL_PROC:-25}" \
  FACETS_MIN_NHET="${FACETS_MIN_NHET:-5}" \
  bash "${ROOT}/scripts/run_facets_sample.sh"
  echo "[$(date -Is)] done ${patient_id} common_snp"
}

run_one \
  ML150006946_L01_137 \
  "${DATA_ROOT}/data/tumor_3yearsAgo/ML150006946_L01_137.align.bam" \
  "${ROOT}/results/chenxiaoliang_facets/ML150006946_L01_137_common_snp_af005_downsample1m"

run_one \
  M1ML150017383 \
  "${DATA_ROOT}/data/liver_0520_WGS_shortReads/seq_liver_26052/M1ML150017383_L01_438.align.bam" \
  "${ROOT}/results/chenxiaoliang_facets/M1ML150017383_common_snp_af005_downsample1m"
