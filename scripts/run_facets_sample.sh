#!/usr/bin/env bash
# Generic FACETS runner for one tumor-normal DNA sample.
#
# Required:
#   PATIENT_ID=sample_id
#   TUMOR_BAM=/path/tumor.bam
#   NORMAL_BAM=/path/normal.bam
#
# Modes:
#   FACETS_MODE=omni2p5  Clean 1000G omni2.5 biallelic SNP-only pileup with stride downsample.
#   FACETS_MODE=common_snp  1000G phase1 common AF>=0.05 biallelic SNP-only pileup.
#   FACETS_MODE=dbsnp       Full dbSNP WGS pileup using snp-pileup -g/-P pseudo-depth sites.
#
# Common examples:
#   PATIENT_ID=S1 TUMOR_BAM=/data/S1.tumor.bam NORMAL_BAM=/data/S1.normal.bam bash scripts/run_facets_sample.sh
#   FACETS_MODE=dbsnp FACETS_CVAL_PRE=150 FACETS_CVAL_PROC=150 PATIENT_ID=S1 TUMOR_BAM=/data/T.bam NORMAL_BAM=/data/N.bam bash scripts/run_facets_sample.sh
#   FACETS_STEP=pileup|downsample|fit|export may be used for omni2p5; FACETS_STEP=pileup|fit|export for dbsnp.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
if [[ -d /home/na/miniforge3 ]]; then
  export NEOAG_CONDA_BASE="${NEOAG_CONDA_BASE:-/home/na/miniforge3}"
  export PATH="${NEOAG_CONDA_BASE}/bin:${PATH}"
fi
# shellcheck source=/dev/null
source "${ROOT}/conf/tools.env.sh"

PATIENT_ID="${PATIENT_ID:?ERROR: set PATIENT_ID=sample_id}"
TUMOR_BAM="${TUMOR_BAM:?ERROR: set TUMOR_BAM=/path/tumor.bam}"
NORMAL_BAM="${NORMAL_BAM:?ERROR: set NORMAL_BAM=/path/normal.bam}"
FACETS_MODE="${FACETS_MODE:-omni2p5}"
COMMON_SNP_VCF="${COMMON_SNP_VCF:-${ROOT}/data/facets/reference/1000G_phase1.common_af0.05.hg38.biallelic.vcf.gz}"

case "${FACETS_MODE}" in
  omni2p5|omni|snp-only|snponly)
    export PATIENT_ID TUMOR_BAM NORMAL_BAM
    export OUTDIR="${OUTDIR:-${ROOT}/results/facets/${PATIENT_ID}/omni2p5_snponly_downsample}"
    export LOG="${LOG:-${OUTDIR}/run.log}"
    exec bash "${ROOT}/scripts/run_facets_omni2p5_snponly_downsample.sh"
    ;;
  common_snp|common|1000g-common)
    export PATIENT_ID TUMOR_BAM NORMAL_BAM
    export FACETS_SNPSET_NAME="common_snp"
    export FACETS_SNP_VCF="${FACETS_SNP_VCF:-${COMMON_SNP_VCF}}"
    export OUTDIR="${OUTDIR:-${ROOT}/results/facets/${PATIENT_ID}/common_snp_snponly_downsample}"
    export LOG="${LOG:-${OUTDIR}/run.log}"
    exec bash "${ROOT}/scripts/run_facets_omni2p5_snponly_downsample.sh"
    ;;
  dbsnp|full-dbsnp|full)
    export PATIENT_ID TUMOR_BAM NORMAL_BAM
    export OUTDIR="${OUTDIR:-${ROOT}/results/facets/${PATIENT_ID}/dbsnp_full}"
    export LOG="${LOG:-${ROOT}/work/run_facets_${PATIENT_ID}.log}"
    exec bash "${ROOT}/scripts/run_facets_chenxiaoliang.sh"
    ;;
  *)
    echo "ERROR: FACETS_MODE must be omni2p5, common_snp, or dbsnp; got ${FACETS_MODE}" >&2
    exit 1
    ;;
esac
