#!/usr/bin/env bash
# FACETS for Chenxiaoliang 3-year-ago tumor WGS (ML150006946_L01_137) vs blood normal.
#
# Usage:
#   source conf/tools.env.sh
#   bash scripts/run_facets_tumor_3yearsAgo.sh
#   FACETS_STEP=pileup|fit|export bash scripts/run_facets_tumor_3yearsAgo.sh
#
# Notes:
#   - Default cval=150 (WGS); cval=25 often fails procSample on this pileup.
#   - Pileup/output live on NAS under chenxiaoliang_data/results/.
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA_ROOT="${CHENXIAOLIANG_DATA_ROOT:-/mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data}"

export PATIENT_ID="${PATIENT_ID:-ML150006946_L01_137}"
export TUMOR_SAMPLE_ID="${TUMOR_SAMPLE_ID:-ML150006946_L01_137}"
export NORMAL_SAMPLE_ID="${NORMAL_SAMPLE_ID:-ML150006927_L01_470}"
export TUMOR_BAM="${TUMOR_BAM:-${DATA_ROOT}/data/tumor_3yearsAgo/ML150006946_L01_137.align.bam}"
export NORMAL_BAM="${NORMAL_BAM:-${DATA_ROOT}/data/blood_0427/retransfer_6927_7362/ML150006927_L01_470/ML150006927_L01_470.markdup.bam}"
export OUTDIR="${OUTDIR:-${DATA_ROOT}/results/ML150006946_L01_137_facets}"
export FACETS_CVAL_PRE="${FACETS_CVAL_PRE:-150}"
export FACETS_CVAL_PROC="${FACETS_CVAL_PROC:-150}"
export FACETS_STEP="${FACETS_STEP:-all}"
export BAM_LINK_DIR="${BAM_LINK_DIR:-${ROOT}/work/bam_links/${PATIENT_ID}}"

exec bash "${ROOT}/scripts/run_facets_chenxiaoliang.sh"
