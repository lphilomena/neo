#!/usr/bin/env bash
# Generic Arriba fusion detection from an aligned RNA BAM.
#
# Prerequisites:
#   bash scripts/install_fusion_tools.sh
#   STAR/chimeric alignment already performed, or a BAM with chimeric supplementary reads
#
# Usage:
#   source conf/tools.env.sh
#   PATIENT_ID=S1 \
#   INPUT_BAM=/path/rna.bam \
#   bash scripts/run_arriba_sample.sh
#
# Optional overrides:
#   REF_FASTA, GTF, BLACKLIST, KNOWN_FUSIONS, PROTEIN_DOMAINS, OUTDIR
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/conf/tools.env.sh"

PATIENT_ID="${PATIENT_ID:?ERROR: set PATIENT_ID}"
INPUT_BAM="${INPUT_BAM:?ERROR: set INPUT_BAM=/path/rna.bam}"
OUT="${OUTDIR:-${ROOT}/results/arriba/${PATIENT_ID}}"
LOG="${LOG:-${ROOT}/work/run_arriba_${PATIENT_ID}.log}"
REF_FASTA="${REF_FASTA:-${NEOAG_REFERENCE_FASTA:-}}"
GTF="${GTF:-${NEOAG_EASYFUSE_REF:-}/Homo_sapiens.GRCh38.110.gtf}"

FUSION_ENV="${NEOAG_CONDA_BASE}/envs/${NEOAG_FUSION_ENV:-neoag-fusion}"
REFERENCE_ROOT="${OPEN_NEO_REFERENCE_ROOT:-${NEOAG_DATA_ROOT:-}}"
ARRIBA_FIXED_SHARE="${NEOAG_ARRIBA_SHARE:-${NEOAG_TOOLS_ROOT:+${NEOAG_TOOLS_ROOT}/data/fusion/arriba}}"
ARRIBA_DEPLOY_SHARE="${OPEN_NEO_DEPLOY_ROOT:+${OPEN_NEO_DEPLOY_ROOT}/refs/data/fusion/arriba}"
if [[ -z "${ARRIBA_ASSET_DIR:-}" && -n "${REFERENCE_ROOT}" ]]; then
  if [[ -d "${REFERENCE_ROOT}/data/fusion/arriba" ]]; then
    ARRIBA_ASSET_DIR="${REFERENCE_ROOT}/data/fusion/arriba"
  elif [[ -d "${REFERENCE_ROOT}/fusion/arriba" ]]; then
    ARRIBA_ASSET_DIR="${REFERENCE_ROOT}/fusion/arriba"
  fi
fi
ARRIBA_ASSET_DIR="${ARRIBA_ASSET_DIR:-}"
ARRIBA_SHARE="${ARRIBA_SHARE:-${FUSION_ENV}/share/arriba}"
ARRIBA_CONDA_DATA="${ARRIBA_CONDA_DATA:-${FUSION_ENV}/var/lib/arriba}"

find_arriba_resource() {
  local override="$1"
  shift
  if [[ -n "${override}" ]]; then
    printf '%s\n' "${override}"
    return 0
  fi
  local candidate
  for candidate in "$@"; do
    if [[ -n "${candidate}" && -f "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  return 1
}

BLACKLIST="$(find_arriba_resource "${BLACKLIST:-}" \
  "${ARRIBA_ASSET_DIR}/blacklist_hg38_GRCh38_v2.5.1.tsv.gz" \
  "${ARRIBA_FIXED_SHARE}/blacklist_hg38_GRCh38_v2.5.1.tsv.gz" \
  "${ARRIBA_DEPLOY_SHARE}/blacklist_hg38_GRCh38_v2.5.1.tsv.gz" \
  "${ARRIBA_CONDA_DATA}/blacklist_hg38_GRCh38_v2.5.1.tsv.gz" \
  "${ARRIBA_SHARE}/blacklist_grch38.tsv.gz" \
  "${ARRIBA_SHARE}/blacklist_grch38.tsv")" || {
    echo "ERROR: Arriba GRCh38 blacklist not found; set BLACKLIST or OPEN_NEO_REFERENCE_ROOT" >&2
    exit 1
  }
KNOWN_FUSIONS="$(find_arriba_resource "${KNOWN_FUSIONS:-}" \
  "${ARRIBA_ASSET_DIR}/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz" \
  "${ARRIBA_FIXED_SHARE}/known_fusions_grch38.tsv.gz" \
  "${ARRIBA_DEPLOY_SHARE}/known_fusions_grch38.tsv.gz" \
  "${ARRIBA_CONDA_DATA}/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz" \
  "${ARRIBA_SHARE}/known_fusions_grch38.tsv.gz" \
  "${ARRIBA_SHARE}/known_fusions_grch38.tsv")" || {
    echo "ERROR: Arriba GRCh38 known-fusions resource not found; set KNOWN_FUSIONS" >&2
    exit 1
  }
PROTEIN_DOMAINS="$(find_arriba_resource "${PROTEIN_DOMAINS:-}" \
  "${ARRIBA_ASSET_DIR}/protein_domains_hg38_GRCh38_v2.5.1.gff3" \
  "${ARRIBA_FIXED_SHARE}/protein_domains_grch38.tsv.gz" \
  "${ARRIBA_DEPLOY_SHARE}/protein_domains_grch38.tsv.gz" \
  "${ARRIBA_CONDA_DATA}/protein_domains_hg38_GRCh38_v2.5.1.gff3" \
  "${ARRIBA_SHARE}/protein_domains_grch38.tsv.gz" \
  "${ARRIBA_SHARE}/protein_domains_grch38.tsv")" || {
    echo "ERROR: Arriba GRCh38 protein-domain resource not found; set PROTEIN_DOMAINS" >&2
    exit 1
  }

mkdir -p "${OUT}" "$(dirname "${LOG}")"
exec > >(tee -a "${LOG}") 2>&1

echo "==> run_arriba_sample $(date -Is)"
echo "    patient=${PATIENT_ID}"
echo "    input_bam=${INPUT_BAM}"
echo "    out=${OUT}"
echo "    arriba_share=${ARRIBA_SHARE}"
echo "    blacklist=${BLACKLIST}"
echo "    known_fusions=${KNOWN_FUSIONS}"
echo "    protein_domains=${PROTEIN_DOMAINS}"

[[ -s "${INPUT_BAM}" ]] || { echo "ERROR: missing BAM: ${INPUT_BAM}" >&2; exit 1; }
[[ -s "${REF_FASTA}" ]] || { echo "ERROR: missing REF_FASTA (set REF_FASTA or NEOAG_REFERENCE_FASTA)" >&2; exit 1; }
[[ -s "${GTF}" ]] || { echo "ERROR: missing GTF (set GTF or NEOAG_EASYFUSE_REF)" >&2; exit 1; }
[[ -s "${BLACKLIST}" ]] || { echo "ERROR: missing Arriba blacklist: ${BLACKLIST}" >&2; exit 1; }
command -v arriba >/dev/null 2>&1 || { echo "ERROR: arriba not on PATH; run bash scripts/install_fusion_tools.sh" >&2; exit 1; }

ARRIBA_EXTRA_ARGS=(-b "${BLACKLIST}")
if [[ -s "${KNOWN_FUSIONS}" ]]; then
  ARRIBA_EXTRA_ARGS+=(-k "${KNOWN_FUSIONS}")
else
  echo "WARN: Arriba known fusions file not found; running without -k: ${KNOWN_FUSIONS}"
fi
if [[ -s "${PROTEIN_DOMAINS}" ]]; then
  ARRIBA_EXTRA_ARGS+=(-p "${PROTEIN_DOMAINS}")
else
  echo "WARN: Arriba protein domains file not found; running without -p: ${PROTEIN_DOMAINS}"
fi

arriba \
  -x "${INPUT_BAM}" \
  -a "${REF_FASTA}" \
  -g "${GTF}" \
  "${ARRIBA_EXTRA_ARGS[@]}" \
  -o "${OUT}/${PATIENT_ID}.fusions.tsv" \
  -O "${OUT}/${PATIENT_ID}.fusions.discarded.tsv"

echo "==> Arriba finished"
ls -lh "${OUT}/${PATIENT_ID}.fusions.tsv" "${OUT}/${PATIENT_ID}.fusions.discarded.tsv"
