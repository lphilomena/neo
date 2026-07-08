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
ARRIBA_SHARE="${ARRIBA_SHARE:-${FUSION_ENV}/share/arriba}"
BLACKLIST="${BLACKLIST:-${ARRIBA_SHARE}/blacklist_grch38.tsv.gz}"
if [[ ! -f "${BLACKLIST}" ]]; then
  BLACKLIST="${ARRIBA_SHARE}/blacklist_grch38.tsv"
fi
KNOWN_FUSIONS="${KNOWN_FUSIONS:-${ARRIBA_SHARE}/known_fusions_grch38.tsv.gz}"
if [[ ! -f "${KNOWN_FUSIONS}" ]]; then
  KNOWN_FUSIONS="${ARRIBA_SHARE}/known_fusions_grch38.tsv"
fi
PROTEIN_DOMAINS="${PROTEIN_DOMAINS:-${ARRIBA_SHARE}/protein_domains_grch38.tsv.gz}"
if [[ ! -f "${PROTEIN_DOMAINS}" ]]; then
  PROTEIN_DOMAINS="${ARRIBA_SHARE}/protein_domains_grch38.tsv"
fi

mkdir -p "${OUT}" "$(dirname "${LOG}")"
exec > >(tee -a "${LOG}") 2>&1

echo "==> run_arriba_sample $(date -Is)"
echo "    patient=${PATIENT_ID}"
echo "    input_bam=${INPUT_BAM}"
echo "    out=${OUT}"

[[ -s "${INPUT_BAM}" ]] || { echo "ERROR: missing BAM: ${INPUT_BAM}" >&2; exit 1; }
[[ -s "${REF_FASTA}" ]] || { echo "ERROR: missing REF_FASTA (set REF_FASTA or NEOAG_REFERENCE_FASTA)" >&2; exit 1; }
[[ -s "${GTF}" ]] || { echo "ERROR: missing GTF (set GTF or NEOAG_EASYFUSE_REF)" >&2; exit 1; }
command -v arriba >/dev/null 2>&1 || { echo "ERROR: arriba not on PATH; run bash scripts/install_fusion_tools.sh" >&2; exit 1; }

arriba \
  -x "${INPUT_BAM}" \
  -a "${REF_FASTA}" \
  -g "${GTF}" \
  -b "${BLACKLIST}" \
  -k "${KNOWN_FUSIONS}" \
  -t "${PROTEIN_DOMAINS}" \
  -o "${OUT}/${PATIENT_ID}.fusions.tsv" \
  -d "${OUT}/${PATIENT_ID}.fusions.discarded.tsv"

echo "==> Arriba finished"
ls -lh "${OUT}/${PATIENT_ID}.fusions.tsv" "${OUT}/${PATIENT_ID}.fusions.discarded.tsv"
