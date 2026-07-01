#!/usr/bin/env bash
# Install VEP offline cache + FASTA (homo_sapiens GRCh38, release 105).
# Large download (~12–15 GB); may take hours. Uses wget/curl resume when available.
#
# Usage:
#   bash scripts/install_vep_cache.sh
#   tail -f work/vep_install.log
#
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "${ROOT}/work"
LOG="${ROOT}/work/vep_install.log"

CACHE_URL="https://ftp.ensembl.org/pub/release-105/variation/indexed_vep_cache/homo_sapiens_vep_105_GRCh38.tar.gz"
CACHE_NAME="homo_sapiens_vep_105_GRCh38.tar.gz"
VEP_DIR="${HOME}/.vep"
TMP_DIR="${VEP_DIR}/tmp"
TARBALL="${TMP_DIR}/${CACHE_NAME}"

CONDA_BASE="$(conda info --base 2>/dev/null || echo /root/miniconda3)"
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate neoag-vep

exec > >(tee -a "${LOG}") 2>&1

echo "==> VEP cache install started $(date -Is)"

mkdir -p "${TMP_DIR}"

if [[ -d "${VEP_DIR}/homo_sapiens" ]] && find "${VEP_DIR}/homo_sapiens" -name 'info.txt' -print -quit | grep -q .; then
  echo "==> Cache already present under ${VEP_DIR}/homo_sapiens — skipping download"
elif [[ -f "${TARBALL}.ok" ]]; then
  echo "==> Complete tarball found at ${TARBALL}"
else
  echo "==> Downloading (resume OK): ${CACHE_URL}"
  if command -v wget >/dev/null 2>&1; then
    wget -c -O "${TARBALL}" "${CACHE_URL}"
  elif command -v curl >/dev/null 2>&1; then
    curl -fL -C - -o "${TARBALL}" "${CACHE_URL}"
  else
    echo "==> wget/curl not found; falling back to vep_install (FTP, may timeout)"
    vep_install -a cf -s homo_sapiens -y GRCh38 -n
    echo "==> Done. Cache dir: ${VEP_DIR}"
    exit 0
  fi
fi

if [[ -f "${TARBALL}" ]]; then
  echo "==> Verifying ${TARBALL} ..."
  gzip -t "${TARBALL}"
  touch "${TARBALL}.ok"
  echo "==> Extracting ${TARBALL} to ${VEP_DIR} ..."
  tar -xzf "${TARBALL}" -C "${VEP_DIR}"
  rm -f "${TARBALL}" "${TARBALL}.ok"
fi

echo "==> Installing/refining with vep_install (plugins, FASTA if needed) ..."
vep_install -a cf -s homo_sapiens -y GRCh38 -n 2>&1 || true

echo "==> Done $(date -Is). Cache dir: ${VEP_DIR}"
echo "Test: bash ${ROOT}/bin/vep-neoag --offline --cache -help 2>&1 | head -5"
