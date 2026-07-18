#!/usr/bin/env bash
set -euo pipefail

CONDA_BASE=~/miniforge3  # 或你的 conda 路径
PACK_DIR=~/Downloads/conda_packs

echo "[1/5] CONDA_BASE=${CONDA_BASE}"
echo "[1/5] PACK_DIR=${PACK_DIR}"

echo "[2/5] sourcing conda.sh ..."
source "${CONDA_BASE}/etc/profile.d/conda.sh"
echo "[2/5] conda.sh loaded"

unpack_one() {
  local env="$1"
  local tgz="${PACK_DIR}/${env}.tar.gz"
  local dest="${CONDA_BASE}/envs/${env}"

  echo "[3/5] env=${env}"
  echo "[3/5] tarball=${tgz}"
  echo "[3/5] dest=${dest}"

  echo "[4/5] mkdir -p ${dest} ..."
  mkdir -p "$dest"
  echo "[4/5] extracting ${tgz} -> ${dest} (may take a while) ..."
  tar -xzf "$tgz" -C "$dest"
  echo "[5/5] running conda-unpack ..."
  "${dest}/bin/conda-unpack"    # 关键：修正硬编码路径
  echo "OK: $env"
}

echo "[3/5] starting unpack loop ..."
for env in neoag-tools; do
  unpack_one "$env"
done
echo "[done] all environments unpacked"