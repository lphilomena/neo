#!/usr/bin/env bash
# Seed Nextflow conda env paths from prebuilt env-* prefixes (same yml hash).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_CACHE="${ROOT}/work/.nextflow_conda"

seed_env_from_prefix() {
  local yml_hash="$1"
  local target_suffix="$2"
  local check_bin="${3:-STAR}"

  local target="${CONDA_CACHE}/env-${yml_hash}-${target_suffix}"
  if [[ -x "${target}/bin/${check_bin}" ]]; then
    echo "    ${target}: ready"
    return 0
  fi

  local src=""
  for candidate in "${CONDA_CACHE}/env-${yml_hash}-"*; do
    [[ -d "${candidate}" ]] || continue
    [[ "${candidate}" == "${target}" ]] && continue
    if [[ -x "${candidate}/bin/${check_bin}" ]]; then
      src="${candidate}"
      break
    fi
  done

  if [[ -z "${src}" ]]; then
    echo "WARN: no source env for hash ${yml_hash} (need ${check_bin})" >&2
    return 1
  fi

  echo "==> seed ${target} <= ${src}"
  rm -rf "${target}"
  cp -a "${src}" "${target}"
}

echo "==> seed_easyfuse_conda_envs $(date -Is)"

# Nextflow 24.10.1 hash for qc.yml (FASTP) when createOptions='-y --yes'.
# Prefer a real mamba create over cp -a (copied prefixes break `conda activate`).
QC_TARGET="${CONDA_CACHE}/env-dd32f47bd0756865-3ca3e4910cc80d7677b2b976f6b1230f"
if [[ ! -x "${QC_TARGET}/bin/fastp" ]]; then
  echo "==> mamba create qc env ${QC_TARGET}"
  QC_YML="${NEOAG_EASYFUSE_HOME:-${ROOT}/../envs/tools/EasyFuse}/environments/qc.yml"
  [[ -s "${QC_YML}" ]] || { echo "ERROR: EasyFuse qc.yml not found: ${QC_YML}" >&2; exit 1; }
  mamba env create -y --prefix "${QC_TARGET}" \
    --file "${QC_YML}"
else
  echo "    ${QC_TARGET}: ready"
fi

# Nextflow 24.10.1 hash for requantification.yml (from work/.nextflow.log).
seed_env_from_prefix \
  "32b8951a86fd0d30" \
  "1ecc45a662409914602fad8296a8299f" \
  "STAR" || true

# Nextflow hash for requantification_wo_easyfuse.yml (STAR_CUSTOM).
seed_env_from_prefix \
  "f48ef6da3d6b5f11" \
  "f3998449ecfa3870b52b2b8e0a3d35a6" \
  "STAR" || true

bash "${ROOT}/scripts/fix_easyfuse_pyeasyfuse_env.sh"
echo "==> seed_easyfuse_conda_envs done"
