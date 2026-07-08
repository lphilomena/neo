#!/usr/bin/env bash
# Sequential prebuild for remaining EasyFuse conda envs (mamba -y; do not rely on Nextflow).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/conf/tools.env.sh"

CONDA_CACHE="${ROOT}/work/.nextflow_conda"
LOG="${ROOT}/work/easyfuse_conda_prebuild.log"

export CONDA_ALWAYS_YES=1

wait_for_mamba_free() {
  while pgrep -f 'mamba env create' >/dev/null 2>&1; do
    sleep 15
  done
  rm -f "${NEOAG_CONDA_BASE}/pkgs/pkgs.lock" 2>/dev/null || true
}

prebuild_conda_env() {
  local env_id="$1"
  local yml="$2"
  local check_bin="$3"
  local prefix="${CONDA_CACHE}/env-${env_id}"

  if [[ -x "${prefix}/bin/${check_bin}" ]]; then
    echo "    ${yml}: already ready (${prefix})"
    return 0
  fi

  wait_for_mamba_free
  echo "==> Pre-building EasyFuse ${yml} ..."
  rm -rf "${prefix}"
  mamba env create -y \
    --prefix "${prefix}" \
    --file "${NEOAG_EASYFUSE_HOME}/environments/${yml}"
  if [[ ! -x "${prefix}/bin/${check_bin}" ]]; then
    echo "ERROR: ${yml} built but ${check_bin} missing under ${prefix}/bin" >&2
    ls -la "${prefix}/bin" >&2 || true
    exit 1
  fi
  bash "${ROOT}/scripts/fix_easyfuse_pyeasyfuse_env.sh" >/dev/null 2>&1 || true
  echo "    ${yml}: done ($(du -sh "${prefix}" | awk '{print $1}'))"
}

exec >> "${LOG}" 2>&1
echo ""
echo "==> easyfuse_prebuild_remaining_envs $(date -Is)"

prebuild_conda_env \
  "6f2b394c864eeaa5-8f88fe4572f59d9bb818f7644ca8f1fa" \
  "alignment.yml" \
  "STAR"

prebuild_conda_env \
  "e4b7dd3f8c4b23e0-ae78dd87a2293bf5178c1752ac54f434" \
  "starfusion.yml" \
  "STAR-Fusion"

prebuild_conda_env \
  "cdd3345c1c7c6ebb-12bec81221d755f4cf369bdee252c72d" \
  "fusioncatcher.yml" \
  "fusioncatcher"

prebuild_conda_env \
  "32b8951a86fd0d30-0412febc37d8c32b8c4b8283292cffe1" \
  "requantification.yml" \
  "STAR"

echo "==> easyfuse_prebuild_remaining_envs done $(date -Is)"
