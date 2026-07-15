#!/usr/bin/env bash
set -euo pipefail
BUNDLE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-${BUNDLE_ROOT}/project/neo-na0707_upload_release-3}"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-${BUNDLE_ROOT}/../test_env_tools}"
MINIFORGE_ROOT="${NEOAG_CONDA_BASE:-${TOOLS_ROOT}/miniforge3}"
NEOAG_CORE_ENV="${NEOAG_CORE_ENV:-neoag-core}"
NEOAG_TOOLS_ENV="${NEOAG_TOOLS_ENV:-neoag-tools}"
log() { printf "[neoag-migrate] %s\n" "$*"; }
warn() { printf "[neoag-migrate] WARN: %s\n" "$*" >&2; }
die() { printf "[neoag-migrate] ERROR: %s\n" "$*" >&2; exit 1; }
ensure_project() {
  if [[ ! -d "$PROJECT_ROOT/src/neoag_v03" ]]; then
    mkdir -p "${BUNDLE_ROOT}/project"
    tar -xzf "${BUNDLE_ROOT}/release/neo-na0707_upload_release-3.tar.gz" -C "${BUNDLE_ROOT}/project"
  fi
}
ensure_conda() {
  if [[ -x "${MINIFORGE_ROOT}/bin/conda" ]]; then
    export PATH="${MINIFORGE_ROOT}/bin:${PATH}"
  elif command -v conda >/dev/null 2>&1; then
    MINIFORGE_ROOT="$(conda info --base)"
  else
    die "conda not found. Install Miniforge first or set NEOAG_CONDA_BASE=/path/to/miniforge3."
  fi
  export NEOAG_CONDA_BASE="$MINIFORGE_ROOT"
  # shellcheck source=/dev/null
  [[ -f "${NEOAG_CONDA_BASE}/etc/profile.d/conda.sh" ]] && source "${NEOAG_CONDA_BASE}/etc/profile.d/conda.sh"
}
activate_paths() {
  ensure_project
  ensure_conda
  export NEOAG_PROJECT_ROOT="$PROJECT_ROOT"
  export NEOAG_TOOLS_ROOT="$TOOLS_ROOT"
  export NXF_HOME="${NXF_HOME:-${TOOLS_ROOT}/nextflow}"
  export PYTHONPATH="${PROJECT_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
  export PATH="${PROJECT_ROOT}/bin:${NEOAG_CONDA_BASE}/envs/${NEOAG_CORE_ENV}/bin:${NEOAG_CONDA_BASE}/envs/${NEOAG_TOOLS_ENV}/bin:${NEOAG_CONDA_BASE}/bin:${PATH}"
  export LD_LIBRARY_PATH="${NEOAG_CONDA_BASE}/envs/${NEOAG_TOOLS_ENV}/lib:${NEOAG_CONDA_BASE}/envs/${NEOAG_CORE_ENV}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
  # MHCflurry 2.0.6 uses legacy Keras session APIs; keep TensorFlow on tf_keras.
  export TF_USE_LEGACY_KERAS="${TF_USE_LEGACY_KERAS:-1}"
  export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:--1}"
  export TF_CPP_MIN_LOG_LEVEL="${TF_CPP_MIN_LOG_LEVEL:-2}"
}
run_in_project() {
  activate_paths
  cd "$PROJECT_ROOT"
  "$@"
}
