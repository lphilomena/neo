#!/usr/bin/env bash
# Install Arriba + Nextflow fusion environment and optional STAR-Fusion / FusionCatcher sources.
#
# Usage:
#   bash scripts/install_fusion_tools.sh
#   source conf/tools.env.sh
#   neoag-v03 check-tools | grep -E 'arriba|star-fusion|fusioncatcher'
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_BASE="${NEOAG_CONDA_BASE:-$(conda info --base)}"
ENV_NAME="${NEOAG_FUSION_ENV:-neoag-fusion}"
FUSION_YML="${ROOT}/conda/env.neoag-fusion.yml"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-${ROOT}}"
STAR_HOME="${NEOAG_STAR_FUSION_HOME:-${TOOLS_ROOT}/tools/STAR-Fusion}"
FC_HOME="${NEOAG_FUSIONCATCHER_HOME:-${TOOLS_ROOT}/tools/fusioncatcher}"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found" >&2
  exit 1
fi
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"

if [[ "${NEOAG_USE_MAMBA:-0}" == "1" ]] && command -v mamba >/dev/null 2>&1; then
  CONDA_RUNNER=mamba
else
  CONDA_RUNNER=conda
fi

env_exists() { conda env list | awk '{print $1}' | grep -qx "$1"; }

fusion_env_has_arriba() {
  env_exists "${ENV_NAME}" && conda run -n "${ENV_NAME}" arriba -h >/dev/null 2>&1
}

if [[ "${NEOAG_FORCE_ENV_UPDATE:-0}" == "1" ]]; then
  if env_exists "${ENV_NAME}"; then
    echo "==> Updating ${ENV_NAME} (NEOAG_FORCE_ENV_UPDATE=1)"
    "${CONDA_RUNNER}" env update -n "${ENV_NAME}" -f "${FUSION_YML}" --prune
  else
    echo "==> Creating ${ENV_NAME} (NEOAG_FORCE_ENV_UPDATE=1)"
    "${CONDA_RUNNER}" env create -n "${ENV_NAME}" -f "${FUSION_YML}" -y
  fi
elif fusion_env_has_arriba; then
  echo "==> ${ENV_NAME} already has arriba; skipping env update (set NEOAG_FORCE_ENV_UPDATE=1 to refresh)"
elif env_exists "${ENV_NAME}"; then
  echo "==> Updating ${ENV_NAME}"
  "${CONDA_RUNNER}" env update -n "${ENV_NAME}" -f "${FUSION_YML}" --prune
else
  echo "==> Creating ${ENV_NAME}"
  "${CONDA_RUNNER}" env create -n "${ENV_NAME}" -f "${FUSION_YML}" -y
fi

if [[ "${NEOAG_SKIP_STAR_FUSION_CLONE:-0}" != "1" && ! -x "${STAR_HOME}/STAR-Fusion" ]]; then
  mkdir -p "$(dirname "${STAR_HOME}")"
  if [[ ! -d "${STAR_HOME}/.git" ]]; then
    git clone --depth 1 https://github.com/STAR-Fusion/STAR-Fusion.git "${STAR_HOME}"
  fi
  chmod +x "${STAR_HOME}/STAR-Fusion" 2>/dev/null || true
fi

if [[ "${NEOAG_SKIP_FUSIONCATCHER_CLONE:-0}" != "1" && ! -d "${FC_HOME}/bin" ]]; then
  mkdir -p "$(dirname "${FC_HOME}")"
  if [[ ! -d "${FC_HOME}/.git" ]]; then
    git clone --depth 1 https://github.com/ndaniel/fusioncatcher.git "${FC_HOME}"
  fi
fi

echo "==> Smoke tests"
if fusion_env_has_arriba; then
  conda run -n "${ENV_NAME}" arriba -h | head -6
elif command -v arriba >/dev/null 2>&1; then
  arriba -h | head -6
else
  echo "ERROR: arriba not available after install" >&2
  exit 1
fi
if [[ -x "${STAR_HOME}/STAR-Fusion" ]]; then
  "${STAR_HOME}/STAR-Fusion" --version | head -1 || true
else
  echo "WARN: STAR-Fusion not installed; set NEOAG_SKIP_STAR_FUSION_CLONE=1 to silence or rerun after network access." >&2
fi


PROJECT_BIN="${ROOT}/bin"
TOOLS_ENV="${ROOT}/conf/tools.env.sh"
ARRIBA_BIN="${CONDA_BASE}/envs/${ENV_NAME}/bin/arriba"
mkdir -p "${PROJECT_BIN}" "${ROOT}/conf"
if [[ -x "${ARRIBA_BIN}" ]]; then
  cat > "${PROJECT_BIN}/arriba" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec "${ARRIBA_BIN}" "\$@"
EOF
  chmod +x "${PROJECT_BIN}/arriba"
fi

fusion_block="export NEOAG_FUSION_ENV=\"${ENV_NAME}\"
export NEOAG_ARRIBA_BIN=\"${ARRIBA_BIN}\"
export PATH=\"${PROJECT_BIN}:${CONDA_BASE}/envs/${ENV_NAME}/bin:\$PATH\""
if [[ -f "${TOOLS_ENV}" ]]; then
  if grep -q 'NEOAG_FUSION_ENV' "${TOOLS_ENV}"; then
    sed -i "s|^export NEOAG_FUSION_ENV=.*|export NEOAG_FUSION_ENV=\"${ENV_NAME}\"|" "${TOOLS_ENV}"
    sed -i "s|^export NEOAG_ARRIBA_BIN=.*|export NEOAG_ARRIBA_BIN=\"${ARRIBA_BIN}\"|" "${TOOLS_ENV}"
  else
    printf '\n# Fusion tools (Arriba/STAR-Fusion/FusionCatcher)\n%s\n' "${fusion_block}" >> "${TOOLS_ENV}"
  fi
else
  printf '%s\n' "${fusion_block}" > "${TOOLS_ENV}"
fi

echo "==> Done. Run: source conf/tools.env.sh && neoag-v03 check-tools | grep -E 'arriba|star-fusion|fusioncatcher'"
