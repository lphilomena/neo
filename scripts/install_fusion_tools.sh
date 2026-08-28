#!/usr/bin/env bash
# Install Arriba + Nextflow fusion environment and optional STAR-Fusion / FusionCatcher sources.
#
# Usage:
#   bash scripts/install_fusion_tools.sh
#   source conf/tools.env.sh
#   neoag check-tools | grep -E 'arriba|star-fusion|fusioncatcher'
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_BASE="${NEOAG_CONDA_BASE:-$(command conda info --base)}"
ENV_NAME="${NEOAG_FUSION_ENV:-neoag-fusion}"
ENV_PREFIX="${CONDA_BASE}/envs/${ENV_NAME}"
FUSION_YML="${ROOT}/conda/env.neoag-fusion.yml"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-${ROOT}}"
STAR_HOME="${NEOAG_STAR_FUSION_HOME:-${TOOLS_ROOT}/tools/STAR-Fusion}"
FC_HOME="${NEOAG_FUSIONCATCHER_HOME:-${TOOLS_ROOT}/tools/fusioncatcher}"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found" >&2
  exit 1
fi
conda_safe() {
  set +u
  conda "$@"
  local rc=$?
  set -u
  return "$rc"
}

# shellcheck disable=SC1091
set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
set -u

if [[ "${NEOAG_USE_MAMBA:-0}" == "1" ]] && command -v mamba >/dev/null 2>&1; then
  CONDA_RUNNER=mamba
else
  CONDA_RUNNER=conda
fi

env_exists() { conda_safe env list | awk '{print $1}' | grep -qx "$1"; }

create_env() {
  if [[ "${CONDA_RUNNER}" == "mamba" ]]; then
    mamba env create -p "${ENV_PREFIX}" -f "${FUSION_YML}" -y \
      --override-channels -c conda-forge -c bioconda
  else
    conda env create -p "${ENV_PREFIX}" -f "${FUSION_YML}" -y
  fi
}

update_env() {
  if [[ "${CONDA_RUNNER}" == "mamba" ]]; then
    mamba env update -p "${ENV_PREFIX}" -f "${FUSION_YML}" --prune
  else
    conda env update -p "${ENV_PREFIX}" -f "${FUSION_YML}" --prune
  fi
}

fusion_env_has_arriba() {
  env_exists "${ENV_NAME}" && conda_safe run -n "${ENV_NAME}" arriba -h >/dev/null 2>&1
}

install_github_snapshot() {
  local label="$1" url="$2" target="$3"
  local work archive
  work="$(mktemp -d)"
  archive="${work}/source.tar.gz"
  echo "==> Downloading ${label} source snapshot"
  if command -v aria2c >/dev/null 2>&1; then
    aria2c -x 8 -s 8 -k 1M --allow-overwrite=true --file-allocation=none \
      -d "${work}" -o source.tar.gz "${url}"
  else
    curl -fL --retry 3 -o "${archive}" "${url}"
  fi
  mkdir -p "${work}/source" "${target}"
  tar -xzf "${archive}" -C "${work}/source" --strip-components=1
  rsync -a --delete "${work}/source/" "${target}/"
  rm -rf "${work}"
}

if [[ "${NEOAG_FORCE_ENV_UPDATE:-0}" == "1" ]]; then
  if env_exists "${ENV_NAME}"; then
    echo "==> Updating ${ENV_NAME} (NEOAG_FORCE_ENV_UPDATE=1)"
    update_env
  else
    echo "==> Creating ${ENV_NAME} (NEOAG_FORCE_ENV_UPDATE=1)"
    create_env
  fi
elif fusion_env_has_arriba; then
  echo "==> ${ENV_NAME} already has arriba; skipping env update (set NEOAG_FORCE_ENV_UPDATE=1 to refresh)"
elif env_exists "${ENV_NAME}"; then
  echo "==> Updating ${ENV_NAME}"
  update_env
else
  echo "==> Creating ${ENV_NAME}"
  create_env
fi

if [[ "${NEOAG_SKIP_STAR_FUSION_CLONE:-0}" != "1" && ! -x "${STAR_HOME}/STAR-Fusion" ]]; then
  mkdir -p "$(dirname "${STAR_HOME}")"
  install_github_snapshot "STAR-Fusion" \
    "${NEOAG_STAR_FUSION_SNAPSHOT_URL:-https://gh-proxy.com/https://github.com/STAR-Fusion/STAR-Fusion/archive/refs/heads/master.tar.gz}" \
    "${STAR_HOME}"
  chmod +x "${STAR_HOME}/STAR-Fusion" 2>/dev/null || true
fi

if [[ "${NEOAG_SKIP_FUSIONCATCHER_CLONE:-0}" != "1" && ! -d "${FC_HOME}/bin" ]]; then
  mkdir -p "$(dirname "${FC_HOME}")"
  install_github_snapshot "FusionCatcher" \
    "${NEOAG_FUSIONCATCHER_SNAPSHOT_URL:-https://gh-proxy.com/https://github.com/ndaniel/fusioncatcher/archive/refs/heads/master.tar.gz}" \
    "${FC_HOME}"
fi

echo "==> Smoke tests"
if fusion_env_has_arriba; then
  conda_safe run -n "${ENV_NAME}" arriba -h | head -6
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
TOOLS_ENV="${NEOAG_TOOLS_ENV:-${ROOT}/conf/tools.env.local.sh}"
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

echo "==> Done. Run: source conf/tools.env.sh && neoag check-tools | grep -E 'arriba|star-fusion|fusioncatcher'"
