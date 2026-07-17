#!/usr/bin/env bash
# Install ASCAT and PyClone-VI into dedicated conda envs.
#
# This script intentionally uses conda by default. Do not `pip install mamba`:
# that package is not the conda-forge mamba solver and will fail on `mamba env create`.
# Set NEOAG_USE_MAMBA=1 only if the real mamba executable is already available.
#
# Usage:
#   bash scripts/install_ascat_pyclone.sh
#   source conf/tools.env.sh
#   neoag-v03 check-tools | grep -E 'ascat|pyclone'
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_BASE="${NEOAG_CONDA_BASE:-$(command conda info --base)}"
ASCAT_ENV="${NEOAG_ASCAT_ENV:-neoag-ascat}"
PYCLONE_ENV="${NEOAG_PYCLONE_ENV:-neoag-pyclone}"
ASCAT_YML="${ROOT}/conda/env.neoag-ascat.yml"
PYCLONE_YML="${ROOT}/conda/env.neoag-pyclone.yml"
TOOLS_ENV="${ROOT}/conf/tools.env.sh"

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

if [[ "${NEOAG_USE_MAMBA:-0}" == "1" ]] && command -v mamba >/dev/null 2>&1 && mamba --version >/dev/null 2>&1; then
  CONDA_RUNNER=mamba
else
  CONDA_RUNNER=conda
fi

env_exists() { conda_safe env list | awk '{print $1}' | grep -qx "$1"; }

create_or_update_env() {
  local env_name="$1" yml="$2"
  echo "==> Installing/updating ${env_name} using ${CONDA_RUNNER}: ${yml}"
  if env_exists "${env_name}"; then
    "${CONDA_RUNNER}" env update -n "${env_name}" -f "${yml}" --prune
  else
    "${CONDA_RUNNER}" env create -n "${env_name}" -f "${yml}" -y
  fi
}

env_has_ascat() {
  conda_safe run -n "$1" Rscript -e 'quit(status=ifelse(requireNamespace("ASCAT", quietly=TRUE),0,1))' >/dev/null 2>&1
}

env_has_pyclone() {
  [[ -x "${CONDA_BASE}/envs/$1/bin/pyclone-vi" ]] || \
    conda_safe run -n "$1" pyclone-vi --version >/dev/null 2>&1
}

ascat_ready() {
  [[ -x "${ROOT}/bin/ascat.R" ]] && "${ROOT}/bin/ascat.R" --version >/dev/null 2>&1
}

pyclone_ready() {
  [[ -x "${ROOT}/bin/pyclone" ]] && "${ROOT}/bin/pyclone" --version >/dev/null 2>&1
}

if ascat_ready && pyclone_ready; then
  echo "==> ASCAT/PyClone wrappers already present; skipping env update"
else
  if [[ "${NEOAG_FORCE_ENV_UPDATE:-0}" == "1" ]] || ! env_exists "${ASCAT_ENV}"; then
    create_or_update_env "${ASCAT_ENV}" "${ASCAT_YML}"
  elif env_has_ascat "${ASCAT_ENV}"; then
    echo "==> ASCAT package present in ${ASCAT_ENV}; skipping env update"
  else
    echo "==> ${ASCAT_ENV} exists but ASCAT R package missing; recreating env"
    conda_safe env remove -n "${ASCAT_ENV}" -y
    create_or_update_env "${ASCAT_ENV}" "${ASCAT_YML}"
  fi

  if [[ "${NEOAG_FORCE_ENV_UPDATE:-0}" == "1" ]] || ! env_exists "${PYCLONE_ENV}"; then
    create_or_update_env "${PYCLONE_ENV}" "${PYCLONE_YML}"
  elif env_has_pyclone "${PYCLONE_ENV}"; then
    echo "==> PyClone-VI present in ${PYCLONE_ENV}; skipping env update"
  else
    echo "==> PyClone-VI missing in ${PYCLONE_ENV}; refreshing env"
    create_or_update_env "${PYCLONE_ENV}" "${PYCLONE_YML}"
  fi
fi

mkdir -p "${ROOT}/bin"
cat > "${ROOT}/bin/ascat.R" <<EOF
#!/usr/bin/env bash
set -euo pipefail
if [[ "\${1:-}" == "--version" || "\${1:-}" == "-v" ]]; then
  "${CONDA_BASE}/bin/conda" run -n "${ASCAT_ENV}" Rscript -e 'cat(as.character(utils::packageVersion("ASCAT")), "\\n")'
  exit 0
fi
if [[ "\$#" -eq 0 ]]; then
  echo "ASCAT wrapper. For custom analyses, run: conda activate ${ASCAT_ENV}; Rscript your_ascat_script.R" >&2
  exit 0
fi
"${CONDA_BASE}/bin/conda" run -n "${ASCAT_ENV}" Rscript "\$@"
EOF
chmod +x "${ROOT}/bin/ascat.R"

cat > "${ROOT}/bin/pyclone" <<EOF
#!/usr/bin/env bash
set -euo pipefail
if command -v "${CONDA_BASE}/envs/${PYCLONE_ENV}/bin/pyclone-vi" >/dev/null 2>&1; then
  exec "${CONDA_BASE}/envs/${PYCLONE_ENV}/bin/pyclone-vi" "\$@"
fi
exec "${CONDA_BASE}/bin/conda" run -n "${PYCLONE_ENV}" pyclone-vi "\$@"
EOF
chmod +x "${ROOT}/bin/pyclone"

echo "==> Smoke tests"
if ! "${ROOT}/bin/ascat.R" --version >/dev/null 2>&1; then
  echo "WARN: ASCAT wrapper version check failed; inspect env ${ASCAT_ENV}" >&2
fi
if ! "${ROOT}/bin/pyclone" --version >/dev/null 2>&1; then
  echo "WARN: PyClone-VI version check failed; inspect env ${PYCLONE_ENV}" >&2
fi

mkdir -p "${ROOT}/conf"
if [[ ! -f "${TOOLS_ENV}" ]]; then
  cat > "${TOOLS_ENV}" <<EOF
export NEOAG_PROJECT_ROOT="${ROOT}"
export NEOAG_TOOLS_ROOT="${ROOT}"
export NEOAG_CONDA_BASE="${CONDA_BASE}"
export NEOAG_CONDA_ENV="neoag-tools"
EOF
fi
if ! grep -q 'ASCAT / PyClone-VI — installed via scripts/install_ascat_pyclone.sh' "${TOOLS_ENV}"; then
  cat >> "${TOOLS_ENV}" <<EOF

# ASCAT / PyClone-VI — installed via scripts/install_ascat_pyclone.sh
export NEOAG_ASCAT_ENV="${ASCAT_ENV}"
export ASCAT_HOME="${ROOT}/bin"
export NEOAG_PYCLONE_ENV="${PYCLONE_ENV}"
export NEOAG_PYCLONE_BIN="${ROOT}/bin/pyclone"
export PATH="${ROOT}/bin:\${PATH}"
EOF
fi

echo "==> Done. Run: source conf/tools.env.sh && neoag-v03 check-tools | grep -E 'ascat|pyclone'"
