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
#   neoag check-tools | grep -E 'ascat|pyclone'
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_BASE="${NEOAG_CONDA_BASE:-$(command conda info --base)}"
ASCAT_ENV="${NEOAG_ASCAT_ENV:-neoag-ascat}"
PYCLONE_ENV="${NEOAG_PYCLONE_ENV:-neoag-pyclone}"
ASCAT_YML="${ROOT}/conda/env.neoag-ascat.yml"
PYCLONE_YML="${ROOT}/conda/env.neoag-pyclone.yml"
TOOLS_ENV="${NEOAG_TOOLS_ENV:-${ROOT}/conf/tools.env.local.sh}"
BIOC_CACHE_HELPER="${ROOT}/.agents/skills/neoag-remote-deploy/scripts/with_bioc_data_cache.sh"

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

detect_env_prefix() {
  local env_name="$1" prefix
  for prefix in \
    "${CONDA_BASE}/envs/${env_name}" \
    "$(dirname "${CONDA_BASE}")/${env_name}" \
    "${NEOAG_TOOLS_ROOT:-}/${env_name}"; do
    [[ -n "${prefix}" && -d "${prefix}" ]] && { printf '%s\n' "${prefix}"; return 0; }
  done
  conda_safe env list | awk -v name="${env_name}" '$1 == name {print $NF; exit}'
}

env_exists() { [[ -n "$(detect_env_prefix "$1")" ]]; }

create_or_update_env() {
  local env_name="$1" yml="$2"
  local -a runner=("${CONDA_RUNNER}")
  echo "==> Installing/updating ${env_name} using ${CONDA_RUNNER}: ${yml}"
  if [[ "${env_name}" == "${ASCAT_ENV}" && -x "${BIOC_CACHE_HELPER}" ]]; then
    runner=("${BIOC_CACHE_HELPER}" --conda-base "${CONDA_BASE}" \
      --cache-root "${NEOAG_TOOLS_ROOT:-${ROOT}}/install_cache" \
      --package-key genomeinfodbdata-1.2.13 -- "${CONDA_RUNNER}")
  fi
  if env_exists "${env_name}"; then
    "${runner[@]}" env update -n "${env_name}" -f "${yml}" --prune \
      --override-channels -c conda-forge -c bioconda
  else
    "${runner[@]}" env create -n "${env_name}" -f "${yml}" -y \
      --override-channels -c conda-forge -c bioconda
  fi
}

env_has_ascat() {
  local prefix
  prefix="$(detect_env_prefix "$1")"
  [[ -x "${prefix}/bin/Rscript" ]] && \
    "${prefix}/bin/Rscript" -e 'quit(status=ifelse(requireNamespace("ASCAT", quietly=TRUE),0,1))' >/dev/null 2>&1
}

env_has_pyclone() {
  local prefix
  prefix="$(detect_env_prefix "$1")"
  [[ -x "${prefix}/bin/pyclone-vi" ]] && \
    "${prefix}/bin/pyclone-vi" --version >/dev/null 2>&1
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

ASCAT_PREFIX="$(detect_env_prefix "${ASCAT_ENV}")"
PYCLONE_PREFIX="$(detect_env_prefix "${PYCLONE_ENV}")"
[[ -n "${ASCAT_PREFIX}" ]] || { echo "ERROR: ASCAT environment prefix not found" >&2; exit 1; }
[[ -n "${PYCLONE_PREFIX}" ]] || { echo "ERROR: PyClone environment prefix not found" >&2; exit 1; }

mkdir -p "${ROOT}/bin"
cat > "${ROOT}/bin/ascat.R" <<EOF
#!/usr/bin/env bash
set -euo pipefail
if [[ "\${1:-}" == "--version" || "\${1:-}" == "-v" ]]; then
  "${ASCAT_PREFIX}/bin/Rscript" -e 'cat(as.character(utils::packageVersion("ASCAT")), "\\n")'
  exit 0
fi
if [[ "\$#" -eq 0 ]]; then
  echo "ASCAT wrapper. For custom analyses, run: conda activate ${ASCAT_ENV}; Rscript your_ascat_script.R" >&2
  exit 0
fi
"${ASCAT_PREFIX}/bin/Rscript" "\$@"
EOF
chmod +x "${ROOT}/bin/ascat.R"

cat > "${ROOT}/bin/pyclone" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec "${PYCLONE_PREFIX}/bin/pyclone-vi" "\$@"
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

echo "==> Done. Run: source conf/tools.env.sh && neoag check-tools | grep -E 'ascat|pyclone'"
