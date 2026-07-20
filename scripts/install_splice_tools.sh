#!/usr/bin/env bash
# Install/register splice-junction neoantigen helper tools.
#
# Installs:
#   - RegTools in a dedicated neoag-splice conda env
#   - pVACsplice wrapper from the existing neoag-tools/pVACtools env
# Optional:
#   - SNAF from local source or approved Git URL
#   - ASNEO / NeoSplice / splice2neo source directories as registered wrappers
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_BASE="${NEOAG_CONDA_BASE:-$(command conda info --base)}"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-${ROOT}}"
ENV_NAME="${NEOAG_SPLICE_ENV:-neoag-splice}"
TOOLS_ENV="${ROOT}/conf/tools.env.sh"
BIN_DIR="${ROOT}/bin"
YML="${ROOT}/conda/env.neoag-splice.yml"
export CONDA_CHANNEL_ALIAS="${NEOAG_CONDA_CHANNEL_ALIAS:-${CONDA_CHANNEL_ALIAS:-https://conda.anaconda.org}}"

INSTALL_SNAF="${NEOAG_INSTALL_SNAF:-0}"
SNAF_SOURCE="${SNAF_SOURCE:-${NEOAG_SNAF_SOURCE:-}}"
SNAF_GIT_URL="${SNAF_GIT_URL:-${NEOAG_SNAF_GIT_URL:-https://github.com/frankligy/SNAF.git}}"
SNAF_GIT_REF="${SNAF_GIT_REF:-${NEOAG_SNAF_GIT_REF:-e23ce39512a1a7f58c74e59b4b7cedc89248b908}}"
ASNEO_SOURCE="${ASNEO_SOURCE:-${NEOAG_ASNEO_SOURCE:-}}"
NEOSPLICE_SOURCE="${NEOSPLICE_SOURCE:-${NEOAG_NEOSPLICE_SOURCE:-}}"
SPLICE2NEO_SOURCE="${SPLICE2NEO_SOURCE:-${NEOAG_SPLICE2NEO_SOURCE:-}}"

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

set +u
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"
set -u

if [[ "${NEOAG_USE_MAMBA:-0}" == "1" ]] && command -v mamba >/dev/null 2>&1; then
  CONDA_RUNNER=mamba
else
  CONDA_RUNNER=conda
fi

env_exists() { conda_safe env list | awk '{print $1}' | grep -qx "$1"; }
env_has_regtools() { env_exists "$1" && conda_safe run -n "$1" regtools junctions extract -h >/dev/null 2>&1; }

if [[ "${NEOAG_FORCE_ENV_UPDATE:-0}" == "1" ]]; then
  if env_exists "${ENV_NAME}"; then
    "${CONDA_RUNNER}" env update -n "${ENV_NAME}" -f "${YML}" --prune
  else
    "${CONDA_RUNNER}" env create -n "${ENV_NAME}" -f "${YML}" -y
  fi
elif env_has_regtools "${ENV_NAME}"; then
  echo "==> ${ENV_NAME} already has RegTools; skipping env update"
elif env_exists "${ENV_NAME}"; then
  "${CONDA_RUNNER}" env update -n "${ENV_NAME}" -f "${YML}" --prune
else
  "${CONDA_RUNNER}" env create -n "${ENV_NAME}" -f "${YML}" -y
fi

mkdir -p "${BIN_DIR}" "${ROOT}/conf" "${TOOLS_ROOT}/tools"

cat > "${BIN_DIR}/regtools-neoag" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec "${CONDA_BASE}/bin/conda" run -n "${ENV_NAME}" regtools "\$@"
EOF
chmod +x "${BIN_DIR}/regtools-neoag"

PVACSPLICE_BIN="${CONDA_BASE}/envs/neoag-tools/bin/pvacsplice"
if [[ -x "${PVACSPLICE_BIN}" ]]; then
  cat > "${BIN_DIR}/pvacsplice-neoag" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec "${PVACSPLICE_BIN}" "\$@"
EOF
  chmod +x "${BIN_DIR}/pvacsplice-neoag"
else
  echo "WARN: pvacsplice not found in neoag-tools; install core env first." >&2
fi

install_python_source() {
  local label="$1" source_dir="$2" target_dir="$3" wrapper="$4" module="$5"
  [[ -n "${source_dir}" ]] || return 0
  [[ -e "${source_dir}" ]] || { echo "ERROR: ${label} source missing: ${source_dir}" >&2; exit 42; }
  mkdir -p "$(dirname "${target_dir}")"
  rsync -a --delete "${source_dir}/" "${target_dir}/"
  if [[ -f "${target_dir}/pyproject.toml" || -f "${target_dir}/setup.py" ]]; then
    "${CONDA_BASE}/bin/conda" run -n "${ENV_NAME}" python -m pip install "${target_dir}"
  fi
  cat > "${BIN_DIR}/${wrapper}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export PYTHONPATH="${target_dir}:\${PYTHONPATH:-}"
exec "${CONDA_BASE}/bin/conda" run -n "${ENV_NAME}" python -m ${module} "\$@"
EOF
  chmod +x "${BIN_DIR}/${wrapper}"
}

if [[ "${INSTALL_SNAF}" == "1" ]]; then
  SNAF_HOME="${TOOLS_ROOT}/tools/SNAF"
  if [[ -n "${SNAF_SOURCE}" ]]; then
    install_python_source "SNAF" "${SNAF_SOURCE}" "${SNAF_HOME}" "snaf-neoag" "snaf"
  else
    echo "==> Installing SNAF from approved Git source ${SNAF_GIT_URL}@${SNAF_GIT_REF}"
    "${CONDA_BASE}/bin/conda" run -n "${ENV_NAME}" python -m pip install "git+${SNAF_GIT_URL}@${SNAF_GIT_REF}"
    cat > "${BIN_DIR}/snaf-neoag" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec "${CONDA_BASE}/bin/conda" run -n "${ENV_NAME}" python -m snaf "\$@"
EOF
    chmod +x "${BIN_DIR}/snaf-neoag"
  fi
fi

install_python_source "ASNEO" "${ASNEO_SOURCE}" "${TOOLS_ROOT}/tools/ASNEO" "asneo-neoag" "ASNEO"
install_python_source "NeoSplice" "${NEOSPLICE_SOURCE}" "${TOOLS_ROOT}/tools/NeoSplice" "neosplice-neoag" "NeoSplice"
install_python_source "splice2neo" "${SPLICE2NEO_SOURCE}" "${TOOLS_ROOT}/tools/splice2neo" "splice2neo-neoag" "splice2neo"

if [[ -f "${TOOLS_ENV}" ]]; then
  if grep -q 'NEOAG_SPLICE_ENV' "${TOOLS_ENV}"; then
    sed -i "s|^export NEOAG_SPLICE_ENV=.*|export NEOAG_SPLICE_ENV="${ENV_NAME}"|" "${TOOLS_ENV}"
    sed -i "s|^export NEOAG_REGTOOLS_BIN=.*|export NEOAG_REGTOOLS_BIN="${BIN_DIR}/regtools-neoag"|" "${TOOLS_ENV}"
    sed -i "s|^export NEOAG_PVACSPLICE_BIN=.*|export NEOAG_PVACSPLICE_BIN="${BIN_DIR}/pvacsplice-neoag"|" "${TOOLS_ENV}"
  else
    cat >> "${TOOLS_ENV}" <<EOF

# Splice neoantigen tools
export NEOAG_SPLICE_ENV="${ENV_NAME}"
export NEOAG_REGTOOLS_BIN="${BIN_DIR}/regtools-neoag"
export NEOAG_PVACSPLICE_BIN="${BIN_DIR}/pvacsplice-neoag"
export PATH="${BIN_DIR}:${CONDA_BASE}/envs/${ENV_NAME}/bin:\$PATH"
EOF
  fi
else
  cat > "${TOOLS_ENV}" <<EOF
export NEOAG_SPLICE_ENV="${ENV_NAME}"
export NEOAG_REGTOOLS_BIN="${BIN_DIR}/regtools-neoag"
export NEOAG_PVACSPLICE_BIN="${BIN_DIR}/pvacsplice-neoag"
export PATH="${BIN_DIR}:${CONDA_BASE}/envs/${ENV_NAME}/bin:\$PATH"
EOF
fi

echo "==> Splice tools smoke"
"${BIN_DIR}/regtools-neoag" junctions extract -h | head -8 || true
if [[ -x "${BIN_DIR}/pvacsplice-neoag" ]]; then
  "${BIN_DIR}/pvacsplice-neoag" --help | head -8 || true
fi
if command -v "${BIN_DIR}/snaf-neoag" >/dev/null 2>&1; then
  "${CONDA_BASE}/bin/conda" run -n "${ENV_NAME}" python - <<'PY'
import snaf
print("SNAF import OK")
PY
fi

echo "==> Done. Run: bash scripts/run_splice_tool_smoke.sh"
