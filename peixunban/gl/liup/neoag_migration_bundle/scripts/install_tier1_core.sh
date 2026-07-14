#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"
ensure_project
ensure_conda
mkdir -p "$TOOLS_ROOT" "${TOOLS_ROOT}/logs" "${TOOLS_ROOT}/nextflow"
if ! conda env list | grep -q "${NEOAG_CONDA_BASE}/envs/${NEOAG_CORE_ENV}"; then
  log "creating ${NEOAG_CORE_ENV} with Python 3.11 and OpenJDK 17"
  conda create -y -n "${NEOAG_CORE_ENV}" -c conda-forge python=3.11 pip setuptools wheel pytest pytest-timeout openjdk=17 > "${TOOLS_ROOT}/logs/create_${NEOAG_CORE_ENV}.log" 2>&1
else
  log "${NEOAG_CORE_ENV} already exists"
fi
activate_paths
python -m pip install -e "$PROJECT_ROOT" > "${TOOLS_ROOT}/logs/pip_install_project_editable.log" 2>&1
find "$PROJECT_ROOT/bin" -maxdepth 1 -type f -exec chmod +x {} \;
python --version
java -version
neoag-v03 --help >/dev/null
log "core entrypoints OK"
