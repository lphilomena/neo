#!/usr/bin/env bash
# Install GATK4 into conda env neoag-gatk and append PATH hints to conf/tools.env.sh.
#
# Usage:
#   bash scripts/install_gatk.sh
#   source conf/tools.env.sh
#   gatk --help
#   neoag-v03 check-tools | grep gatk
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ENV_NAME="${NEOAG_GATK_ENV:-neoag-gatk}"
YML="${ROOT}/conda/env.neoag-gatk.yml"
TOOLS_ENV="${ROOT}/conf/tools.env.sh"
CONDA_BASE="${NEOAG_CONDA_BASE:-$(command conda info --base)}"
GATK_BIN="${CONDA_BASE}/envs/${ENV_NAME}/bin"

echo "==> Installing GATK4 conda env: ${ENV_NAME}"

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

if conda_safe env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "==> Updating existing env ${ENV_NAME} ..."
  conda_safe env update -n "${ENV_NAME}" -f "${YML}" --prune
else
  echo "==> Creating env from ${YML} (may take several minutes) ..."
  conda_safe env create -n "${ENV_NAME}" -f "${YML}"
fi

if [[ ! -x "${GATK_BIN}/gatk" ]]; then
  echo "ERROR: gatk not found at ${GATK_BIN}/gatk" >&2
  exit 1
fi

echo "==> Smoke test:"
"${GATK_BIN}/gatk" --help | head -5

if [[ -f "${TOOLS_ENV}" ]] && ! grep -q 'NEOAG_GATK_ENV' "${TOOLS_ENV}"; then
  cat >> "${TOOLS_ENV}" <<EOF

# GATK4 (Mutect2 / FilterMutectCalls) — installed via scripts/install_gatk.sh
export NEOAG_GATK_ENV="${ENV_NAME}"
if [[ -d "${GATK_BIN}" ]]; then
  export PATH="${GATK_BIN}:\${PATH}"
fi
EOF
  echo "==> Appended NEOAG_GATK_ENV block to conf/tools.env.sh"
else
  echo "==> conf/tools.env.sh already references NEOAG_GATK_ENV (edit PATH manually if needed)"
fi

echo "==> Done. Run: source conf/tools.env.sh && neoag-v03 check-tools"
