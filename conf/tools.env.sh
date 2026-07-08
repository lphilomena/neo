# Source before neoag upstream runs: source conf/tools.env.sh
# Optional site overlay: copy conf/tools.env.local.example.sh -> conf/tools.env.local.sh

_NEOAG_TOOLS_ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export NEOAG_PROJECT_ROOT="$(cd "${_NEOAG_TOOLS_ENV_DIR}/.." && pwd)"
export NEOAG_TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-${NEOAG_PROJECT_ROOT}}"
unset _NEOAG_TOOLS_ENV_DIR

export NEOAG_CONDA_ENV="neoag-tools"
export NEOAG_CONDA_BASE="${NEOAG_CONDA_BASE:-$(conda info --base 2>/dev/null || echo ${HOME}/miniconda3)}"

# neoag-tools stays ahead of gatk/sv/manta python shims, while this checkout
# keeps priority for neoag wrapper scripts such as bin/neoag-v03.
export PATH="${NEOAG_TOOLS_ROOT}/bin:${NEOAG_CONDA_BASE}/envs/neoag-tools/bin:${PATH}"

if [[ -f "${NEOAG_PROJECT_ROOT}/conf/tools.env.local.sh" ]]; then
  # shellcheck source=/dev/null
  source "${NEOAG_PROJECT_ROOT}/conf/tools.env.local.sh"
fi

# VEP — installed via scripts/install_vep.sh
export NEOAG_VEP_ENV="neoag-vep"
export NEOAG_VEP_BIN="/root/miniconda3/envs/neoag-vep/bin/vep"
if [[ -d "/root/miniconda3/envs/neoag-vep/bin" ]]; then
  export PATH="/root/miniconda3/envs/neoag-vep/bin:${PATH}"
fi
