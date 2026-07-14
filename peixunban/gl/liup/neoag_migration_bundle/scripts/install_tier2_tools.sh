#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"
"${SCRIPT_DIR}/install_tier1_core.sh"
activate_paths
if ! conda env list | grep -q "${NEOAG_CONDA_BASE}/envs/${NEOAG_TOOLS_ENV}"; then
  log "creating ${NEOAG_TOOLS_ENV} from bundled lite bioconda environment"
  conda env create -f "${BUNDLE_ROOT}/env/env.neoag-tools-lite.bioconda.yml" > "${TOOLS_ROOT}/logs/create_${NEOAG_TOOLS_ENV}.log" 2>&1
else
  log "${NEOAG_TOOLS_ENV} already exists"
fi
activate_paths
samtools --version | head -1
bcftools --version | head -1
mhcflurry-predict --help >/dev/null
pvacseq --help >/dev/null
log "tier2 lite tools OK"
