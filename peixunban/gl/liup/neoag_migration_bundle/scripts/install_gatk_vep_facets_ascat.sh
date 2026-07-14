#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"
"${SCRIPT_DIR}/install_tier2_tools.sh"
activate_paths
cd "$PROJECT_ROOT"
export NEOAG_CONDA_BASE
export NEOAG_TOOLS_ROOT
log "installing GATK"
bash scripts/install_gatk.sh | tee "${TOOLS_ROOT}/logs/install_gatk.log"
log "installing VEP environment; set INSTALL_VEP_CACHE=1 to also download cache"
if [[ "${INSTALL_VEP_CACHE:-0}" == "1" ]]; then
  bash scripts/install_vep.sh --cache | tee "${TOOLS_ROOT}/logs/install_vep.log"
else
  bash scripts/install_vep.sh | tee "${TOOLS_ROOT}/logs/install_vep.log"
fi
log "installing FACETS"
bash scripts/install_facets.sh | tee "${TOOLS_ROOT}/logs/install_facets.log"
log "installing ASCAT/PyClone-VI"
bash scripts/install_ascat_pyclone.sh | tee "${TOOLS_ROOT}/logs/install_ascat_pyclone.log"
log "GATK/VEP/FACETS/ASCAT installation commands completed"
