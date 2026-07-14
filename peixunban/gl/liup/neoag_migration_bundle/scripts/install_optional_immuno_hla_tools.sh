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
log "installing DeepImmuno"
bash scripts/install_deepimmuno.sh | tee "${TOOLS_ROOT}/logs/install_deepimmuno.log"
log "installing PRIME/MixMHCpred/BigMHC"
bash scripts/install_immunogenicity_tools.sh | tee "${TOOLS_ROOT}/logs/install_immunogenicity_tools.log"
log "installing OptiType"
bash scripts/install_optitype.sh | tee "${TOOLS_ROOT}/logs/install_optitype.log"
log "installing LOHHLA framework"
bash scripts/install_lohhla.sh | tee "${TOOLS_ROOT}/logs/install_lohhla.log"
log "optional immunogenicity/HLA commands completed; licensed components may still require manual staging"
