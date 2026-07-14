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
log "installing fusion tools environment and wrappers"
bash scripts/install_fusion_tools.sh | tee "${TOOLS_ROOT}/logs/install_fusion_tools.log"
log "fusion tool command completed; CTAT/EasyFuse references must still be staged separately"
