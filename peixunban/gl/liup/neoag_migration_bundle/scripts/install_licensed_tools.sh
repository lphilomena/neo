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
if [[ -n "${NETMHCPAN_TARBALL:-}" ]]; then
  log "installing NetMHCpan from NETMHCPAN_TARBALL=${NETMHCPAN_TARBALL}"
  bash scripts/install_netmhcpan.sh "$NETMHCPAN_TARBALL" | tee "${TOOLS_ROOT}/logs/install_netmhcpan.log"
else
  warn "NETMHCPAN_TARBALL not set; skipping NetMHCpan. Set it to an official licensed tarball path."
fi
if [[ -n "${NETMHCSTABPAN_TARBALL:-}" ]]; then
  log "installing NetMHCstabpan from NETMHCSTABPAN_TARBALL=${NETMHCSTABPAN_TARBALL}"
  bash scripts/install_netmhcstabpan.sh "$NETMHCSTABPAN_TARBALL" | tee "${TOOLS_ROOT}/logs/install_netmhcstabpan.log"
elif [[ "${INSTALL_NETMHCSTABPAN_IEDB:-0}" == "1" ]]; then
  log "installing NetMHCstabpan IEDB shim"
  bash scripts/install_netmhcstabpan.sh --iedb | tee "${TOOLS_ROOT}/logs/install_netmhcstabpan.log"
else
  warn "NETMHCSTABPAN_TARBALL not set and INSTALL_NETMHCSTABPAN_IEDB!=1; skipping NetMHCstabpan."
fi
if [[ -n "${POLYSOLVER_PACKAGE:-}" ]]; then
  log "installing Polysolver from POLYSOLVER_PACKAGE=${POLYSOLVER_PACKAGE}"
  bash scripts/install_polysolver.sh "$POLYSOLVER_PACKAGE" | tee "${TOOLS_ROOT}/logs/install_polysolver.log"
else
  warn "POLYSOLVER_PACKAGE not set; skipping Polysolver."
fi
log "licensed tool script finished; verify license compliance before redistributing any binaries"
