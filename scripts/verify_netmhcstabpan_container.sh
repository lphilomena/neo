#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
[[ -f "$REPO_ROOT/conf/tools.env.sh" ]] && source "$REPO_ROOT/conf/tools.env.sh"
NETMHCSTABPAN_HOME=${NETMHCSTABPAN_HOME:-$REPO_ROOT/tools/netMHCstabpan}
if [[ ! -x "$NETMHCSTABPAN_HOME/netMHCstabpan" ]]; then
  echo "WARN: NetMHCstabpan executable is not installed: $NETMHCSTABPAN_HOME/netMHCstabpan"
  exit 0
fi
set +e
"$REPO_ROOT/scripts/run_netmhcstabpan_container.sh" -h >/tmp/netmhcstabpan_container_verify.out 2>&1
rc=$?
set -e
if grep -qi "netMHCstabpan\|usage\|IEDB" /tmp/netmhcstabpan_container_verify.out; then
  echo "PASS: NetMHCstabpan container runtime starts"
else
  echo "WARN: NetMHCstabpan container started with rc=$rc but help text was not recognized"
  head -40 /tmp/netmhcstabpan_container_verify.out
fi
