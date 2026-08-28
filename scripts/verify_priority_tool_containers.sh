#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
cd "$REPO_ROOT"
checks=(
  scripts/verify_netmhcpan_container.sh
  scripts/verify_netmhcstabpan_container.sh
  scripts/verify_hla_la_container.sh
  scripts/verify_spechla_container.sh
  scripts/verify_purple_suite_container.sh
  scripts/verify_easyfuse_container.sh
)
status=0
for check in "${checks[@]}"; do
  echo "===== $check"
  if ! bash "$check"; then
    status=1
  fi
  echo
done
exit "$status"
