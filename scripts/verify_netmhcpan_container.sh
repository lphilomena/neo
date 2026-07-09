#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
TEST_DIR=${NEOAG_NETMHCPAN_TEST_DIR:-$REPO_ROOT/work/netmhcpan_container_verify}
mkdir -p "$TEST_DIR"
PEP_FILE=$(mktemp "$TEST_DIR/netmhcpan_test.XXXXXX.pep")
OUT_FILE=$(mktemp "$TEST_DIR/netmhcpan_container_verify.XXXXXX.out")
cleanup() { rm -f "$PEP_FILE" "$OUT_FILE"; }
trap cleanup EXIT
printf "AAAAAAAAL\n" > "$PEP_FILE"
"$REPO_ROOT/scripts/run_netmhcpan_container.sh" -p "$PEP_FILE" -a HLA-A02:06 >"$OUT_FILE"
if grep -q "Unable to open" "$OUT_FILE"; then
  echo "ERROR: NetMHCpan could not open the peptide file" >&2
  cat "$OUT_FILE" >&2
  exit 1
fi
if grep -q "HLA-A\*02:06\|HLA-A02:06" "$OUT_FILE"; then
  echo "PASS: NetMHCpan container runtime produced HLA-A02:06 output"
else
  echo "ERROR: NetMHCpan ran but expected allele label was not found" >&2
  head -40 "$OUT_FILE" >&2
  exit 1
fi
