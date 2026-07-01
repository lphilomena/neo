#!/usr/bin/env bash
# Fix PRIME temp/test ownership when directories were created as root (e.g. docker install).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PRIME_DIR="${ROOT}/tools/prime"

mkdir -p "${PRIME_DIR}/lib/temp"

if [[ -w "${PRIME_DIR}/temp" ]] 2>/dev/null; then
  echo "OK: ${PRIME_DIR}/temp is already writable"
else
  echo "Attempting chown on legacy temp/test (may need sudo password)..."
  sudo chown -R "$(id -un):$(id -gn)" "${PRIME_DIR}/temp" "${PRIME_DIR}/test" 2>/dev/null || true
fi

if [[ ! -w "${PRIME_DIR}/temp" ]] 2>/dev/null; then
  echo "Note: legacy tools/prime/temp is not writable; neoag uses tools/prime/lib/temp instead."
fi

echo "Smoke test:"
# shellcheck disable=SC1091
source "${ROOT}/conf/tools.env.sh"
"${NEOAG_PRIME_BIN}" -i "${PRIME_DIR}/test/test.txt" -o /tmp/prime_fix_smoke.tsv \
  -a A0101,A2501,B0801,B1801 -mix "${MIXMHCPRED_BIN}"
head -3 /tmp/prime_fix_smoke.tsv
echo "OK"
