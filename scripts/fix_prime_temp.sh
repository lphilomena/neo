#!/usr/bin/env bash
# Create and repair PRIME's real runtime temp directory: $PRIME_HOME/temp.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ -f "${ROOT}/conf/tools.env.sh" ]]; then
  # shellcheck disable=SC1091
  source "${ROOT}/conf/tools.env.sh"
fi

PRIME_DIR="${PRIME_HOME:-${NEOAG_TOOLS_ROOT:-${ROOT}}/tools/prime}"
PRIME_TEMP_DIR="${PRIME_DIR}/temp"
RUNTIME_USER="${NEOAG_PRIME_RUNTIME_USER:-${SUDO_USER:-$(id -un)}}"
RUNTIME_GROUP="$(id -gn "${RUNTIME_USER}" 2>/dev/null || id -gn)"

if [[ ! -d "${PRIME_DIR}" ]]; then
  echo "ERROR: PRIME_HOME does not exist: ${PRIME_DIR}" >&2
  exit 1
fi

if ! mkdir -p "${PRIME_TEMP_DIR}" 2>/dev/null; then
  if [[ "${EUID}" -eq 0 ]]; then
    mkdir -p "${PRIME_TEMP_DIR}"
  elif command -v sudo >/dev/null 2>&1; then
    sudo mkdir -p "${PRIME_TEMP_DIR}"
  else
    echo "ERROR: cannot create ${PRIME_TEMP_DIR}; rerun as its owner or with sudo." >&2
    exit 1
  fi
fi

if [[ "${EUID}" -eq 0 ]]; then
  chown -R "${RUNTIME_USER}:${RUNTIME_GROUP}" "${PRIME_TEMP_DIR}"
elif [[ ! -w "${PRIME_TEMP_DIR}" ]] && command -v sudo >/dev/null 2>&1; then
  sudo chown -R "${RUNTIME_USER}:${RUNTIME_GROUP}" "${PRIME_TEMP_DIR}"
fi
if ! chmod u+rwx "${PRIME_TEMP_DIR}" 2>/dev/null; then
  if command -v sudo >/dev/null 2>&1; then
    sudo chmod u+rwx "${PRIME_TEMP_DIR}"
  else
    echo "ERROR: cannot make ${PRIME_TEMP_DIR} writable; rerun as its owner or with sudo." >&2
    exit 1
  fi
fi

if [[ ! -d "${PRIME_TEMP_DIR}" || ! -w "${PRIME_TEMP_DIR}" ]]; then
  echo "ERROR: PRIME runtime temp directory is not writable: ${PRIME_TEMP_DIR}" >&2
  echo "Set NEOAG_PRIME_RUNTIME_USER to the account that runs neoag and rerun this script." >&2
  exit 1
fi

echo "OK: PRIME runtime temp is writable: ${PRIME_TEMP_DIR}"

echo "Smoke test:"
rm -f /tmp/prime_fix_smoke.tsv
"${NEOAG_PRIME_BIN}" -i "${PRIME_DIR}/test/test.txt" -o /tmp/prime_fix_smoke.tsv \
  -a A0101,A2501,B0801,B1801 -mix "${MIXMHCPRED_BIN}"
test -s /tmp/prime_fix_smoke.tsv
head -3 /tmp/prime_fix_smoke.tsv
echo "OK"
