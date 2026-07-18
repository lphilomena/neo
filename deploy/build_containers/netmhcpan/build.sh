#!/bin/bash
# ============================================================
# Build NetMHCpan 4.2c Docker image
#
# Prerequisites:
#   1. DTU academic license for NetMHCpan
#   2. Tarball placed at vendor/netMHCpan-4.2c.Linux.tar.gz
#      (the tarball is NOT distributed with this repo)
#
# Usage:
#   bash deploy/build_containers/netmhcpan/build.sh
#   bash deploy/build_containers/netmhcpan/build.sh --test   # build + quick smoke test
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
IMAGE="neoag-netmhcpan:4.2c-ubuntu22.04"
TARBALL="${PROJECT_ROOT}/vendor/netMHCpan-4.2c.Linux.tar.gz"

echo "=== Build NetMHCpan 4.2c Docker image ==="
echo "Image:  ${IMAGE}"
echo ""

# --- Check prerequisites ---
if [[ ! -f "${TARBALL}" ]]; then
  echo "ERROR: Tarball not found at ${TARBALL}" >&2
  echo "" >&2
  echo "NetMHCpan requires a DTU academic license.  To obtain the tarball:" >&2
  echo "  1. Register at https://services.healthtech.dtu.dk/cgi-bin/request.cgi?tool_id=NetMHCpan" >&2
  echo "  2. Download netMHCpan-4.2c.Linux.tar.gz from the link in your email" >&2
  echo "  3. Copy it to ${PROJECT_ROOT}/vendor/" >&2
  exit 1
fi

echo "Tarball: ${TARBALL}"
echo "Size:    $(du -h "${TARBALL}" | cut -f1)"
echo ""

# --- Build ---
echo ">>> Building ${IMAGE} ..."
cd "${PROJECT_ROOT}"
docker build \
  -t "${IMAGE}" \
  -f "${SCRIPT_DIR}/Dockerfile" \
  .

echo ""
echo ">>> Image built:"
docker image ls "${IMAGE}"

# --- Optional smoke test ---
if [[ "${1:-}" == "--test" ]]; then
  echo ""
  echo ">>> Running smoke test ..."
  docker run --rm "${IMAGE}" netMHCpan -h 2>&1 | head -10
  echo ""
  echo "=== Build + smoke test OK ==="
else
  echo ""
  echo "=== Build complete ==="
  echo "Run with --test for a quick smoke test after build."
fi
