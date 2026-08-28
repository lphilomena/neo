#!/usr/bin/env bash
# Build and smoke-test the lightweight full-functionality release tarball.
#
# Bundles: source, CLI, Nextflow, tests, fixtures, docs, install scripts.
# Excludes: work/, results/, tools/, large refs, conda envs, patient data.
#
# Full real-data runs still require external NEOAG_TOOLS_ROOT (see conf/tools.env.sh).
#
# Usage:
#   bash scripts/build_lite_release.sh
#   SKIP_SMOKE=1 bash scripts/build_lite_release.sh
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUTDIR="${RELEASE_OUTDIR:-${ROOT}/work/releases}"
SMOKE_DIR="${RELEASE_SMOKE_DIR:-/tmp/neoag_lite_release_smoke}"
SKIP_SMOKE="${SKIP_SMOKE:-0}"

cd "${ROOT}"
mkdir -p "${OUTDIR}"

echo "==> packaging lite release into ${OUTDIR}"
mapfile -t _pack_lines < <(python3 scripts/package_v04_release.py --outdir "${OUTDIR}")
TARBALL="${_pack_lines[0]}"
SHA_FILE="${TARBALL}.sha256"
[[ -f "${TARBALL}" ]] || {
  echo "ERROR: tarball not created: ${TARBALL}" >&2
  exit 1
}

echo "==> tarball: ${TARBALL}"
cat "${SHA_FILE}"

if [[ "${SKIP_SMOKE}" == "1" ]]; then
  echo "==> skip smoke (SKIP_SMOKE=1)"
  exit 0
fi

rm -rf "${SMOKE_DIR}"
mkdir -p "${SMOKE_DIR}"
tar -xzf "${TARBALL}" -C "${SMOKE_DIR}"
PKG_DIR="$(find "${SMOKE_DIR}" -mindepth 1 -maxdepth 1 -type d | head -1)"
[[ -n "${PKG_DIR}" ]] || {
  echo "ERROR: unpack failed under ${SMOKE_DIR}" >&2
  exit 1
}

echo "==> smoke unpack: ${PKG_DIR}"
(
  cd "${PKG_DIR}"
  python3 -m pip install -q -e .
  python3 -m pytest -q
  neoag run-demo --outdir "${SMOKE_DIR}/demo" --sample-id LITE_RELEASE
  test -s "${SMOKE_DIR}/demo/scoring/ranked_peptides.tsv"
  test -s "${SMOKE_DIR}/demo/reports/evidence_report.html"
)

echo "==> lite release OK"
echo "    install: tar -xzf ${TARBALL} && cd <unpack> && pip install -e ."
echo "    demo:    neoag run-demo --outdir work/demo --sample-id DEMO001"
echo "    tools:   export NEOAG_TOOLS_ROOT=/path/to/artifact_bundle && source conf/tools.env.sh"
