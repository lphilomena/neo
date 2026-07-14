#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=common.sh
source "${SCRIPT_DIR}/common.sh"
activate_paths
REF_MANIFEST="${REFERENCE_MANIFEST:-${BUNDLE_ROOT}/refs/reference_manifest.local_template.yaml}"
if [[ ! -f "$REF_MANIFEST" ]]; then
  REF_MANIFEST="${BUNDLE_ROOT}/refs/reference_manifest.test_refs.yaml"
  warn "production reference manifest not found; using test refs manifest: ${REF_MANIFEST}"
fi
TOOLS_MANIFEST="${TOOLS_MANIFEST:-${BUNDLE_ROOT}/configs/local/tools_manifest.yaml}"
SAMPLE_MANIFEST="${SAMPLE_MANIFEST:-${BUNDLE_ROOT}/configs/local/sample_manifest.yaml}"
OUTDIR="${DOCTOR_OUTDIR:-${BUNDLE_ROOT}/doctor_full}"
EXTRA_ARGS=()
if [[ "${SKIP_RELEASE_AUDIT:-1}" == "1" ]]; then
  EXTRA_ARGS+=(--skip-release-audit)
fi
if [[ "${DRY_RUN:-1}" == "1" ]]; then
  EXTRA_ARGS+=(--dry-run)
fi
neoag-doctor \
  --project-root "$PROJECT_ROOT" \
  --tools-manifest "$TOOLS_MANIFEST" \
  --reference-manifest "$REF_MANIFEST" \
  --sample-manifest "$SAMPLE_MANIFEST" \
  --outdir "$OUTDIR" \
  "${EXTRA_ARGS[@]}"
cat "$OUTDIR/blocking_issues.tsv"
