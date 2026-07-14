#!/usr/bin/env bash
set -euo pipefail
BUNDLE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-${BUNDLE_ROOT}/project/neo-na0707_upload_release-3}"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-${BUNDLE_ROOT}/../test_env_tools}"
PYTHON_BIN="${PYTHON_BIN:-${TOOLS_ROOT}/miniforge3/envs/neoag-core/bin/python}"
export NEOAG_BUNDLE_ROOT="$BUNDLE_ROOT"
export NEOAG_PROJECT_ROOT="$PROJECT_ROOT"
export NEOAG_TOOLS_ROOT="$TOOLS_ROOT"
export NEOAG_CONDA_BASE="${NEOAG_CONDA_BASE:-${TOOLS_ROOT}/miniforge3}"
export PYTHONPATH="${PROJECT_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
export PATH="${PROJECT_ROOT}/bin:${NEOAG_CONDA_BASE}/envs/neoag-core/bin:${NEOAG_CONDA_BASE}/envs/neoag-tools/bin:${NEOAG_CONDA_BASE}/bin:${PATH}"
export LD_LIBRARY_PATH="${NEOAG_CONDA_BASE}/envs/neoag-tools/lib:${NEOAG_CONDA_BASE}/envs/neoag-core/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
if [[ ! -d "$PROJECT_ROOT/src/neoag_v03" ]]; then
  mkdir -p "${BUNDLE_ROOT}/project"
  tar -xzf "${BUNDLE_ROOT}/release/neo-na0707_upload_release-3.tar.gz" -C "${BUNDLE_ROOT}/project"
fi
"$PYTHON_BIN" -m pip install -e "$PROJECT_ROOT" >/tmp/neoag_bundle_doctor_pip_install.log 2>&1
neoag-doctor \
  --project-root "$PROJECT_ROOT" \
  --tools-manifest "${BUNDLE_ROOT}/configs/local/tools_manifest.yaml" \
  --reference-manifest "${BUNDLE_ROOT}/refs/reference_manifest.test_refs.yaml" \
  --sample-manifest "${BUNDLE_ROOT}/configs/local/sample_manifest.yaml" \
  --outdir "${BUNDLE_ROOT}/doctor_bundle_test" \
  --skip-release-audit \
  --dry-run
if [[ -s "${BUNDLE_ROOT}/doctor_bundle_test/blocking_issues.tsv" ]] && [[ $(wc -l < "${BUNDLE_ROOT}/doctor_bundle_test/blocking_issues.tsv") -gt 1 ]]; then
  echo "Blocking issues remain:" >&2
  cat "${BUNDLE_ROOT}/doctor_bundle_test/blocking_issues.tsv" >&2
  exit 2
fi
echo "Doctor bundle test passed: no blocking issues. See ${BUNDLE_ROOT}/doctor_bundle_test/doctor_summary.md"
