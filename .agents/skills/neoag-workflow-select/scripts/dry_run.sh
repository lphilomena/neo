#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../../.." && pwd)
OUTDIR=${1:-"${ROOT}/work/skill_smoke/neoag-workflow-select"}
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
PYTHON_BIN=${NEOAG_PYTHON:-}
[[ -n "$PYTHON_BIN" || ! -x "${ROOT}/.venv/bin/python" ]] || PYTHON_BIN="${ROOT}/.venv/bin/python"
[[ -n "$PYTHON_BIN" ]] || PYTHON_BIN=$(command -v python3 || command -v python || true)
[[ -n "$PYTHON_BIN" ]] || { echo "ERROR: Python 3.11+ is required" >&2; exit 127; }
"$PYTHON_BIN" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 11) else 1)' \
  || { echo "ERROR: Python 3.11+ is required; set NEOAG_PYTHON" >&2; exit 2; }

"$PYTHON_BIN" -m neoag.workflow_selection \
  --sample-manifest "${ROOT}/configs/controlled_execution/sample_manifest.example.yaml" \
  --tools-manifest "${ROOT}/configs/controlled_execution/tools_manifest.example.yaml" \
  --reference-manifest "${ROOT}/configs/controlled_execution/reference_manifest.example.yaml" \
  --outdir "$OUTDIR"
