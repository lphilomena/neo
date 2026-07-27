#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../../.." && pwd)
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
PYTHON_BIN=${NEOAG_PYTHON:-}
[[ -n "$PYTHON_BIN" || ! -x "${ROOT}/.venv/bin/python" ]] || PYTHON_BIN="${ROOT}/.venv/bin/python"
[[ -n "$PYTHON_BIN" ]] || PYTHON_BIN=$(command -v python3 || command -v python || true)
[[ -n "$PYTHON_BIN" ]] || { echo "ERROR: Python 3.11+ is required" >&2; exit 127; }
"$PYTHON_BIN" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 11) else 1)' \
  || { echo "ERROR: Python 3.11+ is required; set NEOAG_PYTHON" >&2; exit 2; }
exec "$PYTHON_BIN" -m neoag.workflow_selection "$@"
