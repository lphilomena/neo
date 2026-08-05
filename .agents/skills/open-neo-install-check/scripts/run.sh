#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../../.." && pwd)"

python_is_compatible() {
  "$1" -c 'import sys, tomllib; raise SystemExit(0 if sys.version_info >= (3, 11) else 1)' \
    >/dev/null 2>&1
}

resolve_python() {
  local candidate resolved
  if [[ -n "${NEOAG_PYTHON:-}" ]]; then
    if python_is_compatible "$NEOAG_PYTHON"; then
      printf '%s\n' "$NEOAG_PYTHON"
      return 0
    fi
    echo "OPEN_NEO_PYTHON_UNSUPPORTED: NEOAG_PYTHON must provide Python 3.11+ with tomllib: $NEOAG_PYTHON" >&2
    return 31
  fi

  candidates=()
  if [[ -n "${NEOAG_CONDA_BASE:-}" ]]; then
    candidates+=("$NEOAG_CONDA_BASE/envs/neoag-tools/bin/python" "$NEOAG_CONDA_BASE/bin/python")
  fi
  if [[ -n "${NEOAG_TOOLS_ROOT:-}" ]]; then
    candidates+=("$NEOAG_TOOLS_ROOT/miniforge3/envs/neoag-tools/bin/python" "$NEOAG_TOOLS_ROOT/miniforge3/bin/python")
  fi
  if [[ -n "${NEOAG_DEPLOY_ROOT:-}" ]]; then
    candidates+=("$NEOAG_DEPLOY_ROOT/envs/miniforge3/envs/neoag-tools/bin/python" "$NEOAG_DEPLOY_ROOT/envs/miniforge3/bin/python")
  fi
  candidates+=(
    "$PROJECT_ROOT/.venv/bin/python"
    "$PROJECT_ROOT/../../envs/miniforge3/envs/neoag-tools/bin/python"
    "$PROJECT_ROOT/../../envs/miniforge3/bin/python"
    python3.13 python3.12 python3.11 python3 python
  )
  for candidate in "${candidates[@]}"; do
    if [[ "$candidate" == */* ]]; then
      [[ -x "$candidate" ]] || continue
      resolved="$candidate"
    else
      resolved="$(command -v "$candidate" 2>/dev/null || true)"
      [[ -n "$resolved" ]] || continue
    fi
    if python_is_compatible "$resolved"; then
      printf '%s\n' "$resolved"
      return 0
    fi
  done
  echo "OPEN_NEO_PYTHON_NOT_FOUND: install Python 3.11+ or set NEOAG_PYTHON/NEOAG_CONDA_BASE" >&2
  return 31
}

PYTHON_BIN="$(resolve_python)"
export PYTHONPATH="$PROJECT_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
cd "$PROJECT_ROOT"
exec "$PYTHON_BIN" -m neoag.open_neo.cli install-check "$@"
