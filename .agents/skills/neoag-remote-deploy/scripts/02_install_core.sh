#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="$(pwd)"
PYTHON_BIN="${PYTHON:-python}"
SKIP_INSTALL=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root) PROJECT_ROOT="$2"; shift 2 ;;
    --python) PYTHON_BIN="$2"; shift 2 ;;
    --skip-install) SKIP_INSTALL=1; shift ;;
    -h|--help) echo "Usage: $0 --project-root DIR [--python PATH] [--skip-install]"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done
cd "$PROJECT_ROOT"
"$PYTHON_BIN" - <<'PY'
import sys
if sys.version_info < (3, 11):
    raise SystemExit('PYTHON_ENV_FAILED: Python >=3.11 required')
print(sys.executable, sys.version.split()[0])
PY
if [[ "$SKIP_INSTALL" != "1" ]]; then
  "$PYTHON_BIN" -m pip install -U pip setuptools wheel
  "$PYTHON_BIN" -m pip install -e '.[test]'
fi
echo "core_install=done"
