#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="$(pwd)"
OUTDIR="work/remote_deploy"
PYTHON_BIN="${PYTHON:-python}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root) PROJECT_ROOT="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --python) PYTHON_BIN="$2"; shift 2 ;;
    -h|--help) echo "Usage: $0 --project-root DIR --outdir DIR"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done
cd "$PROJECT_ROOT"
mkdir -p "$OUTDIR"
export PYTHONPATH="$PROJECT_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
LOG="$OUTDIR/runtime_check.log"
{
  "$PYTHON_BIN" -m compileall -q src
  if command -v neoag >/dev/null 2>&1; then neoag --help >/dev/null; else "$PYTHON_BIN" -m neoag.cli --help >/dev/null; fi
  if command -v neoag-skill >/dev/null 2>&1; then neoag-skill --help >/dev/null; else "$PYTHON_BIN" -m neoag.skill_taxonomy.cli --help >/dev/null; fi
} >"$LOG" 2>&1 || { echo "CORE_INSTALL_FAILED: see $LOG" >&2; exit 12; }
echo "runtime_check=$LOG"
