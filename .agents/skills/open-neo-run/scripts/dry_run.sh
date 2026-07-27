#!/usr/bin/env bash
set -euo pipefail
PYTHON_BIN="${NEOAG_PYTHON:-python3}"
OUTDIR="work/open-neo-run"
if [[ $# -gt 0 && "$1" != -* ]]; then
  OUTDIR="$1"
  shift
fi
exec "$PYTHON_BIN" -m neoag.open_neo.cli run --mode dry-run --outdir "$OUTDIR" "$@"
