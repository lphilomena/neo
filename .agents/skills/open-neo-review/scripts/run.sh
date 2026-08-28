#!/usr/bin/env bash
set -euo pipefail
PYTHON_BIN="${NEOAG_PYTHON:-python3}"
exec "$PYTHON_BIN" -m neoag.open_neo.cli review "$@"
