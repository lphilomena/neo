#!/usr/bin/env bash
set -euo pipefail
: "${NEOAG_PYTHON:=python3}"
PYTHONPATH="${PYTHONPATH:-src}" "${NEOAG_PYTHON}" .agents/skills/neoag-sliding-run/scripts/refresh_variant_peptides_annotated.py "$@"
