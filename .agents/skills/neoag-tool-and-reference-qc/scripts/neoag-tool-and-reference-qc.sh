#!/usr/bin/env bash
set -euo pipefail
: "${NEOAG_PYTHON:=python3}"
"${NEOAG_PYTHON}" -m neoag.agent_skills.tool_reference_qc "$@"
