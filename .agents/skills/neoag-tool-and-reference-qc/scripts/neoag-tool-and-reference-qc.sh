#!/usr/bin/env bash
set -euo pipefail
: "${NEOAG_PYTHON:=python3}"
"${NEOAG_PYTHON}" -m neoag_v03.agent_skills.tool_reference_qc "$@"
