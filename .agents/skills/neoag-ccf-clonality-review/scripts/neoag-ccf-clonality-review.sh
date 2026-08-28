#!/usr/bin/env bash
set -euo pipefail
: "${NEOAG_PYTHON:=python3}"
"${NEOAG_PYTHON}" -m neoag.agent_skills.ccf_review "$@"
