#!/usr/bin/env bash
# Environment acceptance: all tools in neoag registry must resolve (exit 0 = OK).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "${ROOT}/conf/tools.env.sh"
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
python -m neoag_v03.cli check-tools "$@"
