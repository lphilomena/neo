#!/usr/bin/env bash
set -euo pipefail
python -m neoag_v03.skill_taxonomy.cli run neoag-technical-report --outdir "${1:-work/neoag-technical-report}" --dry-run
