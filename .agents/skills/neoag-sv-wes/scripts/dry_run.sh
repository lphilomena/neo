#!/usr/bin/env bash
set -euo pipefail
python -m neoag_v03.skill_taxonomy.cli run neoag-sv-wes --outdir "${1:-work/neoag-sv-wes}" --dry-run
