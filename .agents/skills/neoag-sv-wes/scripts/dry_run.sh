#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run neoag-sv-wes --outdir "${1:-work/neoag-sv-wes}" --dry-run
