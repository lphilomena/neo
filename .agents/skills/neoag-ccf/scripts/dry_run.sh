#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run neoag-ccf --outdir "${1:-work/neoag-ccf}" --dry-run
