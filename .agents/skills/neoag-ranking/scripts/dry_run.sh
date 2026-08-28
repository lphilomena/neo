#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run neoag-ranking --outdir "${1:-work/neoag-ranking}" --dry-run
