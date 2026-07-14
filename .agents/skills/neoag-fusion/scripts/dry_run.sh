#!/usr/bin/env bash
set -euo pipefail
python -m neoag_v03.skill_taxonomy.cli run neoag-fusion --outdir "${1:-work/neoag-fusion}" --dry-run
