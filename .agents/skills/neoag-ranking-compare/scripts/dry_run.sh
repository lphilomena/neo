#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run neoag-ranking-compare --outdir "${1:-work/neoag-ranking-compare}" --dry-run
