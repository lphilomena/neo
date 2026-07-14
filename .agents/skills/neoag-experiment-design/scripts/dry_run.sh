#!/usr/bin/env bash
set -euo pipefail
python -m neoag_v03.skill_taxonomy.cli run neoag-experiment-design --outdir "${1:-work/neoag-experiment-design}" --dry-run
