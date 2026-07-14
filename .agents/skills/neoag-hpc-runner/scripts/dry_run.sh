#!/usr/bin/env bash
set -euo pipefail
python -m neoag_v03.skill_taxonomy.cli run neoag-hpc-runner --outdir "${1:-work/neoag-hpc-runner}" --dry-run
