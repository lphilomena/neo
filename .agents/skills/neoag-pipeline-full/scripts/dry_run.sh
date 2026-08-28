#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run neoag-pipeline-full --outdir "${1:-work/neoag-pipeline-full}" --dry-run
