#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run neoag-expression --outdir "${1:-work/neoag-expression}" --dry-run
