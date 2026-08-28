#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run neoag-splice --outdir "${1:-work/neoag-splice}" --dry-run
