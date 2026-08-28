#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run neoag-appm-escape --outdir "${1:-work/neoag-appm-escape}" --dry-run
