#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run open-neo-run --outdir "${1:-work/open-neo-run}" --dry-run
