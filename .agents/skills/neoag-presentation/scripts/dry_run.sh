#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run neoag-presentation --outdir "${1:-work/neoag-presentation}" --dry-run
