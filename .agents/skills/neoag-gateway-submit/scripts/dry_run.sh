#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run neoag-gateway-submit --outdir "${1:-work/neoag-gateway-submit}" --dry-run
