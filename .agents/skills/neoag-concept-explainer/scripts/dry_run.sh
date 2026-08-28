#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run neoag-concept-explainer --outdir "${1:-work/neoag-concept-explainer}" --dry-run
