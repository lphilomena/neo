#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run neoag-rna-evidence --outdir "${1:-work/neoag-rna-evidence}" --dry-run
