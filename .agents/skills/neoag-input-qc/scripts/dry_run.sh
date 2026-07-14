#!/usr/bin/env bash
set -euo pipefail
python -m neoag_v03.skill_taxonomy.cli run neoag-input-qc --outdir "${1:-work/neoag-input-qc}" --dry-run
