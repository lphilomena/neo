#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run neoag-patient-report --outdir "${1:-work/neoag-patient-report}" --dry-run
