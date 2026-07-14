#!/usr/bin/env bash
set -euo pipefail
python -m neoag_v03.skill_taxonomy.cli run neoag-peptide-csv --outdir "${1:-work/neoag-peptide-csv}" --dry-run
