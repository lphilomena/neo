#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run neoag-vcf --outdir "${1:-work/neoag-vcf}" --dry-run
