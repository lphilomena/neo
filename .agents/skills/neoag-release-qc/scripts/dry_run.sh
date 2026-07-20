#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run neoag-release-qc --outdir "${1:-work/neoag-release-qc}" --dry-run
