#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run neoag-doctor --outdir "${1:-work/neoag-doctor}" --dry-run
