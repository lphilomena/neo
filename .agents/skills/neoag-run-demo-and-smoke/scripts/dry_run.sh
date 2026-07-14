#!/usr/bin/env bash
set -euo pipefail
python -m neoag_v03.skill_taxonomy.cli run neoag-run-demo-and-smoke --outdir "${1:-work/neoag-run-demo-and-smoke}" --dry-run
