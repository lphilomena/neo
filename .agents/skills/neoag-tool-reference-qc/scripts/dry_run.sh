#!/usr/bin/env bash
set -euo pipefail
python -m neoag.skill_taxonomy.cli run neoag-tool-reference-qc --outdir "${1:-work/neoag-tool-reference-qc}" --dry-run
