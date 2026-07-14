#!/usr/bin/env bash
set -euo pipefail
python -m neoag_v03.skill_taxonomy.cli run neoag-hla-typing-loh --outdir "${1:-work/neoag-hla-typing-loh}" --dry-run
