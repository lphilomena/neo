#!/usr/bin/env bash
# Pull official images used as optional tool backends.
set -euo pipefail
docker pull griffithlab/pvactools:6.1.0
echo "Set: export NEOAG_PVAC_DOCKER=griffithlab/pvactools:6.1.0"
