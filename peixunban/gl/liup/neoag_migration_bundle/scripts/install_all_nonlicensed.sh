#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${SCRIPT_DIR}/install_tier2_tools.sh"
"${SCRIPT_DIR}/install_gatk_vep_facets_ascat.sh"
"${SCRIPT_DIR}/install_fusion_tools.sh"
"${SCRIPT_DIR}/install_optional_immuno_hla_tools.sh"
"${SCRIPT_DIR}/run_full_doctor.sh"
