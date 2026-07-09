#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
[[ -f "$REPO_ROOT/conf/tools.env.sh" ]] && source "$REPO_ROOT/conf/tools.env.sh"
IMAGE=${NEOAG_SPECHLA_IMAGE:-neoag-spechla:ubuntu22.04}
docker image inspect "$IMAGE" >/dev/null 2>&1 || { echo "WARN: SpecHLA image missing; build with scripts/build_priority_tool_containers.sh spechla"; exit 0; }
docker run --rm "$IMAGE" "python3 --version && samtools --version >/dev/null && bowtie2 --version >/dev/null && echo PASS: SpecHLA container base runtime starts"
SPECHLA_HOME=${SPECHLA_HOME:-$REPO_ROOT/tools/SpecHLA}
[[ -d "$SPECHLA_HOME/script" ]] && echo "PASS: SpecHLA scripts exist: $SPECHLA_HOME/script" || echo "WARN: SpecHLA scripts missing: $SPECHLA_HOME/script"
[[ -d "$SPECHLA_HOME/db" ]] && echo "PASS: SpecHLA db exists: $SPECHLA_HOME/db" || echo "WARN: SpecHLA db missing: $SPECHLA_HOME/db"
