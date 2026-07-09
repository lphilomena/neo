#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
[[ -f "$REPO_ROOT/conf/tools.env.sh" ]] && source "$REPO_ROOT/conf/tools.env.sh"
IMAGE=${NEOAG_EASYFUSE_IMAGE:-neoag-easyfuse:ubuntu22.04}
NEOAG_EASYFUSE_HOME=${NEOAG_EASYFUSE_HOME:-$REPO_ROOT/tools/EasyFuse}
NEOAG_EASYFUSE_REF=${NEOAG_EASYFUSE_REF:-}
docker image inspect "$IMAGE" >/dev/null 2>&1 || { echo "WARN: EasyFuse image missing; build with scripts/build_priority_tool_containers.sh easyfuse"; exit 0; }
docker run --rm "$IMAGE" "java -version >/dev/null && python3 --version && echo PASS: EasyFuse container base runtime starts"
[[ -f "$NEOAG_EASYFUSE_HOME/main.nf" ]] && echo "PASS: EasyFuse main.nf exists: $NEOAG_EASYFUSE_HOME/main.nf" || echo "WARN: NEOAG_EASYFUSE_HOME/main.nf missing: $NEOAG_EASYFUSE_HOME"
[[ -n "$NEOAG_EASYFUSE_REF" && -d "$NEOAG_EASYFUSE_REF" ]] && echo "PASS: EasyFuse ref exists: $NEOAG_EASYFUSE_REF" || echo "WARN: NEOAG_EASYFUSE_REF missing: $NEOAG_EASYFUSE_REF"
