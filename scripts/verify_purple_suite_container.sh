#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
[[ -f "$REPO_ROOT/conf/tools.env.sh" ]] && source "$REPO_ROOT/conf/tools.env.sh"
IMAGE=${NEOAG_PURPLE_IMAGE:-neoag-purple-suite:ubuntu22.04}
HMFTOOLS_HOME=${HMFTOOLS_HOME:-$REPO_ROOT/tools/HMFTOOLS}
PURPLE_JAR=${PURPLE_JAR:-$(find "$HMFTOOLS_HOME" -maxdepth 6 -iname "*purple*.jar" 2>/dev/null | head -1)}
AMBER_JAR=${AMBER_JAR:-$(find "$HMFTOOLS_HOME" -maxdepth 6 -iname "*amber*.jar" 2>/dev/null | head -1)}
COBALT_JAR=${COBALT_JAR:-$(find "$HMFTOOLS_HOME" -maxdepth 6 -iname "*cobalt*.jar" 2>/dev/null | head -1)}
docker image inspect "$IMAGE" >/dev/null 2>&1 || { echo "WARN: PURPLE suite image missing; build with scripts/build_priority_tool_containers.sh purple-suite"; exit 0; }
docker run --rm "$IMAGE" "java -version >/dev/null && echo PASS: PURPLE suite Java runtime starts"
for name in PURPLE_JAR AMBER_JAR COBALT_JAR; do
  value=${!name:-}
  [[ -n "$value" && -f "$value" ]] && echo "PASS: $name=$value" || echo "WARN: $name not set or missing"
done
