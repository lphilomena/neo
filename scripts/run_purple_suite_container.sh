#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
[[ -f "$REPO_ROOT/conf/tools.env.sh" ]] && source "$REPO_ROOT/conf/tools.env.sh"
IMAGE=${NEOAG_PURPLE_IMAGE:-neoag-purple-suite:ubuntu22.04}
HMFTOOLS_HOME=${HMFTOOLS_HOME:-$REPO_ROOT/tools/hmftools}
TOOL=${1:-}
[[ -n "$TOOL" && "$TOOL" != -h && "$TOOL" != --help ]] || { cat <<USAGE
Usage: $0 purple|amber|cobalt [tool args]

Set jar paths when they are outside HMFTOOLS_HOME:
  PURPLE_JAR=/path/purple.jar
  AMBER_JAR=/path/amber.jar
  COBALT_JAR=/path/cobalt.jar
USAGE
exit 0; }
shift
case "$TOOL" in
  purple) JAR=${PURPLE_JAR:-$(find "$HMFTOOLS_HOME" -maxdepth 6 -iname "*purple*.jar" 2>/dev/null | head -1)} ;;
  amber) JAR=${AMBER_JAR:-$(find "$HMFTOOLS_HOME" -maxdepth 6 -iname "*amber*.jar" 2>/dev/null | head -1)} ;;
  cobalt) JAR=${COBALT_JAR:-$(find "$HMFTOOLS_HOME" -maxdepth 6 -iname "*cobalt*.jar" 2>/dev/null | head -1)} ;;
  *) echo "ERROR: tool must be purple, amber, or cobalt" >&2; exit 2 ;;
esac
[[ -n "$JAR" && -f "$JAR" ]] || { echo "ERROR: $TOOL jar missing. Set ${TOOL^^}_JAR or HMFTOOLS_HOME." >&2; exit 2; }
docker image inspect "$IMAGE" >/dev/null 2>&1 || { echo "ERROR: build image first: $REPO_ROOT/scripts/build_priority_tool_containers.sh purple-suite" >&2; exit 127; }
mounts=(-v "$JAR:$JAR:ro" -v "$PWD:$PWD:rw" -v "$REPO_ROOT:$REPO_ROOT:rw")
[[ -d "$HMFTOOLS_HOME" ]] && mounts+=( -v "$HMFTOOLS_HOME:$HMFTOOLS_HOME:ro" )
[[ -d /mnt ]] && mounts+=( -v /mnt:/mnt:rw )
CMD="exec java ${JAVA_OPTS:--Xmx16g} -jar \"$JAR\" \"\$@\""
docker run --rm --user "$(id -u):$(id -g)" --workdir "$PWD" "${mounts[@]}" "$IMAGE" "$CMD" -- "$@"
