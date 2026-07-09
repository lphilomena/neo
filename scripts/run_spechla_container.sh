#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
[[ -f "$REPO_ROOT/conf/tools.env.sh" ]] && source "$REPO_ROOT/conf/tools.env.sh"
SPECHLA_HOME=${SPECHLA_HOME:-$REPO_ROOT/tools/SpecHLA}
IMAGE=${NEOAG_SPECHLA_IMAGE:-neoag-spechla:ubuntu22.04}
MODE=${SPECHLA_MODE:-auto}
[[ ${1:-} == -h || ${1:-} == --help ]] && { cat <<USAGE
Usage: $0 [SpecHLA args]

Default command selects the first available executable:
  spechla from PATH, or $SPECHLA_HOME/script/whole/SpecHLA.sh
Set SPECHLA_MODE=extract to run ExtractHLAread.sh.
Set SPECHLA_CMD=/path/to/custom_command to override.
USAGE
exit 0; }
[[ -d "$SPECHLA_HOME" ]] || { echo "ERROR: SpecHLA home missing: $SPECHLA_HOME" >&2; exit 2; }
docker image inspect "$IMAGE" >/dev/null 2>&1 || { echo "ERROR: build image first: $REPO_ROOT/scripts/build_priority_tool_containers.sh spechla" >&2; exit 127; }
mounts=(-v "$SPECHLA_HOME:$SPECHLA_HOME:rw" -v "$PWD:$PWD:rw" -v "$REPO_ROOT:$REPO_ROOT:rw")
[[ -d /mnt ]] && mounts+=( -v /mnt:/mnt:rw )
if [[ -n ${SPECHLA_CMD:-} ]]; then
  CMD="exec \"$SPECHLA_CMD\" \"\$@\""
elif [[ "$MODE" == extract ]]; then
  CMD="cd \"$SPECHLA_HOME\"; exec bash script/ExtractHLAread.sh \"\$@\""
else
  CMD="cd \"$SPECHLA_HOME\"; if command -v spechla >/dev/null 2>&1; then exec spechla \"\$@\"; elif [[ -f script/whole/SpecHLA.sh ]]; then exec bash script/whole/SpecHLA.sh \"\$@\"; else echo ERROR: SpecHLA command not found >&2; exit 127; fi"
fi
docker run --rm --user "$(id -u):$(id -g)" --workdir "$PWD" -e SPECHLA_HOME="$SPECHLA_HOME" "${mounts[@]}" "$IMAGE" "$CMD" -- "$@"
