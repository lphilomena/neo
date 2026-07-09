#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
[[ -f "$REPO_ROOT/conf/tools.env.sh" ]] && source "$REPO_ROOT/conf/tools.env.sh"
NETMHCSTABPAN_HOME=${NETMHCSTABPAN_HOME:-$REPO_ROOT/tools/netMHCstabpan}
IMAGE=${NEOAG_NETMHCSTABPAN_IMAGE:-neoag-netmhcstabpan:1.0-ubuntu22.04}
TMPDIR_HOST=${NEOAG_NETMHCSTABPAN_TMPDIR:-$REPO_ROOT/work/netmhcstabpan_tmp}
[[ ${1:-} == -h || ${1:-} == --help ]] && { echo "Usage: $0 [netMHCstabpan args]"; exit 0; }
[[ -x "$NETMHCSTABPAN_HOME/netMHCstabpan" ]] || { echo "ERROR: missing $NETMHCSTABPAN_HOME/netMHCstabpan" >&2; exit 2; }
mkdir -p "$TMPDIR_HOST"
ARGS=("$@")
for i in "${!ARGS[@]}"; do [[ ${ARGS[$i]} == HLA-*\** ]] && ARGS[$i]=${ARGS[$i]//\*/}; done
CMD="export TMPDIR=/tmp/netmhcstabpan; exec \"\$NETMHCSTABPAN_HOME/netMHCstabpan\" \"\$@\""
docker image inspect "$IMAGE" >/dev/null 2>&1 || { echo "ERROR: build image first: $REPO_ROOT/scripts/build_priority_tool_containers.sh netmhcstabpan" >&2; exit 127; }
mounts=(-v "$NETMHCSTABPAN_HOME:$NETMHCSTABPAN_HOME:ro" -v "$TMPDIR_HOST:/tmp/netmhcstabpan:rw" -v "$PWD:$PWD:rw" -v "$REPO_ROOT:$REPO_ROOT:rw")
[[ -d /mnt ]] && mounts+=( -v /mnt:/mnt:rw )
docker run --rm --user "$(id -u):$(id -g)" --workdir "$PWD" -e NETMHCSTABPAN_HOME="$NETMHCSTABPAN_HOME" "${mounts[@]}" "$IMAGE" "$CMD" -- "${ARGS[@]}"
