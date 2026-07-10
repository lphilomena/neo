#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
[[ -f "$REPO_ROOT/conf/tools.env.sh" ]] && source "$REPO_ROOT/conf/tools.env.sh"
HLALA_HOME=${HLALA_HOME:-${HLA_LA_HOME:-}}
HLALA_BIN=${HLALA_BIN:-${HLA_LA_BIN:-$HLALA_HOME/bin/HLA-LA.pl}}
[[ -n "$HLALA_HOME" ]] || { echo "ERROR: set HLALA_HOME or HLA_LA_HOME to your HLA-LA installation" >&2; exit 2; }
[[ -x "$HLALA_BIN" ]] || HLALA_BIN=${HLALA_HOME}/HLA-LA.pl
[[ -x "$HLALA_BIN" ]] || HLALA_BIN=${HLALA_HOME}/run_HLA-LA.pl
GRAPH=${HLALA_GRAPH:-${HLA_LA_GRAPH:-$HLALA_HOME/graphs/PRG_MHC_GRCh38_withIMGT}}
IMAGE=${NEOAG_HLALA_IMAGE:-neoag-hla-la:ubuntu22.04}
[[ ${1:-} == -h || ${1:-} == --help ]] && { echo "Usage: $0 [HLA-LA args]"; echo "Set HLALA_HOME, HLALA_BIN, HLALA_GRAPH as needed."; exit 0; }
[[ -e "$HLALA_BIN" ]] || { echo "ERROR: HLA-LA executable missing. Set HLALA_BIN. Tried: $HLALA_BIN" >&2; exit 2; }
docker image inspect "$IMAGE" >/dev/null 2>&1 || { echo "ERROR: build image first: $REPO_ROOT/scripts/build_priority_tool_containers.sh hla-la" >&2; exit 127; }
mounts=(-v "$PWD:$PWD:rw" -v "$REPO_ROOT:$REPO_ROOT:rw")
[[ -d "$HLALA_HOME" ]] && mounts+=( -v "$HLALA_HOME:$HLALA_HOME:ro" )
[[ -d "$GRAPH" ]] && mounts+=( -v "$GRAPH:$GRAPH:ro" )
[[ -d /mnt ]] && mounts+=( -v /mnt:/mnt:rw )
CMD="export HLALA_HOME=\"$HLALA_HOME\"; export HLA_LA_HOME=\"$HLALA_HOME\"; export HLALA_GRAPH=\"$GRAPH\"; exec \"$HLALA_BIN\" \"\$@\""
docker run --rm --user "$(id -u):$(id -g)" --workdir "$PWD" "${mounts[@]}" "$IMAGE" "$CMD" -- "$@"
