#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
[[ -f "$REPO_ROOT/conf/tools.env.sh" ]] && source "$REPO_ROOT/conf/tools.env.sh"
IMAGE=${NEOAG_EASYFUSE_IMAGE:-neoag-easyfuse:ubuntu22.04}
NEOAG_EASYFUSE_HOME=${NEOAG_EASYFUSE_HOME:-$REPO_ROOT/tools/EasyFuse}
NEOAG_EASYFUSE_REF=${NEOAG_EASYFUSE_REF:-}
NEOAG_CONDA_BASE=${NEOAG_CONDA_BASE:-}
[[ ${1:-} == -h || ${1:-} == --help ]] && { cat <<USAGE
Usage: $0 [command]

Default command runs scripts/run_easyfuse_sample.sh inside the container.
Required for real runs: EASYFUSE_FQ1, EASYFUSE_FQ2, SAMPLE_ID/EASYFUSE_SAMPLE_ID,
NEOAG_EASYFUSE_HOME, NEOAG_EASYFUSE_REF, and usually NEOAG_CONDA_BASE.
USAGE
exit 0; }
docker image inspect "$IMAGE" >/dev/null 2>&1 || { echo "ERROR: build image first: $REPO_ROOT/scripts/build_priority_tool_containers.sh easyfuse" >&2; exit 127; }
mounts=(-v "$PWD:$PWD:rw" -v "$REPO_ROOT:$REPO_ROOT:rw")
[[ -d "$NEOAG_EASYFUSE_HOME" ]] && mounts+=( -v "$NEOAG_EASYFUSE_HOME:$NEOAG_EASYFUSE_HOME:rw" )
[[ -n "$NEOAG_EASYFUSE_REF" && -d "$NEOAG_EASYFUSE_REF" ]] && mounts+=( -v "$NEOAG_EASYFUSE_REF:$NEOAG_EASYFUSE_REF:ro" )
[[ -n "$NEOAG_CONDA_BASE" && -d "$NEOAG_CONDA_BASE" ]] && mounts+=( -v "$NEOAG_CONDA_BASE:$NEOAG_CONDA_BASE:rw" )
[[ -d /mnt ]] && mounts+=( -v /mnt:/mnt:rw )
CMD=${1:-bash "$REPO_ROOT/scripts/run_easyfuse_sample.sh"}
if [[ $# -gt 0 ]]; then shift; fi
docker run --rm --user "$(id -u):$(id -g)" --workdir "$PWD" \
  -e NEOAG_EASYFUSE_HOME="$NEOAG_EASYFUSE_HOME" \
  -e NEOAG_EASYFUSE_REF="$NEOAG_EASYFUSE_REF" \
  -e NEOAG_CONDA_BASE="$NEOAG_CONDA_BASE" \
  -e EASYFUSE_FQ1="${EASYFUSE_FQ1:-}" -e EASYFUSE_FQ2="${EASYFUSE_FQ2:-}" \
  -e SAMPLE_ID="${SAMPLE_ID:-}" -e EASYFUSE_SAMPLE_ID="${EASYFUSE_SAMPLE_ID:-}" \
  -e OUTDIR="${OUTDIR:-}" -e LOG="${LOG:-}" \
  "${mounts[@]}" "$IMAGE" "$CMD" -- "$@"
