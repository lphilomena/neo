#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
[[ -f "$REPO_ROOT/conf/tools.env.sh" ]] && source "$REPO_ROOT/conf/tools.env.sh"

HLALA_HOME=${HLALA_HOME:-${HLA_LA_HOME:-${NEOAG_TOOLS_ROOT:-$REPO_ROOT}/tools/HLA-LA}}
HLALA_ENV_PREFIX=${HLALA_ENV_PREFIX:-${HLA_LA_ENV_PREFIX:-$HLALA_HOME/.conda}}
HLALA_BIN=${HLALA_BIN:-${HLA_LA_BIN:-}}
[[ -n "$HLALA_BIN" ]] || HLALA_BIN="$HLALA_ENV_PREFIX/bin/HLA-LA.pl"
[[ -x "$HLALA_BIN" ]] || HLALA_BIN="$HLALA_HOME/bin/HLA-LA.pl"
GRAPH=${HLALA_GRAPH:-${HLA_LA_GRAPH:-$HLALA_HOME/graphs/PRG_MHC_GRCh38_withIMGT}}
BACKEND=${NEOAG_HLALA_BACKEND:-auto}
IMAGE=${NEOAG_HLALA_IMAGE:-neoag-hla-la:ubuntu22.04}

if [[ ${1:-} == --wrapper-help ]]; then
  cat <<USAGE
Usage: $0 [HLA-LA args]

Environment:
  HLALA_HOME             HLA-LA tool home
  HLALA_ENV_PREFIX       Native Bioconda environment (default: HLALA_HOME/.conda)
  HLALA_GRAPH            Prepared PRG graph directory
  NEOAG_HLALA_BACKEND    auto, native, or container (default: auto)
USAGE
  exit 0
fi

[[ -x "$HLALA_BIN" ]] || {
  echo "ERROR: real HLA-LA executable missing: $HLALA_BIN" >&2
  echo "Install Bioconda hla-la or set HLALA_BIN/HLALA_ENV_PREFIX." >&2
  exit 2
}

if [[ "$BACKEND" == auto ]]; then
  if [[ -x "$HLALA_ENV_PREFIX/bin/HLA-LA.pl" ]]; then
    BACKEND=native
  else
    BACKEND=container
  fi
fi

case "$BACKEND" in
  native)
    export PATH="$HLALA_ENV_PREFIX/bin:$PATH"
    export HLALA_HOME HLA_LA_HOME="$HLALA_HOME" HLALA_GRAPH="$GRAPH"
    exec "$HLALA_BIN" "$@"
    ;;
  container)
    docker image inspect "$IMAGE" >/dev/null 2>&1 || {
      echo "ERROR: HLA-LA image missing: $IMAGE" >&2
      exit 127
    }
    mounts=(-v "$PWD:$PWD:rw" -v "$REPO_ROOT:$REPO_ROOT:rw")
    [[ -d "$HLALA_HOME" ]] && mounts+=( -v "$HLALA_HOME:$HLALA_HOME:ro" )
    [[ -d "$HLALA_ENV_PREFIX" ]] && mounts+=( -v "$HLALA_ENV_PREFIX:$HLALA_ENV_PREFIX:ro" )
    [[ -d "$GRAPH" ]] && mounts+=( -v "$GRAPH:$GRAPH:ro" )
    [[ -d /mnt ]] && mounts+=( -v /mnt:/mnt:rw )
    [[ -d /tmp ]] && mounts+=( -v /tmp:/tmp:rw )
    CMD="export PATH=\"$HLALA_ENV_PREFIX/bin:\$PATH\"; export HLALA_GRAPH=\"$GRAPH\"; exec \"$HLALA_BIN\" \"\$@\""
    exec docker run --rm --user "$(id -u):$(id -g)" --workdir "$PWD" "${mounts[@]}" "$IMAGE" "$CMD" -- "$@"
    ;;
  *)
    echo "ERROR: invalid NEOAG_HLALA_BACKEND=$BACKEND; expected auto, native, or container" >&2
    exit 2
    ;;
esac
