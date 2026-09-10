#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
NETMHCPAN_HOME=${NETMHCPAN_HOME:-${NETMHCpan:-$REPO_ROOT/tools/netMHCpan}}
IMAGE=${NEOAG_NETMHCPAN_IMAGE:-neoag-netmhcpan:4.2c-ubuntu22.04}
SIF=${NEOAG_NETMHCPAN_SIF:-$REPO_ROOT/containers/netmhcpan/netmhcpan-4.2c-ubuntu22.04.sif}
TMPDIR_HOST=${NEOAG_NETMHCPAN_TMPDIR:-$REPO_ROOT/work/netmhcpan_tmp}
ENGINE=${NEOAG_NETMHCPAN_ENGINE:-auto}
if [[ -n ${NEOAG_NETMHCPAN_CONTAINER_BIN:-} ]]; then
  CONTAINER_BIN=$NEOAG_NETMHCPAN_CONTAINER_BIN
elif [[ -x $NETMHCPAN_HOME/netMHCpan ]]; then
  CONTAINER_BIN=$NETMHCPAN_HOME/netMHCpan
elif [[ -x $NETMHCPAN_HOME/bin/netMHCpan ]]; then
  CONTAINER_BIN=$NETMHCPAN_HOME/bin/netMHCpan
else
  CONTAINER_BIN=$NETMHCPAN_HOME/bin/netMHCpan-4.2
fi
# Portable frontend (100T/netMHCpan) expects $NETMHCpan = install root (has bin/ wrappers).
# Legacy DTU layout used Linux_$arch as NETMHCpan; only fall back when no root frontend.
if [[ -n ${NEOAG_NETMHCPAN_PLATFORM_HOME:-} ]]; then
  PLATFORM_HOME=$NEOAG_NETMHCPAN_PLATFORM_HOME
elif [[ -x $NETMHCPAN_HOME/bin/netMHCpan-4.2 || -x $NETMHCPAN_HOME/netMHCpan ]]; then
  PLATFORM_HOME=$NETMHCPAN_HOME
else
  PLATFORM_HOME=$NETMHCPAN_HOME/Linux_$(uname -m)
fi
[[ ${1:-} == -h || ${1:-} == --help ]] && { echo "Usage: $0 [netMHCpan args]"; exit 0; }
[[ -x "$CONTAINER_BIN" ]] || { echo "ERROR: missing $CONTAINER_BIN" >&2; exit 2; }
mkdir -p "$TMPDIR_HOST"
NETMHCPAN_ARGS=("$@")
for i in "${!NETMHCPAN_ARGS[@]}"; do
  if [[ ${NETMHCPAN_ARGS[$i]} == HLA-*\** ]]; then
    NETMHCPAN_ARGS[$i]=${NETMHCPAN_ARGS[$i]//\*/}
  fi
done
CONTAINER_CMD="export TMPDIR=/tmp/netmhcpan NEOAG_NETMHCPAN_TMPDIR=/tmp/netmhcpan; exec \"\$NEOAG_NETMHCPAN_CONTAINER_BIN\" \"\$@\""
run_docker() {
  docker image inspect "$IMAGE" >/dev/null 2>&1 || { echo "ERROR: build image first: $REPO_ROOT/scripts/build_netmhcpan_container.sh docker" >&2; return 127; }
  mounts=(-v "$NETMHCPAN_HOME:$NETMHCPAN_HOME:ro" -v "$TMPDIR_HOST:/tmp/netmhcpan:rw")
  if [[ "$PWD" == "$REPO_ROOT" ]]; then
    mounts+=( -v "$REPO_ROOT:$REPO_ROOT:rw" )
  else
    mounts+=( -v "$REPO_ROOT:$REPO_ROOT:ro" -v "$PWD:$PWD:rw" )
  fi
  [[ -n ${NEOAG_NETMHCPAN_EXTRA_MOUNTS:-} ]] && IFS=, read -r -a extra <<< "$NEOAG_NETMHCPAN_EXTRA_MOUNTS" && for m in "${extra[@]}"; do mounts+=( -v "$m" ); done
  docker run --rm --user "$(id -u):$(id -g)" --workdir "$PWD" -e TMPDIR=/tmp/netmhcpan -e NEOAG_NETMHCPAN_TMPDIR=/tmp/netmhcpan -e NMHOME="$NETMHCPAN_HOME" -e NETMHCPAN_HOME="$NETMHCPAN_HOME" -e NETMHCpan="$PLATFORM_HOME" -e NEOAG_NETMHCPAN_CONTAINER_BIN="$CONTAINER_BIN" "${mounts[@]}" "$IMAGE" "$CONTAINER_CMD" -- "$@"
}
run_apptainer() {
  runtime=$1; shift
  [[ -f "$SIF" ]] || { echo "ERROR: build sif first: $REPO_ROOT/scripts/build_netmhcpan_container.sh apptainer" >&2; return 127; }
  binds=(-B "$NETMHCPAN_HOME:$NETMHCPAN_HOME:ro" -B "$TMPDIR_HOST:/tmp/netmhcpan:rw")
  if [[ "$PWD" == "$REPO_ROOT" ]]; then
    binds+=( -B "$REPO_ROOT:$REPO_ROOT:rw" )
  else
    binds+=( -B "$REPO_ROOT:$REPO_ROOT:ro" -B "$PWD:$PWD:rw" )
  fi
  "$runtime" exec --cleanenv "${binds[@]}" --env TMPDIR=/tmp/netmhcpan,NEOAG_NETMHCPAN_TMPDIR=/tmp/netmhcpan,NMHOME="$NETMHCPAN_HOME",NETMHCPAN_HOME="$NETMHCPAN_HOME",NETMHCpan="$PLATFORM_HOME",NEOAG_NETMHCPAN_CONTAINER_BIN="$CONTAINER_BIN" --pwd "$PWD" "$SIF" /bin/bash -lc "$CONTAINER_CMD" -- "$@"
}
case "$ENGINE" in
  docker) run_docker "${NETMHCPAN_ARGS[@]}" ;;
  apptainer) run_apptainer apptainer "${NETMHCPAN_ARGS[@]}" ;;
  singularity) run_apptainer singularity "${NETMHCPAN_ARGS[@]}" ;;
  auto) if command -v docker >/dev/null 2>&1; then run_docker "${NETMHCPAN_ARGS[@]}"; elif command -v apptainer >/dev/null 2>&1; then run_apptainer apptainer "${NETMHCPAN_ARGS[@]}"; elif command -v singularity >/dev/null 2>&1; then run_apptainer singularity "${NETMHCPAN_ARGS[@]}"; else echo "ERROR: no container runtime" >&2; exit 127; fi ;;
  *) echo "ERROR: invalid NEOAG_NETMHCPAN_ENGINE=$ENGINE" >&2; exit 2 ;;
esac
