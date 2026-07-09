#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
NETMHCPAN_HOME=${NETMHCPAN_HOME:-${NETMHCpan:-$REPO_ROOT/tools/netMHCpan}}
IMAGE=${NEOAG_NETMHCPAN_IMAGE:-neoag-netmhcpan:4.2c-ubuntu22.04}
SIF=${NEOAG_NETMHCPAN_SIF:-$REPO_ROOT/containers/netmhcpan/netmhcpan-4.2c-ubuntu22.04.sif}
TMPDIR_HOST=${NEOAG_NETMHCPAN_TMPDIR:-$REPO_ROOT/work/netmhcpan_tmp}
ENGINE=${NEOAG_NETMHCPAN_ENGINE:-auto}
[[ ${1:-} == -h || ${1:-} == --help ]] && { echo "Usage: $0 [netMHCpan args]"; exit 0; }
[[ -x "$NETMHCPAN_HOME/netMHCpan" ]] || { echo "ERROR: missing $NETMHCPAN_HOME/netMHCpan" >&2; exit 2; }
mkdir -p "$TMPDIR_HOST"
NETMHCPAN_ARGS=("$@")
for i in "${!NETMHCPAN_ARGS[@]}"; do
  if [[ ${NETMHCPAN_ARGS[$i]} == HLA-*\** ]]; then
    NETMHCPAN_ARGS[$i]=${NETMHCPAN_ARGS[$i]//\*/}
  fi
done
CONTAINER_CMD="export TMPDIR=/tmp/netmhcpan; exec /bin/tcsh -f \"\$NETMHCPAN_HOME/netMHCpan\" \"\$@\""
run_docker() {
  docker image inspect "$IMAGE" >/dev/null 2>&1 || { echo "ERROR: build image first: $REPO_ROOT/scripts/build_netmhcpan_container.sh docker" >&2; return 127; }
  mounts=(-v "$NETMHCPAN_HOME:$NETMHCPAN_HOME:ro" -v "$TMPDIR_HOST:/tmp/netmhcpan:rw" -v "$PWD:$PWD:rw" -v "$REPO_ROOT:$REPO_ROOT:rw")
  [[ -d /mnt ]] && mounts+=( -v /mnt:/mnt:rw )
  [[ -n ${NEOAG_NETMHCPAN_EXTRA_MOUNTS:-} ]] && IFS=, read -r -a extra <<< "$NEOAG_NETMHCPAN_EXTRA_MOUNTS" && for m in "${extra[@]}"; do mounts+=( -v "$m" ); done
  docker run --rm --user "$(id -u):$(id -g)" --workdir "$PWD" -e TMPDIR=/tmp/netmhcpan -e NETMHCPAN_HOME="$NETMHCPAN_HOME" -e NETMHCpan="$NETMHCPAN_HOME" "${mounts[@]}" "$IMAGE" "$CONTAINER_CMD" -- "$@"
}
run_apptainer() {
  runtime=$1; shift
  [[ -f "$SIF" ]] || { echo "ERROR: build sif first: $REPO_ROOT/scripts/build_netmhcpan_container.sh apptainer" >&2; return 127; }
  binds=(-B "$NETMHCPAN_HOME:$NETMHCPAN_HOME:ro" -B "$TMPDIR_HOST:/tmp/netmhcpan:rw" -B "$PWD:$PWD:rw" -B "$REPO_ROOT:$REPO_ROOT:rw")
  [[ -d /mnt ]] && binds+=( -B /mnt:/mnt:rw )
  "$runtime" exec --cleanenv "${binds[@]}" --env TMPDIR=/tmp/netmhcpan,NETMHCPAN_HOME="$NETMHCPAN_HOME",NETMHCpan="$NETMHCPAN_HOME" --pwd "$PWD" "$SIF" /bin/bash -lc "$CONTAINER_CMD" -- "$@"
}
case "$ENGINE" in
  docker) run_docker "${NETMHCPAN_ARGS[@]}" ;;
  apptainer) run_apptainer apptainer "${NETMHCPAN_ARGS[@]}" ;;
  singularity) run_apptainer singularity "${NETMHCPAN_ARGS[@]}" ;;
  auto) if command -v docker >/dev/null 2>&1; then run_docker "${NETMHCPAN_ARGS[@]}"; elif command -v apptainer >/dev/null 2>&1; then run_apptainer apptainer "${NETMHCPAN_ARGS[@]}"; elif command -v singularity >/dev/null 2>&1; then run_apptainer singularity "${NETMHCPAN_ARGS[@]}"; else echo "ERROR: no container runtime" >&2; exit 127; fi ;;
  *) echo "ERROR: invalid NEOAG_NETMHCPAN_ENGINE=$ENGINE" >&2; exit 2 ;;
esac
