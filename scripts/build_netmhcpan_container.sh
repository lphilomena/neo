#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
IMAGE=${NEOAG_NETMHCPAN_IMAGE:-neoag-netmhcpan:4.2c-ubuntu22.04}
SIF=${NEOAG_NETMHCPAN_SIF:-$REPO_ROOT/containers/netmhcpan/netmhcpan-4.2c-ubuntu22.04.sif}

usage() {
  cat <<USAGE
Usage: $0 [docker|apptainer|all]

Builds a runtime image for the official NetMHCpan 4.2c package.
The image contains only OS dependencies such as tcsh and glibc; it does not
include licensed NetMHCpan files. Mount tools/netMHCpan at runtime.

Environment:
  NEOAG_NETMHCPAN_IMAGE  Docker image tag, default: $IMAGE
  NEOAG_NETMHCPAN_SIF    Apptainer image path, default: $SIF
USAGE
}

mode=${1:-docker}
case "$mode" in
  docker)
    docker build -t "$IMAGE" -f "$REPO_ROOT/containers/netmhcpan/Dockerfile" "$REPO_ROOT/containers/netmhcpan"
    ;;
  apptainer|singularity)
    if command -v apptainer >/dev/null 2>&1; then
      apptainer build "$SIF" "$REPO_ROOT/containers/netmhcpan/netmhcpan-4.2c.apptainer.def"
    elif command -v singularity >/dev/null 2>&1; then
      singularity build "$SIF" "$REPO_ROOT/containers/netmhcpan/netmhcpan-4.2c.apptainer.def"
    else
      echo "ERROR: apptainer/singularity is not available" >&2
      exit 127
    fi
    ;;
  all)
    "$0" docker
    "$0" apptainer
    ;;
  -h|--help|help)
    usage
    ;;
  *)
    usage >&2
    exit 2
    ;;
esac
