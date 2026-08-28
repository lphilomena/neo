#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
OUTDIR=${NEOAG_APPTAINER_OUTDIR:-$REPO_ROOT/containers/sif}
mkdir -p "$OUTDIR"

if command -v apptainer >/dev/null 2>&1; then
  RUNTIME=apptainer
elif command -v singularity >/dev/null 2>&1; then
  RUNTIME=singularity
else
  echo "ERROR: apptainer/singularity is not available" >&2
  exit 127
fi

build_sif() {
  local image=$1
  local sif=$2
  docker image inspect "$image" >/dev/null 2>&1 || {
    echo "ERROR: Docker image missing: $image" >&2
    echo "Build Docker images first: $REPO_ROOT/scripts/build_priority_tool_containers.sh all" >&2
    exit 2
  }
  "$RUNTIME" build "$OUTDIR/$sif" "docker-daemon://$image"
}

mode=${1:-all}
case "$mode" in
  netmhcpan) build_sif neoag-netmhcpan:4.2c-ubuntu22.04 neoag-netmhcpan-4.2c-ubuntu22.04.sif ;;
  netmhcstabpan) build_sif neoag-netmhcstabpan:1.0-ubuntu22.04 neoag-netmhcstabpan-1.0-ubuntu22.04.sif ;;
  hla-la) build_sif neoag-hla-la:ubuntu22.04 neoag-hla-la-ubuntu22.04.sif ;;
  spechla) build_sif neoag-spechla:ubuntu22.04 neoag-spechla-ubuntu22.04.sif ;;
  purple-suite) build_sif neoag-purple-suite:ubuntu22.04 neoag-purple-suite-ubuntu22.04.sif ;;
  easyfuse) build_sif neoag-easyfuse:ubuntu22.04 neoag-easyfuse-ubuntu22.04.sif ;;
  all)
    "$0" netmhcpan
    "$0" netmhcstabpan
    "$0" hla-la
    "$0" spechla
    "$0" purple-suite
    "$0" easyfuse
    ;;
  -h|--help|help)
    echo "Usage: $0 [netmhcpan|netmhcstabpan|hla-la|spechla|purple-suite|easyfuse|all]"
    ;;
  *)
    echo "ERROR: unknown mode: $mode" >&2
    exit 2
    ;;
esac
