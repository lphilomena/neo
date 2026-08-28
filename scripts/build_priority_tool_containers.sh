#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)

usage() {
  cat <<USAGE
Usage: $0 [base|netmhcstabpan|hla-la|spechla|purple-suite|easyfuse|all]

Build Docker runtime images for high-priority external tools.
Licensed software and large reference data are not baked into the images;
mount them at runtime with the corresponding run_*_container.sh wrapper.
USAGE
}

build_one() {
  local name=$1
  local tag=$2
  docker build -t "$tag" -f "$REPO_ROOT/containers/$name/Dockerfile" "$REPO_ROOT/containers/$name"
}

mode=${1:-all}
case "$mode" in
  base)
    build_one base-bioinfo neoag-base-bioinfo:ubuntu22.04
    ;;
  netmhcstabpan)
    "$0" base
    build_one netmhcstabpan neoag-netmhcstabpan:1.0-ubuntu22.04
    ;;
  hla-la)
    "$0" base
    build_one hla-la neoag-hla-la:ubuntu22.04
    ;;
  spechla)
    "$0" base
    build_one spechla neoag-spechla:ubuntu22.04
    ;;
  purple-suite)
    build_one purple-suite neoag-purple-suite:ubuntu22.04
    ;;
  easyfuse)
    "$0" base
    build_one easyfuse neoag-easyfuse:ubuntu22.04
    ;;
  all)
    "$0" base
    build_one netmhcstabpan neoag-netmhcstabpan:1.0-ubuntu22.04
    build_one hla-la neoag-hla-la:ubuntu22.04
    build_one spechla neoag-spechla:ubuntu22.04
    build_one purple-suite neoag-purple-suite:ubuntu22.04
    build_one easyfuse neoag-easyfuse:ubuntu22.04
    ;;
  -h|--help|help)
    usage
    ;;
  *)
    usage >&2
    exit 2
    ;;
esac
