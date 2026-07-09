#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
[[ -f "$REPO_ROOT/conf/tools.env.sh" ]] && source "$REPO_ROOT/conf/tools.env.sh"
IMAGE=${NEOAG_HLALA_IMAGE:-neoag-hla-la:ubuntu22.04}
docker image inspect "$IMAGE" >/dev/null 2>&1 || { echo "WARN: HLA-LA image missing; build with scripts/build_priority_tool_containers.sh hla-la"; exit 0; }
docker run --rm "$IMAGE" "perl -v >/dev/null && samtools --version >/dev/null && echo PASS: HLA-LA container base runtime starts"
HLALA_HOME=${HLALA_HOME:-${HLA_LA_HOME:-/mnt/zjl-bgi-zzb/peixunban/gl/data/tools/hla-la/env/opt/hla-la}}
[[ -d "$HLALA_HOME" ]] && echo "PASS: HLA-LA home exists: $HLALA_HOME" || echo "WARN: HLA-LA home missing: $HLALA_HOME"
