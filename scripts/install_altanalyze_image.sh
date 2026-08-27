#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOURCE="${NEOAG_ALTANALYZE_IMAGE_SOURCE:-frankligy123/altanalyze@sha256:b93cb071af933290daf385a75152349a0f1c75678bd61f72528f372a100e41bd}"
IMAGE="${NEOAG_ALTANALYZE_IMAGE:-neoag-altanalyze:snaf}"
TOOLS_ENV="${NEOAG_TOOLS_ENV:-${ROOT}/conf/tools.env.local.sh}"

command -v docker >/dev/null 2>&1 || {
  echo "ERROR: docker is required for the AltAnalyze image" >&2
  exit 2
}

if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
  docker pull "$SOURCE"
  docker tag "$SOURCE" "$IMAGE"
fi

docker image inspect "$IMAGE" --format 'AltAnalyze image ready: {{.Id}} {{.Size}} bytes'
entrypoint="$(docker image inspect "$IMAGE" --format '{{json .Config.Entrypoint}} {{json .Config.Cmd}}')"
[[ "$entrypoint" == *"AltAnalyze"* || "$entrypoint" == *"identify"* ]] || {
  echo "ERROR: AltAnalyze image does not expose the expected SNAF runtime entrypoint: $entrypoint" >&2
  exit 3
}

mkdir -p "$(dirname "$TOOLS_ENV")"
touch "$TOOLS_ENV"
if grep -q '^export NEOAG_ALTANALYZE_IMAGE=' "$TOOLS_ENV"; then
  sed -i "s|^export NEOAG_ALTANALYZE_IMAGE=.*|export NEOAG_ALTANALYZE_IMAGE=\"${IMAGE}\"|" "$TOOLS_ENV"
else
  printf '\nexport NEOAG_ALTANALYZE_IMAGE="%s"\n' "$IMAGE" >> "$TOOLS_ENV"
fi
