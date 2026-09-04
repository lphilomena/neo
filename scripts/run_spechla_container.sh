#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$SCRIPT_DIR/.." && pwd)
[[ -f "$REPO_ROOT/conf/tools.env.sh" ]] && source "$REPO_ROOT/conf/tools.env.sh"
# shellcheck source=lib/resolve_spechla_db.sh
source "$REPO_ROOT/scripts/lib/resolve_spechla_db.sh"
SPECHLA_HOME=${SPECHLA_HOME:-${NEOAG_SPECHLA_HOME:-${NEOAG_TOOLS_ROOT:-$REPO_ROOT}/tools/SpecHLA}}
SPECHLA_DB="$(neoag_resolve_spechla_db "$REPO_ROOT" || true)"
SPECHLA_ENV=${SPECHLA_ENV:-}
if [[ -z "$SPECHLA_ENV" ]]; then
  for candidate in \
    "${CONDA_BASE:-}/envs/neoag-tools" \
    "${NEOAG_TOOLS_ROOT:-}/miniforge3/envs/neoag-tools" \
    "$SPECHLA_HOME/../../miniforge3/envs/neoag-tools" \
    "$SPECHLA_HOME/spechla_env"; do
    if [[ -x "$candidate/bin/python3" ]]; then SPECHLA_ENV="$candidate"; break; fi
  done
fi
IMAGE=${NEOAG_SPECHLA_IMAGE:-neoag-spechla:ubuntu22.04}
MODE=${SPECHLA_MODE:-auto}
[[ ${1:-} == -h || ${1:-} == --help ]] && { cat <<USAGE
Usage: $0 [SpecHLA args]

Default command selects the first available executable:
  spechla from PATH, or $SPECHLA_HOME/script/whole/SpecHLA.sh
Set SPECHLA_MODE=extract to run ExtractHLAread.sh.
Set SPECHLA_MODE=loh to run the SpecHLA HLA-LOH copy-number module.
Set SPECHLA_CMD=/path/to/custom_command to override.
USAGE
exit 0; }
[[ -d "$SPECHLA_HOME" ]] || { echo "ERROR: SpecHLA home missing: $SPECHLA_HOME" >&2; exit 2; }
[[ -n "$SPECHLA_DB" && -d "$SPECHLA_DB" ]] || {
  echo "ERROR: SpecHLA database missing. Set SPECHLA_DB or stage data/hla/spechla/db (legacy: data/hla/spechla_db) under the reference bundle." >&2
  exit 2
}
docker image inspect "$IMAGE" >/dev/null 2>&1 || { echo "ERROR: build image first: $REPO_ROOT/scripts/build_priority_tool_containers.sh spechla" >&2; exit 127; }
mounts=(-v "$SPECHLA_HOME:$SPECHLA_HOME:rw" -v "$PWD:$PWD:rw" -v "$REPO_ROOT:$REPO_ROOT:rw")
[[ "$SPECHLA_DB" == "$SPECHLA_HOME"/* ]] || mounts+=( -v "$SPECHLA_DB:$SPECHLA_DB:rw" )
if [[ -n "$SPECHLA_ENV" ]]; then
  SPECHLA_ENV=$(cd "$SPECHLA_ENV" && pwd -P)
  mounts+=( -v "$SPECHLA_ENV:$SPECHLA_ENV:ro" )
  CMD_PREFIX="export PATH='$SPECHLA_ENV/bin':\"\$PATH\"; "
else
  CMD_PREFIX=""
fi
[[ -d /mnt ]] && mounts+=( -v /mnt:/mnt:rw )
[[ -d /tmp ]] && mounts+=( -v /tmp:/tmp:rw )

if [[ ${SPECHLA_AUTO_INDEX:-1} == 1 && -x "$SPECHLA_HOME/bin/novoalign" && -f "$SPECHLA_HOME/bin/novoalign.lic" ]]; then
  indexer="$SPECHLA_HOME/bin/novoindex"
  [[ -x "$indexer" ]] || { echo "ERROR: novoalign is licensed but novoindex is missing: $indexer" >&2; exit 3; }
  for prefix in hla_gen.format.filter.extend.DRB.no26789 hla_gen.format.filter.extend.DRB.no26789.v2; do
    fasta="$SPECHLA_DB/ref/$prefix.fasta"
    index="$SPECHLA_DB/ref/$prefix.ndx"
    if [[ -s "$fasta" && ! -s "$index" ]]; then
      echo "Building missing SpecHLA novoalign index: $index" >&2
      tmp="$index.tmp.$$"
      "$indexer" "$tmp" "$fasta"
      mv "$tmp" "$index"
    fi
  done
fi
if [[ -n ${SPECHLA_CMD:-} ]]; then
  CMD="${CMD_PREFIX}exec \"$SPECHLA_CMD\" \"\$@\""
elif [[ "$MODE" == extract ]]; then
  CMD="${CMD_PREFIX}cd \"$SPECHLA_HOME\"; exec bash script/ExtractHLAread.sh \"\$@\""
elif [[ "$MODE" == loh ]]; then
  CMD="${CMD_PREFIX}cd \"$SPECHLA_HOME\"; [[ -f script/cal.hla.copy.pl ]] || { echo ERROR: SpecHLA LOH script missing: $SPECHLA_HOME/script/cal.hla.copy.pl >&2; exit 127; }; exec perl script/cal.hla.copy.pl \"\$@\""
else
  CMD="${CMD_PREFIX}cd \"$SPECHLA_HOME\"; if command -v spechla >/dev/null 2>&1; then exec spechla \"\$@\"; elif [[ -f script/whole/SpecHLA.sh ]]; then exec bash script/whole/SpecHLA.sh \"\$@\"; else echo ERROR: SpecHLA command not found >&2; exit 127; fi"
fi
docker run --rm --user "$(id -u):$(id -g)" --workdir "$PWD" \
  -e SPECHLA_HOME="$SPECHLA_HOME" -e SPECHLA_DB="$SPECHLA_DB" -e SPECHLA_ENV="$SPECHLA_ENV" \
  "${mounts[@]}" --entrypoint /bin/bash "$IMAGE" -lc "$CMD" -- "$@"
