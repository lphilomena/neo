#!/usr/bin/env bash
set -euo pipefail
OUTDIR="work/remote_deploy"
OLD_HOST=""
OLD_ENV_TOOL=""
OLD_REFERENCE_ROOT=""
OLD_LICENSED_ROOT=""
TOOLS_ROOT="/root/neo/env_tool"
REFERENCE_ROOT="/root/neo/neodata4git"
LICENSED_ROOT="/root/neo/licensed_tools"
EXECUTE=0
RSYNC_OPTS=(-aH --info=progress2)
while [[ $# -gt 0 ]]; do
  case "$1" in
    --outdir) OUTDIR="$2"; shift 2 ;;
    --old-host) OLD_HOST="$2"; shift 2 ;;
    --old-env-tool) OLD_ENV_TOOL="$2"; shift 2 ;;
    --old-reference-root) OLD_REFERENCE_ROOT="$2"; shift 2 ;;
    --old-licensed-root) OLD_LICENSED_ROOT="$2"; shift 2 ;;
    --tools-root) TOOLS_ROOT="$2"; shift 2 ;;
    --reference-root) REFERENCE_ROOT="$2"; shift 2 ;;
    --licensed-root) LICENSED_ROOT="$2"; shift 2 ;;
    --execute) EXECUTE=1; shift ;;
    -h|--help) echo "Usage: $0 --old-host user@host --old-env-tool DIR --old-reference-root DIR --old-licensed-root DIR --tools-root DIR --reference-root DIR --licensed-root DIR [--execute]"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done
mkdir -p "$OUTDIR"
LOG="$OUTDIR/asset_sync.log"
REPORT="$OUTDIR/asset_sync_report.md"
if [[ "$EXECUTE" != "1" ]]; then
  RSYNC_OPTS+=(--dry-run)
  MODE="DRY_RUN"
else
  MODE="EXECUTE"
fi
[[ -n "$OLD_HOST" ]] || { echo "ASSET_SOURCE_MISSING: --old-host required" >&2; exit 20; }
run_sync() {
  local name="$1" src="$2" dst="$3"
  if [[ -z "$src" ]]; then
    echo "SKIP $name source unset" | tee -a "$LOG"
    return 0
  fi
  mkdir -p "$(dirname "$dst")"
  echo "[$MODE] $name: $OLD_HOST:$src/ -> $dst/" | tee -a "$LOG"
  rsync "${RSYNC_OPTS[@]}" -e "ssh -o BatchMode=yes" "$OLD_HOST:$src/" "$dst/" | tee -a "$LOG"
}
: > "$LOG"
run_sync env_tool "$OLD_ENV_TOOL" "$TOOLS_ROOT"
run_sync references "$OLD_REFERENCE_ROOT" "$REFERENCE_ROOT"
run_sync licensed_tools "$OLD_LICENSED_ROOT" "$LICENSED_ROOT"
{
  echo "# Asset sync report"
  echo
  echo "Mode: \`$MODE\`"
  echo "Log: \`$LOG\`"
  echo
  if [[ "$EXECUTE" != "1" ]]; then
    echo "This was a dry run. Re-run with \`--execute\` only after user approval and license review."
  else
    echo "Assets were synced. Next run \`10_rewrite_production_activation.sh --write\`."
  fi
} > "$REPORT"
echo "asset_sync_report=$REPORT"
