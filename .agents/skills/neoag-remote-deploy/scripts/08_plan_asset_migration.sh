#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="$(pwd)"
OUTDIR="work/remote_deploy"
TOOLS_ROOT="/root/neo/env_tool"
REFERENCE_ROOT="/root/neo/neodata4git"
LICENSED_ROOT="/root/neo/licensed_tools"
OLD_HOST=""
OLD_ENV_TOOL=""
OLD_REFERENCE_ROOT=""
OLD_LICENSED_ROOT=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root) PROJECT_ROOT="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --tools-root) TOOLS_ROOT="$2"; shift 2 ;;
    --reference-root) REFERENCE_ROOT="$2"; shift 2 ;;
    --licensed-root) LICENSED_ROOT="$2"; shift 2 ;;
    --old-host) OLD_HOST="$2"; shift 2 ;;
    --old-env-tool) OLD_ENV_TOOL="$2"; shift 2 ;;
    --old-reference-root) OLD_REFERENCE_ROOT="$2"; shift 2 ;;
    --old-licensed-root) OLD_LICENSED_ROOT="$2"; shift 2 ;;
    -h|--help) echo "Usage: $0 --project-root DIR --outdir DIR --tools-root DIR --reference-root DIR --licensed-root DIR [--old-host user@host --old-env-tool DIR --old-reference-root DIR --old-licensed-root DIR]"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done
mkdir -p "$OUTDIR"
PLAN="$OUTDIR/asset_migration_plan.tsv"
REPORT="$OUTDIR/asset_migration_report.md"
status_path() { [[ -e "$1" ]] && echo PRESENT || echo MISSING; }
{
  printf 'asset\ttarget_path\ttarget_status\tsource\taction\trisk\n'
  printf 'env_tool\t%s\t%s\t%s\t%s\t%s\n' "$TOOLS_ROOT" "$(status_path "$TOOLS_ROOT")" "${OLD_HOST:+$OLD_HOST:}${OLD_ENV_TOOL:-UNSET}" "sync_or_rebuild_then_rewrite_activation" "HIGH"
  printf 'references\t%s\t%s\t%s\t%s\t%s\n' "$REFERENCE_ROOT" "$(status_path "$REFERENCE_ROOT")" "${OLD_HOST:+$OLD_HOST:}${OLD_REFERENCE_ROOT:-UNSET}" "sync_or_stage_reference_manifest" "HIGH"
  printf 'licensed_tools\t%s\t%s\t%s\t%s\t%s\n' "$LICENSED_ROOT" "$(status_path "$LICENSED_ROOT")" "${OLD_HOST:+$OLD_HOST:}${OLD_LICENSED_ROOT:-UNSET}" "sync_only_if_license_allows" "HIGH"
  printf 'activation\t%s\t%s\tgenerated\t%s\t%s\n' "$TOOLS_ROOT/activate_neoag_production_refs.sh" "$(status_path "$TOOLS_ROOT/activate_neoag_production_refs.sh")" "rewrite_required" "MEDIUM"
  printf 'common_sh\t%s\t%s\tproject\t%s\t%s\n' "$PROJECT_ROOT/scripts/common.sh" "$(status_path "$PROJECT_ROOT/scripts/common.sh")" "verify_default_tools_root" "LOW"
} > "$PLAN"
{
  echo "# Production Asset Migration Plan"
  echo
  echo "Project root: \`$PROJECT_ROOT\`"
  echo "Target tools root: \`$TOOLS_ROOT\`"
  echo "Target reference root: \`$REFERENCE_ROOT\`"
  echo "Target licensed-tool root: \`$LICENSED_ROOT\`"
  echo
  echo "Plan TSV: \`$PLAN\`"
  echo
  echo "## Required approval before execution"
  echo
  echo "No files were copied. Review licensed-tool rights and data volume before running sync with \`--execute\`."
} > "$REPORT"
echo "asset_migration_plan=$PLAN"
echo "asset_migration_report=$REPORT"
