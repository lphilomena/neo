#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="$(pwd)"
OUTDIR="work/remote_deploy"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root) PROJECT_ROOT="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    -h|--help) echo "Usage: $0 --project-root DIR --outdir DIR"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done
mkdir -p "$OUTDIR"
REPORT="$OUTDIR/preflight_report.md"
TSV="$OUTDIR/preflight_status.tsv"
{
  printf 'check\tstatus\tdetail\n'
  for cmd in bash git curl wget tar unzip python java docker apptainer nextflow; do
    if command -v "$cmd" >/dev/null 2>&1; then
      printf '%s\tFOUND\t%s\n' "$cmd" "$(command -v "$cmd")"
    else
      printf '%s\tMISSING\t\n' "$cmd"
    fi
  done
} > "$TSV"
{
  echo "# Preflight report"
  echo
  echo "Project root: \`$PROJECT_ROOT\`"
  echo
  echo '```text'
  uname -a || true
  python --version 2>&1 || true
  df -h "$PROJECT_ROOT" || true
  echo '```'
  echo
  echo "Status table: \`$TSV\`"
} > "$REPORT"
echo "preflight_report=$REPORT"
