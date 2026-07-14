#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="$(pwd)"
MANIFEST_DIR="configs/local"
OUTDIR="work/remote_deploy/doctor"
PYTHON_BIN="${PYTHON:-python}"
MINI_SMOKE=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root) PROJECT_ROOT="$2"; shift 2 ;;
    --manifest-dir) MANIFEST_DIR="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --python) PYTHON_BIN="$2"; shift 2 ;;
    --mini-smoke) MINI_SMOKE=1; shift ;;
    -h|--help) echo "Usage: $0 --project-root DIR --manifest-dir DIR --outdir DIR [--mini-smoke]"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done
cd "$PROJECT_ROOT"
export PYTHONPATH="$PROJECT_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
TOOLS="$MANIFEST_DIR/tools_manifest.yaml"
REFS="$MANIFEST_DIR/reference_manifest.yaml"
SAMPLE="$MANIFEST_DIR/sample_manifest.yaml"
[[ -f "$TOOLS" ]] || { echo "TOOLS_MANIFEST_MISSING: $TOOLS" >&2; exit 14; }
[[ -f "$REFS" ]] || { echo "REFERENCE_MANIFEST_MISSING: $REFS" >&2; exit 15; }
ARGS=(--project-root . --tools-manifest "$TOOLS" --reference-manifest "$REFS" --outdir "$OUTDIR" --dry-run)
[[ -f "$SAMPLE" ]] && ARGS+=(--sample-manifest "$SAMPLE")
[[ "$MINI_SMOKE" == "1" ]] && ARGS+=(--mini-smoke)
if command -v neoag-doctor >/dev/null 2>&1; then
  neoag-doctor "${ARGS[@]}" || true
else
  "$PYTHON_BIN" -m neoag_v03.controlled_execution.doctor "${ARGS[@]}" || true
fi
echo "doctor_outdir=$OUTDIR"
