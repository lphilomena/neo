#!/usr/bin/env bash
# Skill-first bootstrap for deploying Project B on a new machine.
# This script is safe-by-default: it installs Python entry points, creates local
# manifest templates, validates skills, runs read-only Doctor, and generates a
# pipeline-full dry-run plan. It does not download references, install licensed
# tools, submit HPC jobs, delete files, or run production workflows.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PREFIX="conf"
OUTDIR="work/agent_deploy"
PYTHON_BIN="${PYTHON:-python}"
DO_INSTALL=1
MINI_SMOKE=0
STRICT=0

usage() {
  cat <<'USAGE'
Usage: bash scripts/bootstrap_agent_deploy.sh [options]

Options:
  --prefix DIR       Local manifest directory (default: conf)
  --outdir DIR       Output directory for checks (default: work/agent_deploy)
  --python PATH      Python executable to use (default: $PYTHON or python)
  --skip-install     Skip python -m pip install -e '.[test]'
  --mini-smoke       Ask Doctor to run optional mini smoke checks
  --strict           Exit non-zero when Doctor or pipeline dry-run reports issues
  -h, --help         Show this help

Safe-by-default boundary:
  No external reference downloads, licensed-tool installation, HPC submission,
  deletion, overwrite, or heavy pipeline execution is performed.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix)
      PREFIX="$2"; shift 2 ;;
    --outdir)
      OUTDIR="$2"; shift 2 ;;
    --python)
      PYTHON_BIN="$2"; shift 2 ;;
    --skip-install)
      DO_INSTALL=0; shift ;;
    --mini-smoke)
      MINI_SMOKE=1; shift ;;
    --strict)
      STRICT=1; shift ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2 ;;
  esac
done

cd "$ROOT"
mkdir -p "$PREFIX" "$OUTDIR"

log() { printf '[agent-deploy] %s\n' "$*"; }
copy_if_missing() {
  local src="$1"
  local dst="$2"
  if [[ -f "$dst" ]]; then
    log "keep existing $dst"
  else
    cp "$src" "$dst"
    log "created $dst from $src"
  fi
}

log "project root: $ROOT"
log "read first: .agents/config/skills_registry.abcd.json"
log "read first: docs/SKILL_FIRST_MIGRATION.md"
log "read first: docs/SKILLS_TAXONOMY_ABCD.md"
log "read first: docs/CONTROLLED_EXECUTION_PHASE0_2.md"

"$PYTHON_BIN" - <<'PY'
import sys
if sys.version_info < (3, 11):
    raise SystemExit("Python >=3.11 is required by pyproject.toml")
print(f"python={sys.executable} version={sys.version.split()[0]}")
PY

if [[ "$DO_INSTALL" == "1" ]]; then
  log "install editable package entry points"
  "$PYTHON_BIN" -m pip install -e '.[test]'
else
  log "skip install by request"
fi

copy_if_missing "configs/controlled_execution/tools_manifest.example.yaml" "$PREFIX/tools_manifest.yaml"
copy_if_missing "configs/controlled_execution/reference_manifest.example.yaml" "$PREFIX/reference_manifest.yaml"
copy_if_missing "configs/controlled_execution/sample_manifest.example.yaml" "$PREFIX/sample_manifest.yaml"

export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

run_neoag_skill() {
  if command -v neoag-skill >/dev/null 2>&1; then
    neoag-skill "$@"
  else
    "$PYTHON_BIN" -m neoag.skill_taxonomy.cli "$@"
  fi
}

run_neoag_doctor() {
  if command -v neoag-doctor >/dev/null 2>&1; then
    neoag-doctor "$@"
  else
    "$PYTHON_BIN" -m neoag.controlled_execution.doctor "$@"
  fi
}

run_neoag_pipeline_full() {
  if command -v neoag-pipeline-full >/dev/null 2>&1; then
    neoag-pipeline-full "$@"
  else
    "$PYTHON_BIN" -m neoag.controlled_execution.pipeline_runner "$@"
  fi
}

log "validate skill package"
run_neoag_skill validate --root . --outdir "$OUTDIR/skill_validate"

DOCTOR_ARGS=(
  --project-root .
  --tools-manifest "$PREFIX/tools_manifest.yaml"
  --reference-manifest "$PREFIX/reference_manifest.yaml"
  --sample-manifest "$PREFIX/sample_manifest.yaml"
  --outdir "$OUTDIR/doctor"
  --dry-run
)
if [[ "$MINI_SMOKE" == "1" ]]; then
  DOCTOR_ARGS+=(--mini-smoke)
fi

log "run read-only Doctor"
DOCTOR_EXIT=0
run_neoag_doctor "${DOCTOR_ARGS[@]}" || DOCTOR_EXIT=$?
if [[ "$DOCTOR_EXIT" != "0" ]]; then
  log "Doctor reported non-ready status or failed with exit code $DOCTOR_EXIT; continuing to dry-run planning"
fi

log "run pipeline-full dry-run"
PIPELINE_EXIT=0
run_neoag_pipeline_full \
  --sample-manifest "$PREFIX/sample_manifest.yaml" \
  --tools-manifest "$PREFIX/tools_manifest.yaml" \
  --reference-manifest "$PREFIX/reference_manifest.yaml" \
  --outdir "$OUTDIR/pipeline_full" \
  --profile local || PIPELINE_EXIT=$?
if [[ "$PIPELINE_EXIT" != "0" ]]; then
  log "pipeline-full dry-run returned exit code $PIPELINE_EXIT"
fi

SUMMARY="$OUTDIR/agent_deploy_summary.md"
cat > "$SUMMARY" <<EOF
# Skill-first deployment summary

Project root: \`$ROOT\`

Generated or reused manifests:

- \`$PREFIX/tools_manifest.yaml\`
- \`$PREFIX/reference_manifest.yaml\`
- \`$PREFIX/sample_manifest.yaml\`

Check outputs:

- Doctor exit code: \`$DOCTOR_EXIT\`
- Pipeline dry-run exit code: \`$PIPELINE_EXIT\`

- Skill validation: \`$OUTDIR/skill_validate\`
- Doctor: \`$OUTDIR/doctor\`
- Pipeline dry-run: \`$OUTDIR/pipeline_full\`

Next files to inspect:

- \`$OUTDIR/doctor/doctor_summary.md\`
- \`$OUTDIR/doctor/blocking_issues.tsv\`
- \`$OUTDIR/doctor/recommended_fixes.md\`
- \`$OUTDIR/pipeline_full/pipeline_plan.md\`
- \`$OUTDIR/pipeline_full/pipeline_status.tsv\`

Boundary: this bootstrap did not install external bioinformatics tools, download
large references, submit HPC jobs, delete files, overwrite results, or execute a
production pipeline.
EOF

log "done: $SUMMARY"

if [[ "$STRICT" == "1" ]]; then
  if [[ "$DOCTOR_EXIT" != "0" || "$PIPELINE_EXIT" != "0" ]]; then
    exit 1
  fi
fi

