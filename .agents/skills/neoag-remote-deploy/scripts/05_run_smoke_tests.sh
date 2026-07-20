#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="$(pwd)"
OUTDIR="work/remote_deploy"
PYTHON_BIN="${PYTHON:-python}"
TIER="tier0"
RUN_NEXTFLOW=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root) PROJECT_ROOT="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --python) PYTHON_BIN="$2"; shift 2 ;;
    --tier) TIER="$2"; shift 2 ;;
    --nextflow) RUN_NEXTFLOW=1; shift ;;
    -h|--help) echo "Usage: $0 --project-root DIR --outdir DIR [--tier tier0|tier1|tier2|tier3] [--nextflow]"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done
cd "$PROJECT_ROOT"
mkdir -p "$OUTDIR"
export PYTHONPATH="$PROJECT_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
REPORT="$OUTDIR/smoke_test_report.md"
STATUS="$OUTDIR/smoke_status.tsv"
printf 'test\tstatus\tlog\n' > "$STATUS"
run_step() {
  local name="$1"; shift
  local log="$OUTDIR/${name}.log"
  if "$@" >"$log" 2>&1; then
    printf '%s\tPASS\t%s\n' "$name" "$log" >> "$STATUS"
  else
    printf '%s\tFAIL\t%s\n' "$name" "$log" >> "$STATUS"
    return 1
  fi
}
FAIL=0
run_step compileall "$PYTHON_BIN" -m compileall -q src || FAIL=1
run_step skill_validate "$PYTHON_BIN" -m neoag.skill_taxonomy.cli validate --root . --outdir "$OUTDIR/skill_validate" || FAIL=1
run_step run_demo "$PYTHON_BIN" -m neoag.cli run-demo --outdir "$OUTDIR/demo_v043" --sample-id DEMO001 || FAIL=1
if [[ "$TIER" != "tier0" ]]; then
  run_step pytest_skills "$PYTHON_BIN" -m pytest -q tests/test_skills_taxonomy_abcd.py || FAIL=1
fi
if [[ "$RUN_NEXTFLOW" == "1" ]]; then
  find bin -maxdepth 1 -type f -exec chmod +x {} \;
  run_step nextflow_fixture bin/neoag-nextflow run workflows/main.nf --pvac_files data/fixtures/pvacseq_aggregated.tsv --outdir "$OUTDIR/demo_nf" --sample_id NF_DEMO || FAIL=1
fi
{
  echo "# Smoke test report"
  echo
  echo "Tier: \`$TIER\`"
  echo
  echo "Status: \`$STATUS\`"
} > "$REPORT"
[[ "$FAIL" == "0" ]] || { echo "SMOKE_FAILED: see $STATUS" >&2; exit 13; }
echo "smoke_report=$REPORT"
