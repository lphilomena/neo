#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="$(pwd)"
TOOLS_ROOT="/root/neo/env_tool"
OUTDIR="work/remote_deploy/production_runtime"
MINI_PRIME=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root) PROJECT_ROOT="$2"; shift 2 ;;
    --tools-root) TOOLS_ROOT="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --mini-prime) MINI_PRIME=1; shift ;;
    -h|--help) echo "Usage: $0 --project-root DIR --tools-root DIR --outdir DIR [--mini-prime]"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done
mkdir -p "$OUTDIR/logs"
STATUS="$OUTDIR/production_runtime_status.tsv"
REPORT="$OUTDIR/production_runtime_report.md"
printf 'check\tstatus\tdetail\n' > "$STATUS"
check() {
  local name="$1"; shift
  local log="$OUTDIR/logs/$name.log"
  if "$@" >"$log" 2>&1; then
    printf '%s\tPASS\t%s\n' "$name" "$log" >> "$STATUS"
  else
    printf '%s\tFAIL\t%s\n' "$name" "$log" >> "$STATUS"
  fi
}
ACT="$TOOLS_ROOT/activate_neoag_production_refs.sh"
if [[ -f "$ACT" ]]; then
  source "$ACT"
  printf 'activation\tPASS\t%s\n' "$ACT" >> "$STATUS"
else
  printf 'activation\tFAIL\t%s\n' "$ACT" >> "$STATUS"
fi
export PYTHONPATH="$PROJECT_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
PY="${NEOAG_CONDA_BASE:-$TOOLS_ROOT/miniforge3}/envs/neoag-core/bin/python"
[[ -x "$PY" ]] || PY="$(command -v python3 || command -v python)"
check python_deps "$PY" -c "import torch,numpy,pandas,scipy,sklearn,psutil"
check neoag_doctor "$PY" -m neoag.controlled_execution.doctor --help
check vep "${NEOAG_VEP_BIN:-vep}" --help
check reference_fasta test -f "${NEOAG_REFERENCE_FASTA:-}"
check reference_fasta_fai test -f "${NEOAG_REFERENCE_FASTA:-}.fai"
check vep_cache test -d "${NEOAG_VEP_CACHE:-}/homo_sapiens/${NEOAG_VEP_CACHE_VERSION:-105}_GRCh38"
check netmhcpan netMHCpan -h
check mixmhcpred "${MIXMHCPRED_BIN:-$TOOLS_ROOT/wrappers/mixMHCpred_install/MixMHCpred}" -h
check prime_help "${NEOAG_PRIME_BIN:-$TOOLS_ROOT/tools/prime/PRIME}" -h
if [[ "$MINI_PRIME" == "1" ]]; then
  PEPS="$OUTDIR/prime_smoke_peptides.txt"
  OUT="$OUTDIR/prime_smoke.tsv"
  printf 'SYFPEITHI\nLLFGYPVYV\n' > "$PEPS"
  check prime_smoke "${NEOAG_PRIME_BIN:-$TOOLS_ROOT/tools/prime/PRIME}" -i "$PEPS" -o "$OUT" -a A0201 -mix "${MIXMHCPRED_BIN:-$TOOLS_ROOT/wrappers/mixMHCpred_install/MixMHCpred}"
  if [[ -s "$OUT" && $(wc -c < "$OUT") -gt 1 ]]; then
    printf 'prime_smoke_output\tPASS\t%s\n' "$OUT" >> "$STATUS"
  else
    printf 'prime_smoke_output\tFAIL\t%s\n' "$OUT" >> "$STATUS"
  fi
fi
FAILS=$(awk -F '\t' 'NR>1 && $2=="FAIL" {n++} END {print n+0}' "$STATUS")
{
  echo "# Production runtime report"
  echo
  echo "Project root: \`$PROJECT_ROOT\`"
  echo "Tools root: \`$TOOLS_ROOT\`"
  echo "Status TSV: \`$STATUS\`"
  echo "Failures: \`$FAILS\`"
  echo
  echo "Real data execution is allowed only when required checks are PASS or explicitly waived."
} > "$REPORT"
echo "production_runtime_status=$STATUS"
echo "production_runtime_report=$REPORT"
[[ "$FAILS" == "0" ]] || exit 21
