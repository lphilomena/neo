#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
usage() {
  cat >&2 <<'USAGE'
Usage: run_hla_loh_multi_tool.sh --sample-id ID --tumor-id ID --normal-id ID \
  --tumor-bam BAM --normal-bam BAM --hla-file FILE --purity-tsv FILE \
  --outdir DIR [--tools lohhla,spechla] [--threads N]
USAGE
}

SAMPLE_ID=""; TUMOR_ID=""; NORMAL_ID=""; TUMOR_BAM=""; NORMAL_BAM=""
HLA_FILE=""; PURITY_TSV=""; OUTDIR=""; TOOLS="lohhla,spechla"; THREADS=8
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --tumor-id) TUMOR_ID="$2"; shift 2 ;;
    --normal-id) NORMAL_ID="$2"; shift 2 ;;
    --tumor-bam) TUMOR_BAM="$2"; shift 2 ;;
    --normal-bam) NORMAL_BAM="$2"; shift 2 ;;
    --hla-file) HLA_FILE="$2"; shift 2 ;;
    --purity-tsv) PURITY_TSV="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --tools) TOOLS="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage; exit 2 ;;
  esac
done

for value in SAMPLE_ID TUMOR_ID NORMAL_ID TUMOR_BAM NORMAL_BAM HLA_FILE PURITY_TSV OUTDIR; do
  [[ -n "${!value}" ]] || { usage; exit 2; }
done
for path in "$TUMOR_BAM" "$NORMAL_BAM" "$HLA_FILE" "$PURITY_TSV"; do
  [[ -s "$path" ]] || { echo "ERROR: required input missing: $path" >&2; exit 3; }
done

read_recommendation() {
  awk -F '\t' '
    NR == 1 { for (i=1; i<=NF; i++) h[$i]=i; next }
    NR == 2 { print $(h["purity"]), $(h["ploidy"]); exit }
  ' "$PURITY_TSV"
}
read -r PURITY PLOIDY < <(read_recommendation)
awk -v p="$PURITY" -v q="$PLOIDY" 'BEGIN { exit !(p > 0 && p <= 1 && q > 0) }' || {
  echo "ERROR: recommended purity/ploidy are invalid: purity=$PURITY ploidy=$PLOIDY" >&2
  exit 4
}

require_multi_tool_purity_consensus() {
  local purity_dir
  local consensus_tsv
  local tool_summary_tsv

  purity_dir="$(cd "$(dirname "$PURITY_TSV")" && pwd)"
  consensus_tsv="$purity_dir/purity_cnv_consensus.tsv"
  tool_summary_tsv="$purity_dir/purity_cnv_tool_summary.tsv"

  if [[ "${ALLOW_SINGLE_TOOL_PURITY_FOR_HLA_LOH:-0}" == "1" ]]; then
    return 0
  fi

  [[ -s "$consensus_tsv" && -s "$tool_summary_tsv" ]] || {
    echo "ERROR: HLA LOH requires purity_cnv_consensus.tsv and purity_cnv_tool_summary.tsv next to the purity TSV" >&2
    return 1
  }

  awk -F '\t' '
    NR == 1 {
      for (i=1; i<=NF; i++) h[$i]=i
      next
    }
    NR == 2 {
      status = $(h["consensus_status"])
      exit !(status != "" && status != "SINGLE_TOOL")
    }
  ' "$consensus_tsv" || {
    echo "ERROR: HLA LOH requires non-SINGLE_TOOL purity consensus" >&2
    return 1
  }

  awk -F '\t' '
    NR == 1 {
      for (i=1; i<=NF; i++) h[$i]=i
      next
    }
    NR > 1 {
      if ($(h["status"]) == "FOUND") found++
    }
    END {
      exit !(found >= 2)
    }
  ' "$tool_summary_tsv" || {
    echo "ERROR: HLA LOH requires at least 2 purity tools with FOUND results" >&2
    return 1
  }
}
require_multi_tool_purity_consensus || exit 4

wants() { [[ ",${TOOLS}," == *",$1,"* ]]; }
mkdir -p "$OUTDIR"
LOHHLA_TSV=""; SPECHLA_TSV=""

if [[ "${ALLOW_SINGLE_HLA_LOH_TOOL:-0}" != "1" ]]; then
  wants lohhla && wants spechla || {
    echo "ERROR: HLA LOH requires both LOHHLA and SpecHLA before downstream stages" >&2
    exit 6
  }
fi

run_lohhla() {
  LOHHLA_DIR="$OUTDIR/lohhla"
  mkdir -p "$LOHHLA_DIR"
  COPYNUM="$LOHHLA_DIR/lohhla_copy_number_input.tsv"
  {
    printf '\ttumorPurity\ttumorPloidy\n'
    printf '%s\t%s\t%s\n' "$TUMOR_ID" "$PURITY" "$PLOIDY"
  } > "$COPYNUM"
  PATIENT_ID="$SAMPLE_ID" TUMOR_SAMPLE_ID="$TUMOR_ID" NORMAL_SAMPLE_ID="$NORMAL_ID" \
    TUMOR_BAM="$TUMOR_BAM" NORMAL_BAM="$NORMAL_BAM" HLA_FILE="$HLA_FILE" \
    OUTDIR="$LOHHLA_DIR" LOHHLA_NAS_ROOT="$LOHHLA_DIR/work" COPYNUM_LOC="$COPYNUM" \
    POLYSOLVER_THREADS="$THREADS" bash "$ROOT/scripts/run_lohhla_sample.sh"
  prediction="$(find "$LOHHLA_DIR" -type f -name '*HLAlossPrediction_CI*' -size +0c -print -quit)"
  [[ -n "$prediction" ]] || { echo "ERROR: LOHHLA prediction output missing" >&2; exit 5; }
  LOHHLA_TSV="$LOHHLA_DIR/hla_loh.tsv"
  "$ROOT/bin/neoag" convert-lohhla -i "$prediction" -o "$LOHHLA_TSV"
}

run_spechla() {
  TYPING_DIR="$OUTDIR/spechla_typing"
  bash "$ROOT/scripts/run_spechla_sample.sh" --bam "$TUMOR_BAM" --sample-id "$TUMOR_ID" \
    --threads "$THREADS" --outdir "$TYPING_DIR"
  typing_result="$(find "$TYPING_DIR/typing" -type f -name 'hla.result.txt' -size +0c -print -quit)"
  [[ -n "$typing_result" ]] || { echo "ERROR: tumor SpecHLA typing output missing" >&2; exit 5; }
  spechla_args=(--sample-id "$SAMPLE_ID" --typing-dir "$(dirname "$typing_result")" \
    --purity "$PURITY" --ploidy "$PLOIDY" --outdir "$OUTDIR/spechla" --force)
  bash "$ROOT/scripts/run_spechla_loh.sh" "${spechla_args[@]}"
  SPECHLA_TSV="$OUTDIR/spechla/hla_loh.tsv"
}

LOG_DIR="$OUTDIR/logs"
mkdir -p "$LOG_DIR"
LOHHLA_PID=""; SPECHLA_PID=""
LOHHLA_STATUS=0; SPECHLA_STATUS=0

if wants lohhla; then
  (
    run_lohhla
  ) > "$LOG_DIR/lohhla.log" 2>&1 & LOHHLA_PID=$!
  echo "$LOHHLA_PID" > "$OUTDIR/lohhla.pid"
fi

if wants spechla; then
  (
    run_spechla
  ) > "$LOG_DIR/spechla.log" 2>&1 & SPECHLA_PID=$!
  echo "$SPECHLA_PID" > "$OUTDIR/spechla.pid"
fi

if [[ -n "$LOHHLA_PID" ]]; then
  wait "$LOHHLA_PID" || LOHHLA_STATUS=$?
  LOHHLA_TSV="$OUTDIR/lohhla/hla_loh.tsv"
fi
if [[ -n "$SPECHLA_PID" ]]; then
  wait "$SPECHLA_PID" || SPECHLA_STATUS=$?
  SPECHLA_TSV="$OUTDIR/spechla/hla_loh.tsv"
fi

if [[ "$LOHHLA_STATUS" -ne 0 || "$SPECHLA_STATUS" -ne 0 ]]; then
  echo "ERROR: HLA LOH tool failure: LOHHLA=$LOHHLA_STATUS SpecHLA=$SPECHLA_STATUS" >&2
  exit 6
fi

if [[ "${ALLOW_SINGLE_HLA_LOH_TOOL:-0}" != "1" ]]; then
  [[ -s "$LOHHLA_TSV" ]] || {
    echo "ERROR: LOHHLA output is required before building HLA LOH consensus" >&2
    exit 6
  }
  [[ -s "$SPECHLA_TSV" ]] || {
    echo "ERROR: SpecHLA LOH output is required before building HLA LOH consensus" >&2
    exit 6
  }
fi

consensus_args=(--sample-id "$SAMPLE_ID" --outdir "$OUTDIR")
[[ -n "$LOHHLA_TSV" ]] && consensus_args+=(--lohhla "$LOHHLA_TSV")
[[ -n "$SPECHLA_TSV" ]] && consensus_args+=(--spechla "$SPECHLA_TSV")
"${PYTHON:-python3}" "$ROOT/scripts/build_hla_loh_consensus.py" "${consensus_args[@]}"

if [[ "${ALLOW_SINGLE_HLA_LOH_TOOL:-0}" != "1" ]]; then
  [[ -s "$OUTDIR/hla_loh_consensus.tsv" ]] || {
    echo "ERROR: HLA LOH consensus was not created" >&2
    exit 6
  }
  awk -F '\t' '
    NR == 1 {
      for (i=1; i<=NF; i++) h[$i]=i
      next
    }
    NR > 1 {
      rows++
      if ($(h["source_tools"]) == "lohhla;spechla") both++
    }
    END {
      exit !(rows > 0 && both > 0)
    }
  ' "$OUTDIR/hla_loh_consensus.tsv" || {
    echo "ERROR: HLA LOH consensus requires evidence from both LOHHLA and SpecHLA" >&2
    exit 6
  }
fi

{
  printf 'key\tvalue\n'
  printf 'sample_id\t%s\n' "$SAMPLE_ID"
  printf 'tumor_id\t%s\n' "$TUMOR_ID"
  printf 'normal_id\t%s\n' "$NORMAL_ID"
  printf 'purity\t%s\n' "$PURITY"
  printf 'ploidy\t%s\n' "$PLOIDY"
  printf 'tools\t%s\n' "$TOOLS"
} > "$OUTDIR/run_metadata.tsv"
date -Is > "$OUTDIR/.complete"
echo "HLA LOH multi-tool workflow completed: $OUTDIR"
