#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 --rscript R --form-script SCRIPT --split-script SCRIPT --introns FILE --txdb FILE --bsgenome-package PKG --bsgenome-object OBJ --functions FILE --outdir DIR --chunks N" >&2
}

RSCRIPT=""; FORM_SCRIPT=""; SPLIT_SCRIPT=""; INTRONS=""; TXDB=""
BSGENOME_PACKAGE=""; BSGENOME_OBJECT=""; FUNCTIONS=""; OUTDIR=""; CHUNKS=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --rscript) RSCRIPT="$2"; shift 2 ;;
    --form-script) FORM_SCRIPT="$2"; shift 2 ;;
    --split-script) SPLIT_SCRIPT="$2"; shift 2 ;;
    --introns) INTRONS="$2"; shift 2 ;;
    --txdb) TXDB="$2"; shift 2 ;;
    --bsgenome-package) BSGENOME_PACKAGE="$2"; shift 2 ;;
    --bsgenome-object) BSGENOME_OBJECT="$2"; shift 2 ;;
    --functions) FUNCTIONS="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --chunks) CHUNKS="$2"; shift 2 ;;
    *) usage; exit 2 ;;
  esac
done
for value in "$RSCRIPT" "$FORM_SCRIPT" "$SPLIT_SCRIPT" "$INTRONS" "$TXDB" \
             "$BSGENOME_PACKAGE" "$BSGENOME_OBJECT" "$FUNCTIONS" "$OUTDIR"; do
  [[ -n "$value" ]] || { usage; exit 2; }
done

mkdir -p "$OUTDIR/chunks" "$OUTDIR/logs"
"$RSCRIPT" "$SPLIT_SCRIPT" --input "$INTRONS" --outdir "$OUTDIR/chunks" --chunks "$CHUNKS"

pids=()
cleanup_children() {
  if [[ ${#pids[@]} -gt 0 ]]; then
    kill "${pids[@]}" 2>/dev/null || true
    wait "${pids[@]}" 2>/dev/null || true
  fi
}
trap cleanup_children INT TERM EXIT

for chunk in "$OUTDIR"/chunks/introns.chunk_*.rds; do
  name="$(basename "$chunk" .rds)"
  prefix="$OUTDIR/$name"
  "$RSCRIPT" "$FORM_SCRIPT" -o "$prefix" -t "$TXDB" -j "$chunk" \
    -b "$BSGENOME_OBJECT" -p "$BSGENOME_PACKAGE" -f "$FUNCTIONS" \
    > "$OUTDIR/logs/$name.log" 2>&1 &
  pids+=("$!")
done

status=0
for pid in "${pids[@]}"; do
  wait "$pid" || status=1
done
pids=()
trap - INT TERM EXIT
[[ $status -eq 0 ]] || { echo "ERROR: one or more SpliceMutr transcript chunks failed" >&2; exit 1; }

find "$OUTDIR" -maxdepth 1 -type f -name '*_data_splicemutr.rds' -print | sort > "$OUTDIR/formed_metadata_files.txt"
[[ -s "$OUTDIR/formed_metadata_files.txt" ]] || { echo "ERROR: no formed transcript metadata was produced" >&2; exit 1; }
touch "$OUTDIR/form_transcripts.done"
