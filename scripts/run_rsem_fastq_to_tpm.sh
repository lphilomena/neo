#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_rsem_fastq_to_tpm.sh --fastq1 R1.fq.gz --fastq2 R2.fq.gz --outdir OUTDIR [options]

Options:
  --sample-id SAMPLE        Sample ID used as output prefix, default: sample
  --rsem-reference PREFIX   RSEM reference prefix; default: $RSEM_REFERENCE
  --threads N              Threads; default: $RSEM_THREADS or 8

Environment:
  RSEM_BIN                 rsem-calculate-expression executable, default: command on PATH
  RSEM_REFERENCE           RSEM reference prefix built from matching transcriptome/GTF
  RSEM_THREADS             thread count
  RSEM_ALIGNER             auto, bowtie, bowtie2, or star; default: auto
  RSEM_RUNTIME_LIBDIR      optional library directory for the selected aligner
USAGE
}

FASTQ1=""; FASTQ2=""; OUTDIR=""; SAMPLE_ID="sample"; RSEM_REFERENCE="${RSEM_REFERENCE:-}"; THREADS="${RSEM_THREADS:-8}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq1) FASTQ1="$2"; shift 2 ;;
    --fastq2) FASTQ2="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --rsem-reference) RSEM_REFERENCE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$FASTQ1" && -f "$FASTQ1" ]] || { echo "ERROR: --fastq1 missing or not a file: ${FASTQ1:-unset}" >&2; exit 3; }
[[ -n "$FASTQ2" && -f "$FASTQ2" ]] || { echo "ERROR: --fastq2 missing or not a file: ${FASTQ2:-unset}" >&2; exit 3; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required" >&2; exit 3; }
[[ -n "$RSEM_REFERENCE" ]] || { echo "ERROR: RSEM_REFERENCE/--rsem-reference is required" >&2; exit 3; }

select_rsem_bin() {
  local requested="${RSEM_BIN:-}" root candidate
  candidates=()
  [[ -n "$requested" ]] && candidates+=("$requested")
  candidate="$(command -v rsem-calculate-expression || true)"
  [[ -n "$candidate" ]] && candidates+=("$candidate")
  if [[ -n "${CONDA_EXE:-}" ]]; then
    root="$(cd "$(dirname "$CONDA_EXE")/.." 2>/dev/null && pwd -P || true)"
    [[ -n "$root" ]] && candidates+=("$root"/envs/*/bin/rsem-calculate-expression)
  fi
  [[ -n "${MAMBA_ROOT_PREFIX:-}" ]] && candidates+=("$MAMBA_ROOT_PREFIX"/envs/*/bin/rsem-calculate-expression)
  for candidate in "${candidates[@]}"; do
    [[ -x "$candidate" ]] || continue
    if "$candidate" --version >/dev/null 2>&1; then
      RSEM_BIN="$candidate"
      export PATH="$(dirname "$candidate"):$PATH"
      echo "==> RSEM executable: $RSEM_BIN" >&2
      return 0
    fi
  done
  echo "ERROR: no functional rsem-calculate-expression found; set RSEM_BIN to a relocated Conda environment executable" >&2
  return 1
}
select_rsem_bin

mkdir -p "$OUTDIR"
LOG="$OUTDIR/rsem_quant.log"
OUT_PREFIX="$OUTDIR/$SAMPLE_ID"

probe_aligner() {
  local executable="$1"
  "$executable" --version >/dev/null 2>&1
}

configure_aligner_runtime() {
  local executable="$1"
  local prefix root candidate
  probe_aligner "$executable" && return 0

  # A cached Conda executable can otherwise load the host libstdc++ and fail
  # with a GLIBCXX version error. Prefer an explicit override, then inspect
  # installed Conda environments without assuming a particular environment name.
  candidates=()
  [[ -n "${RSEM_RUNTIME_LIBDIR:-}" ]] && candidates+=("$RSEM_RUNTIME_LIBDIR")
  prefix="$(cd "$(dirname "$executable")/.." 2>/dev/null && pwd -P || true)"
  [[ -n "$prefix" ]] && candidates+=("$prefix/lib")
  prefix="$(cd "$(dirname "$RSEM_BIN")/.." 2>/dev/null && pwd -P || true)"
  [[ -n "$prefix" ]] && candidates+=("$prefix/lib")
  [[ -n "${CONDA_PREFIX:-}" ]] && candidates+=("$CONDA_PREFIX/lib")
  if [[ -n "${CONDA_EXE:-}" ]]; then
    root="$(cd "$(dirname "$CONDA_EXE")/.." 2>/dev/null && pwd -P || true)"
    [[ -n "$root" ]] && candidates+=("$root/lib" "$root"/envs/*/lib)
  fi
  [[ -n "${MAMBA_ROOT_PREFIX:-}" ]] && candidates+=("$MAMBA_ROOT_PREFIX/lib" "$MAMBA_ROOT_PREFIX"/envs/*/lib)

  for candidate in "${candidates[@]}"; do
    [[ -r "$candidate/libstdc++.so.6" ]] || continue
    if LD_LIBRARY_PATH="$candidate${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
        "$executable" --version >/dev/null 2>&1; then
      export LD_LIBRARY_PATH="$candidate${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
      echo "==> RSEM aligner runtime: $candidate" >&2
      return 0
    fi
  done
  echo "ERROR: aligner exists but is not functional: $executable" >&2
  "$executable" --version >&2 || true
  return 1
}

ALIGNER_ARGS=()
case "${RSEM_ALIGNER:-auto}" in
  auto)
    if command -v bowtie >/dev/null 2>&1; then
      :
    elif command -v bowtie2 >/dev/null 2>&1; then
      ALIGNER_ARGS+=(--bowtie2)
    elif command -v STAR >/dev/null 2>&1; then
      ALIGNER_ARGS+=(--star)
    else
      echo "ERROR: no supported RSEM aligner found (bowtie, bowtie2, or STAR)" >&2
      exit 3
    fi
    ;;
  bowtie) command -v bowtie >/dev/null 2>&1 || { echo "ERROR: bowtie not found" >&2; exit 3; } ;;
  bowtie2) command -v bowtie2 >/dev/null 2>&1 || { echo "ERROR: bowtie2 not found" >&2; exit 3; }; ALIGNER_ARGS+=(--bowtie2) ;;
  star) command -v STAR >/dev/null 2>&1 || { echo "ERROR: STAR not found" >&2; exit 3; }; ALIGNER_ARGS+=(--star) ;;
  *) echo "ERROR: unsupported RSEM_ALIGNER=${RSEM_ALIGNER}" >&2; exit 3 ;;
esac

case "${RSEM_ALIGNER:-auto}" in
  star) SELECTED_ALIGNER="$(command -v STAR)" ;;
  bowtie) SELECTED_ALIGNER="$(command -v bowtie)" ;;
  bowtie2) SELECTED_ALIGNER="$(command -v bowtie2)" ;;
  auto)
    if command -v bowtie >/dev/null 2>&1; then
      SELECTED_ALIGNER="$(command -v bowtie)"
    elif command -v bowtie2 >/dev/null 2>&1; then
      SELECTED_ALIGNER="$(command -v bowtie2)"
    else
      SELECTED_ALIGNER="$(command -v STAR)"
    fi
    ;;
esac
configure_aligner_runtime "$SELECTED_ALIGNER"

"$RSEM_BIN" --paired-end -p "$THREADS" --estimate-rspd --append-names "${ALIGNER_ARGS[@]}" \
  "$FASTQ1" "$FASTQ2" "$RSEM_REFERENCE" "$OUT_PREFIX" >"$LOG" 2>&1

GENES_RESULTS="$OUT_PREFIX.genes.results"
[[ -f "$GENES_RESULTS" ]] || { echo "ERROR: expected RSEM genes result missing: $GENES_RESULTS" >&2; exit 4; }

awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++){if($i=="gene_id") gid=i; if($i=="TPM") tpm=i} print "gene_id","tpm"; next} gid && tpm {print $gid,$tpm}' \
  "$GENES_RESULTS" > "$OUTDIR/gene_tpm.tsv"

ISOFORM_RESULTS="$OUT_PREFIX.isoforms.results"
[[ -f "$ISOFORM_RESULTS" ]] || { echo "ERROR: expected RSEM isoform result missing: $ISOFORM_RESULTS" >&2; exit 4; }
awk 'BEGIN{FS=OFS="\t"} NR==1{for(i=1;i<=NF;i++){if($i=="transcript_id") tid=i; if($i=="TPM") tpm=i} print "transcript_id","tpm"; next} tid && tpm {print $tid,$tpm}' \
  "$ISOFORM_RESULTS" > "$OUTDIR/transcript_tpm.tsv"

cat > "$OUTDIR/rna_fastq_to_tpm.summary.json" <<JSON
{
  "method": "rsem",
  "sample_id": "$SAMPLE_ID",
  "fastq1": "$FASTQ1",
  "fastq2": "$FASTQ2",
  "rsem_reference": "$RSEM_REFERENCE",
  "genes_results": "$GENES_RESULTS",
  "gene_tpm": "$OUTDIR/gene_tpm.tsv",
  "transcript_tpm": "$OUTDIR/transcript_tpm.tsv",
  "log": "$LOG"
}
JSON
