#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CONFIG="$ROOT/configs/workflows/longrna_splice_profile.yaml"
SAMPLE_ID=""; INPUT_DIR=""; WORKDIR=""; DRY_RUN=0
usage() { echo "Usage: $0 --sample-id ID --input-dir FASTQ_DIR --workdir DIR [--config FILE] [--dry-run]"; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG=$2; shift 2;;
    --sample-id) SAMPLE_ID=$2; shift 2;;
    --input-dir) INPUT_DIR=$2; shift 2;;
    --workdir) WORKDIR=$2; shift 2;;
    --dry-run) DRY_RUN=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2;;
  esac
done
[[ -s "$CONFIG" ]] || { echo "ERROR: missing profile: $CONFIG" >&2; exit 1; }
[[ -n "$SAMPLE_ID" && -n "$WORKDIR" ]] || { usage >&2; exit 2; }
[[ -d "$INPUT_DIR" ]] || { echo "ERROR: input directory not found: $INPUT_DIR" >&2; exit 1; }
mapfile -t FASTQS < <(find "$INPUT_DIR" -maxdepth 1 -type f \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \) -print | sort)
(( ${#FASTQS[@]} > 0 )) || { echo "ERROR: no FASTQ files found" >&2; exit 1; }
mkdir -p "$WORKDIR"/{logs,isoquant,sqanti3,translation,splice_candidates,cross_validation,presentation}
STATUS="$WORKDIR/longrna_splice_status.tsv"
printf 'stage\tstatus\tupdated_at\tdetail\n' > "$STATUS"
record() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$(date -Iseconds)" "$3" >> "$STATUS"; }
record fastq_validation PASS "files=${#FASTQS[@]}"
if (( DRY_RUN )); then
  printf 'profile=%s\nsample_id=%s\ninput_dir=%s\nworkdir=%s\n' "$CONFIG" "$SAMPLE_ID" "$INPUT_DIR" "$WORKDIR"
  printf '%s\n' fastq_validation isoquant_transcript_reconstruction sqanti3_structure_orf_annotation td2_protein_translation snaf_splicemutr_cross_validation junction_peptide_generation mhc_presentation_prediction comprehensive_evidence
  exit 0
fi
cat > "$WORKDIR/README.workflow.txt" <<EOF
sample_id=$SAMPLE_ID
profile=$CONFIG
fastq_files=${#FASTQS[@]}
status=$STATUS
EOF
record workflow_initialized PASS "canonical output root created"
for stage in isoquant_transcript_reconstruction sqanti3_structure_orf_annotation td2_protein_translation snaf_splicemutr_cross_validation junction_peptide_generation mhc_presentation_prediction comprehensive_evidence; do
  record "$stage" PENDING "run the installed stage driver"
done
echo "$STATUS"
