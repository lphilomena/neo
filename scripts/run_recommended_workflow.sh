#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
usage() {
  cat <<USAGE
Usage: $0 --input-dir DIR --sample-id ID --tumor-id ID --normal-id ID \
  --tumor-bam BAM --normal-bam BAM --somatic-vcf VCF --outdir DIR [options]

Options:
  --steps LIST              Comma-separated 0,1,2,3,4 (default: all)
  --threads N               General threads (default: 8)
  --expression FILE         Gene TPM table
  --transcript-expression FILE
  --rna-vaf FILE            RNA ALT/depth/VAF table
  --min-free-gb N           Required output filesystem space (default: 300)
  --execute                 Execute; otherwise print the complete plan
USAGE
}

INPUT_DIR=""; SAMPLE_ID=""; TUMOR_ID=""; NORMAL_ID=""; TUMOR_BAM=""; NORMAL_BAM=""; SOMATIC_VCF=""; OUTDIR=""
STEPS=0,1,2,3,4; THREADS=8; EXPRESSION=""; TRANSCRIPT_EXPRESSION=""; RNA_VAF=""; MIN_FREE_GB=300; EXECUTE=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir) INPUT_DIR="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --tumor-id) TUMOR_ID="$2"; shift 2 ;;
    --normal-id) NORMAL_ID="$2"; shift 2 ;;
    --tumor-bam) TUMOR_BAM="$2"; shift 2 ;;
    --normal-bam) NORMAL_BAM="$2"; shift 2 ;;
    --somatic-vcf) SOMATIC_VCF="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --steps) STEPS="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --expression) EXPRESSION="$2"; shift 2 ;;
    --transcript-expression) TRANSCRIPT_EXPRESSION="$2"; shift 2 ;;
    --rna-vaf) RNA_VAF="$2"; shift 2 ;;
    --min-free-gb) MIN_FREE_GB="$2"; shift 2 ;;
    --execute) EXECUTE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage; exit 2 ;;
  esac
done
for value in "$INPUT_DIR" "$SAMPLE_ID" "$TUMOR_ID" "$NORMAL_ID" "$TUMOR_BAM" "$NORMAL_BAM" "$SOMATIC_VCF" "$OUTDIR"; do
  [[ -n "$value" ]] || { usage; exit 2; }
done
[[ "$OUTDIR" != /root/* ]] || { echo "ERROR: WGS workflow output must not be placed under /root" >&2; exit 3; }
mkdir -p "$OUTDIR/logs" "$OUTDIR/tmp"

source /root/neo/envs/activate_neoag_production_refs.sh
[[ -f "$ROOT/conf/tools.env.local.sh" ]] && source "$ROOT/conf/tools.env.local.sh"
export NEOAG_PROJECT_ROOT="$ROOT" PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}" TMPDIR="$OUTDIR/tmp"
PYTHON="${NEOAG_CONDA_BASE}/envs/neoag-tools/bin/python"
[[ -x "$PYTHON" ]] || PYTHON=python3
FACETS_VCF="${FACETS_SNP_VCF:-${NEOAG_REF_BUNDLE}/data/facets/reference/1000G_omni2.5.hg38.biallelic.vcf.gz}"
BAM_MATCH_REF="${BAM_MATCHER_REFERENCE:-$NEOAG_REFERENCE_FASTA}"

available_kb="$(df -Pk "$OUTDIR" | awk 'NR==2 {print $4}')"
required_kb=$((MIN_FREE_GB * 1024 * 1024))
(( available_kb >= required_kb )) || { echo "ERROR: output filesystem has less than ${MIN_FREE_GB} GiB free: $OUTDIR" >&2; exit 3; }

wants() { [[ ",${STEPS}," == *",$1,"* ]]; }
run() {
  echo "+ $*"
  [[ "$EXECUTE" == 1 ]] || return 0
  "$@"
}

QC="$OUTDIR/input_qc"
HLA="$OUTDIR/hla_typing"
CNV="$OUTDIR/purity_cnv"
LOH="$OUTDIR/hla_loh"
RANK="$OUTDIR/neoantigen_ranking"

if wants 0; then
  run "$PYTHON" "$ROOT/scripts/run_case_input_qc.py" --input-dir "$INPUT_DIR" --tumor-bam "$TUMOR_BAM" --normal-bam "$NORMAL_BAM" \
    --somatic-vcf "$SOMATIC_VCF" --reference "$BAM_MATCH_REF" --bam-matcher-loci "$FACETS_VCF" --threads "$THREADS" --outdir "$QC"
fi

if wants 1; then
  run bash "$ROOT/scripts/run_spechla_sample.sh" --bam "$NORMAL_BAM" --sample-id "$NORMAL_ID" --threads "$THREADS" --outdir "$HLA/spechla"
  run bash "$ROOT/scripts/run_hla_la_sample.sh" --bam "$NORMAL_BAM" --sample-id "$NORMAL_ID" --threads "$THREADS" --outdir "$HLA/hla-la"
  run bash "$ROOT/scripts/run_optitype_sample.sh" --bam "$NORMAL_BAM" --sample-id "$NORMAL_ID" --threads "$THREADS" --outdir "$HLA/optitype"
  run "$PYTHON" -m neoag.agent_skills.hla_typing_compare --result-dir "$HLA/spechla" --result-dir "$HLA/hla-la" --result-dir "$HLA/optitype" --outdir "$HLA"
fi

if wants 2; then
  if [[ "$EXECUTE" == 1 ]]; then
    PATIENT_ID="$SAMPLE_ID" TUMOR_NAME="$TUMOR_ID" NORMAL_NAME="$NORMAL_ID" TUMOR_BAM="$TUMOR_BAM" NORMAL_BAM="$NORMAL_BAM" \
      OUTDIR="$CNV/facets" FACETS_SNP_VCF="$FACETS_VCF" FACETS_TARGET_ROWS=1000000 FACETS_CVAL_PRE=50 FACETS_CVAL_PROC=300 FACETS_MIN_NHET=10 \
      bash "$ROOT/scripts/run_facets_omni2p5_snponly_downsample.sh"
    date -Is > "$CNV/facets/.complete"
    SAMPLE_ID="$SAMPLE_ID" TUMOR_BAM="$TUMOR_BAM" NORMAL_BAM="$NORMAL_BAM" OUTDIR="$CNV/sequenza" \
      REF_FASTA="$SEQUENZA_FASTA" GC_WIGGLE="$SEQUENZA_GC_WIG" CHUNK_JOBS=3 bash "$ROOT/scripts/run_sequenza_sample_by_chrom.sh"
    date -Is > "$CNV/sequenza/.complete"
  else
    echo "+ FACETS omni2p5 CVAL_PRE=50 CVAL_PROC=300 MIN_NHET=10 TARGET_ROWS=1000000 -> $CNV/facets"
    echo "+ Sequenza tumor-normal -> $CNV/sequenza"
  fi
  run bash "$ROOT/scripts/run_purple_sample.sh" --sample-id "$SAMPLE_ID" --tumor-id "$TUMOR_ID" --normal-id "$NORMAL_ID" \
    --tumor-bam "$TUMOR_BAM" --normal-bam "$NORMAL_BAM" --threads "$THREADS" --outdir "$CNV/purple"
  run "$PYTHON" -m neoag.agent_skills.purity_cnv_review --sample-id "$TUMOR_ID" --result-dir "$CNV/facets" --result-dir "$CNV/sequenza" --result-dir "$CNV/purple" --outdir "$CNV"
fi

if wants 3; then
  if [[ "$EXECUTE" == 1 ]]; then
    purity="$(awk -F '\t' 'NR==2 {print $2}' "$CNV/recommended_purity.tsv")"
    ploidy="$(awk -F '\t' 'NR==2 {print $3}' "$CNV/recommended_purity.tsv")"
    [[ "$purity" =~ ^[0-9]+([.][0-9]+)?$ && "$ploidy" =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "ERROR: invalid recommended purity/ploidy" >&2; exit 4; }
    mkdir -p "$LOH/lohhla"
    { printf '\ttumorPurity\ttumorPloidy\n'; printf '%s\t%s\t%s\n' "$TUMOR_ID" "$purity" "$ploidy"; } > "$LOH/lohhla/lohhla_copy_number_input.tsv"
    PATIENT_ID="$SAMPLE_ID" TUMOR_SAMPLE_ID="$TUMOR_ID" NORMAL_SAMPLE_ID="$NORMAL_ID" TUMOR_BAM="$TUMOR_BAM" NORMAL_BAM="$NORMAL_BAM" \
      HLA_FILE="$HLA/recommended_hla.txt" OUTDIR="$LOH/lohhla" LOHHLA_NAS_ROOT="$LOH/lohhla/work" COPYNUM_LOC="$LOH/lohhla/lohhla_copy_number_input.tsv" \
      POLYSOLVER_THREADS="$THREADS" bash "$ROOT/scripts/run_lohhla_sample.sh"
    prediction="$(find "$LOH/lohhla" -type f -name '*HLAlossPrediction_CI*' -size +0c -print -quit)"
    "$ROOT/bin/neoag" convert-lohhla -i "$prediction" -o "$LOH/lohhla/hla_loh.tsv"
    bash "$ROOT/scripts/run_spechla_sample.sh" --bam "$TUMOR_BAM" --sample-id "$TUMOR_ID" --threads "$THREADS" --outdir "$LOH/spechla_typing"
    typing_result="$(find "$LOH/spechla_typing/typing" -type f -name 'hla.result.txt' -size +0c -print -quit)"
    typing_dir="$(dirname "$typing_result")"
    bash "$ROOT/scripts/run_spechla_loh.sh" --sample-id "$SAMPLE_ID" --typing-dir "$typing_dir" --purity "$purity" --ploidy "$ploidy" \
      --lohhla-hla-loh "$LOH/lohhla/hla_loh.tsv" --outdir "$LOH/spechla"
    "$PYTHON" "$ROOT/scripts/build_hla_loh_consensus.py" --sample-id "$SAMPLE_ID" --lohhla "$LOH/lohhla/hla_loh.tsv" \
      --spechla "$LOH/spechla/hla_loh.tsv" --outdir "$LOH"
  else
    echo "+ LOHHLA with recommended HLA and purity/ploidy -> $LOH/lohhla"
    echo "+ tumor SpecHLA typing + SpecHLA LOH -> $LOH/spechla"
    echo "+ explicit LOST/RETAINED/UNASSESSED/CONFLICT consensus -> $LOH"
  fi
fi

if wants 4; then
  config="$RANK/run_full.toml"
  config_args=(--sample-id "$SAMPLE_ID" --outdir "$RANK/run-full" --vcf "$SOMATIC_VCF" --tumor-sample-name "$TUMOR_ID" --normal-sample-name "$NORMAL_ID" \
    --hla-file "$HLA/recommended_hla.txt" --purity "$CNV/recommended_purity.tsv" --cnv "$CNV/recommended_cnv_segments.tsv" \
    --hla-loh "$LOH/recommended_hla_loh.tsv" --output "$config")
  [[ -n "$EXPRESSION" ]] && config_args+=(--expression "$EXPRESSION")
  [[ -n "$TRANSCRIPT_EXPRESSION" ]] && config_args+=(--transcript-expression "$TRANSCRIPT_EXPRESSION")
  [[ -n "$RNA_VAF" ]] && config_args+=(--rna-vaf "$RNA_VAF")
  run "$PYTHON" "$ROOT/scripts/generate_recommended_run_config.py" "${config_args[@]}"
  run "$PYTHON" -m neoag.cli run-full --config "$config" --outdir "$RANK/run-full"
fi

if [[ "$EXECUTE" == 1 ]]; then
  date -Is > "$OUTDIR/.complete"
  echo "Recommended workflow execution finished: $OUTDIR"
else
  echo "Recommended workflow dry-run plan finished: $OUTDIR"
fi
