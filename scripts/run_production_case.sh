#!/usr/bin/env bash
set -euo pipefail

# Generic production wrapper.
# Aligned with the validated production launcher:
# - Do NOT pass --skip-netmhcstabpan
# - Require local DTU NetMHCstabpan
# - IEDB fallback behavior is controlled by profile:
#   [immunogenicity] use_iedb_fallback=true

usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_production_case.sh \
    --sample-id <sample_id> \
    --case-root <case_root> \
    --outdir <outdir> \
    --somatic-vcf <somatic.vcf.gz> \
    [--project-root <neo_repo_root>] \
    [--profile <profiles/*.toml>] \
    [--evidence-consensus-rules <configs/ranking/*.toml>] \
    [--asset-root <liup_neodata4git>] \
    [--reference-fasta <GRCh38.fasta>] \
    [--gencode-gtf <gencode.gtf>] \
    [--sequenza <result_file_or_dir>] \
    [--purple <result_file_or_dir>] \
    [--rna-fastq1 <R1.fq.gz[,lane2_R1.fq.gz]>] \
    [--rna-fastq2 <R2.fq.gz[,lane2_R2.fq.gz]>] \
    [--rna-bam <sorted_rna.bam> | --rna-vaf <rna_alt_vaf.tsv>] \
    [--star-index <GRCh38_STAR_index>] \
    [--star-executable <STAR>] \
    [--samtools-executable <samtools>] \
    [--rna-threads <N>] \
    [--pred-deps <predictor_deps_dir>] \
    [--netmhcpan-home <netMHCpan_home>] \
    [--netmhcstabpan-home <netMHCstabpan_home>] \
    [--python <python_bin>]

Environment defaults may be supplied through conf/tools.env.local.sh or:
  NEOAG_PYTHON, NEOAG_ASSET_ROOT, NEOAG_PREDICTOR_DEPS,
  NETMHCPAN_HOME, and NETMHCSTABPAN_HOME.
EOF
}

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${NEOAG_PYTHON:-$(command -v python3 || command -v python || true)}"
PROFILE="profiles/sarcoma_rna_supported_v2_provisional.toml"
EVIDENCE_CONSENSUS_RULES="configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml"
ASSET="${NEOAG_ASSET_ROOT:-${NEOAG_TOOLS_ROOT:-}}"
PRED_DEPS="${NEOAG_PREDICTOR_DEPS:-${NEOAG_TOOL_QUARANTINE:-}}"
NETMHCPAN_HOME_DEFAULT="${NETMHCPAN_HOME:-}"
NETMHCSTABPAN_HOME="${NETMHCSTABPAN_HOME:-}"
CLI_ASSET=""
CLI_PRED_DEPS=""
CLI_NETMHCPAN_HOME=""
CLI_NETMHCSTABPAN_HOME=""

SAMPLE_ID=""
CASE_ROOT=""
OUTDIR=""
SOMATIC_VCF=""
REFERENCE_FASTA=""
GENCODE_GTF=""
SEQUENZA=""
PURPLE=""
RNA_FASTQ1=""
RNA_FASTQ2=""
RNA_BAM=""
RNA_VAF=""
STAR_INDEX=""
STAR_EXECUTABLE=""
SAMTOOLS_EXECUTABLE="samtools"
RNA_THREADS=16

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --case-root) CASE_ROOT="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --somatic-vcf) SOMATIC_VCF="$2"; shift 2 ;;
    --project-root) PROJECT_ROOT="$2"; shift 2 ;;
    --profile) PROFILE="$2"; shift 2 ;;
    --evidence-consensus-rules) EVIDENCE_CONSENSUS_RULES="$2"; shift 2 ;;
    --asset-root) ASSET="$2"; CLI_ASSET="$2"; shift 2 ;;
    --reference-fasta) REFERENCE_FASTA="$2"; shift 2 ;;
    --gencode-gtf) GENCODE_GTF="$2"; shift 2 ;;
    --sequenza) SEQUENZA="$2"; shift 2 ;;
    --purple) PURPLE="$2"; shift 2 ;;
    --rna-fastq1) RNA_FASTQ1="$2"; shift 2 ;;
    --rna-fastq2) RNA_FASTQ2="$2"; shift 2 ;;
    --rna-bam) RNA_BAM="$2"; shift 2 ;;
    --rna-vaf) RNA_VAF="$2"; shift 2 ;;
    --star-index) STAR_INDEX="$2"; shift 2 ;;
    --star-executable) STAR_EXECUTABLE="$2"; shift 2 ;;
    --samtools-executable) SAMTOOLS_EXECUTABLE="$2"; shift 2 ;;
    --rna-threads) RNA_THREADS="$2"; shift 2 ;;
    --pred-deps) PRED_DEPS="$2"; CLI_PRED_DEPS="$2"; shift 2 ;;
    --netmhcpan-home) NETMHCPAN_HOME_DEFAULT="$2"; CLI_NETMHCPAN_HOME="$2"; shift 2 ;;
    --netmhcstabpan-home) NETMHCSTABPAN_HOME="$2"; CLI_NETMHCSTABPAN_HOME="$2"; shift 2 ;;
    --python) PY="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2 ;;
  esac
done

[[ -n "$SAMPLE_ID" ]] || { echo "missing --sample-id" >&2; exit 2; }
[[ -n "$CASE_ROOT" ]] || { echo "missing --case-root" >&2; exit 2; }
[[ -n "$OUTDIR" ]] || { echo "missing --outdir" >&2; exit 2; }
[[ -n "$SOMATIC_VCF" ]] || { echo "missing --somatic-vcf" >&2; exit 2; }
[[ -x "$PY" ]] || { echo "python not executable: $PY" >&2; exit 2; }
[[ -d "$PROJECT_ROOT" ]] || { echo "project root missing: $PROJECT_ROOT" >&2; exit 2; }
[[ "$RNA_THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "--rna-threads must be a positive integer" >&2; exit 2; }
[[ -z "$RNA_FASTQ1" && -z "$RNA_FASTQ2" ]] || {
  [[ -n "$RNA_FASTQ1" && -n "$RNA_FASTQ2" ]] || {
    echo "--rna-fastq1 and --rna-fastq2 must be supplied together" >&2
    exit 2
  }
}
rna_input_modes=0
[[ -n "$RNA_FASTQ1" ]] && ((rna_input_modes+=1))
[[ -n "$RNA_BAM" ]] && ((rna_input_modes+=1))
[[ -n "$RNA_VAF" ]] && ((rna_input_modes+=1))
((rna_input_modes <= 1)) || {
  echo "Use only one RNA allele-evidence input mode: FASTQ pair, BAM, or existing RNA VAF" >&2
  exit 2
}
if [[ -n "$RNA_FASTQ1" ]]; then
  [[ -n "$STAR_INDEX" ]] || { echo "RNA FASTQ mode requires --star-index" >&2; exit 2; }
  [[ -n "$GENCODE_GTF" ]] || { echo "RNA FASTQ mode requires --gencode-gtf" >&2; exit 2; }
fi

cd "$PROJECT_ROOT"
mkdir -p "$OUTDIR/manifest" "$OUTDIR/logs" "$OUTDIR/tools"

if [[ -f "$PROJECT_ROOT/conf/tools.env.sh" ]]; then
  # shellcheck disable=SC1091
  source "$PROJECT_ROOT/conf/tools.env.sh"
fi
if [[ -f "$PROJECT_ROOT/conf/tools.env.local.sh" ]]; then
  # shellcheck disable=SC1091
  source "$PROJECT_ROOT/conf/tools.env.local.sh"
fi

# Command-line values take precedence; otherwise accept values loaded from the
# project tool configuration.
ASSET="${ASSET:-${NEOAG_ASSET_ROOT:-${NEOAG_TOOLS_ROOT:-}}}"
PRED_DEPS="${PRED_DEPS:-${NEOAG_PREDICTOR_DEPS:-${NEOAG_TOOL_QUARANTINE:-}}}"
NETMHCPAN_HOME_DEFAULT="${NETMHCPAN_HOME_DEFAULT:-${NETMHCPAN_HOME:-}}"
NETMHCSTABPAN_HOME="${NETMHCSTABPAN_HOME:-}"
[[ -n "$CLI_ASSET" ]] && ASSET="$CLI_ASSET"
[[ -n "$CLI_PRED_DEPS" ]] && PRED_DEPS="$CLI_PRED_DEPS"
[[ -n "$CLI_NETMHCPAN_HOME" ]] && NETMHCPAN_HOME_DEFAULT="$CLI_NETMHCPAN_HOME"
[[ -n "$CLI_NETMHCSTABPAN_HOME" ]] && NETMHCSTABPAN_HOME="$CLI_NETMHCSTABPAN_HOME"

PROFILE_PATH="$PROFILE"
[[ "$PROFILE_PATH" = /* ]] || PROFILE_PATH="$PROJECT_ROOT/$PROFILE"
[[ -f "$PROFILE_PATH" ]] || { echo "profile missing: $PROFILE_PATH" >&2; exit 2; }
CONSENSUS_RULES_PATH="$EVIDENCE_CONSENSUS_RULES"
[[ "$CONSENSUS_RULES_PATH" = /* ]] || CONSENSUS_RULES_PATH="$PROJECT_ROOT/$EVIDENCE_CONSENSUS_RULES"
[[ -f "$CONSENSUS_RULES_PATH" ]] || { echo "evidence-consensus rules missing: $CONSENSUS_RULES_PATH" >&2; exit 2; }

STABPAN_BIN="${NETMHCSTABPAN_HOME}/Linux_x86_64/bin/netMHCstabpan"
[[ -x "$STABPAN_BIN" && -d "${NETMHCSTABPAN_HOME}/data" ]] || {
  echo "NetMHCstabpan DTU tree required: ${NETMHCSTABPAN_HOME}" >&2
  echo "Need Linux_x86_64/bin/netMHCstabpan + data/" >&2
  exit 2
}

add_if() {
  local flag="$1" path="$2"
  if [[ -e "$path" ]]; then
    GEN_ARGS+=("$flag" "$path")
  fi
}

GEN_ARGS=(
  --project-root "$PROJECT_ROOT"
  --sample-id "$SAMPLE_ID"
  --profile "$PROFILE_PATH"
  --evidence-consensus-rules "$CONSENSUS_RULES_PATH"
  --outdir "$OUTDIR"
  --output "$OUTDIR/manifest/production.results.toml"
  --somatic-vcf "$SOMATIC_VCF"
)

add_if --optitype "$CASE_ROOT/hla/optitype/${SAMPLE_ID}_blood_result.tsv"
add_if --optitype "$CASE_ROOT/hla/optitype/${SAMPLE_ID}_result.tsv"
add_if --spechla-typing "$CASE_ROOT/hla/spechla/typing/normal/${SAMPLE_ID}_blood/hla.result.txt"
add_if --hla-la "$CASE_ROOT/hla/hla_la/working/${SAMPLE_ID}_blood/hla/R1_bestguess_G.txt"
add_if --facets "$CASE_ROOT/facets/omni2p5_snponly_downsample"
add_if --ascat "$CASE_ROOT/ascat"
if [[ -n "$SEQUENZA" ]]; then GEN_ARGS+=(--sequenza "$SEQUENZA"); else add_if --sequenza "$CASE_ROOT/sequenza"; fi
if [[ -n "$PURPLE" ]]; then GEN_ARGS+=(--purple "$PURPLE"); else add_if --purple "$CASE_ROOT/purple"; fi
add_if --purity "$CASE_ROOT/evidence/purity.tsv"
add_if --cnv "$CASE_ROOT/evidence/cnv_segments.tsv"
add_if --lohhla "$CASE_ROOT/evidence/hla_loh.tsv"
add_if --spechla-loh "$CASE_ROOT/evidence/hla_loh.spechla.tsv"
add_if --expression "$CASE_ROOT/short-rna/evidence/gene_expression.tsv"
add_if --transcript-expression "$CASE_ROOT/short-rna/evidence/transcript_quant.sf"
add_if --easyfuse "$CASE_ROOT/short-rna/evidence/easyfuse.fusions.pass.csv"
add_if --star-fusion "$CASE_ROOT/short-rna/evidence/star-fusion.fusion_predictions.tsv"
add_if --arriba "$CASE_ROOT/short-rna/evidence/arriba.fusions.tsv"
add_if --junctions "$CASE_ROOT/short-rna/evidence/regtools_junctions.tsv"
add_if --snaf "$CASE_ROOT/short-rna/snaf/snaf_candidates.tsv"
add_if --splicemutr "$CASE_ROOT/short-rna/splicemutr"

if [[ -n "$ASSET" ]]; then
  add_if --normal-junctions "$ASSET/data/normal/junctions/normal_junctions.GRCh38.tsv.gz"
  add_if --normal-expression "$ASSET/data/normal/expression/normal_expression.gtex_v11_hpa_hspc.tsv"
  add_if --normal-hla-ligands "$ASSET/data/normal/ligandome/normal_ms_ligands.tsv"
  add_if --reference-proteome "$ASSET/data/normal/proteome/gencode.v49.pc_translations.clean.fa"
  add_if --netchop-executable "$ASSET/data/predictors/netchop/netchop-3.1/Linux_x86_64/bin/netChop"
  add_if --netchop-home "$ASSET/data/predictors/netchop/netchop-3.1"
fi

# Explicit reference and RNA inputs are passed through unchanged. They are not
# inferred from patient-specific paths because build/version mismatches are unsafe.
[[ -n "$REFERENCE_FASTA" ]] && GEN_ARGS+=(--reference-fasta "$REFERENCE_FASTA")
[[ -n "$GENCODE_GTF" ]] && GEN_ARGS+=(--gencode-gtf "$GENCODE_GTF")
[[ -n "$RNA_FASTQ1" ]] && GEN_ARGS+=(--rna-fastq1 "$RNA_FASTQ1" --rna-fastq2 "$RNA_FASTQ2")
[[ -n "$RNA_BAM" ]] && GEN_ARGS+=(--rna-bam "$RNA_BAM")
[[ -n "$RNA_VAF" ]] && GEN_ARGS+=(--rna-vaf "$RNA_VAF")
[[ -n "$STAR_INDEX" ]] && GEN_ARGS+=(--star-index "$STAR_INDEX")
[[ -n "$STAR_EXECUTABLE" ]] && GEN_ARGS+=(--star-executable "$STAR_EXECUTABLE")
GEN_ARGS+=(--samtools-executable "$SAMTOOLS_EXECUTABLE" --rna-threads "$RNA_THREADS")

echo "[INFO] generate manifest: $OUTDIR/manifest/production.results.toml"
"$PY" scripts/generate_production_from_results_manifest.py "${GEN_ARGS[@]}"

[[ -n "$ASSET" ]] && export NEOAG_TOOLS_ROOT="$ASSET"
[[ -n "$PRED_DEPS" ]] && export NEOAG_TOOL_QUARANTINE="$PRED_DEPS"
export NEOAG_FORCE_CPU=1
export NEOAG_PRIME_JOBS=4
export NEOAG_NETMHCPAN_LOCAL_CHUNK_SIZE=5000

if [[ -z "$NETMHCPAN_HOME_DEFAULT" && -n "$ASSET" ]]; then
  NETMHCPAN_HOME_DEFAULT="$ASSET/data/predictors/netMHCpan"
fi
[[ -n "$NETMHCPAN_HOME_DEFAULT" ]] && export NETMHCPAN_HOME="$NETMHCPAN_HOME_DEFAULT"
[[ -n "${NETMHCPAN_HOME:-}" ]] && export NETMHCpan="$NETMHCPAN_HOME"
export NEOAG_NETMHCPAN_BIN="${NEOAG_NETMHCPAN_BIN:-$OUTDIR/tools/netMHCpan-local}"
export NETMHCSTABPAN_HOME="$NETMHCSTABPAN_HOME"
export NETMHCSTABPAN_BIN="${NETMHCSTABPAN_BIN:-$STABPAN_BIN}"

pick_tool_dir() {
  local current="$1" deps_name="$2" asset_name="$3" marker="$4"
  if [[ -n "$current" && -e "$current/$marker" ]]; then printf '%s' "$current"; return; fi
  if [[ -n "$PRED_DEPS" && -e "$PRED_DEPS/$deps_name/$marker" ]]; then printf '%s' "$PRED_DEPS/$deps_name"; return; fi
  if [[ -n "$ASSET" && -e "$ASSET/data/predictors/$asset_name/$marker" ]]; then printf '%s' "$ASSET/data/predictors/$asset_name"; return; fi
  printf '%s' "$current"
}

PRIME_HOME="$(pick_tool_dir "${PRIME_HOME:-}" prime prime PRIME)"
MIXMHCPRED_HOME="$(pick_tool_dir "${MIXMHCPRED_HOME:-}" mixMHCpred_install mixMHCpred_install MixMHCpred)"
BIGMHC_DIR="$(pick_tool_dir "${BIGMHC_DIR:-}" bigmhc bigmhc src/predict.py)"
DEEPIMMUNO_DIR="$(pick_tool_dir "${DEEPIMMUNO_DIR:-}" DeepImmuno DeepImmuno deepimmuno-cnn.py)"
[[ -n "$PRIME_HOME" ]] && export PRIME_HOME NEOAG_PRIME_BIN="${NEOAG_PRIME_BIN:-$PRIME_HOME/PRIME}"
[[ -n "$MIXMHCPRED_HOME" ]] && export MIXMHCPRED_HOME MIXMHCPRED_BIN="${MIXMHCPRED_BIN:-$MIXMHCPRED_HOME/MixMHCpred}"
[[ -n "$BIGMHC_DIR" ]] && export BIGMHC_DIR
[[ -n "$DEEPIMMUNO_DIR" ]] && export DEEPIMMUNO_DIR
export NEOAG_PRIME_PYTHON="${NEOAG_PRIME_PYTHON:-$PY}"
export BIGMHC_PYTHON="${BIGMHC_PYTHON:-$PY}"

if [[ -n "$ASSET" ]]; then
  export NEOAG_NETCHOP_BIN="${NEOAG_NETCHOP_BIN:-$ASSET/data/predictors/netchop/netchop-3.1/Linux_x86_64/bin/netChop}"
  export NETCHOP_HOME="${NETCHOP_HOME:-$ASSET/data/predictors/netchop/netchop-3.1}"
fi

tool_path_entries=("${NETMHCPAN_HOME:-}" "$NETMHCSTABPAN_HOME" "$PRIME_HOME" "$MIXMHCPRED_HOME")
for tool_path in "${tool_path_entries[@]}"; do
  [[ -n "$tool_path" ]] && PATH="$tool_path:$PATH"
done
export PATH

echo "[INFO] run production_runner --execute"
PYTHONPATH="$PROJECT_ROOT/src" "$PY" -m neoag.production_runner \
  --manifest "$OUTDIR/manifest/production.results.toml" \
  --project-root "$PROJECT_ROOT" \
  --outdir "$OUTDIR" \
  --execute

echo "[OK] done: $OUTDIR/final/reports/"
