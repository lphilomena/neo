#!/usr/bin/env bash
# Two-pass production driver for NeoAg v0.5.0 Splice Provenance Layer.
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/run_splice_provenance_v050.sh \
    --sample-id SAMPLE001 --outdir results/SAMPLE001/splice_v050 \
    --junctions regtools.tsv \
    --spladder-gff3 merge_graphs_exon_skip_C3.confirmed.gff3 \
    --irfinder IRFinder-IR-nondir.txt \
    --immunopepper-meta ref_sample_peptides_meta.tsv \
    [--immunopepper-kmers graph_kmer_JuncExpr.tsv] \
    [--normal-junctions normal_junctions.tsv] \
    [--normal-coverage normal_coverage.tsv] \
    [--hla HLA-A*02:01,HLA-B*07:02 --pvacbind-algorithms MHCflurry]

The driver first builds the formal junction→event→transcript→ORF→peptide layer,
then optionally runs pVACbind on the exact generated FASTA and rebuilds the final
layer with exact FASTA-index presentation provenance.
USAGE
}

SAMPLE_ID=""; OUTDIR=""; GENOME_BUILD="GRCh38"; DISEASE_PROFILE="default"
JUNCTIONS=""; JUNCTION_COORD="auto"; STAR_JUNCTIONS=""
JUNCTION_ASSAY=""; STAR_JUNCTION_ASSAY=""
SPLADDER_GFF3=(); SPLADDER_TXT=(); IRFINDER=(); IMMUNO_META=(); IMMUNO_KMERS=(); HIGH_ORDER=()
NORMAL_JUNCTIONS=(); NORMAL_COVERAGE=(); TOOL_VERSIONS=()
IR_COORD=""; NORMAL_COORD="auto"
HLA=""; HLA_FILE=""; PVAC_ALGORITHMS="MHCflurry"; PVAC_THREADS="4"; PVAC_REF_PROTEOME=""
STRICT=0; SKIP_PVAC=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --genome-build) GENOME_BUILD="$2"; shift 2 ;;
    --disease-profile) DISEASE_PROFILE="$2"; shift 2 ;;
    --junctions) JUNCTIONS="$2"; shift 2 ;;
    --junction-coordinate-system) JUNCTION_COORD="$2"; shift 2 ;;
    --junction-source-assay-id) JUNCTION_ASSAY="$2"; shift 2 ;;
    --star-junctions) STAR_JUNCTIONS="$2"; shift 2 ;;
    --star-junction-source-assay-id) STAR_JUNCTION_ASSAY="$2"; shift 2 ;;
    --spladder-gff3) SPLADDER_GFF3+=("$2"); shift 2 ;;
    --spladder-txt) SPLADDER_TXT+=("$2"); shift 2 ;;
    --irfinder) IRFINDER+=("$2"); shift 2 ;;
    --irfinder-coordinate-system) IR_COORD="$2"; shift 2 ;;
    --immunopepper-meta) IMMUNO_META+=("$2"); shift 2 ;;
    --immunopepper-kmers) IMMUNO_KMERS+=("$2"); shift 2 ;;
    --normal-junctions) NORMAL_JUNCTIONS+=("$2"); shift 2 ;;
    --normal-coordinate-system) NORMAL_COORD="$2"; shift 2 ;;
    --normal-coverage) NORMAL_COVERAGE+=("$2"); shift 2 ;;
    --high-order-evidence) HIGH_ORDER+=("$2"); shift 2 ;;
    --tool-version) TOOL_VERSIONS+=("$2"); shift 2 ;;
    --hla) HLA="$2"; shift 2 ;;
    --hla-file) HLA_FILE="$2"; shift 2 ;;
    --pvacbind-algorithms) PVAC_ALGORITHMS="$2"; shift 2 ;;
    --pvacbind-threads) PVAC_THREADS="$2"; shift 2 ;;
    --reference-proteome) PVAC_REF_PROTEOME="$2"; shift 2 ;;
    --skip-pvacbind) SKIP_PVAC=1; shift ;;
    --strict) STRICT=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done
[[ -n "$SAMPLE_ID" ]] || { echo "ERROR: --sample-id required" >&2; exit 2; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir required" >&2; exit 2; }
if [[ ${#IRFINDER[@]} -gt 0 && -z "$IR_COORD" ]]; then
  echo "ERROR: --irfinder-coordinate-system is required when --irfinder is used" >&2
  exit 2
fi
mkdir -p "$OUTDIR"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${NEOAG_PYTHON:-python}"
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
PRE="$OUTDIR/pre_pvacbind"
FINAL="$OUTDIR/splice_provenance"
PVAC_OUT="$OUTDIR/pvacbind"
rm -rf "$PRE" "$FINAL"

build_args=(
  -m neoag.splice.cli build --sample-id "$SAMPLE_ID" --outdir "$PRE"
  --genome-build "$GENOME_BUILD" --disease-profile "$DISEASE_PROFILE"
  --junction-coordinate-system "$JUNCTION_COORD"
  --normal-coordinate-system "$NORMAL_COORD"
)
[[ -n "$IR_COORD" ]] && build_args+=(--irfinder-coordinate-system "$IR_COORD")
[[ -n "$JUNCTIONS" ]] && build_args+=(--junctions "$JUNCTIONS")
[[ -n "$JUNCTION_ASSAY" ]] && build_args+=(--junction-source-assay-id "$JUNCTION_ASSAY")
[[ -n "$STAR_JUNCTIONS" ]] && build_args+=(--star-junctions "$STAR_JUNCTIONS")
[[ -n "$STAR_JUNCTION_ASSAY" ]] && build_args+=(--star-junction-source-assay-id "$STAR_JUNCTION_ASSAY")
for value in "${SPLADDER_GFF3[@]}"; do build_args+=(--spladder-gff3 "$value"); done
for value in "${SPLADDER_TXT[@]}"; do build_args+=(--spladder-txt "$value"); done
for value in "${IRFINDER[@]}"; do build_args+=(--irfinder "$value"); done
for value in "${IMMUNO_META[@]}"; do build_args+=(--immunopepper-meta "$value"); done
for value in "${IMMUNO_KMERS[@]}"; do build_args+=(--immunopepper-kmers "$value"); done
for value in "${NORMAL_JUNCTIONS[@]}"; do build_args+=(--normal-junctions "$value"); done
for value in "${NORMAL_COVERAGE[@]}"; do build_args+=(--normal-coverage "$value"); done
for value in "${HIGH_ORDER[@]}"; do build_args+=(--high-order-evidence "$value"); done
for value in "${TOOL_VERSIONS[@]}"; do build_args+=(--tool-version "$value"); done
[[ "$STRICT" == 1 ]] && build_args+=(--strict)
"$PYTHON_BIN" "${build_args[@]}" > "$OUTDIR/pre_pvacbind.outputs.json"

PVAC_LIST=""
if [[ "$SKIP_PVAC" == 0 && ( -n "$HLA" || -n "$HLA_FILE" ) && -s "$PRE/splice_pvacbind_input.fasta" ]]; then
  pvac_args=(bash "$ROOT/scripts/run_pvacbind_sample.sh" --fasta "$PRE/splice_pvacbind_input.fasta" --sample-id "$SAMPLE_ID" --outdir "$PVAC_OUT" --algorithms "$PVAC_ALGORITHMS" --threads "$PVAC_THREADS")
  [[ -n "$HLA" ]] && pvac_args+=(--hla "$HLA")
  [[ -n "$HLA_FILE" ]] && pvac_args+=(--hla-file "$HLA_FILE")
  [[ -n "$PVAC_REF_PROTEOME" ]] && pvac_args+=(--reference-proteome "$PVAC_REF_PROTEOME")
  "${pvac_args[@]}" > "$OUTDIR/pvacbind.driver.out"
  PVAC_LIST="$PVAC_OUT/pvacbind_all_epitopes.list"
fi

if [[ -n "$PVAC_LIST" && -s "$PVAC_LIST" ]]; then
  final_args=("${build_args[@]}")
  # Replace the preliminary output directory with the final output directory.
  for ((i=0; i<${#final_args[@]}; i++)); do
    if [[ "${final_args[$i]}" == "$PRE" ]]; then final_args[$i]="$FINAL"; fi
  done
  final_args+=(--pvacbind-fasta-map "$PRE/splice_pvacbind_fasta_map.tsv")
  while IFS= read -r path; do [[ -n "$path" ]] && final_args+=(--pvacbind "$path"); done < "$PVAC_LIST"
  "$PYTHON_BIN" "${final_args[@]}" > "$OUTDIR/final.outputs.json"
else
  mkdir -p "$FINAL"
  cp -a "$PRE/." "$FINAL/"
  cp "$OUTDIR/pre_pvacbind.outputs.json" "$OUTDIR/final.outputs.json"
fi

validate_args=(-m neoag.splice.cli validate --layer-dir "$FINAL" --report "$OUTDIR/validation.json")
[[ "$STRICT" == 1 ]] && validate_args+=(--strict)
"$PYTHON_BIN" "${validate_args[@]}"
printf '%s\n' "$FINAL"
