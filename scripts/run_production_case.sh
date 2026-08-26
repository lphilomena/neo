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
    [--event-top-n <N; default 20>] \
    [--candidate-top-n <N; default 100>] \
    [--asset-root <liup_neodata4git>] \
    [--reference-fasta <GRCh38.fasta>] \
    [--gencode-gtf <matching.gtf; also used for exact splice strand/origin reconstruction>] \
    [--sequenza <result_file_or_dir>] \
    [--purple <result_file_or_dir>] \
    [--expression <gene_tpm.tsv>] \
    [--transcript-expression <transcript_tpm.tsv|quant.sf>] \
    [--rna-fastq1 <R1.fq.gz[,lane2_R1.fq.gz]>] \
    [--rna-fastq2 <R2.fq.gz[,lane2_R2.fq.gz]>] \
    [--rna-bam <sorted_rna.bam> | --rna-vaf <rna_alt_vaf.tsv>] \
    [--star-index <GRCh38_STAR_index>] \
    [--easyfuse-star-index <EasyFuse_STAR_index>] \
    [--star-index-build-dir <new_STAR_index_dir>] \
    [--star-sjdb-overhang <N>] \
    [--star-executable <STAR>] \
    [--samtools-executable <samtools>] \
    [--rna-threads <N>] \
    [--fusion-caller-root <completed_caller_results_dir>] \
    [--star-chimeric <STAR/Chimeric.out.junction; repeatable>] \
    [--normal-readthrough <normal_readthrough.tsv>] \
    [--prime-evidence <prime_evidence.tsv>] \
    [--bigmhc-evidence <bigmhc_im_evidence.tsv>] \
    [--deepimmuno-evidence <deepimmuno_evidence.tsv>] \
    [--pred-deps <predictor_deps_dir>] \
    [--netmhcpan-home <netMHCpan_home>] \
    [--netmhcstabpan-home <netMHCstabpan_home>] \
    [--python <python_bin>]

Environment defaults may be supplied through conf/tools.env.local.sh or:
  NEOAG_PYTHON, NEOAG_ASSET_ROOT, NEOAG_PREDICTOR_DEPS,
  NETMHCPAN_HOME, NETMHCSTABPAN_HOME, EASYFUSE_STAR_INDEX,
  and NEOAG_EASYFUSE_REF.
EOF
}

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${NEOAG_PYTHON:-$(command -v python3 || command -v python || true)}"
PROFILE="profiles/sarcoma_rna_supported_v2_provisional.toml"
EVIDENCE_CONSENSUS_RULES="configs/ranking/sarcoma_evidence_consensus_v3_source_chain.toml"
EVENT_TOP_N=20
CANDIDATE_TOP_N=100
ASSET="${NEOAG_ASSET_ROOT:-${NEOAG_TOOLS_ROOT:-}}"
PRED_DEPS="${NEOAG_PREDICTOR_DEPS:-${NEOAG_TOOL_QUARANTINE:-}}"
NETMHCPAN_HOME_DEFAULT="${NETMHCPAN_HOME:-}"
NETMHCSTABPAN_HOME="${NETMHCSTABPAN_HOME:-}"
LICENSED_TOOLS_ROOT="${NEOAG_LICENSED_TOOLS_ROOT:-}"
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
EXPRESSION=""
TRANSCRIPT_EXPRESSION=""
RNA_FASTQ1=""
RNA_FASTQ2=""
RNA_BAM=""
RNA_VAF=""
STAR_INDEX=""
EASYFUSE_STAR_INDEX="${EASYFUSE_STAR_INDEX:-}"
STAR_INDEX_BUILD_DIR=""
STAR_SJDB_OVERHANG=149
STAR_EXECUTABLE=""
SAMTOOLS_EXECUTABLE="samtools"
RNA_THREADS=16
FUSION_CALLER_ROOTS=()
STAR_CHIMERIC_FILES=()
NORMAL_READTHROUGH=""
PRIME_EVIDENCE=""
BIGMHC_EVIDENCE=""
DEEPIMMUNO_EVIDENCE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --case-root) CASE_ROOT="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --somatic-vcf) SOMATIC_VCF="$2"; shift 2 ;;
    --project-root) PROJECT_ROOT="$2"; shift 2 ;;
    --profile) PROFILE="$2"; shift 2 ;;
    --evidence-consensus-rules) EVIDENCE_CONSENSUS_RULES="$2"; shift 2 ;;
    --event-top-n) EVENT_TOP_N="$2"; shift 2 ;;
    --candidate-top-n) CANDIDATE_TOP_N="$2"; shift 2 ;;
    --asset-root) ASSET="$2"; CLI_ASSET="$2"; shift 2 ;;
    --reference-fasta) REFERENCE_FASTA="$2"; shift 2 ;;
    --gencode-gtf) GENCODE_GTF="$2"; shift 2 ;;
    --sequenza) SEQUENZA="$2"; shift 2 ;;
    --purple) PURPLE="$2"; shift 2 ;;
    --expression) EXPRESSION="$2"; shift 2 ;;
    --transcript-expression) TRANSCRIPT_EXPRESSION="$2"; shift 2 ;;
    --rna-fastq1) RNA_FASTQ1="$2"; shift 2 ;;
    --rna-fastq2) RNA_FASTQ2="$2"; shift 2 ;;
    --rna-bam) RNA_BAM="$2"; shift 2 ;;
    --rna-vaf) RNA_VAF="$2"; shift 2 ;;
    --star-index) STAR_INDEX="$2"; shift 2 ;;
    --easyfuse-star-index) EASYFUSE_STAR_INDEX="$2"; shift 2 ;;
    --star-index-build-dir) STAR_INDEX_BUILD_DIR="$2"; shift 2 ;;
    --star-sjdb-overhang) STAR_SJDB_OVERHANG="$2"; shift 2 ;;
    --star-executable) STAR_EXECUTABLE="$2"; shift 2 ;;
    --samtools-executable) SAMTOOLS_EXECUTABLE="$2"; shift 2 ;;
    --rna-threads) RNA_THREADS="$2"; shift 2 ;;
    --fusion-caller-root) FUSION_CALLER_ROOTS+=("$2"); shift 2 ;;
    --star-chimeric) STAR_CHIMERIC_FILES+=("$2"); shift 2 ;;
    --normal-readthrough) NORMAL_READTHROUGH="$2"; shift 2 ;;
    --prime-evidence) PRIME_EVIDENCE="$2"; shift 2 ;;
    --bigmhc-evidence) BIGMHC_EVIDENCE="$2"; shift 2 ;;
    --deepimmuno-evidence) DEEPIMMUNO_EVIDENCE="$2"; shift 2 ;;
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
[[ "$STAR_SJDB_OVERHANG" =~ ^[1-9][0-9]*$ ]] || { echo "--star-sjdb-overhang must be a positive integer" >&2; exit 2; }
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
  [[ -n "$GENCODE_GTF" ]] || { echo "RNA FASTQ mode requires --gencode-gtf" >&2; exit 2; }
fi

verify_event_track_precedence() {
  # Fail before heavy execution when a deployed code version can still turn
  # a DNA SNV/InDel into an RNA splice event solely because VEP reports a
  # splice-related peptide consequence.
  local check='from neoag.evidence_states import event_track
from neoag.source_chain import source_chain_track
from neoag.reports_dual import _patient_track
from neoag.open_neo.review import _event_kind
snv = {"event_type": "SNV", "mutation_source": "SNV_INDEL", "peptide_consequence": "splice_junction"}
indel = {"event_type": "InDel", "mutation_source": "INDEL", "peptide_consequence": "splice_junction"}
splice = {"event_type": "Splice", "mutation_source": "SPLICE", "peptide_consequence": "splice_junction"}
assert event_track(snv) == "MISSENSE"
assert source_chain_track(snv) == "SNV"
assert _patient_track(snv) == "SNV"
assert _event_kind("SNV", "splice_junction") == "MISSENSE"
assert event_track(indel) != "SPLICE"
assert source_chain_track(indel) == "INDEL"
assert _patient_track(indel) == "InDel"
assert _event_kind("InDel", "splice_junction") != "SPLICE"
assert event_track(splice) == "SPLICE"
assert source_chain_track(splice) == "SPLICE"
assert _patient_track(splice) == "Splice"
assert _event_kind("Splice", "splice_junction") == "SPLICE"'
  if ! PYTHONPATH="$PROJECT_ROOT/src" "$PY" -c "$check"; then
    echo "event-track precedence preflight failed; update NeoAg before production execution" >&2
    return 1
  fi
}

verify_mtwt_interpretation_rules() {
  # Confirm that this checkout applies the same cautious, structure-aware
  # MT/WT policy used by Skill2 and Skill3 before starting expensive work.
  local check='from neoag.mutant_specificity import evaluate_mutant_specificity
profile = {"mutant_specificity": {"near_equal_el_rank_difference": 0.01, "positive_agretopicity_ratio": 2.0, "positive_el_rank_difference": 0.10}}
row = evaluate_mutant_specificity(
    {"peptide": "ABCXEFGHI", "wildtype_peptide": "ABCDEFGHI", "mhc_class": "I"},
    {"netmhcpan_mt_rank_el": "0.4", "netmhcpan_wt_rank_el": "0.7", "netmhcpan_wt_rank_ba": "0.9", "netmhcpan_wt_ic50": "42"},
    profile,
)
assert row["mutation_position_role"] == "PUTATIVE_TCR_FACING"
assert row["wt_self_reactivity_risk_status"] == "WT_STRONG_BINDING_REVIEW"
assert row["mutant_specificity_gate_status"] == "CAUTION"
assert row["mutant_specificity_priority_cap"]
mhc2 = evaluate_mutant_specificity(
    {"peptide": "ABCDEFGHIJKLMNO", "wildtype_peptide": "ABCXEFGHIJKLMNO", "hla_allele": "HLA-DRB1*04:01"},
    {"netmhcpan_mt_rank_el": "0.2", "netmhcpan_wt_rank_el": "3.0"},
    profile,
)
assert mhc2["mutation_position_role"] == "STRUCTURAL_ROLE_UNCERTAIN"'
  if ! PYTHONPATH="$PROJECT_ROOT/src" "$PY" -c "$check"; then
    echo "MT/WT structure/risk preflight failed; update NeoAg before production execution" >&2
    return 1
  fi
}

verify_mtwt_output_fields() {
  local ranked="$1"
  [[ -s "$ranked" ]] || return 0
  local check='import csv, sys
path = sys.argv[1]
missing = []
with open(path, encoding="utf-8", newline="") as handle:
    for line_no, row in enumerate(csv.DictReader(handle, delimiter="\t"), start=2):
        event_type = str(row.get("event_type", "")).strip().lower()
        if event_type not in {"snv", "missense", "substitution", "indel", "insertion", "deletion", "frameshift"}:
            continue
        if not str(row.get("wildtype_peptide", "")).strip():
            continue
        required = ("mutation_position_role", "mutation_position_interpretation", "wt_self_reactivity_risk_status", "wt_self_reactivity_risk_reason", "mt_wt_interpretation_caution")
        absent = [field for field in required if not str(row.get(field, "")).strip()]
        wt_prediction = any(str(row.get(field, "")).strip() for field in ("netmhcpan_wt_rank_el", "netmhcpan_wt_rank_ba", "netmhcpan_wt_ic50"))
        if wt_prediction and str(row.get("wt_self_reactivity_risk_status", "")).upper() in {"", "UNASSESSED"}:
            absent.append("assessed_wt_self_reactivity_risk_status")
        if absent:
            missing.append("line {}: {}".format(line_no, ",".join(dict.fromkeys(absent))))
if missing:
    raise SystemExit("MT/WT structure/risk fields missing from final evidence ranking: " + "; ".join(missing[:20]))'
  PYTHONPATH="$PROJECT_ROOT/src" "$PY" -c "$check" "$ranked"
}

verify_splice_prefilter_outputs() {
  local funnel="$1" decisions="$2"
  [[ -s "$funnel" ]] || { echo "splice prefilter funnel missing: $funnel" >&2; return 1; }
  [[ -s "$decisions" ]] || { echo "splice prefilter decisions missing: $decisions" >&2; return 1; }
  local required=(
    ALIGNMENT_COORDINATE_QC UNIQUE_JUNCTION_READS TOTAL_JUNCTION_COVERAGE PSI
    MATCHED_NORMAL_JUNCTION NORMAL_COHORT_JUNCTION ANNOTATED_NORMAL_ISOFORM
    CREDIBLE_ORF NMD JUNCTION_SPANNING_PEPTIDE NORMAL_PROTEOME_EXCLUSION
    SELECTED_FOR_PRESENTATION
  )
  local stage
  for stage in "${required[@]}"; do
    grep -q "^${stage}[[:space:]]" "$funnel" || {
      echo "splice prefilter funnel missing stage: $stage" >&2
      return 1
    }
  done
}

cd "$PROJECT_ROOT"
verify_event_track_precedence
verify_mtwt_interpretation_rules
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

if [[ -z "$EASYFUSE_STAR_INDEX" && -n "${NEOAG_EASYFUSE_REF:-}" ]]; then
  for candidate in \
    "$NEOAG_EASYFUSE_REF/starfusion_index/ref_genome.fa.star.idx" \
    "$NEOAG_EASYFUSE_REF/star_index"; do
    if [[ -d "$candidate" ]]; then EASYFUSE_STAR_INDEX="$candidate"; break; fi
  done
fi
if [[ -z "$EASYFUSE_STAR_INDEX" && -n "$ASSET" ]]; then
  for candidate in \
    "$ASSET/data/easyfuse/easyfuse_ref_v4/starfusion_index/ref_genome.fa.star.idx" \
    "$ASSET/data/ref/ctat/current/ctat_genome_lib_build_dir/ref_genome.fa.star.idx"; do
    if [[ -d "$candidate" ]]; then EASYFUSE_STAR_INDEX="$candidate"; break; fi
  done
fi
if [[ -z "$LICENSED_TOOLS_ROOT" && -n "$ASSET" ]]; then
  LICENSED_TOOLS_ROOT="$(cd "$ASSET/.." && pwd)/licensed_tools"
fi

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

add_file_if() {
  local flag="$1" path="$2"
  if [[ -f "$path" && -s "$path" ]]; then
    GEN_ARGS+=("$flag" "$path")
  fi
}

ensure_normal_junction_index() {
  local normal_junctions="$1"
  local index_path="${normal_junctions}.sqlite"
  if [[ ! -s "$normal_junctions" || -s "$index_path" ]]; then
    return 0
  fi
  echo "[INFO] build normal junction sqlite index: $index_path"
  PYTHONPATH="$PROJECT_ROOT/src" "$PY" scripts/build_normal_junction_index.py \
    --input "$normal_junctions" \
    --output "$index_path"
}

latest_matching_file() {
  local root="$1"
  shift
  [[ -d "$root" ]] || return 1
  find "$root" -maxdepth 7 -type f \( "$@" \) -size +0c -printf '%T@ %p\n' 2>/dev/null \
    | sort -nr \
    | awk 'NR == 1 {sub(/^[^ ]+ /, ""); print; exit}'
}

add_first_existing() {
  local flag="$1"
  shift
  local path=""
  for path in "$@"; do
    if [[ -s "$path" ]]; then
      GEN_ARGS+=("$flag" "$path")
      return 0
    fi
  done
  return 1
}

GEN_ARGS=(
  --project-root "$PROJECT_ROOT"
  --sample-id "$SAMPLE_ID"
  --profile "$PROFILE_PATH"
  --evidence-consensus-rules "$CONSENSUS_RULES_PATH"
  --event-top-n "$EVENT_TOP_N"
  --candidate-top-n "$CANDIDATE_TOP_N"
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
if [[ -n "$EXPRESSION" ]]; then
  add_if --expression "$EXPRESSION"
else
  discovered_expression="$(latest_matching_file "$CASE_ROOT" -name gene_tpm.tsv -o -name '*.genes.results')"
  [[ -z "$discovered_expression" && -d "$OUTDIR" ]] && discovered_expression="$(latest_matching_file "$OUTDIR" -name gene_tpm.tsv -o -name '*.genes.results')"
  add_first_existing --expression \
    "$CASE_ROOT/short-rna/evidence/gene_expression.tsv" \
    "$CASE_ROOT/short-rna/evidence/gene_tpm.tsv" \
    "$discovered_expression" || true
fi
if [[ -n "$TRANSCRIPT_EXPRESSION" ]]; then
  add_if --transcript-expression "$TRANSCRIPT_EXPRESSION"
else
  discovered_transcript_expression="$(latest_matching_file "$CASE_ROOT" -name transcript_tpm.tsv -o -name quant.sf -o -name '*.isoforms.results')"
  [[ -z "$discovered_transcript_expression" && -d "$OUTDIR" ]] && discovered_transcript_expression="$(latest_matching_file "$OUTDIR" -name transcript_tpm.tsv -o -name quant.sf -o -name '*.isoforms.results')"
  add_first_existing --transcript-expression \
    "$CASE_ROOT/short-rna/evidence/transcript_tpm.tsv" \
    "$CASE_ROOT/short-rna/evidence/transcript_quant.sf" \
    "$discovered_transcript_expression" || true
fi
if [[ -z "$RNA_FASTQ1" && -z "$RNA_BAM" && -z "$RNA_VAF" ]]; then
  discovered_rna_vaf="$(latest_matching_file "$CASE_ROOT" -name rna_alt_vaf.tsv -o -name '*rna*vaf*.tsv' -o -name '*rna*alt*.tsv')"
  [[ -z "$discovered_rna_vaf" && -d "$OUTDIR" ]] && discovered_rna_vaf="$(latest_matching_file "$OUTDIR" -name rna_alt_vaf.tsv -o -name '*rna*vaf*.tsv' -o -name '*rna*alt*.tsv')"
  if [[ -n "$discovered_rna_vaf" ]]; then
    RNA_VAF="$discovered_rna_vaf"
  else
    discovered_rna_bam="$(latest_matching_file "$CASE_ROOT" -name Aligned.sortedByCoord.out.bam -o -name '*.Aligned.sortedByCoord.out.bam' -o -name '*.rna.bam')"
    [[ -z "$discovered_rna_bam" && -d "$OUTDIR" ]] && discovered_rna_bam="$(latest_matching_file "$OUTDIR" -name Aligned.sortedByCoord.out.bam -o -name '*.Aligned.sortedByCoord.out.bam' -o -name '*.rna.bam')"
    if [[ -z "$discovered_rna_bam" ]]; then
      for candidate in \
        "$CASE_ROOT/pipeline/production/rna/star/Aligned.sortedByCoord.out.bam" \
        "$CASE_ROOT/pipeline/production/fusion_fixed_20260811_171125/star/Aligned.sortedByCoord.out.bam" \
        "$CASE_ROOT/pipeline/production/rna_parallel_core_20260811_095246/star/Aligned.sortedByCoord.out.bam" \
        "$OUTDIR/pipeline/production/rna/star/Aligned.sortedByCoord.out.bam" \
        "$OUTDIR/rna/star/Aligned.sortedByCoord.out.bam"; do
        if [[ -s "$candidate" ]]; then
          discovered_rna_bam="$candidate"
          break
        fi
      done
    fi
    [[ -n "$discovered_rna_bam" ]] && RNA_BAM="$discovered_rna_bam"
  fi
fi
add_if --easyfuse "$CASE_ROOT/short-rna/evidence/easyfuse.fusions.pass.csv"
add_if --star-fusion "$CASE_ROOT/short-rna/evidence/star-fusion.fusion_predictions.tsv"
add_if --arriba "$CASE_ROOT/short-rna/evidence/arriba.fusions.tsv"
add_if --fusioncatcher "$CASE_ROOT/short-rna/evidence/fusioncatcher.final-list.txt"
add_if --jaffal "$CASE_ROOT/long-rna/evidence/jaffa_results.csv"
add_if --jaffal "$CASE_ROOT/long-rna/jaffal/output/jaffa_results.csv"
for fusion_root in "${FUSION_CALLER_ROOTS[@]}"; do
  [[ -e "$fusion_root" ]] && GEN_ARGS+=(--fusion-caller-root "$fusion_root")
done
for fusion_root in "$CASE_ROOT/short-rna/fusion" "$CASE_ROOT/long-rna/fusion"; do
  [[ -d "$fusion_root" ]] && GEN_ARGS+=(--fusion-caller-root "$fusion_root")
done
for chimeric in "${STAR_CHIMERIC_FILES[@]}" \
  "$CASE_ROOT/short-rna/star/Chimeric.out.junction" \
  "$CASE_ROOT/short-rna/fusion/Chimeric.out.junction" \
  "$OUTDIR/pipeline/production/rna/star/Chimeric.out.junction" \
  "$OUTDIR/rna/star/Chimeric.out.junction"; do
  [[ -s "$chimeric" ]] && GEN_ARGS+=(--star-chimeric "$chimeric")
done
[[ -n "$NORMAL_READTHROUGH" ]] && GEN_ARGS+=(--normal-readthrough "$NORMAL_READTHROUGH")
[[ -n "$PRIME_EVIDENCE" ]] && GEN_ARGS+=(--prime-evidence "$PRIME_EVIDENCE")
[[ -n "$BIGMHC_EVIDENCE" ]] && GEN_ARGS+=(--bigmhc-evidence "$BIGMHC_EVIDENCE")
[[ -n "$DEEPIMMUNO_EVIDENCE" ]] && GEN_ARGS+=(--deepimmuno-evidence "$DEEPIMMUNO_EVIDENCE")
add_if --junctions "$CASE_ROOT/short-rna/evidence/regtools_junctions.tsv"
add_if --snaf "$CASE_ROOT/short-rna/snaf/snaf_candidates.tsv"
add_if --splicemutr "$CASE_ROOT/short-rna/splicemutr"

if [[ -n "$ASSET" ]]; then
  NORMAL_JUNCTIONS="$ASSET/data/normal/junctions/normal_junctions.GRCh38.tsv.gz"
  ensure_normal_junction_index "$NORMAL_JUNCTIONS"
  add_file_if --normal-junctions "$NORMAL_JUNCTIONS"
  add_file_if --normal-expression "$ASSET/data/normal/expression/normal_expression.gtex_v11_hpa_hspc.tsv"
  add_file_if --normal-hla-ligands "$ASSET/data/normal/ligandome/normal_ms_ligands.tsv"
  discovered_reference_proteome=""
  if [[ -d "$ASSET/data/normal/proteome" ]]; then
    discovered_reference_proteome="$(latest_matching_file "$ASSET/data/normal/proteome" -name '*.fa' -o -name '*.fasta' -o -name '*.faa')"
  fi
  add_first_existing --reference-proteome \
    "${NEOAG_NORMAL_PROTEOME_FASTA:-}" \
    "${NEOAG_NORMAL_PROTEOME:-}" \
    "$ASSET/data/normal/proteome/gencode.v49.pc_translations.clean.fa" \
    "$ASSET/data/normal/proteome/Homo_sapiens.GRCh38.pep.all.fa" \
    "$ASSET/data/normal/proteome/Homo_sapiens.GRCh38.pep.all.fa.gz" \
    "$discovered_reference_proteome" || true
  add_first_existing --netchop-executable \
    "${NEOAG_NETCHOP_BIN:-}" \
    "$LICENSED_TOOLS_ROOT/netchop/netchop-3.1/Linux_x86_64/bin/netChop" \
    "$ASSET/data/predictors/netchop/netchop-3.1/Linux_x86_64/bin/netChop" || true
  for netchop_home_candidate in \
    "${NETCHOP_HOME:-}" \
    "$LICENSED_TOOLS_ROOT/netchop/netchop-3.1" \
    "$ASSET/data/predictors/netchop/netchop-3.1"; do
    if [[ -d "$netchop_home_candidate" ]]; then
      GEN_ARGS+=(--netchop-home "$netchop_home_candidate")
      break
    fi
  done
fi

# Explicit reference and RNA inputs are passed through unchanged. They are not
# inferred from patient-specific paths because build/version mismatches are unsafe.
[[ -n "$REFERENCE_FASTA" ]] && GEN_ARGS+=(--reference-fasta "$REFERENCE_FASTA")
[[ -n "$GENCODE_GTF" ]] && GEN_ARGS+=(--gencode-gtf "$GENCODE_GTF")
[[ -n "$RNA_FASTQ1" ]] && GEN_ARGS+=(--rna-fastq1 "$RNA_FASTQ1" --rna-fastq2 "$RNA_FASTQ2")
[[ -n "$RNA_BAM" ]] && GEN_ARGS+=(--rna-bam "$RNA_BAM")
[[ -n "$RNA_VAF" ]] && GEN_ARGS+=(--rna-vaf "$RNA_VAF")
[[ -n "$STAR_INDEX" ]] && GEN_ARGS+=(--star-index "$STAR_INDEX")
[[ -n "$EASYFUSE_STAR_INDEX" ]] && GEN_ARGS+=(--easyfuse-star-index "$EASYFUSE_STAR_INDEX")
[[ -n "$STAR_INDEX_BUILD_DIR" ]] && GEN_ARGS+=(--star-index-build-dir "$STAR_INDEX_BUILD_DIR")
GEN_ARGS+=(--star-sjdb-overhang "$STAR_SJDB_OVERHANG")
[[ -n "$STAR_EXECUTABLE" ]] && GEN_ARGS+=(--star-executable "$STAR_EXECUTABLE")
GEN_ARGS+=(--samtools-executable "$SAMTOOLS_EXECUTABLE" --rna-threads "$RNA_THREADS")

echo "[INFO] generate manifest: $OUTDIR/manifest/production.results.toml"
"$PY" scripts/generate_production_from_results_manifest.py "${GEN_ARGS[@]}"

[[ -n "$ASSET" ]] && export NEOAG_TOOLS_ROOT="$ASSET"
[[ -n "$PRED_DEPS" ]] && export NEOAG_TOOL_QUARANTINE="$PRED_DEPS"
if [[ -n "$ASSET" ]]; then
  export NEOAG_VEP_CACHE="${NEOAG_VEP_CACHE:-$ASSET/data/vep}"
fi
export NEOAG_VEP_CACHE_VERSION="${NEOAG_VEP_CACHE_VERSION:-105}"
PY_PREFIX="$(cd "$(dirname "$PY")/.." && pwd)"
for lib_dir in "$PY_PREFIX/lib" "${NEOAG_ENV_TOOL_ROOT:-}/lib"; do
  if [[ -d "$lib_dir" ]]; then
    export LD_LIBRARY_PATH="$lib_dir:${LD_LIBRARY_PATH:-}"
  fi
done
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
  for netchop_bin_candidate in \
    "${NEOAG_NETCHOP_BIN:-}" \
    "$LICENSED_TOOLS_ROOT/netchop/netchop-3.1/Linux_x86_64/bin/netChop" \
    "$ASSET/data/predictors/netchop/netchop-3.1/Linux_x86_64/bin/netChop"; do
    if [[ -x "$netchop_bin_candidate" ]]; then
      export NEOAG_NETCHOP_BIN="$netchop_bin_candidate"
      break
    fi
  done
  for netchop_home_candidate in \
    "${NETCHOP_HOME:-}" \
    "$LICENSED_TOOLS_ROOT/netchop/netchop-3.1" \
    "$ASSET/data/predictors/netchop/netchop-3.1"; do
    if [[ -d "$netchop_home_candidate" ]]; then
      export NETCHOP_HOME="$netchop_home_candidate"
      break
    fi
  done
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

verify_mtwt_output_fields "$OUTDIR/final/scoring/ranked_peptides.evidence_consensus.tsv"
if [[ -d "$OUTDIR/final" ]]; then
  verify_splice_prefilter_outputs \
    "$OUTDIR/final/parsed/splice_prefilter_funnel.tsv" \
    "$OUTDIR/final/parsed/splice_prefilter_decisions.tsv"
fi

echo "[OK] done: $OUTDIR/final/reports/"
