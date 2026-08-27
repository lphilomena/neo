#!/usr/bin/env bash
# Run Sequenza for one tumor-normal WGS sample by chromosome chunks.
# Resume order: binned seqz -> merged seqz -> chromosome seqz -> bam2seqz.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/conf/tools.env.sh"
: "${NEOAG_CONDA_BASE:?ERROR: set NEOAG_CONDA_BASE to your conda/mamba installation root}"
ENV="${SEQUENZA_ENV:-neoag-tools}"
R_ENV="${SEQUENZA_R_ENV:-neoag-r}"
CONDA_SH="${NEOAG_CONDA_BASE}/etc/profile.d/conda.sh"
source "${CONDA_SH}"
if [[ ! -d "${NEOAG_CONDA_BASE}/envs/${ENV}" ]]; then
  if [[ -d "${NEOAG_CONDA_BASE}/envs/neoag-tools" ]]; then
    echo "WARN Sequenza command env ${ENV} not found; falling back to neoag-tools" >&2
    ENV="neoag-tools"
  else
    echo "ERROR Sequenza command env not found: ${ENV}" >&2
    exit 1
  fi
fi
r_env_has_sequenza() {
  local candidate="$1"
  [[ -x "${NEOAG_CONDA_BASE}/envs/${candidate}/bin/Rscript" ]] || return 1
  env PATH="${NEOAG_CONDA_BASE}/envs/${candidate}/bin:${PATH}" \
    LD_LIBRARY_PATH="${NEOAG_CONDA_BASE}/envs/${candidate}/lib:${LD_LIBRARY_PATH:-}" \
    Rscript -e 'quit(status = ifelse(requireNamespace("sequenza", quietly = TRUE), 0, 1))' \
    >/dev/null 2>&1
}
if ! r_env_has_sequenza "${R_ENV}"; then
  for candidate in neoag-r "${ENV}" neoag-tools neoag-fusion; do
    if r_env_has_sequenza "${candidate}"; then
      echo "WARN Sequenza R env ${R_ENV} cannot load package sequenza; falling back to ${candidate}" >&2
      R_ENV="${candidate}"
      break
    fi
  done
fi
if ! r_env_has_sequenza "${R_ENV}"; then
  echo "ERROR no conda R env with package sequenza found. Set SEQUENZA_R_ENV to a working env." >&2
  exit 1
fi

SAMPLE_ID="${SAMPLE_ID:?ERROR: set SAMPLE_ID}"
TUMOR_BAM="${TUMOR_BAM:-}"
NORMAL_BAM="${NORMAL_BAM:-}"
REF="${REF_FASTA:-${SEQUENZA_FASTA:-${NEOAG_TOOLS_ROOT:-${ROOT}}/data/sequenza/reference/GRCh38.primary_assembly.chr.fa}}"
OUTDIR="${OUTDIR:-${ROOT}/results/sequenza/${SAMPLE_ID}}"
GC="${GC_WIGGLE:-${SEQUENZA_GC_WIG:-${ROOT}/work/sequenza/reference/gc${GC_WINDOW:-50}.wig.gz}}"
CHROMS="${CHROMS:-chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY}"
CHUNK_JOBS="${CHUNK_JOBS:-3}"
BIN_WINDOW="${BIN_WINDOW:-${SEQUENZA_BIN_WINDOW:-500}}"
GC_WINDOW="${GC_WINDOW:-50}"
QLIMIT="${QLIMIT:-20}"
MIN_DEPTH_N="${MIN_DEPTH_N:-20}"
HOM="${HOM:-0.9}"
HET="${HET:-0.25}"
FORCE="${FORCE:-0}"
REUSE_MERGED="${REUSE_MERGED:-1}"
REUSE_BINNED="${REUSE_BINNED:-1}"
MERGED_SEQZ="${MERGED_SEQZ:-}"
BINNED_SEQZ="${BINNED_SEQZ:-}"

SAMTOOLS="${NEOAG_CONDA_BASE}/envs/${ENV}/bin/samtools"
TABIX="${NEOAG_CONDA_BASE}/envs/${ENV}/bin/tabix"
LOG="${LOG:-${OUTDIR}/run.log}"
mkdir -p "${OUTDIR}/chrom" "${OUTDIR}/sequenza_fit" "$(dirname "${LOG}")"
exec > >(tee -a "${LOG}") 2>&1

run_env() {
  PATH="${NEOAG_CONDA_BASE}/envs/${ENV}/bin:${PATH}" LD_LIBRARY_PATH="${NEOAG_CONDA_BASE}/envs/${ENV}/lib:${LD_LIBRARY_PATH:-}" "$@"
}
run_r_env() {
  PATH="${NEOAG_CONDA_BASE}/envs/${R_ENV}/bin:${PATH}" LD_LIBRARY_PATH="${NEOAG_CONDA_BASE}/envs/${R_ENV}/lib:${LD_LIBRARY_PATH:-}" "$@"
}

validate_seqz() {
  local seqz="$1"
  [[ -s "${seqz}" ]] || return 1
  gzip -cd "${seqz}" 2>/dev/null | awk -F '\t' '
    NR == 1 {
      if (NF != 14 || $1 != "chromosome" || $2 != "position") exit 1
      next
    }
    NF != 14 { exit 1 }
    END { if (NR < 2) exit 1 }
  '
}

echo "[$(date -Is)] run_sequenza_sample_by_chrom sample=${SAMPLE_ID}"
echo "    tumor=${TUMOR_BAM}"
echo "    normal=${NORMAL_BAM}"
echo "    ref=${REF}"
echo "    gc=${GC}"
echo "    outdir=${OUTDIR}"
echo "    chroms=${CHROMS}"
echo "    chunk_jobs=${CHUNK_JOBS}"
echo "    bin_window=${BIN_WINDOW}"
echo "    merged_seqz=${MERGED_SEQZ}"
echo "    binned_seqz=${BINNED_SEQZ}"
echo "    sequenza_cmd_env=${ENV}"
echo "    sequenza_r_env=${R_ENV}"

require_bam_inputs() {
  [[ -n "${TUMOR_BAM}" ]] || { echo "ERROR set TUMOR_BAM unless BINNED_SEQZ or MERGED_SEQZ is provided" >&2; exit 1; }
  [[ -n "${NORMAL_BAM}" ]] || { echo "ERROR set NORMAL_BAM unless BINNED_SEQZ or MERGED_SEQZ is provided" >&2; exit 1; }
  for f in "${TUMOR_BAM}" "${NORMAL_BAM}" "${REF}"; do [[ -s "$f" ]] || { echo "ERROR missing $f" >&2; exit 1; }; done
  for bai in "${TUMOR_BAM}.bai" "${NORMAL_BAM}.bai"; do [[ -s "$bai" ]] || echo "WARN missing BAI by .bam.bai convention: $bai"; done
}

if [[ -z "${BINNED_SEQZ}" && -z "${MERGED_SEQZ}" ]]; then
  require_bam_inputs
fi

if [[ -z "${BINNED_SEQZ}" && -z "${MERGED_SEQZ}" && ! -s "${GC}" ]]; then
  mkdir -p "$(dirname "${GC}")"
  echo "[$(date -Is)] generating GC wiggle"
  run_env sequenza-utils gc_wiggle -f "${REF}" -w "${GC_WINDOW}" -o "${GC}"
fi

cat > "${OUTDIR}/run_parameters.tsv" <<EOF
sample_id	tumor_bam	normal_bam	reference	gc_wiggle	chroms	chunk_jobs	qlimit	min_depth_N	hom	het	gc_window	bin_window	merged_seqz	binned_seqz	force
${SAMPLE_ID}	${TUMOR_BAM}	${NORMAL_BAM}	${REF}	${GC}	${CHROMS}	${CHUNK_JOBS}	${QLIMIT}	${MIN_DEPTH_N}	${HOM}	${HET}	${GC_WINDOW}	${BIN_WINDOW}	${MERGED_SEQZ}	${BINNED_SEQZ}	${FORCE}
EOF

run_chrom() {
  local chrom="$1"
  local seqz="${OUTDIR}/chrom/${SAMPLE_ID}.${chrom}.seqz.gz"
  local seqz_tmp="${seqz}.tmp.$$"
  if validate_seqz "${seqz}"; then
    echo "[$(date -Is)] reuse ${chrom} ${seqz}"
    return 0
  fi
  if [[ -e "${seqz}" ]]; then
    mv "${seqz}" "${seqz}.partial.$(date +%Y%m%d_%H%M%S).$$"
  fi
  echo "[$(date -Is)] bam2seqz ${SAMPLE_ID} ${chrom}"
  rm -f "${seqz_tmp}" "${seqz_tmp}.tbi"
  env PATH="${NEOAG_CONDA_BASE}/envs/${ENV}/bin:${PATH}" \
    LD_LIBRARY_PATH="${NEOAG_CONDA_BASE}/envs/${ENV}/lib:${LD_LIBRARY_PATH:-}" \
    sequenza-utils bam2seqz \
    -n "${NORMAL_BAM}" \
    -t "${TUMOR_BAM}" \
    -gc "${GC}" \
    -F "${REF}" \
    -S "${SAMTOOLS}" \
    -T "${TABIX}" \
    -q "${QLIMIT}" \
    -N "${MIN_DEPTH_N}" \
    --hom "${HOM}" \
    --het "${HET}" \
    -C "${chrom}" \
    -o "${seqz_tmp}"
  validate_seqz "${seqz_tmp}"
  mv "${seqz_tmp}" "${seqz}"
  if [[ -s "${seqz_tmp}.tbi" ]]; then
    mv "${seqz_tmp}.tbi" "${seqz}.tbi"
  fi
}
export -f run_chrom validate_seqz
export SAMPLE_ID TUMOR_BAM NORMAL_BAM REF GC OUTDIR ENV NEOAG_CONDA_BASE SAMTOOLS TABIX QLIMIT MIN_DEPTH_N HOM HET PATH LD_LIBRARY_PATH

if [[ -n "${BINNED_SEQZ}" ]]; then
  validate_seqz "${BINNED_SEQZ}" || { echo "ERROR invalid or truncated BINNED_SEQZ ${BINNED_SEQZ}" >&2; exit 1; }
  BINNED="${BINNED_SEQZ}"
  echo "[$(date -Is)] reuse provided binned seqz ${BINNED}"
else
  MERGED="${MERGED_SEQZ:-${OUTDIR}/${SAMPLE_ID}.merged.seqz.gz}"
  if [[ -n "${MERGED_SEQZ}" ]]; then
    validate_seqz "${MERGED}" || { echo "ERROR invalid or truncated MERGED_SEQZ ${MERGED}" >&2; exit 1; }
    echo "[$(date -Is)] reuse provided merged seqz ${MERGED}"
  elif [[ "${FORCE}" != 1 && "${REUSE_MERGED}" == 1 ]] && validate_seqz "${MERGED}"; then
    echo "[$(date -Is)] reuse merged seqz ${MERGED}"
  else
    require_bam_inputs
    printf "%s\n" ${CHROMS} | xargs -I{} -P "${CHUNK_JOBS}" bash -c "run_chrom \"{}\""

    echo "[$(date -Is)] merge chrom seqz"
    {
      first=1
      for chrom in ${CHROMS}; do
        f="${OUTDIR}/chrom/${SAMPLE_ID}.${chrom}.seqz.gz"
        [[ -s "$f" ]] || { echo "ERROR missing chrom seqz $f" >&2; exit 1; }
        if ! validate_seqz "$f"; then
          echo "ERROR invalid or truncated chrom seqz $f" >&2
          exit 1
        fi
        if [[ "$first" == 1 ]]; then
          gzip -cd "$f"
          first=0
        else
          gzip -cd "$f" | tail -n +2
        fi
      done
    } | gzip -c > "${MERGED}.tmp"
    mv "${MERGED}.tmp" "${MERGED}"
  fi

  BINNED="${OUTDIR}/${SAMPLE_ID}.w${BIN_WINDOW}.seqz.gz"
  if [[ "${FORCE}" != 1 && "${REUSE_BINNED}" == 1 ]] && validate_seqz "${BINNED}"; then
    echo "[$(date -Is)] reuse binned seqz ${BINNED}"
  else
    echo "[$(date -Is)] seqz_binning window=${BIN_WINDOW}"
    TMP_BINNED="${BINNED}.tmp.gz"
    run_env sequenza-utils seqz_binning -s "${MERGED}" -w "${BIN_WINDOW}" -T "${TABIX}" -o "${TMP_BINNED}"
    mv "${TMP_BINNED}" "${BINNED}"
    if [[ -s "${TMP_BINNED}.tbi" ]]; then
      mv "${TMP_BINNED}.tbi" "${BINNED}.tbi"
    fi
  fi
  ln -sfn "$(basename "${BINNED}")" "${OUTDIR}/${SAMPLE_ID}.small.seqz.gz"
fi

SUMMARY="${OUTDIR}/sequenza_fit/${SAMPLE_ID}.sequenza_summary.tsv"
if [[ "${FORCE}" != 1 && -s "${SUMMARY}" ]]; then
  echo "[$(date -Is)] reuse Sequenza R fit ${SUMMARY}"
else
  echo "[$(date -Is)] R fit"
  CHECK_R="${OUTDIR}/sequenza_fit/check_sequenza_package.R"
  cat > "${CHECK_R}" <<'EOF_R'
suppressPackageStartupMessages(library(sequenza))
cat("sequenza_ok\n")
EOF_R
  if ! timeout "${SEQUENZA_R_LOAD_TIMEOUT:-60s}" env PATH="${NEOAG_CONDA_BASE}/envs/${R_ENV}/bin:${PATH}" LD_LIBRARY_PATH="${NEOAG_CONDA_BASE}/envs/${R_ENV}/lib:${LD_LIBRARY_PATH:-}" Rscript "${CHECK_R}" >/dev/null 2>&1; then
    echo "ERROR R package sequenza is not available or did not load within ${SEQUENZA_R_LOAD_TIMEOUT:-60s} in conda env ${R_ENV}. Set SEQUENZA_R_ENV to a working env or reuse existing purity consensus." >&2
    exit 1
  fi
  run_r_env Rscript "${ROOT}/scripts/run_sequenza_fit.R" "${BINNED}" "${OUTDIR}/sequenza_fit" "${SAMPLE_ID}"
fi
echo "[$(date -Is)] finished ${SAMPLE_ID}"
