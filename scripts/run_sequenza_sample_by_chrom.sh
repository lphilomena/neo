#!/usr/bin/env bash
# Run Sequenza for one tumor-normal WGS sample by chromosome chunks.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export NEOAG_CONDA_BASE="${NEOAG_CONDA_BASE:-/home/na/miniforge3}"
ENV="${SEQUENZA_ENV:-neoag-sequenza}"
CONDA_SH="${NEOAG_CONDA_BASE}/etc/profile.d/conda.sh"
source "${CONDA_SH}"

SAMPLE_ID="${SAMPLE_ID:?ERROR: set SAMPLE_ID}"
TUMOR_BAM="${TUMOR_BAM:?ERROR: set TUMOR_BAM}"
NORMAL_BAM="${NORMAL_BAM:?ERROR: set NORMAL_BAM}"
REF="${REF_FASTA:-/mnt/zjl-bgi-zzb/peixunban/gl/data/reference/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa}"
OUTDIR="${OUTDIR:-${ROOT}/results/sequenza/${SAMPLE_ID}}"
GC="${GC_WIGGLE:-${ROOT}/results/chenxiaoliang_sequenza/reference/Homo_sapiens.GRCh38.dna.primary_assembly.chr.gc50.wig.gz}"
CHROMS="${CHROMS:-chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY}"
CHUNK_JOBS="${CHUNK_JOBS:-3}"
BIN_WINDOW="${BIN_WINDOW:-50}"
GC_WINDOW="${GC_WINDOW:-50}"
QLIMIT="${QLIMIT:-20}"
MIN_DEPTH_N="${MIN_DEPTH_N:-20}"
HOM="${HOM:-0.9}"
HET="${HET:-0.25}"

SAMTOOLS="${NEOAG_CONDA_BASE}/envs/${ENV}/bin/samtools"
TABIX="${NEOAG_CONDA_BASE}/envs/${ENV}/bin/tabix"
LOG="${LOG:-${OUTDIR}/run.log}"
mkdir -p "${OUTDIR}/chrom" "${OUTDIR}/sequenza_fit" "$(dirname "${LOG}")"
exec > >(tee -a "${LOG}") 2>&1

run_env() { mamba run -n "${ENV}" "$@"; }

echo "[$(date -Is)] run_sequenza_sample_by_chrom sample=${SAMPLE_ID}"
echo "    tumor=${TUMOR_BAM}"
echo "    normal=${NORMAL_BAM}"
echo "    ref=${REF}"
echo "    gc=${GC}"
echo "    outdir=${OUTDIR}"
echo "    chroms=${CHROMS}"
echo "    chunk_jobs=${CHUNK_JOBS}"

for f in "${TUMOR_BAM}" "${NORMAL_BAM}" "${REF}"; do [[ -s "$f" ]] || { echo "ERROR missing $f" >&2; exit 1; }; done
for bai in "${TUMOR_BAM}.bai" "${NORMAL_BAM}.bai"; do [[ -s "$bai" ]] || echo "WARN missing BAI by .bam.bai convention: $bai"; done
if [[ ! -s "${GC}" ]]; then
  mkdir -p "$(dirname "${GC}")"
  echo "[$(date -Is)] generating GC wiggle"
  run_env sequenza-utils gc_wiggle -f "${REF}" -w "${GC_WINDOW}" -o "${GC}"
fi

cat > "${OUTDIR}/run_parameters.tsv" <<EOF
sample_id	tumor_bam	normal_bam	reference	gc_wiggle	chroms	chunk_jobs	qlimit	min_depth_N	hom	het	gc_window	bin_window
${SAMPLE_ID}	${TUMOR_BAM}	${NORMAL_BAM}	${REF}	${GC}	${CHROMS}	${CHUNK_JOBS}	${QLIMIT}	${MIN_DEPTH_N}	${HOM}	${HET}	${GC_WINDOW}	${BIN_WINDOW}
EOF

run_chrom() {
  local chrom="$1"
  local seqz="${OUTDIR}/chrom/${SAMPLE_ID}.${chrom}.seqz.gz"
  if [[ -s "${seqz}" ]]; then
    echo "[$(date -Is)] reuse ${chrom} ${seqz}"
    return 0
  fi
  echo "[$(date -Is)] bam2seqz ${SAMPLE_ID} ${chrom}"
  "${NEOAG_CONDA_BASE}/bin/mamba" run -n "${ENV}" sequenza-utils bam2seqz \
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
    -o "${seqz}"
}
export -f run_chrom
export SAMPLE_ID TUMOR_BAM NORMAL_BAM REF GC OUTDIR ENV NEOAG_CONDA_BASE SAMTOOLS TABIX QLIMIT MIN_DEPTH_N HOM HET
printf "%s\n" ${CHROMS} | xargs -I{} -P "${CHUNK_JOBS}" bash -c "run_chrom \"{}\""

echo "[$(date -Is)] merge chrom seqz"
MERGED="${OUTDIR}/${SAMPLE_ID}.merged.seqz.gz"
{
  first=1
  for chrom in ${CHROMS}; do
    f="${OUTDIR}/chrom/${SAMPLE_ID}.${chrom}.seqz.gz"
    [[ -s "$f" ]] || { echo "ERROR missing chrom seqz $f" >&2; exit 1; }
    if [[ "$first" == 1 ]]; then
      zcat "$f"
      first=0
    else
      zcat "$f" | tail -n +2
    fi
  done
} | gzip -c > "${MERGED}.tmp"
mv "${MERGED}.tmp" "${MERGED}"

BINNED="${OUTDIR}/${SAMPLE_ID}.small.seqz.gz"
echo "[$(date -Is)] seqz_binning"
run_env sequenza-utils seqz_binning -s "${MERGED}" -w "${BIN_WINDOW}" -T "${TABIX}" -o "${BINNED}"

echo "[$(date -Is)] R fit"
run_env Rscript "${ROOT}/scripts/run_sequenza_fit.R" "${BINNED}" "${OUTDIR}/sequenza_fit" "${SAMPLE_ID}"
echo "[$(date -Is)] finished ${SAMPLE_ID}"
