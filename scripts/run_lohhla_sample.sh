#!/usr/bin/env bash
# Generic WGS LOHHLA runner: HLA typing or HLA file -> patient HLA FASTA -> LOHHLA.
#
# Prerequisites:
#   bash scripts/install_polysolver.sh
#   LOHHLA tools (LOHHLAscript.R + data/hla.dat + picard.jar)
#   neoag-fusion conda env (R + LOHHLA R packages)
#
# Usage:
#   source conf/tools.env.sh
#   PATIENT_ID=sample TUMOR_BAM=/path/tumor.bam NORMAL_BAM=/path/normal.bam bash scripts/run_lohhla_sample.sh
#   LOHHLA_STEP=polysolver bash scripts/run_lohhla_sample.sh
#   LOHHLA_STEP=fasta      bash scripts/run_lohhla_sample.sh
#   LOHHLA_STEP=lohhla     bash scripts/run_lohhla_sample.sh
#
# Environment overrides:
#   TUMOR_BAM, NORMAL_BAM, PATIENT_ID, OUTDIR, LOHHLA_STEP
#   LOHHLA_NAS_ROOT  - scratch for LOHHLA intermediates (default under project work/lohhla)
#   LOHHLA_OUT       — LOHHLA --outputDir (defaults to ${LOHHLA_NAS_ROOT}/lohhla)
#   POLYSOLVER_RACE=Unknown  POLYSOLVER_BUILD=hg19  POLYSOLVER_FORMAT=STDFQ
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
_LOHHLA_GATK_DIR_OVERRIDE="${LOHHLA_GATK_DIR:-}"
# shellcheck source=/dev/null
source "${ROOT}/conf/tools.env.sh"
if [[ -n "${_LOHHLA_GATK_DIR_OVERRIDE}" ]]; then
  LOHHLA_GATK_DIR="${_LOHHLA_GATK_DIR_OVERRIDE}"
fi
unset _LOHHLA_GATK_DIR_OVERRIDE
: "${NEOAG_CONDA_BASE:?ERROR: set NEOAG_CONDA_BASE to your conda/mamba installation root}"
export PATH="${NEOAG_CONDA_BASE}/bin:${PATH}"
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"

PATIENT_ID="${PATIENT_ID:-sample}"
TUMOR_SAMPLE="${TUMOR_SAMPLE_ID:-${PATIENT_ID}}"
NORMAL_SAMPLE="${NORMAL_SAMPLE_ID:-normal}"

TUMOR_BAM="${TUMOR_BAM:?ERROR: set TUMOR_BAM=/path/tumor.bam}"
NORMAL_BAM="${NORMAL_BAM:?ERROR: set NORMAL_BAM=/path/normal.bam}"
HLA_FILE="${HLA_FILE:-${ROOT}/work/${PATIENT_ID}.hla.txt}"
OUT="${OUTDIR:-${ROOT}/results/lohhla/${PATIENT_ID}}"
LOHHLA_NAS_ROOT="${LOHHLA_NAS_ROOT:-${ROOT}/work/lohhla/${PATIENT_ID}}"
LOHHLA_OUT="${LOHHLA_OUT:-${LOHHLA_NAS_ROOT}/lohhla}"
LOHHLA_STEP="${LOHHLA_STEP:-all}"
LOG="${LOG:-${ROOT}/work/run_lohhla_${PATIENT_ID}.log}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-${NEOAG_CONDA_BASE}/envs/${NEOAG_GATK_ENV:-neoag-gatk}/bin/samtools}"
BAM_LINK_DIR="${BAM_LINK_DIR:-${ROOT}/work/bam_links/${PATIENT_ID}}"

resolve_polysolver_home() {
  local candidate
  for candidate in \
    "${POLYSOLVER_HOME:-}" \
    "${NEOAG_LICENSED_ROOT:+${NEOAG_LICENSED_ROOT}/polysolver}" \
    "${NEOAG_ASSET_ROOT:+${NEOAG_ASSET_ROOT}/data/lohhla/polysolver}" \
    "${NEOAG_ASSET_ROOT:+${NEOAG_ASSET_ROOT}/data/polysolver}" \
    "${NEOAG_PUBLIC_ASSET_ROOT:+${NEOAG_PUBLIC_ASSET_ROOT}/data/lohhla/polysolver}" \
    "${NEOAG_TOOLS_ROOT:-${ROOT}}/tools/polysolver"; do
    if [[ -n "${candidate}" \
      && -x "${candidate}/scripts/shell_call_hla_type" \
      && -s "${candidate}/data/abc_complete.fasta" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  return 1
}

PSHOME="$(resolve_polysolver_home || true)"
PSHOME="${PSHOME:-${NEOAG_TOOLS_ROOT:-${ROOT}}/tools/polysolver}"
LOHHLA_HOME="${LOHHLA_HOME:-${NEOAG_TOOLS_ROOT}/tools/lohhla}"
GATK_ENV="${NEOAG_CONDA_BASE}/envs/${NEOAG_GATK_ENV:-neoag-gatk}"

resolve_lohhla_env() {
  local candidate
  local -a candidates=()
  [[ -n "${LOHHLA_ENV_PREFIX:-}" ]] && candidates+=("${LOHHLA_ENV_PREFIX}")
  candidates+=(
    "${NEOAG_CONDA_BASE}/envs/neoag-lohhla"
    "${NEOAG_CONDA_BASE}/envs/${NEOAG_FUSION_ENV:-neoag-fusion}"
    "${NEOAG_CONDA_BASE}/envs/neoag-r"
    "${NEOAG_CONDA_BASE}/envs/neoag-splicemutr"
  )
  for candidate in "${candidates[@]}"; do
    [[ -x "${candidate}/bin/Rscript" ]] || continue
    if "${candidate}/bin/Rscript" -e \
      'quit(status=ifelse(requireNamespace("optparse",quietly=TRUE)&&requireNamespace("Rsamtools",quietly=TRUE)&&requireNamespace("Biostrings",quietly=TRUE)&&requireNamespace("seqinr",quietly=TRUE),0,1))' \
      >/dev/null 2>&1; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  echo "ERROR: no LOHHLA R environment provides optparse, Rsamtools, Biostrings, and seqinr" >&2
  return 1
}

FUSION_ENV="$(resolve_lohhla_env)"

# Some deployments split the Bioconductor stack and the small CRAN plotting
# dependencies across R environments. Reuse a compatible installed library
# instead of letting LOHHLA attempt an interactive install into a read-only env.
if ! "${FUSION_ENV}/bin/Rscript" -e \
  'quit(status=ifelse(requireNamespace("beeswarm", quietly=TRUE), 0, 1))' \
  >/dev/null 2>&1; then
  for candidate in \
    "${NEOAG_CONDA_BASE}/envs/neoag-r/lib/R/library" \
    "${NEOAG_CONDA_BASE}/envs/neoag-lohhla-r36/lib/R/library"; do
    [[ -d "${candidate}/beeswarm" ]] || continue
    export R_LIBS_USER="${candidate}${R_LIBS_USER:+:${R_LIBS_USER}}"
    break
  done
fi

resolve_bedtools_dir() {
  local candidate
  if command -v bedtools >/dev/null 2>&1; then
    dirname "$(command -v bedtools)"
    return 0
  fi
  for candidate in \
    "${NEOAG_TOOLS_ROOT:-${ROOT}}/bin" \
    "${NEOAG_CONDA_BASE}/envs/neoag-tools/bin" \
    "${NEOAG_CONDA_BASE}/envs/neoag-splice/bin" \
    "${FUSION_ENV}/bin"; do
    if [[ -x "${candidate}/bedtools" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  echo "ERROR: bedtools is required by LOHHLA but was not found" >&2
  return 1
}

LOHHLA_BEDTOOLS_DIR="$(resolve_bedtools_dir)"

POLYSOLVER_RACE="${POLYSOLVER_RACE:-Unknown}"
POLYSOLVER_INCLUDE_FREQ="${POLYSOLVER_INCLUDE_FREQ:-1}"
POLYSOLVER_BUILD="${POLYSOLVER_BUILD:-hg19}"
POLYSOLVER_FORMAT="${POLYSOLVER_FORMAT:-STDFQ}"
POLYSOLVER_INSERT_CALC="${POLYSOLVER_INSERT_CALC:-1}"
POLYSOLVER_THREADS="${POLYSOLVER_THREADS:-8}"
COPYNUM_LOC="${COPYNUM_LOC:-FALSE}"
MIN_COVERAGE="${MIN_COVERAGE_FILTER:-10}"
LOHHLA_MAPPING_STEP="${LOHHLA_MAPPING_STEP:-TRUE}"
LOHHLA_FISHING_STEP="${LOHHLA_FISHING_STEP:-TRUE}"
LOHHLA_GENOME_ASSEMBLY="${LOHHLA_GENOME_ASSEMBLY:-grch38}"

PS_OUT="${OUT}/polysolver"
HLA_DIR="${OUT}/hla"
WINNERS="${PS_OUT}/winners.hla.txt"
PATIENT_HLA_FASTA="${HLA_DIR}/${PATIENT_ID}.patient.hlaFasta.fa"
LOHHLA_SCRIPT="${LOHHLA_HOME}/LOHHLAscript.R"
HLA_EXON_LOC="${LOHHLA_HOME}/data/hla.dat"
if [[ -n "${LOHHLA_GATK_DIR:-}" ]]; then
  LOHHLA_GATK_RUNTIME_DIR="${LOHHLA_GATK_DIR}"
elif [[ -f "${PSHOME}/binaries/SamToFastq.jar" \
  && -f "${PSHOME}/binaries/SortSam.jar" \
  && -f "${PSHOME}/binaries/FilterSamReads.jar" ]]; then
  LOHHLA_GATK_RUNTIME_DIR="${PSHOME}/binaries"
else
  LOHHLA_GATK_RUNTIME_DIR="${LOHHLA_HOME}"
fi
NOVO_DIR="${PSHOME}/binaries"
NOVOINDEX="${PSHOME}/scripts/novoindex"
ABC_FASTA="${PSHOME}/data/abc_complete.fasta"
ABC_COMPLETE_DIR="${PSHOME}/data/complete"

mkdir -p "${OUT}" "${PS_OUT}" "${HLA_DIR}" "${LOHHLA_NAS_ROOT}" "${LOHHLA_OUT}" "$(dirname "${LOG}")"
LOCAL_LOHHLA_LINK="${OUT}/lohhla"
if [[ "${LOHHLA_OUT}" != "${LOCAL_LOHHLA_LINK}" ]]; then
  rm -rf "${LOCAL_LOHHLA_LINK}"
  ln -sfn "${LOHHLA_OUT}" "${LOCAL_LOHHLA_LINK}"
fi
exec > >(tee -a "${LOG}") 2>&1

echo "==> run_lohhla_sample $(date -Is)"
echo "    step=${LOHHLA_STEP}"
echo "    patient=${PATIENT_ID}"
echo "    tumor_bam=${TUMOR_BAM}"
echo "    normal_bam=${NORMAL_BAM}"
echo "    out=${OUT}"
echo "    lohhla_nas_root=${LOHHLA_NAS_ROOT}"
echo "    lohhla_out=${LOHHLA_OUT}"

resolve_lohhla_home() {
  if [[ -f "${LOHHLA_HOME}/LOHHLAscript.R" ]]; then
    return 0
  fi
  echo "ERROR: LOHHLAscript.R not found. Set LOHHLA_HOME to a tree with LOHHLAscript.R, data/hla.dat, picard.jar" >&2
  exit 1
}

step_wanted() {
  [[ "${LOHHLA_STEP}" == "all" || "${LOHHLA_STEP}" == "$1" ]]
}

bam_ready() {
  local bam="$1"
  [[ -s "${bam}" && -f "${bam}.bai" ]]
}

# When BAM lives on a root-owned NFS dir, index to work/ and expose via symlink pair.
resolve_bam_for_tools() {
  local var_name="$1"
  local bam="$2"
  mkdir -p "${BAM_LINK_DIR}"
  local base
  base="$(basename "${bam}")"
  local role="${var_name%_bam}"
  local link_bam="${BAM_LINK_DIR}/${PATIENT_ID}_${role}.bam"
  local link_bai="${link_bam}.bai"
  local bai_src="${ROOT}/work/${base}.bai"

  if [[ ! -L "${link_bam}" ]]; then
    ln -sf "${bam}" "${link_bam}"
  fi
  if [[ -f "${bai_src}" && ! -f "${link_bai}" ]]; then
    ln -sf "${bai_src}" "${link_bai}"
  elif [[ -f "${bam}.bai" && ! -f "${link_bai}" ]]; then
    ln -sf "${bam}.bai" "${link_bai}"
  fi

  if bam_ready "${link_bam}"; then
    printf -v "${var_name}" '%s' "${link_bam}"
  else
    printf -v "${var_name}" '%s' "${bam}"
  fi
}

ensure_bam_index() {
  local bam="$1"
  [[ -s "${bam}" ]] || {
    echo "ERROR: BAM missing: ${bam}" >&2
    exit 1
  }
  if [[ -f "${bam}.bai" ]]; then
    return 0
  fi
  local base
  base="$(basename "${bam}")"
  local bai_out="${ROOT}/work/${base}.bai"
  if [[ -f "${bai_out}" ]]; then
  mkdir -p "${BAM_LINK_DIR}"
    ln -sf "${bai_out}" "${BAM_LINK_DIR}/${base}.bai"
    return 0
  fi
  [[ -x "${SAMTOOLS_BIN}" ]] || {
    echo "ERROR: samtools not found: ${SAMTOOLS_BIN}" >&2
    exit 1
  }
  echo "==> indexing BAM with ${SAMTOOLS_BIN} (output ${bai_out}; may take 1–2 h on WGS)"
  "${SAMTOOLS_BIN}" index -@ "${POLYSOLVER_THREADS}" -o "${bai_out}" "${bam}"
  [[ -f "${bai_out}" ]] || {
    echo "ERROR: failed to create ${bai_out}" >&2
    exit 1
  }
  mkdir -p "${BAM_LINK_DIR}"
  ln -sf "${bai_out}" "${BAM_LINK_DIR}/${base}.bai"
}

# Map HLA-A*02:06 → polysolver abc_complete ID (hla_a_02_06_01).
hla_to_polysolver_id() {
  local allele="$1"
  local locus="${allele#HLA-}"
  locus="${locus%%\**}"
  local nums="${allele#*\*}"
  nums="${nums//:/_}"
  case "${locus}" in
    A) echo "hla_a_${nums}_01" ;;
    B) echo "hla_b_${nums}_01" ;;
    C)
      if [[ "${nums}" == "06_02" ]]; then
        echo "hla_c_06_02_01_01"
      else
        echo "hla_c_${nums}_01"
      fi
      ;;
    *) echo "hla_${locus,,}_${nums}_01" ;;
  esac
}

write_winners_from_hla_file() {
  [[ -f "${HLA_FILE}" ]] || {
    echo "ERROR: HLA_FILE not found: ${HLA_FILE}" >&2
    return 1
  }
  mapfile -t _alleles < <(grep -v '^#' "${HLA_FILE}" | tr ', ' '\n' | sed '/^$/d')
  [[ ${#_alleles[@]} -ge 2 ]] || {
    echo "ERROR: need at least 2 HLA alleles in ${HLA_FILE}" >&2
    return 1
  }
  declare -A by_locus=()
  for a in "${_alleles[@]}"; do
    locus="${a%%\**}"
    by_locus["${locus}"]+="${a} "
  done
  mkdir -p "${PS_OUT}"
  {
    for locus in HLA-A HLA-B HLA-C; do
      read -r -a pair <<< "${by_locus[$locus]:-}"
      [[ ${#pair[@]} -ge 1 ]] || continue
      id1="$(hla_to_polysolver_id "${pair[0]}")"
      # Homozygous consensus (one allele/locus): duplicate for Polysolver/LOHHLA pair format.
      if [[ ${#pair[@]} -ge 2 ]]; then
        id2="$(hla_to_polysolver_id "${pair[1]}")"
      else
        id2="${id1}"
      fi
      printf '%s\t%s\t%s\n' "${locus}" "${id1}" "${id2}"
    done
  } > "${WINNERS}"
  echo "==> wrote winners from ${HLA_FILE}:"
  cat "${WINNERS}"
}

ensure_polysolver() {
  local resolved_pshome="${PSHOME}"
  [[ -f "${PSHOME}/scripts/config.local.bash" ]] || {
    echo "ERROR: Polysolver not configured. Run: bash scripts/install_polysolver.sh" >&2
    exit 1
  }
  # shellcheck source=/dev/null
  source "${PSHOME}/scripts/config.local.bash"
  # Licensed asset bundles are commonly copied between machines. Their
  # generated config.local.bash may still contain the source machine's
  # absolute PSHOME, so the explicitly resolved deployment path must win.
  PSHOME="${resolved_pshome}"
  SAMTOOLS_DIR="${PSHOME}/binaries"
  NOVOALIGN_DIR="${PSHOME}/binaries"
  GATK_DIR="${PSHOME}/binaries"
  MUTECT_DIR="${PSHOME}/binaries"
  STRELKA_DIR="${PSHOME}/binaries"
  TMP_DIR="${PSHOME}/sachet"
  export PSHOME SAMTOOLS_DIR NOVOALIGN_DIR GATK_DIR MUTECT_DIR STRELKA_DIR TMP_DIR
  export PATH="${NEOAG_CONDA_BASE}/envs/${NEOAG_TOOLS_ENV:-neoag-tools}/bin:${PSHOME}/binaries:${PSHOME}/scripts:${PATH}"
  if [[ -f "${NOVOALIGN_LICENSE_FILE:-}" ]]; then
    if [[ "$(readlink -f "${NOVOALIGN_LICENSE_FILE}")" != "$(readlink -f "${NOVO_DIR}/novoalign.lic" 2>/dev/null || true)" ]]; then
      cp -f "${NOVOALIGN_LICENSE_FILE}" "${NOVO_DIR}/novoalign.lic"
    fi
  fi
}

run_polysolver() {
  if [[ -s "${WINNERS}" ]]; then
    echo "==> skip polysolver: ${WINNERS} exists"
    return 0
  fi

  if [[ "${USE_HLA_FILE_FOR_WINNERS:-1}" == "1" && -f "${HLA_FILE}" ]]; then
    echo "==> using HLA_FILE for winners (skip Polysolver WGS typing): ${HLA_FILE}"
    write_winners_from_hla_file
    return 0
  fi

  ensure_polysolver
  ensure_bam_index "${NORMAL_BAM}"
  bam_ready "${NORMAL_BAM}" || {
    echo "ERROR: normal BAM+BAI required for Polysolver: ${NORMAL_BAM}" >&2
    exit 1
  }

  export NUM_THREADS="${POLYSOLVER_THREADS}"
  echo "==> Polysolver HLA typing on normal BAM (heavy; may take hours on WGS) ..."
  bash "${PSHOME}/scripts/shell_call_hla_type" \
    "${NORMAL_BAM}" \
    "${POLYSOLVER_RACE}" \
    "${POLYSOLVER_INCLUDE_FREQ}" \
    "${POLYSOLVER_BUILD}" \
    "${POLYSOLVER_FORMAT}" \
    "${POLYSOLVER_INSERT_CALC}" \
    "${PS_OUT}"

  [[ -s "${WINNERS}" ]] || {
    echo "ERROR: Polysolver did not produce ${WINNERS}" >&2
    exit 1
  }
  echo "==> Polysolver winners:"
  cat "${WINNERS}"
}

run_build_hla_fasta() {
  ensure_polysolver
  [[ -s "${WINNERS}" ]] || {
    echo "ERROR: ${WINNERS} missing. Run LOHHLA_STEP=polysolver first." >&2
    exit 1
  }

  if [[ -s "${PATIENT_HLA_FASTA}" && -f "${PATIENT_HLA_FASTA}.fai" ]]; then
    echo "==> skip patient HLA fasta: ${PATIENT_HLA_FASTA} exists"
    return 0
  fi

  export PATH="${GATK_ENV}/bin:${PATH}"
  echo "==> building patient HLAfastaLoc from ${WINNERS} ..."
  python3 "${ROOT}/scripts/build_patient_hla_fasta.py" \
    --winners "${WINNERS}" \
    --ref-fasta "${ABC_FASTA}" \
    --complete-dir "${ABC_COMPLETE_DIR}" \
    --out-fasta "${PATIENT_HLA_FASTA}" \
    --samtools samtools \
    --novoindex "${NOVOINDEX}" || {
      [[ -s "${PATIENT_HLA_FASTA}" ]] || exit 1
      echo "WARN: novoindex failed; continuing with FASTA only" >&2
    }

  samtools faidx "${PATIENT_HLA_FASTA}"
  echo "==> HLAfastaLoc ready: ${PATIENT_HLA_FASTA}"
}

run_lohhla() {
  resolve_lohhla_home
  ensure_bam_index "${TUMOR_BAM}"
  ensure_bam_index "${NORMAL_BAM}"
  local tumor_bam="${TUMOR_BAM}"
  local normal_bam="${NORMAL_BAM}"
  local flagstat_dir="${LOHHLA_OUT}/flagstat"
  local override_dir="FALSE"
  resolve_bam_for_tools tumor_bam "${TUMOR_BAM}"
  resolve_bam_for_tools normal_bam "${NORMAL_BAM}"
  bam_ready "${tumor_bam}" && bam_ready "${normal_bam}" || {
    echo "ERROR: tumor/normal BAM+BAI required for LOHHLA." >&2
    echo "  tumor:  ${tumor_bam} (src ${TUMOR_BAM})" >&2
    echo "  normal: ${normal_bam} (src ${NORMAL_BAM})" >&2
    exit 1
  }
  [[ -s "${PATIENT_HLA_FASTA}" ]] || {
    echo "ERROR: patient HLA FASTA missing. Run LOHHLA_STEP=fasta first." >&2
    exit 1
  }
  [[ -f "${LOHHLA_SCRIPT}" && -f "${HLA_EXON_LOC}" ]] || {
    echo "ERROR: LOHHLA install incomplete under ${LOHHLA_HOME}" >&2
    exit 1
  }
  [[ -f "${LOHHLA_GATK_RUNTIME_DIR}/picard.jar" ]] || {
    echo "ERROR: picard.jar missing in ${LOHHLA_GATK_RUNTIME_DIR}" >&2
    exit 1
  }

  for _jhome in "${GATK_ENV}/lib/jvm" "${FUSION_ENV}/lib/jvm"; do
    if [[ -x "${_jhome}/bin/java" && -f "${_jhome}/lib/libjli.so" ]]; then
      export JAVA_HOME="${_jhome}"
      export LD_LIBRARY_PATH="${JAVA_HOME}/lib:${JAVA_HOME}/lib/server:${LD_LIBRARY_PATH:-}"
      break
    fi
  done
  unset _jhome
  export PATH="${FUSION_ENV}/bin:${LOHHLA_BEDTOOLS_DIR}:${NOVO_DIR}:${PSHOME}/scripts:${GATK_ENV}/bin:${PATH}"
  if [[ -n "${JAVA_HOME:-}" ]]; then
    export PATH="${JAVA_HOME}/bin:${PATH}"
  fi
  if [[ -f "${NOVOALIGN_LICENSE_FILE:-}" ]]; then
    if [[ "$(readlink -f "${NOVOALIGN_LICENSE_FILE}")" != "$(readlink -f "${NOVO_DIR}/novoalign.lic" 2>/dev/null || true)" ]]; then
      cp -f "${NOVOALIGN_LICENSE_FILE}" "${NOVO_DIR}/novoalign.lic"
    fi
  fi

  echo "==> LOHHLA (tumor vs normal HLA LOH); intermediates on NAS: ${LOHHLA_OUT}"
  [[ -n "${JAVA_HOME:-}" ]] && echo "    JAVA_HOME=${JAVA_HOME}"
  mkdir -p "${LOHHLA_OUT}"
  if [[ -s "${flagstat_dir}/$(basename "${tumor_bam}").proc.flagstat" && -s "${flagstat_dir}/$(basename "${normal_bam}").proc.flagstat" ]]; then
    # LOHHLA reuses files in its default flagstat directory automatically.
    # Passing that same directory through --overrideDir triggers an upstream
    # character-to-logical conversion bug in some LOHHLA revisions.
    echo "    reusing flagstat from ${flagstat_dir}"
  fi
  local bam_dir_abs normal_bam_dir_abs tumor_bam_name normal_bam_name
  bam_dir_abs="$(cd "$(dirname "${tumor_bam}")" && pwd -P)"
  normal_bam_dir_abs="$(cd "$(dirname "${normal_bam}")" && pwd -P)"
  [[ "${bam_dir_abs}" == "${normal_bam_dir_abs}" ]] || {
    echo "ERROR: LOHHLA requires tumor and normal BAM links in one directory" >&2
    exit 1
  }
  tumor_bam_name="$(basename "${tumor_bam}")"
  normal_bam_name="$(basename "${normal_bam}")"
  local lohhla_out_abs winners_abs hla_fasta_abs copynum_abs
  lohhla_out_abs="$(mkdir -p "${LOHHLA_OUT}" && cd "${LOHHLA_OUT}" && pwd -P)"
  winners_abs="$(cd "$(dirname "${WINNERS}")" && pwd -P)/$(basename "${WINNERS}")"
  hla_fasta_abs="$(cd "$(dirname "${PATIENT_HLA_FASTA}")" && pwd -P)/$(basename "${PATIENT_HLA_FASTA}")"
  if [[ "${COPYNUM_LOC}" == "FALSE" ]]; then
    copynum_abs="FALSE"
  else
    copynum_abs="$(cd "$(dirname "${COPYNUM_LOC}")" && pwd -P)/$(basename "${COPYNUM_LOC}")"
  fi
  Rscript "${LOHHLA_SCRIPT}" \
    --patientId "${PATIENT_ID}" \
    --outputDir "${lohhla_out_abs}" \
    --BAMDir "${bam_dir_abs}" \
    --normalBAMfile "${normal_bam_name}" \
    --tumorBAMfile "${tumor_bam_name}" \
    --hlaPath "${winners_abs}" \
    --HLAfastaLoc "${hla_fasta_abs}" \
    --CopyNumLoc "${copynum_abs}" \
    --overrideDir "${override_dir}" \
    --mappingStep "${LOHHLA_MAPPING_STEP}" \
    --fishingStep "${LOHHLA_FISHING_STEP}" \
    --plottingStep FALSE \
    --coverageStep TRUE \
    --minCoverageFilter "${MIN_COVERAGE}" \
    --cleanUp FALSE \
    --gatkDir "${LOHHLA_GATK_RUNTIME_DIR}" \
    --LOHHLA_loc "${LOHHLA_HOME}" \
    --novoDir "${NOVO_DIR}" \
    --HLAexonLoc "${HLA_EXON_LOC}"

  echo ""
  echo "==> Done. Key outputs:"
  echo "    Polysolver:  ${WINNERS}"
  echo "    HLAfastaLoc: ${PATIENT_HLA_FASTA}"
  echo "    LOHHLA dir:  ${LOHHLA_OUT}"
  echo "    Convert:     neoag convert-lohhla -i <*HLAlossPrediction_CI*> -o hla_loh.tsv"
}

if step_wanted polysolver; then
  run_polysolver
fi
if step_wanted fasta; then
  run_build_hla_fasta
fi
if step_wanted lohhla; then
  run_lohhla
fi

echo "==> run_lohhla_sample finished $(date -Is)"
