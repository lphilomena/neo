#!/usr/bin/env bash
# Closed-loop deployment for external neoag-v03 tools.
#
# Designed for a fresh machine/release checkout:
#   - finds conda even when it is not yet on PATH
#   - creates missing conf/bin/tools directories
#   - runs each installer only when the tool is not already detected
#   - re-sources the generated environment before verification
#
# Usage:
#   bash scripts/deploy_external_tools.sh
#   bash scripts/deploy_external_tools.sh --smoke
#   bash scripts/deploy_external_tools.sh --preflight-only
#   NEOAG_REF_BUNDLE=/path/to/neodata4git bash scripts/deploy_external_tools.sh --smoke
#   NEOAG_CONDA_BASE=/path/to/miniforge3 bash scripts/deploy_external_tools.sh
#   NEOAG_FORCE_INSTALL=1 bash scripts/deploy_external_tools.sh
#   SKIP_LOHHLA=1 SKIP_PRIME=1 SKIP_OPTITYPE=1 bash scripts/deploy_external_tools.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SMOKE=0
PREFLIGHT_ONLY=0
SKIP_VERIFY="${NEOAG_SKIP_VERIFY:-0}"

usage() {
  cat <<USAGE
Usage: bash scripts/deploy_external_tools.sh [--smoke] [--preflight-only] [--skip-verify]

Environment switches:
  NEOAG_CONDA_BASE=/path/to/miniforge3  Conda installation root if conda is not on PATH
  NEOAG_TOOLS_ROOT=/path/to/tools_root   External tool/data root; defaults to this checkout
  NEOAG_REF_BUNDLE=/path/to/neodata4git  Portable reference bundle; auto-sources neodata4git.env.sh
  NEOAG_FORCE_INSTALL=1                  Re-run installers even when wrappers are detected
  NEOAG_FORCE_ENV_UPDATE=1               Force conda env update inside installers
  SKIP_LOHHLA=1                          Skip LOHHLA
  SKIP_FACETS=1                          Skip FACETS
  SKIP_ASCAT=1                           Skip ASCAT / PyClone-VI
  SKIP_ARRIBA=1                          Skip Arriba / fusion env
  SKIP_PRIME=1                           Skip PRIME / MixMHCpred / BigMHC
  SKIP_OPTITYPE=1                        Skip OptiType
USAGE
}

for arg in "$@"; do
  case "${arg}" in
    --smoke) SMOKE=1 ;;
    --preflight-only) PREFLIGHT_ONLY=1 ;;
    --skip-verify) SKIP_VERIFY=1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: ${arg}" >&2; usage >&2; exit 2 ;;
  esac
done

info() { echo "==> $*"; }
warn() { echo "WARN: $*" >&2; }
die() { echo "ERROR: $*" >&2; exit 1; }

find_conda() {
  if command -v conda >/dev/null 2>&1; then
    command -v conda
    return 0
  fi

  local base="${NEOAG_CONDA_BASE:-}"
  local candidates=()
  [[ -n "${base}" ]] && candidates+=("${base}/bin/conda")
  candidates+=(
    "${HOME}/miniforge3/bin/conda"
    "${HOME}/mambaforge/bin/conda"
    "${HOME}/miniconda3/bin/conda"
    "${HOME}/anaconda3/bin/conda"
    "/opt/conda/bin/conda"
  )

  local c
  for c in "${candidates[@]}"; do
    if [[ -x "${c}" ]]; then
      echo "${c}"
      return 0
    fi
  done
  return 1
}

bootstrap_conda() {
  local conda_bin
  conda_bin="$(find_conda)" || die "conda not found. Install Miniforge/Mambaforge, or set NEOAG_CONDA_BASE=/path/to/conda_root."
  export PATH="$(dirname "${conda_bin}"):${PATH}"
  export NEOAG_CONDA_BASE="${NEOAG_CONDA_BASE:-$("${conda_bin}" info --base 2>/dev/null)}"
  [[ -n "${NEOAG_CONDA_BASE}" && -d "${NEOAG_CONDA_BASE}" ]] || die "Cannot determine conda base from ${conda_bin}"
  if [[ -f "${NEOAG_CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
    # shellcheck source=/dev/null
    source "${NEOAG_CONDA_BASE}/etc/profile.d/conda.sh"
  fi
  info "Conda: ${NEOAG_CONDA_BASE}"
}

ensure_layout() {
  mkdir -p "${ROOT}/bin" "${ROOT}/conf" "${ROOT}/tools" "${ROOT}/data" "${ROOT}/work"
  export NEOAG_PROJECT_ROOT="${ROOT}"
  export NEOAG_TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-${ROOT}}"
  export PATH="${ROOT}/bin:${PATH}"

  if [[ ! -f "${ROOT}/conf/tools.env.sh" ]]; then
    cat > "${ROOT}/conf/tools.env.sh" <<EOF_ENV
# Source before neoag-v03 runs: source conf/tools.env.sh
_NEOAG_TOOLS_ENV_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
export NEOAG_PROJECT_ROOT="\$(cd "\${_NEOAG_TOOLS_ENV_DIR}/.." && pwd)"
export NEOAG_TOOLS_ROOT="\${NEOAG_TOOLS_ROOT:-\${NEOAG_PROJECT_ROOT}}"
export NEOAG_CONDA_BASE="\${NEOAG_CONDA_BASE:-${NEOAG_CONDA_BASE}}"
export NEOAG_CONDA_ENV="\${NEOAG_CONDA_ENV:-neoag-tools}"
export PATH="\${NEOAG_PROJECT_ROOT}/bin:\${NEOAG_CONDA_BASE}/envs/\${NEOAG_CONDA_ENV}/bin:\${PATH}"
if [[ -f "\${NEOAG_PROJECT_ROOT}/conf/tools.env.local.sh" ]]; then
  source "\${NEOAG_PROJECT_ROOT}/conf/tools.env.local.sh"
fi
unset _NEOAG_TOOLS_ENV_DIR
EOF_ENV
    info "Created conf/tools.env.sh"
  fi
}

source_tools_env_optional() {
  # tools.env.sh may reference optional tools that are not present on a fresh machine.
  # Source it only after conda bootstrap, and tolerate failures during early deployment.
  if [[ -f "${ROOT}/conf/tools.env.sh" ]]; then
    # shellcheck source=/dev/null
    source "${ROOT}/conf/tools.env.sh" || warn "conf/tools.env.sh returned non-zero; continuing deployment"
  fi

  local ref_env=""
  if [[ -n "${NEOAG_REF_BUNDLE:-}" && -f "${NEOAG_REF_BUNDLE}/neodata4git.env.sh" ]]; then
    ref_env="${NEOAG_REF_BUNDLE}/neodata4git.env.sh"
  elif [[ -f "${NEOAG_TOOLS_ROOT:-${ROOT}}/neodata4git.env.sh" ]]; then
    ref_env="${NEOAG_TOOLS_ROOT:-${ROOT}}/neodata4git.env.sh"
    export NEOAG_REF_BUNDLE="${NEOAG_TOOLS_ROOT:-${ROOT}}"
  fi
  if [[ -n "${ref_env}" ]]; then
    # shellcheck source=/dev/null
    source "${ref_env}" || warn "${ref_env} returned non-zero; continuing deployment"
    export NEOAG_REF_BUNDLE="$(cd "$(dirname "${ref_env}")" && pwd)"
  fi
  if [[ -f "${ROOT}/conf/tools.env.local.sh" ]]; then
    # shellcheck source=/dev/null
    source "${ROOT}/conf/tools.env.local.sh" || warn "conf/tools.env.local.sh returned non-zero; continuing deployment"
  fi

  export NEOAG_CONDA_BASE="${NEOAG_CONDA_BASE:-$(conda info --base 2>/dev/null || true)}"
  export NEOAG_TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-${ROOT}}"
  export PATH="${ROOT}/bin:${PATH}"
}

preflight() {
  bootstrap_conda
  ensure_layout
  source_tools_env_optional

  command -v conda >/dev/null 2>&1 || die "conda still not available after bootstrap"
  command -v git >/dev/null 2>&1 || die "git not found. Install git and retry."

  local missing=0
  local required_scripts=(
    scripts/install_lohhla.sh
    scripts/install_facets.sh
    scripts/install_ascat_pyclone.sh
    scripts/install_fusion_tools.sh
    scripts/install_immunogenicity_tools.sh
    scripts/install_optitype.sh
    scripts/verify_external_tools.sh
    scripts/verify_reference_bundle.sh
    scripts/verify_all_tools_and_refs.sh
  )
  local s
  for s in "${required_scripts[@]}"; do
    if [[ ! -f "${ROOT}/${s}" ]]; then
      warn "Missing ${s}"
      missing=1
    fi
  done
  [[ "${missing}" == "0" ]] || die "One or more deployment scripts are missing from this release."

  if ! command -v curl >/dev/null 2>&1; then
    warn "curl not found; PRIME download fallback may fail. Install curl if PRIME is enabled."
  fi
  if ! command -v g++ >/dev/null 2>&1; then
    warn "g++ not found; PRIME compilation may fail. Install a C++ compiler if PRIME is enabled."
  fi
}

run_step() {
  local label="$1"
  shift
  echo
  info "${label}"
  "$@"
  source_tools_env_optional
}

skip_step() { info "Skip ${1} (already installed)"; }

lohhla_installed() {
  command -v LOHHLA >/dev/null 2>&1 && [[ -f "${LOHHLA_HOME:-${NEOAG_TOOLS_ROOT:-${ROOT}}/tools/lohhla}/LOHHLAscript.R" ]]
}

facets_installed() {
  command -v runFACETS.R >/dev/null 2>&1 && runFACETS.R --version >/dev/null 2>&1
}

ascat_installed() {
  command -v ascat.R >/dev/null 2>&1 && ascat.R --version >/dev/null 2>&1
}

arriba_installed() {
  local fusion_env="${NEOAG_FUSION_ENV:-neoag-fusion}"
  if [[ -n "${NEOAG_CONDA_BASE:-}" && -x "${NEOAG_CONDA_BASE}/envs/${fusion_env}/bin/arriba" ]]; then
    return 0
  fi
  command -v arriba >/dev/null 2>&1
}

prime_installed() {
  local prime_bin="${NEOAG_PRIME_BIN:-${NEOAG_TOOLS_ROOT:-${ROOT}}/tools/prime/PRIME}"
  [[ -x "${prime_bin}" ]]
}

optitype_installed() {
  command -v optitype >/dev/null 2>&1 && optitype check-deps >/dev/null 2>&1
}

should_install() {
  local skip_flag="$1" installed_fn="$2"
  [[ "${skip_flag}" != "1" ]] || return 1
  [[ "${NEOAG_FORCE_INSTALL:-0}" == "1" ]] && return 0
  if "${installed_fn}"; then
    return 1
  fi
  return 0
}

preflight
if [[ "${PREFLIGHT_ONLY}" == "1" ]]; then
  info "Preflight passed. ROOT=${ROOT}; NEOAG_CONDA_BASE=${NEOAG_CONDA_BASE}; NEOAG_TOOLS_ROOT=${NEOAG_TOOLS_ROOT}"
  exit 0
fi

if should_install "${SKIP_LOHHLA:-0}" lohhla_installed; then
  run_step "Install LOHHLA" bash "${ROOT}/scripts/install_lohhla.sh"
elif [[ "${SKIP_LOHHLA:-0}" == "1" ]]; then
  info "Skip LOHHLA (SKIP_LOHHLA=1)"
else
  skip_step "LOHHLA"
fi

if should_install "${SKIP_FACETS:-0}" facets_installed; then
  run_step "Install FACETS" bash "${ROOT}/scripts/install_facets.sh"
elif [[ "${SKIP_FACETS:-0}" == "1" ]]; then
  info "Skip FACETS (SKIP_FACETS=1)"
else
  skip_step "FACETS"
fi

if should_install "${SKIP_ASCAT:-0}" ascat_installed; then
  run_step "Install ASCAT / PyClone-VI" bash "${ROOT}/scripts/install_ascat_pyclone.sh"
elif [[ "${SKIP_ASCAT:-0}" == "1" ]]; then
  info "Skip ASCAT (SKIP_ASCAT=1)"
else
  skip_step "ASCAT"
fi

if should_install "${SKIP_ARRIBA:-0}" arriba_installed; then
  run_step "Install Arriba / fusion env" bash "${ROOT}/scripts/install_fusion_tools.sh"
elif [[ "${SKIP_ARRIBA:-0}" == "1" ]]; then
  info "Skip Arriba (SKIP_ARRIBA=1)"
else
  skip_step "Arriba"
fi

if should_install "${SKIP_PRIME:-0}" prime_installed; then
  run_step "Install PRIME / MixMHCpred / BigMHC" bash "${ROOT}/scripts/install_immunogenicity_tools.sh"
elif [[ "${SKIP_PRIME:-0}" == "1" ]]; then
  info "Skip PRIME (SKIP_PRIME=1)"
else
  skip_step "PRIME"
fi

if should_install "${SKIP_OPTITYPE:-0}" optitype_installed; then
  run_step "Install OptiType" bash "${ROOT}/scripts/install_optitype.sh"
elif [[ "${SKIP_OPTITYPE:-0}" == "1" ]]; then
  info "Skip OptiType (SKIP_OPTITYPE=1)"
else
  skip_step "OptiType"
fi

source_tools_env_optional

if [[ -n "${NEOAG_REF_BUNDLE:-}" && -d "${NEOAG_REF_BUNDLE}" ]]; then
  run_step "Verify reference bundle" bash "${ROOT}/scripts/verify_reference_bundle.sh" "${NEOAG_REF_BUNDLE}"
fi

VERIFY_ARGS=()
[[ "${SMOKE}" == "1" ]] && VERIFY_ARGS+=(--smoke)
if [[ "${SKIP_VERIFY}" == "1" ]]; then
  info "Skip verification (NEOAG_SKIP_VERIFY=1 or --skip-verify)"
else
  run_step "Verify external tools" bash "${ROOT}/scripts/verify_external_tools.sh" "${VERIFY_ARGS[@]}"
fi

cat <<EOF_DONE

==> Deployment closed loop complete.

Environment:
  source conf/tools.env.sh

Next steps by tool:
  LOHHLA : configure POLYSOLVER_HOME + NOVOALIGN_LICENSE_FILE in conf/tools.env.local.sh, then bash scripts/run_lohhla_example.sh
  FACETS : stage bin/snp-pileup + SNP VCF, then PATIENT_ID=S1 TUMOR_BAM=... NORMAL_BAM=... bash scripts/run_facets_sample.sh
  ASCAT  : run FACETS pileup first, then PILEUP=... PATIENT_ID=S1 bash scripts/run_ascat_sample.sh
  Arriba : RNA BAM + references, then PATIENT_ID=S1 INPUT_BAM=... bash scripts/run_arriba_sample.sh
  PRIME  : neoag-v03 peptide-predict -i peptides.tsv -o results/sample --profile default
  OptiType: optitype run -i tumor_R1.fastq.gz -i tumor_R2.fastq.gz --dna -o results/optitype_sample --solver cbc --threads 8

Reference bundle:
  NEOAG_REF_BUNDLE=/path/to/neodata4git bash scripts/verify_reference_bundle.sh /path/to/neodata4git

Full release acceptance:
  NEOAG_REF_BUNDLE=/path/to/neodata4git bash scripts/verify_all_tools_and_refs.sh --smoke
  NEOAG_REF_BUNDLE=/path/to/neodata4git bash scripts/verify_all_tools_and_refs.sh --strict

Convert tool outputs into neoag evidence:
  neoag-v03 convert-lohhla -i <HLAlossPrediction_CI*> -o hla_loh.tsv
  neoag-v03 convert-facets --purity-input facets_purity.txt --purity-output purity.tsv
  neoag-v03 convert-ascat --summary-input ascat_summary.tsv --purity-output purity.tsv
EOF_DONE
