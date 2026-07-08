#!/usr/bin/env bash
# Configure an existing Polysolver distribution and verify dependencies.
#
# Polysolver is licensed/external and is not downloaded by this project.
# Stage the Polysolver tree on the target machine first, then run:
#
#   POLYSOLVER_HOME=/path/to/polysolver \
#   NOVOALIGN_LICENSE_FILE=/path/to/novoalign.lic \
#   bash scripts/install_polysolver.sh
#
# The Polysolver tree is expected to contain scripts/, binaries/, and data/.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TOOLS_ENV="${ROOT}/conf/tools.env.sh"

if [[ -f "${TOOLS_ENV}" ]]; then
  # shellcheck source=/dev/null
  source "${TOOLS_ENV}"
fi

find_conda_base() {
  if [[ -n "${NEOAG_CONDA_BASE:-}" && -d "${NEOAG_CONDA_BASE}" ]]; then
    echo "${NEOAG_CONDA_BASE}"
    return 0
  fi
  if command -v conda >/dev/null 2>&1; then
    conda info --base
    return 0
  fi
  for base in "${HOME}/miniforge3" "${HOME}/mambaforge" "${HOME}/miniconda3" "${HOME}/anaconda3" /opt/conda; do
    if [[ -x "${base}/bin/conda" ]]; then
      echo "${base}"
      return 0
    fi
  done
  return 1
}

PSHOME="${POLYSOLVER_HOME:-}"
if [[ -z "${PSHOME}" ]]; then
  cat >&2 <<ERR
ERROR: POLYSOLVER_HOME is required.
Set it to an existing Polysolver distribution, for example:
  POLYSOLVER_HOME=/path/to/polysolver NOVOALIGN_LICENSE_FILE=/path/to/novoalign.lic bash scripts/install_polysolver.sh
ERR
  exit 2
fi

NOVO_LICENSE_SRC="${NOVOALIGN_LICENSE_FILE:-}"
CONDA_ENV="${NEOAG_CONDA_ENV:-${NEOAG_TOOLS_ENV:-neoag-tools}}"
CONDA_BASE="$(find_conda_base || true)"

[[ -d "${PSHOME}/scripts" && -d "${PSHOME}/binaries" && -d "${PSHOME}/data" ]] || {
  echo "ERROR: Polysolver tree missing at ${PSHOME}" >&2
  echo "  Expected scripts/, binaries/, data/ from the Polysolver distribution." >&2
  exit 1
}

echo "==> Polysolver home: ${PSHOME}"
mkdir -p "${PSHOME}/sachet" "${ROOT}/work/polysolver" "${ROOT}/bin" "${ROOT}/conf"

if [[ -n "${NOVO_LICENSE_SRC}" ]]; then
  if [[ -f "${NOVO_LICENSE_SRC}" ]]; then
    cp -f "${NOVO_LICENSE_SRC}" "${PSHOME}/binaries/novoalign.lic"
    echo "==> Novoalign license: ${PSHOME}/binaries/novoalign.lic"
  else
    echo "WARN: NOVOALIGN_LICENSE_FILE does not exist: ${NOVO_LICENSE_SRC}" >&2
  fi
else
  echo "WARN: NOVOALIGN_LICENSE_FILE unset; Polysolver alignment may fail without a Novoalign license." >&2
fi

chmod +x \
  "${PSHOME}/binaries/novoalign" \
  "${PSHOME}/binaries/samtools" \
  "${PSHOME}/binaries/java" \
  "${PSHOME}/scripts/novoindex" \
  "${PSHOME}/scripts/shell_call_hla_type" \
  "${PSHOME}/scripts/shell_call_hla_type_test" 2>/dev/null || true

cat > "${PSHOME}/scripts/config.local.bash" <<EOF_CFG
#!/usr/bin/env bash
PSHOME="${PSHOME}"
SAMTOOLS_DIR="\${PSHOME}/binaries"
JAVA_DIR="\${PSHOME}/binaries"
NOVOALIGN_DIR="\${PSHOME}/binaries"
GATK_DIR="\${PSHOME}/binaries"
MUTECT_DIR="\${PSHOME}/binaries"
STRELKA_DIR="\${PSHOME}/binaries"
TMP_DIR="\${PSHOME}/sachet"
NUM_THREADS="\${NUM_THREADS:-8}"
export PSHOME SAMTOOLS_DIR JAVA_DIR NOVOALIGN_DIR GATK_DIR MUTECT_DIR STRELKA_DIR TMP_DIR NUM_THREADS
EOF_CFG
chmod +x "${PSHOME}/scripts/config.local.bash"

if [[ -n "${CONDA_BASE}" && -x "${CONDA_BASE}/bin/mamba" ]]; then
  echo "==> Installing Perl deps in conda env: ${CONDA_ENV}"
  "${CONDA_BASE}/bin/mamba" install -y -n "${CONDA_ENV}" -c bioconda -c conda-forge \
    perl-list-moreutils perl-bioperl perl-parallel-forkmanager \
    >/tmp/install_polysolver_mamba.log 2>&1 || {
      echo "WARN: mamba install had issues; see /tmp/install_polysolver_mamba.log" >&2
    }
elif [[ -n "${CONDA_BASE}" && -x "${CONDA_BASE}/bin/conda" ]]; then
  echo "==> Installing Perl deps in conda env: ${CONDA_ENV}"
  "${CONDA_BASE}/bin/conda" install -y -n "${CONDA_ENV}" -c bioconda -c conda-forge \
    perl-list-moreutils perl-bioperl perl-parallel-forkmanager \
    >/tmp/install_polysolver_conda.log 2>&1 || {
      echo "WARN: conda install had issues; see /tmp/install_polysolver_conda.log" >&2
    }
else
  echo "WARN: conda not found; skipping Perl dependency install." >&2
fi

if [[ -n "${CONDA_BASE}" && -d "${CONDA_BASE}/envs/${CONDA_ENV}/bin" ]]; then
  export PATH="${CONDA_BASE}/envs/${CONDA_ENV}/bin:${PATH}"
fi
export PATH="${PSHOME}/binaries:${PSHOME}/scripts:${PATH}"

echo "==> Verifying binaries"
"${PSHOME}/binaries/novoalign" 2>&1 | head -1 || true
"${PSHOME}/binaries/samtools" 2>&1 | head -1 || true
java -version 2>&1 | head -1 || "${PSHOME}/binaries/java" -version 2>&1 | head -1 || true

echo "==> Verifying Perl modules"
perl -MList::MoreUtils -e 'print "List::MoreUtils OK\n"' || echo "WARN: List::MoreUtils missing" >&2
perl -MBio::SeqIO -e 'print "BioPerl OK\n"' || echo "WARN: BioPerl missing" >&2
perl -MParallel::ForkManager -e 'print "Parallel::ForkManager OK\n"' || echo "WARN: Parallel::ForkManager missing" >&2

cat > "${ROOT}/bin/run-polysolver" <<'WRAP'
#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/conf/tools.env.sh"
PSHOME="${POLYSOLVER_HOME:?Set POLYSOLVER_HOME in conf/tools.env.local.sh or the shell}"
# shellcheck source=/dev/null
source "${PSHOME}/scripts/config.local.bash"
if [[ -n "${NEOAG_CONDA_BASE:-}" && -d "${NEOAG_CONDA_BASE}/envs/${NEOAG_CONDA_ENV:-neoag-tools}/bin" ]]; then
  export PATH="${NEOAG_CONDA_BASE}/envs/${NEOAG_CONDA_ENV:-neoag-tools}/bin:${PATH}"
fi
export PATH="${PSHOME}/binaries:${PSHOME}/scripts:${PATH}"
exec bash "${PSHOME}/scripts/shell_call_hla_type" "$@"
WRAP
chmod +x "${ROOT}/bin/run-polysolver"

if [[ ! -f "${TOOLS_ENV}" ]]; then
  cat > "${TOOLS_ENV}" <<EOF_ENV
export NEOAG_PROJECT_ROOT="${ROOT}"
export NEOAG_TOOLS_ROOT="${ROOT}"
export NEOAG_CONDA_ENV="neoag-tools"
EOF_ENV
fi
if ! grep -q 'Polysolver - configured via scripts/install_polysolver.sh' "${TOOLS_ENV}"; then
  cat >> "${TOOLS_ENV}" <<EOF_ENV

# Polysolver - configured via scripts/install_polysolver.sh
export POLYSOLVER_HOME="\${POLYSOLVER_HOME:-${PSHOME}}"
export NOVOALIGN_LICENSE_FILE="\${NOVOALIGN_LICENSE_FILE:-${NOVO_LICENSE_SRC}}"
export PATH="${ROOT}/bin:\${POLYSOLVER_HOME}/binaries:\${POLYSOLVER_HOME}/scripts:\${PATH}"
EOF_ENV
fi

cat <<EOF_DONE

==> Polysolver ready.
    POLYSOLVER_HOME=${PSHOME}
    source ${PSHOME}/scripts/config.local.bash
    ${ROOT}/bin/run-polysolver <bam> <race> <includeFreq> <build> <format> <insertCalc> <outDir>

Example:
    POLYSOLVER_HOME=${PSHOME} ${ROOT}/bin/run-polysolver tumor.bam Unknown 1 hg19 STDFQ 1 results/ps_out
EOF_DONE
