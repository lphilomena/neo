#!/usr/bin/env bash
# Configure local Polysolver install and verify dependencies.
#
# Polysolver tree (external to this repo):
#   /home/na/project/neoantigen/software/polysolver
#
# Usage:
#   bash scripts/install_polysolver.sh
#   source /home/na/project/neoantigen/software/polysolver/scripts/config.local.bash
#   bash /home/na/project/neoantigen/software/polysolver/scripts/shell_call_hla_type ...
#
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/conf/tools.env.sh"

PSHOME="${POLYSOLVER_HOME:-/home/na/project/neoantigen/software/polysolver}"
NOVO_LICENSE_SRC="${NOVOALIGN_LICENSE_FILE:-/home/na/project/neoantigen/wbscript/novoalign.lic}"
CONDA_ENV="${NEOAG_TOOLS_ENV:-neoag-tools}"
CONDA_BASE="${NEOAG_CONDA_BASE:-${HOME}/miniforge3}"

[[ -d "${PSHOME}/scripts" && -d "${PSHOME}/binaries" && -d "${PSHOME}/data" ]] || {
  echo "ERROR: Polysolver tree missing at ${PSHOME}" >&2
  echo "  Expected scripts/, binaries/, data/ from the Polysolver distribution." >&2
  exit 1
}

echo "==> Polysolver home: ${PSHOME}"

mkdir -p "${PSHOME}/sachet" "${ROOT}/work/polysolver"

# Novoalign license (required for alignment step)
if [[ -f "${NOVO_LICENSE_SRC}" ]]; then
  cp -f "${NOVO_LICENSE_SRC}" "${PSHOME}/binaries/novoalign.lic"
  echo "==> Novoalign license: ${PSHOME}/binaries/novoalign.lic"
else
  echo "WARN: Novoalign license not found at ${NOVO_LICENSE_SRC}" >&2
fi

chmod +x "${PSHOME}/binaries/novoalign" \
  "${PSHOME}/binaries/samtools" \
  "${PSHOME}/binaries/java" \
  "${PSHOME}/scripts/novoindex" \
  "${PSHOME}/scripts/shell_call_hla_type" \
  "${PSHOME}/scripts/shell_call_hla_type_test" 2>/dev/null || true

# Write local config used by wrapper scripts
cat > "${PSHOME}/scripts/config.local.bash" <<EOF
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
EOF

# Perl modules for Polysolver scripts
if [[ -x "${CONDA_BASE}/bin/mamba" ]]; then
  echo "==> installing Perl deps in conda env: ${CONDA_ENV}"
  "${CONDA_BASE}/bin/mamba" install -y -n "${CONDA_ENV}" -c bioconda -c conda-forge \
    perl-list-moreutils perl-bioperl perl-parallel-forkmanager \
    >/tmp/install_polysolver_mamba.log 2>&1 || {
      echo "WARN: mamba install had issues; see /tmp/install_polysolver_mamba.log" >&2
    }
fi

export PATH="${CONDA_BASE}/envs/${CONDA_ENV}/bin:${PSHOME}/binaries:${PSHOME}/scripts:${PATH}"

echo "==> verifying binaries"
"${PSHOME}/binaries/novoalign" 2>&1 | head -1 || true
"${PSHOME}/binaries/samtools" 2>&1 | head -1 || true
java -version 2>&1 | head -1 || "${PSHOME}/binaries/java" -version 2>&1 | head -1

echo "==> verifying Perl modules"
perl -MList::MoreUtils -e 'print "List::MoreUtils OK\n"'
perl -MBio::SeqIO -e 'print "BioPerl OK\n"'
perl -MParallel::ForkManager -e 'print "Parallel::ForkManager OK\n"'

cat > "${ROOT}/bin/run-polysolver" <<'WRAP'
#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck source=/dev/null
source "${ROOT}/conf/tools.env.sh"
PSHOME="${POLYSOLVER_HOME:-/home/na/project/neoantigen/software/polysolver}"
# shellcheck source=/dev/null
source "${PSHOME}/scripts/config.local.bash"
export PATH="${NEOAG_CONDA_BASE}/envs/${NEOAG_TOOLS_ENV:-neoag-tools}/bin:${PSHOME}/binaries:${PSHOME}/scripts:${PATH}"
exec bash "${PSHOME}/scripts/shell_call_hla_type" "$@"
WRAP
chmod +x "${ROOT}/bin/run-polysolver"

echo ""
echo "==> Polysolver ready."
echo "    POLYSOLVER_HOME=${PSHOME}"
echo "    source ${PSHOME}/scripts/config.local.bash"
echo "    ${ROOT}/bin/run-polysolver <bam> <race> <includeFreq> <build> <format> <insertCalc> <outDir>"
echo ""
echo "Example (hg38 BAM with chr prefix — use build hg19 for Polysolver coordinates):"
echo "  ${ROOT}/bin/run-polysolver tumor.bam Unknown 1 hg19 STDFQ 1 results/ps_out"
