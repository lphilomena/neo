#!/usr/bin/env bash
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${NEOAG_CONDA_BASE:-${CONDA_PREFIX:-${HOME}/miniforge3}}"
ENV_ROOT="${NEOAG_ENV_ROOT:-${CONDA_BASE}/envs}"
SPLADDER_ENV="${NEOAG_SPLADDER_ENV:-${ENV_ROOT}/neoag-spladder}"
IMMUNOPEPPER_ENV="${NEOAG_IMMUNOPEPPER_ENV:-${ENV_ROOT}/neoag-immunopepper}"
BIN_DIR="${NEOAG_BIN_DIR:-${ROOT}/bin}"
TOOLS_ROOT="${NEOAG_SPLICE_TOOLS_ROOT:-${NEOAG_TOOLS_ROOT:-${ROOT}}/tools/splice}"
PVACTOOLS_BASE_ENV="${NEOAG_PVACTOOLS_BASE_ENV:-${CONDA_BASE}/envs/neoag-tools}"
IRFINDER_IMAGE="${NEOAG_IRFINDER_IMAGE:-cloxd/irfinder:2.0.1}"
failed=0

check() {
  local tool="$1" command="$2" expected="$3" output
  if output=$(bash -c "${command}" 2>&1) && grep -Eiq "${expected}" <<<"${output}"; then
    printf '%s\tREADY\t%s\n' "${tool}" "$(head -1 <<<"${output}")"
  else
    printf '%s\tBLOCKED\t%s\n' "${tool}" "$(tail -1 <<<"${output}")"
    failed=1
  fi
}

printf 'tool\tstatus\tdetail\n'
check SplAdder "'${SPLADDER_ENV}/bin/python' -c 'import spladder; print(spladder.__version__)'" '^3\.1\.1$'
check ImmunoPepper "'${IMMUNOPEPPER_ENV}/bin/python' -c 'import immunopepper; print(getattr(immunopepper, \"__version__\", \"unknown\"))'" '^2\.0\.0$'
check pVACbind "PYTHONPATH='${TOOLS_ROOT}/pvactools-7.1.1-overlay' '${PVACTOOLS_BASE_ENV}/bin/python' -c 'import importlib.metadata as m; print(m.version(\"pvactools\"))'" '^7\.1\.1$'
if docker image inspect "${IRFINDER_IMAGE}" >/dev/null 2>&1; then
  printf 'IRFinder-S\tREADY\t%s\n' "${IRFINDER_IMAGE}"
else
  printf 'IRFinder-S\tBLOCKED\timage missing: %s\n' "${IRFINDER_IMAGE}"
  failed=1
fi

for wrapper in spladder-neoag immunopepper-neoag pvacbind-neoag7 irfinder-s-neoag; do
  [[ -x "${BIN_DIR}/${wrapper}" ]] || { printf '%s\tBLOCKED\twrapper missing\n' "${wrapper}"; failed=1; }
done
exit "${failed}"
