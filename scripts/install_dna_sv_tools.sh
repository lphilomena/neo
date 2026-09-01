#!/usr/bin/env bash
# Install the required-all production DNA-SV caller group in an isolated env.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONDA_BASE="${NEOAG_CONDA_BASE:-}"
ENV_PREFIX="${NEOAG_SV_ENV_PREFIX:-}"
GRIDSS_ENV_PREFIX="${NEOAG_GRIDSS_ENV_PREFIX:-}"
YML="${NEOAG_SV_ENV_YML:-${ROOT}/conda/env.neoag-sv.yml}"
GRIDSS_YML="${NEOAG_GRIDSS_ENV_YML:-${ROOT}/conda/env.neoag-gridss.yml}"
BIOC_CACHE_HELPER="${ROOT}/.agents/skills/neoag-remote-deploy/scripts/with_bioc_data_cache.sh"
CACHE_ROOT="${NEOAG_INSTALL_CACHE_ROOT:-${NEOAG_TOOLS_ROOT:-${ROOT}}/install_cache}"

if [[ -z "${CONDA_BASE}" ]]; then
  CONDA_BIN="$(command -v conda || true)"
  [[ -n "${CONDA_BIN}" ]] || { echo "ERROR: conda not found" >&2; exit 2; }
  CONDA_BASE="$("${CONDA_BIN}" info --base)"
fi
CONDA_BIN="${CONDA_BASE}/bin/conda"
[[ -x "${CONDA_BIN}" ]] || { echo "ERROR: conda not found at ${CONDA_BIN}" >&2; exit 2; }
[[ -f "${YML}" ]] || { echo "ERROR: DNA-SV environment file missing: ${YML}" >&2; exit 2; }
[[ -f "${GRIDSS_YML}" ]] || { echo "ERROR: GRIDSS environment file missing: ${GRIDSS_YML}" >&2; exit 2; }
ENV_PREFIX="${ENV_PREFIX:-${CONDA_BASE}/envs/neoag-sv}"
GRIDSS_ENV_PREFIX="${GRIDSS_ENV_PREFIX:-${CONDA_BASE}/envs/neoag-gridss}"

conda_safe() {
  set +u
  "${CONDA_BIN}" "$@"
  local rc=$?
  set -u
  return "${rc}"
}

GRIDSS_CONDA_RUNNER=("${CONDA_BIN}")
if [[ -x "${BIOC_CACHE_HELPER}" ]]; then
  GRIDSS_CONDA_RUNNER=(
    "${BIOC_CACHE_HELPER}"
    --conda-base "${CONDA_BASE}"
    --cache-root "${CACHE_ROOT}"
    --package-key genomeinfodbdata-1.2.13
    -- "${CONDA_BIN}"
  )
else
  echo "ERROR: Bioconductor data cache helper missing: ${BIOC_CACHE_HELPER}" >&2
  exit 2
fi

if [[ -x "${ENV_PREFIX}/bin/svaba" ]]; then
  echo "==> Updating existing DNA-SV environment: ${ENV_PREFIX}"
  conda_safe env update -p "${ENV_PREFIX}" -f "${YML}" --prune
else
  echo "==> Creating DNA-SV environment: ${ENV_PREFIX}"
  conda_safe env create -p "${ENV_PREFIX}" -f "${YML}"
fi

if [[ -x "${GRIDSS_ENV_PREFIX}/bin/gridss" || -x "${GRIDSS_ENV_PREFIX}/bin/gridss.sh" ]]; then
  echo "==> Updating existing GRIDSS environment: ${GRIDSS_ENV_PREFIX}"
  "${GRIDSS_CONDA_RUNNER[@]}" env update -p "${GRIDSS_ENV_PREFIX}" -f "${GRIDSS_YML}" --prune
else
  echo "==> Creating isolated GRIDSS environment: ${GRIDSS_ENV_PREFIX}"
  "${GRIDSS_CONDA_RUNNER[@]}" env create -p "${GRIDSS_ENV_PREFIX}" -f "${GRIDSS_YML}"
fi

MANTA_BIN="${ENV_PREFIX}/bin/configManta.py"
SVABA_BIN="${ENV_PREFIX}/bin/svaba"
GRIDSS_BIN=""
for candidate in "${GRIDSS_ENV_PREFIX}/bin/gridss" "${GRIDSS_ENV_PREFIX}/bin/gridss.sh"; do
  if [[ -x "${candidate}" ]]; then GRIDSS_BIN="${candidate}"; break; fi
done
[[ -x "${MANTA_BIN}" ]] || { echo "ERROR: Manta configManta.py missing" >&2; exit 3; }
[[ -x "${SVABA_BIN}" ]] || { echo "ERROR: SvABA executable missing" >&2; exit 3; }
[[ -n "${GRIDSS_BIN}" ]] || { echo "ERROR: GRIDSS executable missing" >&2; exit 3; }

mkdir -p "${ROOT}/bin" "${ROOT}/conf"
write_wrapper() {
  local target="$1" real="$2" env_bin="$3"
  printf '#!/usr/bin/env bash\nexport PATH=%q:"$PATH"\nexec %q "$@"\n' "${env_bin}" "${real}" > "${ROOT}/bin/${target}"
  chmod +x "${ROOT}/bin/${target}"
}
write_wrapper configManta.py "${MANTA_BIN}" "${ENV_PREFIX}/bin"
write_wrapper svaba "${SVABA_BIN}" "${ENV_PREFIX}/bin"
write_wrapper gridss "${GRIDSS_BIN}" "${GRIDSS_ENV_PREFIX}/bin"

LOCAL_ENV="${ROOT}/conf/tools.env.local.sh"
touch "${LOCAL_ENV}"
python3 - "${LOCAL_ENV}" "${ENV_PREFIX}" "${GRIDSS_ENV_PREFIX}" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
prefix = sys.argv[2]
gridss_prefix = sys.argv[3]
start = "# >>> OPEN_NEO_DNA_SV >>>"
end = "# <<< OPEN_NEO_DNA_SV <<<"
text = path.read_text(encoding="utf-8")
block = (
    f"{start}\n"
    f"export NEOAG_SV_ENV_PREFIX={prefix!r}\n"
    f"export NEOAG_GRIDSS_ENV_PREFIX={gridss_prefix!r}\n"
    f"export PATH={str(Path(prefix) / 'bin')!r}:\"$PATH\"\n"
    f"{end}\n"
)
if start in text and end in text:
    before, rest = text.split(start, 1)
    _, after = rest.split(end, 1)
    text = before.rstrip() + "\n" + block + after.lstrip("\n")
else:
    text = text.rstrip() + "\n\n" + block
path.write_text(text, encoding="utf-8")
PY

"${ROOT}/bin/configManta.py" --help >/dev/null
"${ROOT}/bin/svaba" --help >/dev/null 2>&1 || "${ROOT}/bin/svaba" run --help >/dev/null 2>&1
"${ROOT}/bin/gridss" --help >/dev/null 2>&1 || true
echo "==> DNA-SV capability group READY: manta, svaba, gridss"
