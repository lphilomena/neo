#!/usr/bin/env bash
# Install OptiType in a portable conda environment.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_NAME="${NEOAG_OPTITYPE_ENV:-neoag-optitype}"
CONDA_BASE="${NEOAG_CONDA_BASE:-}"
FORCE="${NEOAG_FORCE_ENV_UPDATE:-0}"

usage() {
  cat <<USAGE
Usage: bash scripts/install_optitype.sh [--force]

Environment:
  NEOAG_CONDA_BASE=/path/to/miniforge3   Conda root; auto-detected if conda is on PATH
  NEOAG_OPTITYPE_ENV=neoag-optitype      Conda env name
  NEOAG_FORCE_ENV_UPDATE=1               Reinstall/update packages even if env exists

Installs:
  optitype, razers3, glpk, coincbc

Verify:
  source conf/tools.env.sh
  optitype check-deps
USAGE
}

for arg in "$@"; do
  case "$arg" in
    --force) FORCE=1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $arg" >&2; usage >&2; exit 2 ;;
  esac
done

find_conda() {
  if command -v conda >/dev/null 2>&1; then command -v conda; return 0; fi
  local c
  for c in \
    "${CONDA_BASE:+${CONDA_BASE}/bin/conda}" \
    "${HOME}/miniforge3/bin/conda" \
    "${HOME}/mambaforge/bin/conda" \
    "${HOME}/miniconda3/bin/conda" \
    "${HOME}/anaconda3/bin/conda" \
    "/opt/conda/bin/conda"; do
    [[ -n "$c" && -x "$c" ]] && echo "$c" && return 0
  done
  return 1
}

CONDA_BIN="$(find_conda)" || { echo "ERROR: conda not found; set NEOAG_CONDA_BASE" >&2; exit 1; }
CONDA_BASE="${NEOAG_CONDA_BASE:-$("$CONDA_BIN" info --base)}"
export PATH="${CONDA_BASE}/bin:${PATH}"
# shellcheck source=/dev/null
source "${CONDA_BASE}/etc/profile.d/conda.sh"

MAMBA_BIN="${CONDA_BASE}/bin/mamba"
CREATE_CMD=(conda create -y -n "$ENV_NAME" -c conda-forge -c bioconda optitype glpk coincbc razers3)
UPDATE_CMD=(conda install -y -n "$ENV_NAME" -c conda-forge -c bioconda optitype glpk coincbc razers3)
if [[ -x "$MAMBA_BIN" ]]; then
  CREATE_CMD=("$MAMBA_BIN" create -y -n "$ENV_NAME" -c conda-forge -c bioconda optitype glpk coincbc razers3)
  UPDATE_CMD=("$MAMBA_BIN" install -y -n "$ENV_NAME" -c conda-forge -c bioconda optitype glpk coincbc razers3)
fi

if conda env list | awk "{print \$1}" | grep -qx "$ENV_NAME"; then
  if [[ "$FORCE" == "1" ]]; then
    echo "==> Updating existing env: $ENV_NAME"
    "${UPDATE_CMD[@]}"
  else
    echo "==> Env exists: $ENV_NAME"
  fi
else
  echo "==> Creating env: $ENV_NAME"
  "${CREATE_CMD[@]}"
fi

ENV_PREFIX="$(conda env list | awk -v n="$ENV_NAME" "\$1==n {print \$NF}")"
[[ -n "$ENV_PREFIX" && -x "$ENV_PREFIX/bin/optitype" ]] || { echo "ERROR: optitype missing after install" >&2; exit 1; }
# Some wrappers/tools look for this historical capitalization.
if [[ -x "$ENV_PREFIX/bin/razers3" && ! -e "$ENV_PREFIX/bin/RazerS3" ]]; then
  ln -s "$ENV_PREFIX/bin/razers3" "$ENV_PREFIX/bin/RazerS3"
fi

TOOLS_ENV="${ROOT}/conf/tools.env.sh"
if [[ -f "$TOOLS_ENV" ]] && ! grep -q "NEOAG_OPTITYPE_ENV" "$TOOLS_ENV"; then
  cat >> "$TOOLS_ENV" <<EOF

# OptiType (HLA-I typing from DNA/RNA FASTQ or BAM)
export NEOAG_OPTITYPE_ENV="\${NEOAG_OPTITYPE_ENV:-$ENV_NAME}"
if [[ -z "\${OPTITYPE_ENV:-}" ]]; then
  if [[ -d "\${NEOAG_CONDA_BASE}/envs/\${NEOAG_OPTITYPE_ENV}" ]]; then
    export OPTITYPE_ENV="\${NEOAG_CONDA_BASE}/envs/\${NEOAG_OPTITYPE_ENV}"
  fi
fi
if [[ -n "\${OPTITYPE_ENV:-}" && -x "\${OPTITYPE_ENV}/bin/optitype" ]]; then
  export OPTITYPE_BIN="\${OPTITYPE_ENV}/bin/optitype"
  export OPTITYPE_REFERENCE="\${OPTITYPE_ENV}/share/optitype/data"
  export PATH="\${OPTITYPE_ENV}/bin:\${PATH}"
fi
EOF
fi

export PATH="$ENV_PREFIX/bin:$PATH"
echo "==> OptiType version"
optitype --version
echo "==> OptiType dependency check"
optitype_status=0
optitype check-deps || optitype_status=$?
if [[ "$optitype_status" != "0" ]]; then
  echo "ERROR: optitype check-deps failed" >&2
  exit "$optitype_status"
fi

echo "==> Installed OptiType env: $ENV_PREFIX"
echo "==> Reference data: $ENV_PREFIX/share/optitype/data"
