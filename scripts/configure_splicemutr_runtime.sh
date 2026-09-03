#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
HOME_DIR="${NEOAG_SPLICEMUTR_HOME:-}"
ENV_PREFIX="${NEOAG_SPLICEMUTR_ENV_PREFIX:-}"
LEAFCUTTER_PREFIX="${NEOAG_SPLICEMUTR_LEAFCUTTER_ENV_PREFIX:-}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --home) HOME_DIR="$2"; shift 2;;
    --env-prefix) ENV_PREFIX="$2"; shift 2;;
    --leafcutter-env-prefix) LEAFCUTTER_PREFIX="$2"; shift 2;;
    *) echo "ERROR: unknown option $1" >&2; exit 2;;
  esac
done
[[ -d "$HOME_DIR/Rscripts" && -d "$HOME_DIR/python_scripts" ]] || { echo "ERROR: invalid SpliceMutr home: $HOME_DIR" >&2; exit 2; }
[[ -x "$ENV_PREFIX/bin/Rscript" && -x "$ENV_PREFIX/bin/snakemake" ]] || { echo "ERROR: incomplete SpliceMutr environment: $ENV_PREFIX" >&2; exit 2; }
[[ -x "$LEAFCUTTER_PREFIX/bin/python" ]] || { echo "ERROR: incomplete LeafCutter environment: $LEAFCUTTER_PREFIX" >&2; exit 2; }

mkdir -p "$ROOT/bin"
cat > "$ROOT/bin/splicemutr-neoag" <<EOF
#!/usr/bin/env bash
set -euo pipefail
HOME_DIR=$(printf '%q' "$HOME_DIR")
ENV_PREFIX=$(printf '%q' "$ENV_PREFIX")
LEAFCUTTER_PREFIX=$(printf '%q' "$LEAFCUTTER_PREFIX")
run_isolated() (
  unset VIRTUAL_ENV PYTHONHOME PYTHONPATH PYTHONSTARTUP
  unset CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_PROMPT_MODIFIER CONDA_SHLVL
  unset CONDA_PYTHON_EXE CONDA_EXE CONDA_ROOT
  unset R_LIBS R_LIBS_USER R_LIBS_SITE R_PROFILE R_PROFILE_USER
  unset BASH_ENV ENV
  export PATH="\${ENV_PREFIX}/bin:/usr/bin:/bin"
  export LD_LIBRARY_PATH="\${ENV_PREFIX}/lib"
  export CONDA_PREFIX="\${ENV_PREFIX}"
  exec "\$@"
)
case "\${1:-doctor}" in
  doctor)
    test -d "\$HOME_DIR/Rscripts"
    run_isolated "\$ENV_PREFIX/bin/python" -c 'import Bio, numpy, pandas, sklearn; print("SpliceMutr Python runtime OK")'
    run_isolated "\$ENV_PREFIX/bin/Rscript" -e 'suppressPackageStartupMessages({library(BSgenome); library(GenomicFeatures); library(optparse)}); cat("SpliceMutr R runtime OK\\n")'
    run_isolated "\$ENV_PREFIX/bin/snakemake" --version
    env -u VIRTUAL_ENV -u PYTHONHOME -u PYTHONPATH PATH="\$LEAFCUTTER_PREFIX/bin:/usr/bin:/bin" LD_LIBRARY_PATH="\$LEAFCUTTER_PREFIX/lib" "\$LEAFCUTTER_PREFIX/bin/python" -c 'import numpy, scipy; print("LeafCutter runtime OK")'
    ;;
  workflow)
    shift; workflow="\${1:?workflow path required}"; shift
    run_isolated "\$ENV_PREFIX/bin/snakemake" -s "\$workflow" "\$@"
    ;;
  r)
    shift; script="\${1:?R script required}"; shift
    [[ "\$script" = /* ]] || script="\$HOME_DIR/Rscripts/\$script"
    run_isolated "\$ENV_PREFIX/bin/Rscript" "\$script" "\$@"
    ;;
  python)
    shift; script="\${1:?Python script required}"; shift
    [[ "\$script" = /* ]] || script="\$HOME_DIR/python_scripts/\$script"
    run_isolated "\$ENV_PREFIX/bin/python" "\$script" "\$@"
    ;;
  leafcutter)
    shift
    env -u VIRTUAL_ENV -u PYTHONHOME -u PYTHONPATH PATH="\$LEAFCUTTER_PREFIX/bin:/usr/bin:/bin" LD_LIBRARY_PATH="\$LEAFCUTTER_PREFIX/lib" exec "\$LEAFCUTTER_PREFIX/bin/python" "\$LEAFCUTTER_PREFIX/bin/leafcutter_cluster_regtools.py" "\$@"
    ;;
  *) echo "ERROR: unknown command: \$1" >&2; exit 2;;
esac
EOF
chmod +x "$ROOT/bin/splicemutr-neoag"
"$ROOT/bin/splicemutr-neoag" doctor
echo "Configured $ROOT/bin/splicemutr-neoag"
