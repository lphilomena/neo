#!/usr/bin/env bash
# Reproducible installer skeleton for v0.5.1 external tools.
# Network access and large reference assets are intentionally operator-controlled.
set -euo pipefail
usage() {
 cat <<'USAGE'
Usage: bash scripts/install_splice_v051_tools.sh --env-dir ENV [options]
Options:
  --mopepgen-version VERSION     Default stable: 1.4.6
  --easyquant-version VERSION    Default: 0.6.0
  --pvactools-version VERSION    Required to install pVACtools; site-controlled
  --install-splice2neo           Install R package from GitHub ref
  --splice2neo-ref REF           Default: v0.6.14
  --rscript PATH                Rscript used for splice2neo (R >= 4.4)
  --install-k4neo                Install k4neo source after license review
  --k4neo-ref REF                Required with --install-k4neo
  --k4neo-license-accepted       Explicit acknowledgement
  --k4neo-env-dir DIR            Dedicated Conda env; default: ENV-k4neo
  --k4neo-source-dir DIR         Optional complete recursive source checkout
  --conda PATH                   Conda executable; default: conda
  --python PYTHON                Default: python3

The script records exact versions. splice2neo is installed from its tagged
source after removing its unused liftOver workflow import; coordinate lifting
continues to use rtracklayer. The script does not download GTEx/TCGA/k4neo indices,
IEDB assets, licensed NetMHC software, or patient data.
USAGE
}
ENV_DIR=""; PY="python3"; MO="1.4.6"; EQ="0.6.0"; PV=""; S2N=0; S2N_REF="v0.6.14"; RSCRIPT="Rscript"; K4=0; K4_REF=""; K4_ACCEPT=0; K4_ENV=""; K4_SOURCE=""; CONDA="conda"
while [[ $# -gt 0 ]]; do
 case "$1" in
  --env-dir) ENV_DIR="$2"; shift 2;; --python) PY="$2"; shift 2;;
  --mopepgen-version) MO="$2"; shift 2;; --easyquant-version) EQ="$2"; shift 2;;
  --pvactools-version) PV="$2"; shift 2;; --install-splice2neo) S2N=1; shift;;
  --splice2neo-ref) S2N_REF="$2"; shift 2;; --rscript) RSCRIPT="$2"; shift 2;;
  --install-k4neo) K4=1; shift;;
  --k4neo-ref) K4_REF="$2"; shift 2;; --k4neo-license-accepted) K4_ACCEPT=1; shift;;
  --k4neo-env-dir) K4_ENV="$2"; shift 2;; --k4neo-source-dir) K4_SOURCE="$2"; shift 2;;
  --conda) CONDA="$2"; shift 2;;
  -h|--help) usage; exit 0;; *) echo "ERROR: unknown argument $1" >&2; exit 2;;
 esac
done
[[ -n "$ENV_DIR" ]] || { echo "ERROR: --env-dir required" >&2; exit 2; }
if [[ "$K4" == 1 ]]; then [[ "$K4_ACCEPT" == 1 && -n "$K4_REF" ]] || { echo "ERROR: k4neo requires --k4neo-ref and --k4neo-license-accepted" >&2; exit 2; }; fi
[[ -n "$K4_ENV" ]] || K4_ENV="${ENV_DIR}-k4neo"
"$PY" -c 'import sys; assert (3, 10) <= sys.version_info[:2] <= (3, 12), "use Python 3.10-3.12 for pinned moPepGen/EasyQuant wheels"' || exit 3
"$PY" -m venv "$ENV_DIR"
# shellcheck disable=SC1091
source "$ENV_DIR/bin/activate"
python -m pip install --upgrade pip wheel "setuptools<81"
# moPepGen/matplotlib wheels currently require the NumPy 1.x ABI.
python -m pip install "numpy<2" "mopepgen==$MO" "bp-quant==$EQ"
[[ -n "$PV" ]] && python -m pip install "pvactools==$PV"
if [[ "$S2N" == 1 ]]; then
 [[ -x "$RSCRIPT" ]] || command -v "$RSCRIPT" >/dev/null || { echo "ERROR: Rscript required for splice2neo" >&2; exit 3; }
 "$RSCRIPT" -e "if (getRversion() < '4.4.0') stop('splice2neo requires R >= 4.4')"
 command -v git >/dev/null || { echo "ERROR: git required for splice2neo" >&2; exit 3; }
 s2n_src="$(mktemp -d)"
 trap 'rm -rf "${s2n_src:-}"' EXIT
 git clone --quiet --depth 1 --branch "$S2N_REF" https://github.com/TRON-Bioinformatics/splice2neo.git "$s2n_src"
 # liftOver is a Bioconductor workflow imported only by an upstream no-op
 # placeholder. rtracklayer supplies the actual liftOver implementation.
 sed -i '/^[[:space:]]*liftOver[,]\{0,1\}[[:space:]]*$/d' "$s2n_src/DESCRIPTION"
 sed -i '/^[[:space:]]*liftOver::cur[[:space:]]*$/d' "$s2n_src/R/globals.R"
 "$RSCRIPT" -e "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='https://cloud.r-project.org'); remotes::install_local('$s2n_src', dependencies=TRUE, upgrade='never')"
 rm -rf "$s2n_src"
 trap - EXIT
fi
if [[ "$K4" == 1 ]]; then
 command -v "$CONDA" >/dev/null || [[ -x "$CONDA" ]] || { echo "ERROR: Conda is required for k4neo" >&2; exit 3; }
 if [[ -n "$K4_SOURCE" ]]; then
  [[ -f "$K4_SOURCE/k4neo.yaml" && -f "$K4_SOURCE/pyproject.toml" ]] || { echo "ERROR: incomplete --k4neo-source-dir" >&2; exit 3; }
  k4_src="$K4_SOURCE"
 else
  command -v git >/dev/null || { echo "ERROR: git is required for k4neo" >&2; exit 3; }
  k4_src="$(mktemp -d)"
  git -C "$k4_src" init -q
  git -C "$k4_src" remote add origin https://github.com/TRON-Bioinformatics/k4neo.git
  git -C "$k4_src" config extensions.partialClone origin
  git -C "$k4_src" config remote.origin.promisor true
  git -C "$k4_src" config remote.origin.partialclonefilter blob:none
  git -C "$k4_src" -c protocol.version=2 fetch --depth=1 --filter=blob:none origin "$K4_REF"
  git -C "$k4_src" checkout -q --detach FETCH_HEAD
  git -C "$k4_src" submodule update --init --recursive --depth=1
 fi
 k4_commit="$(git -C "$k4_src" rev-parse HEAD 2>/dev/null || cat "$k4_src/.k4neo_source_commit" 2>/dev/null || printf '%s' "$K4_REF")"
 k4_submodules="$(git -C "$k4_src" submodule status --recursive 2>/dev/null | tr '\n' ';' || true)"
 "$CONDA" env create -p "$K4_ENV" -f "$k4_src/k4neo.yaml"
 (
  cd "$k4_src"
  POETRY_VIRTUALENVS_CREATE=false "$K4_ENV/bin/poetry" build --format wheel
  k4_wheel="$(ls -1t dist/k4neo-*-py3-none-any.whl | head -1)"
  "$K4_ENV/bin/python" -m pip install "$k4_wheel"
 )
 # This deployment profile is local-only. Upstream's unpinned Slurm plugin
 # currently conflicts with k4neo's pinned pandas/PyYAML requirements.
 "$K4_ENV/bin/python" -m pip uninstall -y snakemake-executor-plugin-slurm snakemake-executor-plugin-slurm-jobstep || true
 "$K4_ENV/bin/python" -m pip install 'packaging>=24,<26'
 "$K4_ENV/bin/python" -m pip check
 PATH="$K4_ENV/bin:$PATH" k4neo-annotator --help >/dev/null
 PATH="$K4_ENV/bin:$PATH" raptor --version >/dev/null
 [[ -n "$K4_SOURCE" ]] || rm -rf "$k4_src"
fi
{
 echo "created_at=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
 echo "python=$(python --version 2>&1)"
 echo "mopepgen=$(python -c 'import moPepGen; print(moPepGen.__version__)')"
 echo "easyquant=$(bp_quant --help 2>&1 | sed -n 's/.*bpquant version \([^ ]*\).*/\1/p' | head -1)"
 echo "pvactools=$(pvacsplice --version 2>&1 | head -1 || true)"
 echo "splice2neo_ref=$S2N_REF"
 echo "splice2neo_rscript=$RSCRIPT"
 echo "k4neo_ref=$K4_REF"
 echo "k4neo_commit=${k4_commit:-}"
 echo "k4neo_submodules=${k4_submodules:-}"
 echo "k4neo_env=${K4_ENV:-}"
 echo "k4neo_license_acknowledged=$K4_ACCEPT"
} > "$ENV_DIR/neoag_splice_v051_versions.txt"
cat "$ENV_DIR/neoag_splice_v051_versions.txt"
