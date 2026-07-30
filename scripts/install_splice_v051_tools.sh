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
  --install-k4neo                Install k4neo source after license review
  --k4neo-ref REF                Required with --install-k4neo
  --k4neo-license-accepted       Explicit acknowledgement
  --python PYTHON                Default: python3

The script records exact versions. It does not download GTEx/TCGA/k4neo indices,
IEDB assets, licensed NetMHC software, or patient data.
USAGE
}
ENV_DIR=""; PY="python3"; MO="1.4.6"; EQ="0.6.0"; PV=""; S2N=0; S2N_REF="v0.6.14"; K4=0; K4_REF=""; K4_ACCEPT=0
while [[ $# -gt 0 ]]; do
 case "$1" in
  --env-dir) ENV_DIR="$2"; shift 2;; --python) PY="$2"; shift 2;;
  --mopepgen-version) MO="$2"; shift 2;; --easyquant-version) EQ="$2"; shift 2;;
  --pvactools-version) PV="$2"; shift 2;; --install-splice2neo) S2N=1; shift;;
  --splice2neo-ref) S2N_REF="$2"; shift 2;; --install-k4neo) K4=1; shift;;
  --k4neo-ref) K4_REF="$2"; shift 2;; --k4neo-license-accepted) K4_ACCEPT=1; shift;;
  -h|--help) usage; exit 0;; *) echo "ERROR: unknown argument $1" >&2; exit 2;;
 esac
done
[[ -n "$ENV_DIR" ]] || { echo "ERROR: --env-dir required" >&2; exit 2; }
if [[ "$K4" == 1 ]]; then [[ "$K4_ACCEPT" == 1 && -n "$K4_REF" ]] || { echo "ERROR: k4neo requires --k4neo-ref and --k4neo-license-accepted" >&2; exit 2; }; fi
"$PY" -m venv "$ENV_DIR"
# shellcheck disable=SC1091
source "$ENV_DIR/bin/activate"
python -m pip install --upgrade pip setuptools wheel
python -m pip install "mopepgen==$MO" "bp-quant==$EQ"
[[ -n "$PV" ]] && python -m pip install "pvactools==$PV"
if [[ "$S2N" == 1 ]]; then
 command -v Rscript >/dev/null || { echo "ERROR: Rscript required for splice2neo" >&2; exit 3; }
 Rscript -e "if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='https://cloud.r-project.org'); remotes::install_github('TRON-Bioinformatics/splice2neo@${S2N_REF}')"
fi
if [[ "$K4" == 1 ]]; then
 python -m pip install "git+https://github.com/TRON-Bioinformatics/k4neo.git@${K4_REF}"
fi
{
 echo "created_at=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
 echo "python=$(python --version 2>&1)"
 echo "mopepgen=$(moPepGen --version 2>&1 | head -1 || true)"
 echo "easyquant=$(bp_quant --version 2>&1 | head -1 || true)"
 echo "pvactools=$(pvacsplice --version 2>&1 | head -1 || true)"
 echo "splice2neo_ref=$S2N_REF"
 echo "k4neo_ref=$K4_REF"
 echo "k4neo_license_acknowledged=$K4_ACCEPT"
} > "$ENV_DIR/neoag_splice_v051_versions.txt"
cat "$ENV_DIR/neoag_splice_v051_versions.txt"
