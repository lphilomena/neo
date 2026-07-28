#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
[[ -f "$ROOT/conf/tools.env.sh" ]] && source "$ROOT/conf/tools.env.sh"
: "${NEOAG_CONDA_BASE:?Set NEOAG_CONDA_BASE to a Miniforge/Mambaforge root}"
ENV_NAME="${NEOAG_BAM_MATCHER_ENV:-neoag-bam-matcher}"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-$ROOT}"
INSTALL_DIR="${BAM_MATCHER_HOME:-$TOOLS_ROOT/tools/bam-matcher}"
REPO="${BAM_MATCHER_REPOSITORY:-https://bitbucket.org/sacgf/bam-matcher.git}"
REVISION="${BAM_MATCHER_REVISION:-e0c50af498d35a487cbb0dcff6a3a4506dba695a}"
MAMBA="$NEOAG_CONDA_BASE/bin/mamba"; [[ -x "$MAMBA" ]] || MAMBA="$NEOAG_CONDA_BASE/bin/conda"

if [[ ! -x "$NEOAG_CONDA_BASE/envs/$ENV_NAME/bin/python2" ]]; then
  "$MAMBA" create -y -n "$ENV_NAME" -c conda-forge -c bioconda \
    python=2.7 pysam pyvcf cheetah3 fisher freebayes samtools
fi
if [[ ! -d "$INSTALL_DIR/.git" ]]; then
  mkdir -p "$(dirname "$INSTALL_DIR")"
  git clone "$REPO" "$INSTALL_DIR"
fi
git -C "$INSTALL_DIR" fetch --tags origin
git -C "$INSTALL_DIR" checkout --detach "$REVISION"
mkdir -p "$ROOT/bin"
cat > "$ROOT/bin/bam-matcher" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec "$NEOAG_CONDA_BASE/bin/conda" run -n "$ENV_NAME" --no-capture-output python2 "$INSTALL_DIR/bam-matcher.py" "\$@"
EOF
chmod +x "$ROOT/bin/bam-matcher"
"$ROOT/bin/bam-matcher" --help >/dev/null
echo "BAM-matcher installed at revision $REVISION"
echo "Provide a GRCh38-compatible identity SNP VCF through bam_matcher_loci in the reference manifest."
