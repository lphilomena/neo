#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
[[ -f "$ROOT/conf/tools.env.sh" ]] && source "$ROOT/conf/tools.env.sh"
: "${NEOAG_CONDA_BASE:?Set NEOAG_CONDA_BASE to a Miniforge/Mambaforge root}"
ENV_NAME="${NEOAG_BAM_MATCHER_ENV:-neoag-bam-matcher}"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-$ROOT}"
INSTALL_DIR="${BAM_MATCHER_HOME:-$TOOLS_ROOT/tools/bam-matcher}"
ENV_PREFIX="${NEOAG_BAM_MATCHER_ENV_PREFIX:-$TOOLS_ROOT/conda_envs/$ENV_NAME}"
WRAPPER_DIR="${NEOAG_WRAPPER_DIR:-$TOOLS_ROOT/bin}"
REPO="${BAM_MATCHER_REPOSITORY:-https://bitbucket.org/sacgf/bam-matcher.git}"
REVISION="${BAM_MATCHER_REVISION:-e0c50af498d35a487cbb0dcff6a3a4506dba695a}"
MAMBA="$NEOAG_CONDA_BASE/bin/mamba"; [[ -x "$MAMBA" ]] || MAMBA="$NEOAG_CONDA_BASE/bin/conda"

if [[ ! -x "$ENV_PREFIX/bin/python2" ]]; then
  mkdir -p "$(dirname "$ENV_PREFIX")"
  "$MAMBA" create -y -p "$ENV_PREFIX" --override-channels -c conda-forge -c bioconda \
    python=2.7 pysam=0.15 pyvcf fisher freebayes samtools pip
fi
if ! "$ENV_PREFIX/bin/python2" -c 'from Cheetah.Template import Template' >/dev/null 2>&1; then
  # Conda's current Cheetah3 builds no longer support Python 2, while the pinned
  # BAM-matcher source imports the historical Cheetah.Template API.
  "$ENV_PREFIX/bin/pip" install --no-cache-dir 'Cheetah3==3.2.6.post2'
fi
"$ENV_PREFIX/bin/python2" -c 'import pysam, vcf; from fisher import pvalue; from Cheetah.Template import Template'
if [[ ! -d "$INSTALL_DIR/.git" ]]; then
  mkdir -p "$(dirname "$INSTALL_DIR")"
  git clone "$REPO" "$INSTALL_DIR"
fi
git -C "$INSTALL_DIR" fetch --tags origin
git -C "$INSTALL_DIR" checkout --detach "$REVISION"
if ! "$ENV_PREFIX/bin/freebayes" --help 2>&1 | grep -q -- '--no-indels'; then
  "$ENV_PREFIX/bin/python2" - "$INSTALL_DIR/bam-matcher.py" <<'PY'
import sys

path = sys.argv[1]
with open(path) as handle:
    text = handle.read()
old = '"--no-indels", "--min-coverage"'
new = '"-i", "-X", "-u", "--min-coverage"'
count = text.count(old)
if count not in (0, 2):
    raise SystemExit("unexpected BAM-matcher FreeBayes command layout: %d matches" % count)
if count:
    with open(path, "w") as handle:
        handle.write(text.replace(old, new))
PY
fi
mkdir -p "$WRAPPER_DIR"
cat > "$WRAPPER_DIR/bam-matcher" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export PATH="$ENV_PREFIX/bin:\${PATH}"
exec "$ENV_PREFIX/bin/python2" "$INSTALL_DIR/bam-matcher.py" "\$@"
EOF
chmod +x "$WRAPPER_DIR/bam-matcher"
"$WRAPPER_DIR/bam-matcher" --help >/dev/null
REF_BUNDLE="${NEOAG_REF_BUNDLE:-}"
BCFTOOLS_BIN="$(command -v bcftools 2>/dev/null || true)"
[[ -x "$BCFTOOLS_BIN" ]] || BCFTOOLS_BIN="${NEOAG_CONDA_BASE}/envs/neoag-tools/bin/bcftools"
if [[ -n "$REF_BUNDLE" && -x "$BCFTOOLS_BIN" ]]; then
  ref_dir="$REF_BUNDLE/data/facets/reference"
  omni="$ref_dir/1000G_omni2.5.hg38.biallelic.vcf.gz"
  common="$ref_dir/common_snp.hg38.vcf.gz"
  identity="$ref_dir/bam_matcher.identity.hg38.vcf.gz"
  contamination="$ref_dir/contamination.common.hg38.vcf.gz"
  primary="$(printf 'chr%s,' {1..22})chrX,chrY"
  if [[ -s "$omni" && ! -s "$identity" ]]; then
    "$BCFTOOLS_BIN" view -r "$primary" -v snps -m2 -M2 -i 'POS%1000000<5000' -Oz -o "$identity" "$omni"
    "$BCFTOOLS_BIN" index -f -t "$identity"
  fi
  if [[ -s "$common" && ! -s "$contamination" ]]; then
    "$BCFTOOLS_BIN" view -r "$primary" -v snps -m2 -M2 -i 'AF>=0.01 && AF<=0.99 && POS%1000000<2000' -Oz -o "$contamination" "$common"
    "$BCFTOOLS_BIN" index -f -t "$contamination"
  fi
fi
echo "BAM-matcher installed at revision $REVISION"
echo "BAM-matcher environment: $ENV_PREFIX"
echo "Provide a GRCh38-compatible identity SNP VCF through bam_matcher_loci in the reference manifest."
