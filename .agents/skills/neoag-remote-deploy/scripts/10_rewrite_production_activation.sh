#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="$(pwd)"
TOOLS_ROOT="/root/neo/env_tool"
REFERENCE_ROOT="/root/neo/neodata4git"
LICENSED_ROOT="/root/neo/licensed_tools"
WRITE=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root) PROJECT_ROOT="$2"; shift 2 ;;
    --tools-root) TOOLS_ROOT="$2"; shift 2 ;;
    --reference-root) REFERENCE_ROOT="$2"; shift 2 ;;
    --licensed-root) LICENSED_ROOT="$2"; shift 2 ;;
    --write) WRITE=1; shift ;;
    -h|--help) echo "Usage: $0 --project-root DIR --tools-root DIR --reference-root DIR --licensed-root DIR [--write]"; exit 0 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done
ACT="$TOOLS_ROOT/activate_neoag_production_refs.sh"
if [[ "$WRITE" != "1" ]]; then
  echo "ACTIVATION_REWRITE_REQUIRED: dry run only. Re-run with --write after approval."
  echo "would_write=$ACT"
  exit 0
fi
mkdir -p "$TOOLS_ROOT/bin" "$TOOLS_ROOT/wrappers/mixMHCpred_install"
[[ -f "$ACT" ]] && cp "$ACT" "$ACT.bak_$(date +%Y%m%d_%H%M%S)"
cat > "$TOOLS_ROOT/bin/vep" <<EOF
#!/usr/bin/env bash
set -euo pipefail
VEP_ENV="\${NEOAG_VEP_ENV_PATH:-$TOOLS_ROOT/miniforge3/envs/neoag-vep}"
PERL_BIN="\$VEP_ENV/bin/perl"
VEP_SCRIPT="\$(readlink -f "\$VEP_ENV/bin/vep" 2>/dev/null || echo "\$VEP_ENV/share/ensembl-vep-105.0-0/vep")"
[[ -x "\$PERL_BIN" && -f "\$VEP_SCRIPT" ]] || { echo "ERROR: VEP not found under \$VEP_ENV" >&2; exit 127; }
unset PERL5LIB PERLLIB
export PATH="\$VEP_ENV/bin:/usr/bin:/bin"
exec "\$PERL_BIN" "\$VEP_SCRIPT" "\$@"
EOF
chmod +x "$TOOLS_ROOT/bin/vep"
cat > "$TOOLS_ROOT/bin/netMHCpan" <<EOF
#!/usr/bin/env bash
set -euo pipefail
NMHOME="\${NEOAG_NETMHCPAN_HOME:-$LICENSED_ROOT/netMHCpan}"
PLATFORM="Linux_\$(uname -m)"
BIN="\$NMHOME/\$PLATFORM/bin/netMHCpan-4.2"
[[ -x "\$BIN" ]] || BIN="\$NMHOME/netMHCpan"
[[ -x "\$BIN" ]] || { echo "ERROR: netMHCpan not found under \$NMHOME" >&2; exit 127; }
export NMHOME NETMHCpan="\$NMHOME/\$PLATFORM" TMPDIR="\${NEOAG_NETMHCPAN_TMPDIR:-/tmp}"
exec "\$BIN" "\$@"
EOF
chmod +x "$TOOLS_ROOT/bin/netMHCpan"
cp "$TOOLS_ROOT/bin/netMHCpan" "$TOOLS_ROOT/bin/netmhcpan"
cat > "$TOOLS_ROOT/wrappers/mixMHCpred_install/MixMHCpred" <<EOF
#!/usr/bin/env bash
set -euo pipefail
CONDA_BASE="\${NEOAG_CONDA_BASE:-$TOOLS_ROOT/miniforge3}"
export PATH="\$CONDA_BASE/envs/neoag-core/bin:\$CONDA_BASE/envs/neoag-tools/bin:\${PATH}"
REAL_MIX="\${MIXMHCPRED_REAL_BIN:-$LICENSED_ROOT/mixMHCpred_install/MixMHCpred}"
[[ -x "\$REAL_MIX" ]] || REAL_MIX="$TOOLS_ROOT/tools/mixMHCpred_install/MixMHCpred"
[[ -x "\$REAL_MIX" ]] || { echo "ERROR: set MIXMHCPRED_REAL_BIN to real MixMHCpred" >&2; exit 127; }
exec "\$REAL_MIX" "\$@"
EOF
chmod +x "$TOOLS_ROOT/wrappers/mixMHCpred_install/MixMHCpred"
cat > "$ACT" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export NEOAG_TOOLS_ROOT="$TOOLS_ROOT"
export NEOAG_CONDA_BASE="\${NEOAG_CONDA_BASE:-$TOOLS_ROOT/miniforge3}"
export NEOAG_PROJECT_ROOT="\${NEOAG_PROJECT_ROOT:-$PROJECT_ROOT}"
export PATH="$TOOLS_ROOT/bin:$TOOLS_ROOT/tools/prime:\${NEOAG_PROJECT_ROOT}/bin:\${NEOAG_CONDA_BASE}/envs/neoag-sequenza/bin:\${NEOAG_CONDA_BASE}/envs/neoag-core/bin:\${NEOAG_CONDA_BASE}/envs/neoag-tools/bin:\${NEOAG_CONDA_BASE}/bin:\${PATH}"
export LD_LIBRARY_PATH="\${NEOAG_CONDA_BASE}/envs/neoag-sequenza/lib:\${NEOAG_CONDA_BASE}/envs/neoag-core/lib:\${NEOAG_CONDA_BASE}/envs/neoag-tools/lib\${LD_LIBRARY_PATH:+:\${LD_LIBRARY_PATH}}"
export NEOAG_REFERENCE_FASTA="\${NEOAG_REFERENCE_FASTA:-$REFERENCE_ROOT/data/ref/hg38/Homo_sapiens_assembly38.fasta}"
export NEOAG_GENCODE_GTF="\${NEOAG_GENCODE_GTF:-$REFERENCE_ROOT/data/ref/hg38/gencode.gtf}"
export NEOAG_VEP_CACHE="\${NEOAG_VEP_CACHE:-$REFERENCE_ROOT/data/vep}"
export NEOAG_VEP_CACHE_VERSION="\${NEOAG_VEP_CACHE_VERSION:-105}"
export NEOAG_VEP_BIN="\${NEOAG_VEP_BIN:-$TOOLS_ROOT/bin/vep}"
export PRIME_HOME="\${PRIME_HOME:-$TOOLS_ROOT/tools/prime}"
export NEOAG_PRIME_BIN="\${NEOAG_PRIME_BIN:-\${PRIME_HOME}/PRIME}"
export MIXMHCPRED_REAL_BIN="\${MIXMHCPRED_REAL_BIN:-$LICENSED_ROOT/mixMHCpred_install/MixMHCpred}"
export MIXMHCPRED_BIN="\${MIXMHCPRED_BIN:-$TOOLS_ROOT/wrappers/mixMHCpred_install/MixMHCpred}"
export BIGMHC_DIR="\${BIGMHC_DIR:-$TOOLS_ROOT/tools/bigmhc}"
export BIGMHC_PYTHON="\${BIGMHC_PYTHON:-\${NEOAG_CONDA_BASE}/envs/neoag-tools/bin/python}"
export NEOAG_NETMHCPAN_HOME="\${NEOAG_NETMHCPAN_HOME:-$LICENSED_ROOT/netMHCpan}"
export NEOAG_NETMHCPAN_TMPDIR="\${NEOAG_NETMHCPAN_TMPDIR:-/tmp}"
export NETMHCSTABPAN_HOME="\${NETMHCSTABPAN_HOME:-$LICENSED_ROOT/netMHCstabpan}"
export POLYSOLVER_HOME="\${POLYSOLVER_HOME:-$LICENSED_ROOT/polysolver}"
export NOVOALIGN_LICENSE_FILE="\${NOVOALIGN_LICENSE_FILE:-$LICENSED_ROOT/novoalign/novoalign.lic}"
export SPECHLA_HOME="\${SPECHLA_HOME:-$TOOLS_ROOT/tools/SpecHLA}"
export NEOAG_SPECHLA_HOME="\${NEOAG_SPECHLA_HOME:-\${SPECHLA_HOME}}"
export NEOAG_SPECHLA_IMAGE="\${NEOAG_SPECHLA_IMAGE:-neoag-spechla:ubuntu22.04}"
export HLALA_HOME="\${HLALA_HOME:-$TOOLS_ROOT/tools/HLA-LA}"
export HLA_LA_HOME="\${HLA_LA_HOME:-\${HLALA_HOME}}"
export HLALA_BIN="\${HLALA_BIN:-\${HLALA_HOME}/bin/HLA-LA.pl}"
export HLA_LA_BIN="\${HLA_LA_BIN:-\${HLALA_BIN}}"
export HLALA_GRAPH="\${HLALA_GRAPH:-$REFERENCE_ROOT/data/hla/PRG_MHC_GRCh38_withIMGT}"
export HLA_LA_GRAPH="\${HLA_LA_GRAPH:-\${HLALA_GRAPH}}"
export NEOAG_HLALA_IMAGE="\${NEOAG_HLALA_IMAGE:-neoag-hla-la:ubuntu22.04}"
export HMFTOOLS_HOME="\${HMFTOOLS_HOME:-$TOOLS_ROOT/tools/HMFTOOLS}"
export NEOAG_HMFTOOLS_HOME="\${NEOAG_HMFTOOLS_HOME:-\${HMFTOOLS_HOME}}"
export NEOAG_PURPLE_IMAGE="\${NEOAG_PURPLE_IMAGE:-neoag-purple-suite:ubuntu22.04}"
export HMFTOOLS_AMBER_LOCI="\${HMFTOOLS_AMBER_LOCI:-$REFERENCE_ROOT/data/hmf/purple_reference/amber/GermlineHetPon.38.vcf.gz}"
export HMFTOOLS_GC_PROFILE="\${HMFTOOLS_GC_PROFILE:-$REFERENCE_ROOT/data/hmf/purple_reference/cobalt/GC_profile.1000bp.38.cnp}"
export HMFTOOLS_ENSEMBL_DATA_DIR="\${HMFTOOLS_ENSEMBL_DATA_DIR:-$REFERENCE_ROOT/data/hmf/purple_reference/ensembl_data_cache_38}"
export SEQUENZA_FASTA="\${SEQUENZA_FASTA:-$REFERENCE_ROOT/data/sequenza/reference/GRCh38.primary_assembly.chr.fa}"
export SEQUENZA_GC_WIG="\${SEQUENZA_GC_WIG:-$REFERENCE_ROOT/data/sequenza/reference/gc.wig.gz}"
export TF_USE_LEGACY_KERAS="\${TF_USE_LEGACY_KERAS:-1}"
export CUDA_VISIBLE_DEVICES="\${CUDA_VISIBLE_DEVICES:--1}"
export TF_CPP_MIN_LOG_LEVEL="\${TF_CPP_MIN_LOG_LEVEL:-2}"
EOF
chmod +x "$ACT"
COMMON="$PROJECT_ROOT/scripts/common.sh"
if [[ -f "$COMMON" ]]; then
  cp "$COMMON" "$COMMON.bak_$(date +%Y%m%d_%H%M%S)"
  python3 - "$COMMON" "$TOOLS_ROOT" <<'PY'
import sys
from pathlib import Path
p=Path(sys.argv[1]); tools=sys.argv[2]
text=p.read_text()
lines=[]
for line in text.splitlines():
    if line.startswith('TOOLS_ROOT='):
        lines.append(f'TOOLS_ROOT="${{NEOAG_TOOLS_ROOT:-{tools}}}"')
    else:
        lines.append(line)
p.write_text('\n'.join(lines)+'\n')
PY
fi
echo "activation=$ACT"
