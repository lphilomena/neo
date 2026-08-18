#!/usr/bin/env bash
set -euo pipefail
PROJECT_ROOT="$(pwd)"
DEPLOY_ROOT="${NEOAG_DEPLOY_ROOT:-/opt/neoag}"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-$DEPLOY_ROOT/env_tool}"
REFERENCE_ROOT="${NEOAG_REFERENCE_ROOT:-$DEPLOY_ROOT/refs}"
LICENSED_ROOT="${NEOAG_LICENSED_ROOT:-$DEPLOY_ROOT/licensed_tools}"
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
if [[ -x "$TOOLS_ROOT/conda_pkgs/bedtools-2.31.1-h13024bc_3/bin/bedtools" ]]; then
  mkdir -p "$PROJECT_ROOT/bin"
  cat > "$PROJECT_ROOT/bin/bedtools" <<EOF
#!/usr/bin/env bash
set -euo pipefail
BEDTOOLS_BIN="\${BEDTOOLS_BIN:-$TOOLS_ROOT/conda_pkgs/bedtools-2.31.1-h13024bc_3/bin/bedtools}"
LIB_ROOT="\${BEDTOOLS_LIB_ROOT:-\${NEOAG_CONDA_BASE:-$TOOLS_ROOT/miniforge3}/envs/neoag-tools/lib}"
export LD_LIBRARY_PATH="\${LIB_ROOT}:\${LD_LIBRARY_PATH:-}"
exec "\${BEDTOOLS_BIN}" "\$@"
EOF
  chmod +x "$PROJECT_ROOT/bin/bedtools"
fi
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
export NEOAG_REF_BUNDLE="\${NEOAG_REF_BUNDLE:-$REFERENCE_ROOT}"
export NEOAG_CONDA_BASE="\${NEOAG_CONDA_BASE:-$TOOLS_ROOT/miniforge3}"
export NEOAG_PROJECT_ROOT="\${NEOAG_PROJECT_ROOT:-$PROJECT_ROOT}"
export OPTITYPE_ENV_PREFIX="\${OPTITYPE_ENV_PREFIX:-$TOOLS_ROOT/conda_envs/neoag-optitype}"
export NEOAG_BAM_MATCHER_ENV_PREFIX="\${NEOAG_BAM_MATCHER_ENV_PREFIX:-$TOOLS_ROOT/conda_envs/neoag-bam-matcher}"
export PATH="$TOOLS_ROOT/bin:$TOOLS_ROOT/tools/prime:\${OPTITYPE_ENV_PREFIX}/bin:\${NEOAG_BAM_MATCHER_ENV_PREFIX}/bin:\${NEOAG_PROJECT_ROOT}/bin:\${NEOAG_CONDA_BASE}/envs/neoag-tools/bin:\${NEOAG_CONDA_BASE}/envs/neoag-core/bin:\${NEOAG_CONDA_BASE}/envs/neoag-sequenza/bin:\${NEOAG_CONDA_BASE}/envs/neoag-gatk/bin:\${NEOAG_CONDA_BASE}/bin:\${PATH}"
export LD_LIBRARY_PATH="\${NEOAG_CONDA_BASE}/envs/neoag-tools/lib:\${NEOAG_CONDA_BASE}/envs/neoag-core/lib:\${NEOAG_CONDA_BASE}/envs/neoag-sequenza/lib\${LD_LIBRARY_PATH:+:\${LD_LIBRARY_PATH}}"
export NEOAG_REFERENCE_FASTA="\${NEOAG_REFERENCE_FASTA:-$REFERENCE_ROOT/data/ref/hg38/Homo_sapiens_assembly38.fasta}"
export NEOAG_GENCODE_GTF="\${NEOAG_GENCODE_GTF:-$REFERENCE_ROOT/data/ref/hg38/gencode.gtf}"
if [[ -z "\${NEOAG_NORMAL_JUNCTIONS:-}" ]]; then
  if [[ -f "$REFERENCE_ROOT/data/normal/junctions/normal_junctions.GRCh38.tsv.gz" ]]; then
    export NEOAG_NORMAL_JUNCTIONS="$REFERENCE_ROOT/data/normal/junctions/normal_junctions.GRCh38.tsv.gz"
  elif [[ -f "$REFERENCE_ROOT/data/normal/junctions/gtex_v8_liver.GRCh38.tsv.gz" ]]; then
    export NEOAG_NORMAL_JUNCTIONS="$REFERENCE_ROOT/data/normal/junctions/gtex_v8_liver.GRCh38.tsv.gz"
  fi
fi
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
export NETMHCSTABPAN_HOME="$LICENSED_ROOT/netMHCstabpan"
export NEOAG_NETMHCSTABPAN_BIN="\${NETMHCSTABPAN_HOME}/netMHCstabpan"
export NEOAG_NETMHCSTABPAN_IMAGE="\${NEOAG_NETMHCSTABPAN_IMAGE:-neoag-netmhcstabpan:1.0-ubuntu22.04}"
if [[ -x "\${NEOAG_NETMHCSTABPAN_BIN}" ]]; then
  export PATH="\${NETMHCSTABPAN_HOME}:\${PATH}"
fi
export NETCHOP_HOME="$LICENSED_ROOT/netchop/netchop-3.1"
export NETCHOP="\${NETCHOP_HOME}/Linux_x86_64"
export NEOAG_NETCHOP_BIN="$TOOLS_ROOT/bin/netChop"
export NETCHOP_BIN="\${NEOAG_NETCHOP_BIN}"
if [[ -x "\${NEOAG_NETCHOP_BIN}" ]]; then
  export PATH="$TOOLS_ROOT/bin:\${PATH}"
fi
export POLYSOLVER_HOME="\${POLYSOLVER_HOME:-$LICENSED_ROOT/polysolver}"
export POLYSOLVER_CONDA_ENV="\${POLYSOLVER_CONDA_ENV:-neoag-polysolver}"
export NOVOALIGN_LICENSE_FILE="\${NOVOALIGN_LICENSE_FILE:-$LICENSED_ROOT/novoalign/novoalign.lic}"
export OPTITYPE_BIN="\${OPTITYPE_BIN:-\${OPTITYPE_ENV_PREFIX}/bin/optitype}"
export OPTITYPE_REFERENCE="\${OPTITYPE_REFERENCE:-\${OPTITYPE_ENV_PREFIX}/share/optitype/data}"
export BAM_MATCHER_HOME="\${BAM_MATCHER_HOME:-$TOOLS_ROOT/tools/bam-matcher}"
export BAM_MATCHER_BIN="\${BAM_MATCHER_BIN:-$TOOLS_ROOT/bin/bam-matcher}"
export BAM_MATCHER_PYTHON="\${BAM_MATCHER_PYTHON:-\${NEOAG_BAM_MATCHER_ENV_PREFIX}/bin/python}"
export BAM_MATCHER_REFERENCE="\${BAM_MATCHER_REFERENCE:-$REFERENCE_ROOT/data/sequenza/reference/GRCh38.primary_assembly.chr.fa}"
export BAM_MATCHER_LOCI="\${BAM_MATCHER_LOCI:-$REFERENCE_ROOT/data/facets/reference/bam_matcher.identity.hg38.vcf.gz}"
export FACETS_SNP_VCF="\${FACETS_SNP_VCF:-$REFERENCE_ROOT/data/facets/reference/1000G_omni2.5.hg38.biallelic.vcf.gz}"
export FACETS_R_ENV_PREFIX="\${FACETS_R_ENV_PREFIX:-\${NEOAG_CONDA_BASE}/envs/neoag-fusion}"
export SNP_PILEUP_BIN="\${SNP_PILEUP_BIN:-\${NEOAG_CONDA_BASE}/envs/neoag-tools/bin/snp-pileup}"
export NEOAG_CONTAMINATION_SITES="\${NEOAG_CONTAMINATION_SITES:-$REFERENCE_ROOT/data/facets/reference/contamination.common.hg38.vcf.gz}"
export SPECHLA_HOME="\${SPECHLA_HOME:-$TOOLS_ROOT/tools/SpecHLA}"
export NEOAG_SPECHLA_HOME="\${NEOAG_SPECHLA_HOME:-\${SPECHLA_HOME}}"
export SPECHLA_DB="\${SPECHLA_DB:-$REFERENCE_ROOT/data/hla/spechla/db}"
export SPECHLA_ENV="\${SPECHLA_ENV:-\${NEOAG_CONDA_BASE}/envs/neoag-tools}"
export NEOAG_SPECHLA_IMAGE="\${NEOAG_SPECHLA_IMAGE:-neoag-spechla:ubuntu22.04}"
export HLALA_HOME="\${HLALA_HOME:-$TOOLS_ROOT/tools/HLA-LA}"
export HLA_LA_HOME="\${HLA_LA_HOME:-\${HLALA_HOME}}"
export HLALA_ENV_PREFIX="\${HLALA_ENV_PREFIX:-\${HLALA_HOME}/.conda}"
export HLA_LA_ENV_PREFIX="\${HLA_LA_ENV_PREFIX:-\${HLALA_ENV_PREFIX}}"
export HLALA_BIN="\${HLALA_BIN:-\${HLALA_ENV_PREFIX}/bin/HLA-LA.pl}"
export HLA_LA_BIN="\${HLA_LA_BIN:-\${HLALA_BIN}}"
export HLALA_GRAPH="\${HLALA_GRAPH:-$REFERENCE_ROOT/data/hla/PRG_MHC_GRCh38_withIMGT}"
export HLA_LA_GRAPH="\${HLA_LA_GRAPH:-\${HLALA_GRAPH}}"
export NEOAG_HLALA_BACKEND="\${NEOAG_HLALA_BACKEND:-auto}"
export NEOAG_HLALA_IMAGE="\${NEOAG_HLALA_IMAGE:-neoag-hla-la:ubuntu22.04}"
export HMFTOOLS_HOME="\${HMFTOOLS_HOME:-$TOOLS_ROOT/tools/HMFTOOLS}"
export NEOAG_HMFTOOLS_HOME="\${NEOAG_HMFTOOLS_HOME:-\${HMFTOOLS_HOME}}"
export HMF_ENV="\${HMF_ENV:-\${HMFTOOLS_HOME}/.conda}"
export NEOAG_PURPLE_IMAGE="\${NEOAG_PURPLE_IMAGE:-neoag-purple-suite:ubuntu22.04}"
export HMFTOOLS_REFERENCE_ROOT="\${HMFTOOLS_REFERENCE_ROOT:-$REFERENCE_ROOT/data/hmf/purple_reference}"
export HMFTOOLS_REFERENCE_FASTA="\${HMFTOOLS_REFERENCE_FASTA:-$REFERENCE_ROOT/data/sequenza/reference/GRCh38.primary_assembly.chr.fa}"
export HMFTOOLS_AMBER_LOCI="\${HMFTOOLS_AMBER_LOCI:-$REFERENCE_ROOT/data/hmf/purple_reference/amber/GermlineHetPon.38.vcf.gz}"
export HMFTOOLS_GC_PROFILE="\${HMFTOOLS_GC_PROFILE:-$REFERENCE_ROOT/data/hmf/purple_reference/cobalt/GC_profile.1000bp.38.cnp}"
export HMFTOOLS_ENSEMBL_DATA_DIR="\${HMFTOOLS_ENSEMBL_DATA_DIR:-$REFERENCE_ROOT/data/hmf/purple_reference/ensembl_data_cache_38}"
export SEQUENZA_FASTA="\${SEQUENZA_FASTA:-$REFERENCE_ROOT/data/sequenza/reference/GRCh38.primary_assembly.chr.fa}"
export SEQUENZA_GC_WIG="\${SEQUENZA_GC_WIG:-$REFERENCE_ROOT/data/sequenza/reference/Homo_sapiens.GRCh38.dna.primary_assembly.chr.gc50.wig.gz}"
export SAMTOOLS="\${SAMTOOLS:-\${NEOAG_CONDA_BASE}/envs/neoag-tools/bin/samtools}"
export TF_USE_LEGACY_KERAS="\${TF_USE_LEGACY_KERAS:-1}"
export CUDA_VISIBLE_DEVICES="\${CUDA_VISIBLE_DEVICES:--1}"
export TF_CPP_MIN_LOG_LEVEL="\${TF_CPP_MIN_LOG_LEVEL:-2}"
EOF
chmod +x "$ACT"
LOCAL_ENV="$PROJECT_ROOT/conf/tools.env.local.sh"
python3 - "$LOCAL_ENV" "$TOOLS_ROOT" "$REFERENCE_ROOT" "$PROJECT_ROOT" <<'PY'
import sys
from pathlib import Path

path = Path(sys.argv[1])
tools_root, reference_root, project_root = sys.argv[2:]
begin = "# BEGIN NEOAG PRODUCTION OVERRIDES"
end = "# END NEOAG PRODUCTION OVERRIDES"
text = path.read_text() if path.exists() else ""
if begin in text and end in text:
    prefix, remainder = text.split(begin, 1)
    _, suffix = remainder.split(end, 1)
    text = prefix.rstrip() + "\n\n" + suffix.lstrip()
block = f'''{begin}
# Unconditional site paths override portable bundle paths from another host.
export OPTITYPE_ENV="{tools_root}/conda_envs/neoag-optitype"
export OPTITYPE_ENV_PREFIX="${{OPTITYPE_ENV}}"
export OPTITYPE_BIN="${{OPTITYPE_ENV}}/bin/optitype"
export OPTITYPE_REFERENCE="${{OPTITYPE_ENV}}/share/optitype/data"
export NEOAG_BAM_MATCHER_ENV_PREFIX="{tools_root}/conda_envs/neoag-bam-matcher"
export BAM_MATCHER_HOME="{tools_root}/tools/bam-matcher"
export BAM_MATCHER_BIN="{tools_root}/bin/bam-matcher"
export BAM_MATCHER_PYTHON="${{NEOAG_BAM_MATCHER_ENV_PREFIX}}/bin/python"
export BAM_MATCHER_REFERENCE="{reference_root}/data/sequenza/reference/GRCh38.primary_assembly.chr.fa"
export BAM_MATCHER_LOCI="{reference_root}/data/facets/reference/bam_matcher.identity.hg38.vcf.gz"
export FACETS_SNP_VCF="{reference_root}/data/facets/reference/1000G_omni2.5.hg38.biallelic.vcf.gz"
export FACETS_R_ENV_PREFIX="${{NEOAG_CONDA_BASE}}/envs/neoag-fusion"
export SNP_PILEUP_BIN="${{NEOAG_CONDA_BASE}}/envs/neoag-tools/bin/snp-pileup"
export NEOAG_CONTAMINATION_SITES="{reference_root}/data/facets/reference/contamination.common.hg38.vcf.gz"
export PRIME_HOME="{reference_root}/data/predictors/prime"
export NEOAG_PRIME_BIN="${{PRIME_HOME}}/PRIME"
export MIXMHCPRED_HOME="{reference_root}/data/predictors/mixMHCpred_install"
export MIXMHCPRED_BIN="${{MIXMHCPRED_HOME}}/MixMHCpred"
export HMFTOOLS_HOME="{tools_root}/tools/HMFTOOLS"
export NEOAG_HMFTOOLS_HOME="${{HMFTOOLS_HOME}}"
export HMF_ENV="${{HMFTOOLS_HOME}}/.conda"
export HMFTOOLS_REFERENCE_ROOT="{reference_root}/data/hmf/purple_reference"
export HMFTOOLS_REFERENCE_FASTA="{reference_root}/data/sequenza/reference/GRCh38.primary_assembly.chr.fa"
export HMFTOOLS_AMBER_LOCI="{reference_root}/data/hmf/purple_reference/amber/GermlineHetPon.38.vcf.gz"
export HMFTOOLS_GC_PROFILE="{reference_root}/data/hmf/purple_reference/cobalt/GC_profile.1000bp.38.cnp"
export HMFTOOLS_ENSEMBL_DATA_DIR="{reference_root}/data/hmf/purple_reference/ensembl_data_cache_38"
export SEQUENZA_FASTA="{reference_root}/data/sequenza/reference/GRCh38.primary_assembly.chr.fa"
export SEQUENZA_GC_WIG="{reference_root}/data/sequenza/reference/Homo_sapiens.GRCh38.dna.primary_assembly.chr.gc50.wig.gz"
export SAMTOOLS="${{NEOAG_CONDA_BASE}}/envs/neoag-tools/bin/samtools"
export PATH="{project_root}/bin:${{OPTITYPE_ENV}}/bin:${{NEOAG_BAM_MATCHER_ENV_PREFIX}}/bin:${{PRIME_HOME}}:${{MIXMHCPRED_HOME}}:${{NEOAG_CONDA_BASE}}/envs/neoag-tools/bin:${{PATH}}"
{end}
'''
path.parent.mkdir(parents=True, exist_ok=True)
path.write_text(text.rstrip() + "\n\n" + block)
PY
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
