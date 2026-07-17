#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(pwd)"
TOOLS_ROOT="/root/neo/env_tool"
LICENSED_ROOT="/root/neo/licensed_tools"
REFERENCE_ROOT="/root/neo/neodata4git"
CONDA_BASE=""
OUTDIR="work/remote_deploy"
EXECUTE=0
ALLOW_DOWNLOAD=0
INSTALL_MINIFORGE=1
MINIFORGE_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"

INSTALL_CORE_ENV=0
INSTALL_VEP=0
INSTALL_VEP_CACHE=0
VEP_VERSION="105"
INSTALL_GATK=0
INSTALL_IMMUNOGENICITY=0
INSTALL_DEEPIMMUNO=0
INSTALL_NETMHCSTABPAN=0
INSTALL_LOHHLA=0
INSTALL_POLYSOLVER=0
INSTALL_OPTITYPE=0
INSTALL_FACETS=0
INSTALL_ASCAT_PYCLONE=0
INSTALL_FUSION=0
RUN_VERIFY=0
STRICT_VERIFY=0
RUN_REAL_VCF_SMOKE=0
REAL_VCF_SMOKE_TOP_N=50
REAL_VCF_SMOKE_SKIP_MHCFLURRY=0
REAL_VCF_RAW=""
REAL_VCF_ANNOTATED=""
REAL_VCF_HLA_ALLELES=""
REAL_VCF_HLA_FILE=""
BIGMHC_MODELS_DIR=""
BIGMHC_MODELS_HOST=""
ASSET_MANIFEST="configs/assets/production_assets.tsv"
REFERENCE_MANIFEST="configs/references/reference_manifest.yaml"
SYNC_ASSETS=0
ASSET_SOURCE_HOST=""
CORE_ENV_LITE=1
SKIP_TORCH_INSTALL=1

NETMHCPAN_TAR=""
NETMHCPAN_DIR=""
NETMHCPAN_URL=""
MIXMHCPRED_DIR=""
MIXMHCPRED_ARCHIVE=""
MIXMHCPRED_URL=""
NETMHCSTABPAN_DIR=""
NETMHCSTABPAN_ARCHIVE=""
NETMHCSTABPAN_URL=""
POLYSOLVER_HOME_ARG=""
NOVOALIGN_LICENSE_FILE_ARG=""
DEEPIMMUNO_SOURCE=""

usage() {
  cat <<'USAGE'
Usage: 13_install_readme_tools.sh [options]

Install the external tools listed in README.md by orchestrating the repository's
existing installation scripts. The default mode is dry-run; add --execute to
actually install. Licensed tools still require local archives/directories or an
explicit approved URL plus --allow-download.

Common options:
  --project-root DIR          Project checkout (default: current directory)
  --tools-root DIR            Tool/env root (default: /root/neo/env_tool)
  --licensed-root DIR         Licensed tool root (default: /root/neo/licensed_tools)
  --reference-root DIR        Reference root (default: /root/neo/neodata4git)
  --conda-base DIR            Miniforge/conda base (default: tools-root/miniforge3)
  --install-miniforge         Install Miniforge3 if conda is missing (default enabled)
  --no-install-miniforge      Do not install Miniforge automatically if conda is missing
  --miniforge-url URL         Miniforge installer URL
  --outdir DIR                Log/report directory (default: work/remote_deploy)
  --allow-download            Permit downloads from official/user-approved URLs
  --execute                   Actually run installation commands
  --full-core-env             Use full core env instead of default lite env
  --install-torch             Let immunogenicity installer install torch if missing

Tool groups:
  --core-env                  pVACtools/MHCflurry core conda env via scripts/setup_tools_env.sh
  --vep                      VEP conda env via scripts/install_vep.sh
  --vep-cache                VEP cache via scripts/install_vep_cache.sh (large download)
  --vep-version VERSION      Ensembl VEP/cache release to install/use (default: 105)
  --gatk                     GATK4 / Mutect2 via scripts/install_gatk.sh
  --immunogenicity           PRIME + MixMHCpred + BigMHC via scripts/install_immunogenicity_tools.sh
  --deepimmuno               DeepImmuno via scripts/install_deepimmuno.sh
  --netmhcstabpan            NetMHCstabpan or IEDB shim via scripts/install_netmhcstabpan.sh
  --lohhla                   LOHHLA source wrapper via scripts/install_lohhla.sh
  --polysolver               Configure existing Polysolver; requires --polysolver-home
  --optitype                 OptiType via scripts/install_optitype.sh
  --facets                   FACETS via scripts/install_facets.sh
  --ascat-pyclone            ASCAT + PyClone-VI via scripts/install_ascat_pyclone.sh
  --fusion                   Arriba/Nextflow fusion env plus STAR-Fusion/FusionCatcher clones
  --all-open                 Install open/conda/git tools except very large VEP cache and licensed packages
  --all                      Install all supported groups, including VEP cache; licensed packages still need sources/URLs
  --verify                   Run scripts/verify_all_tools_and_refs.sh after installs
  --strict-verify            Treat optional missing tools as failure during verify
  --run-real-vcf-smoke       Run default M1ML150017383 real VCF top-N smoke test
  --real-vcf-smoke-top-n N   Number of unique peptides for real VCF smoke (default: 50)
  --skip-real-vcf-mhcflurry Skip MHCflurry only for the real VCF smoke fallback
  --real-vcf FILE           Override raw VCF path for real VCF smoke
  --real-annotated-vcf FILE Override VEP-annotated VCF path for real VCF smoke
  --real-vcf-hla-alleles L  Override comma-separated HLA alleles for real VCF smoke
  --real-vcf-hla-file FILE  Override HLA file for real VCF smoke
  --bigmhc-models-dir DIR   Copy BigMHC models from local/source directory into tools-root
  --bigmhc-models-host HOST Optional source host for --bigmhc-models-dir, e.g. na@10.200.50.134
  --asset-manifest FILE    TSV manifest for large assets (default: configs/assets/production_assets.tsv)
  --reference-manifest FILE
                          YAML reference manifest verified after asset sync
  --sync-assets            Sync large assets from manifest (dry-run unless --execute)
  --asset-source-host HOST Default source host for manifest source_path values

Licensed/restricted source options:
  --netmhcpan-tar FILE       Local NetMHCpan archive
  --netmhcpan-dir DIR        Existing NetMHCpan directory to copy
  --netmhcpan-url URL        Approved NetMHCpan archive URL
  --mixmhcpred-dir DIR       Existing MixMHCpred directory to copy
  --mixmhcpred-archive FILE  Local MixMHCpred archive
  --mixmhcpred-url URL       Approved MixMHCpred archive URL
  --netmhcstabpan-dir DIR    Existing NetMHCstabpan directory to copy
  --netmhcstabpan-archive FILE
  --netmhcstabpan-url URL    Approved NetMHCstabpan archive URL
  --polysolver-home DIR      Existing Polysolver distribution
  --novoalign-license-file FILE
  --deepimmuno-source DIR    Existing DeepImmuno checkout; otherwise script clones official repo

Examples:
  bash .agents/skills/neoag-remote-deploy/scripts/13_install_readme_tools.sh \
    --project-root /root/neo/src/na0707_upload_release \
    --tools-root /root/neo/env_tool \
    --conda-base /root/neo/env_tool/miniforge3 \
    --core-env --vep --gatk --optitype --allow-download --execute

  bash .agents/skills/neoag-remote-deploy/scripts/13_install_readme_tools.sh \
    --all-open --verify --execute

Notes:
  - Use only official or user-approved URLs.
  - This script does not bypass registration, login, license, or institutional access controls.
  - README tools that require external databases or licensed resources may install wrappers first and still need references configured before production use.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root) PROJECT_ROOT="$2"; shift 2 ;;
    --tools-root) TOOLS_ROOT="$2"; shift 2 ;;
    --licensed-root) LICENSED_ROOT="$2"; shift 2 ;;
    --reference-root) REFERENCE_ROOT="$2"; shift 2 ;;
    --conda-base) CONDA_BASE="$2"; shift 2 ;;
    --install-miniforge) INSTALL_MINIFORGE=1; shift ;;
    --no-install-miniforge) INSTALL_MINIFORGE=0; shift ;;
    --miniforge-url) MINIFORGE_URL="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --allow-download) ALLOW_DOWNLOAD=1; shift ;;
    --execute) EXECUTE=1; shift ;;
    --full-core-env) CORE_ENV_LITE=0; shift ;;
    --install-torch) SKIP_TORCH_INSTALL=0; shift ;;
    --core-env) INSTALL_CORE_ENV=1; shift ;;
    --vep) INSTALL_VEP=1; shift ;;
    --vep-cache) INSTALL_VEP=1; INSTALL_VEP_CACHE=1; shift ;;
    --vep-version) VEP_VERSION="$2"; shift 2 ;;
    --gatk) INSTALL_GATK=1; shift ;;
    --immunogenicity) INSTALL_IMMUNOGENICITY=1; shift ;;
    --deepimmuno) INSTALL_DEEPIMMUNO=1; shift ;;
    --netmhcstabpan) INSTALL_NETMHCSTABPAN=1; shift ;;
    --lohhla) INSTALL_LOHHLA=1; shift ;;
    --polysolver) INSTALL_POLYSOLVER=1; shift ;;
    --optitype) INSTALL_OPTITYPE=1; shift ;;
    --facets) INSTALL_FACETS=1; shift ;;
    --ascat-pyclone) INSTALL_ASCAT_PYCLONE=1; shift ;;
    --fusion) INSTALL_FUSION=1; shift ;;
    --all-open)
      INSTALL_CORE_ENV=1; INSTALL_VEP=1; INSTALL_GATK=1; INSTALL_IMMUNOGENICITY=1
      INSTALL_DEEPIMMUNO=1; INSTALL_NETMHCSTABPAN=1; INSTALL_LOHHLA=1
      INSTALL_OPTITYPE=1; INSTALL_FACETS=1; INSTALL_ASCAT_PYCLONE=1; INSTALL_FUSION=1
      shift ;;
    --all)
      INSTALL_CORE_ENV=1; INSTALL_VEP=1; INSTALL_VEP_CACHE=1; INSTALL_GATK=1; INSTALL_IMMUNOGENICITY=1
      INSTALL_DEEPIMMUNO=1; INSTALL_NETMHCSTABPAN=1; INSTALL_LOHHLA=1
      INSTALL_OPTITYPE=1; INSTALL_FACETS=1; INSTALL_ASCAT_PYCLONE=1; INSTALL_FUSION=1
      shift ;;
    --verify) RUN_VERIFY=1; shift ;;
    --strict-verify) RUN_VERIFY=1; STRICT_VERIFY=1; shift ;;
    --run-real-vcf-smoke) RUN_REAL_VCF_SMOKE=1; shift ;;
    --real-vcf-smoke-top-n) REAL_VCF_SMOKE_TOP_N="$2"; shift 2 ;;
    --skip-real-vcf-mhcflurry) REAL_VCF_SMOKE_SKIP_MHCFLURRY=1; shift ;;
    --real-vcf) REAL_VCF_RAW="$2"; shift 2 ;;
    --real-annotated-vcf) REAL_VCF_ANNOTATED="$2"; shift 2 ;;
    --real-vcf-hla-alleles) REAL_VCF_HLA_ALLELES="$2"; shift 2 ;;
    --real-vcf-hla-file) REAL_VCF_HLA_FILE="$2"; shift 2 ;;
    --bigmhc-models-dir) BIGMHC_MODELS_DIR="$2"; shift 2 ;;
    --bigmhc-models-host) BIGMHC_MODELS_HOST="$2"; shift 2 ;;
    --asset-manifest) ASSET_MANIFEST="$2"; shift 2 ;;
    --reference-manifest) REFERENCE_MANIFEST="$2"; shift 2 ;;
    --sync-assets) SYNC_ASSETS=1; shift ;;
    --asset-source-host) ASSET_SOURCE_HOST="$2"; shift 2 ;;
    --netmhcpan-tar) NETMHCPAN_TAR="$2"; shift 2 ;;
    --netmhcpan-dir) NETMHCPAN_DIR="$2"; shift 2 ;;
    --netmhcpan-url) NETMHCPAN_URL="$2"; shift 2 ;;
    --mixmhcpred-dir) MIXMHCPRED_DIR="$2"; shift 2 ;;
    --mixmhcpred-archive) MIXMHCPRED_ARCHIVE="$2"; shift 2 ;;
    --mixmhcpred-url) MIXMHCPRED_URL="$2"; shift 2 ;;
    --netmhcstabpan-dir) NETMHCSTABPAN_DIR="$2"; shift 2 ;;
    --netmhcstabpan-archive) NETMHCSTABPAN_ARCHIVE="$2"; shift 2 ;;
    --netmhcstabpan-url) NETMHCSTABPAN_URL="$2"; shift 2 ;;
    --polysolver-home) POLYSOLVER_HOME_ARG="$2"; shift 2 ;;
    --novoalign-license-file) NOVOALIGN_LICENSE_FILE_ARG="$2"; shift 2 ;;
    --deepimmuno-source) DEEPIMMUNO_SOURCE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

mkdir -p "$OUTDIR"
LOG="$OUTDIR/readme_tools_install.log"
REPORT="$OUTDIR/readme_tools_install_report.md"
: > "$LOG"
MODE="DRY_RUN"
[[ "$EXECUTE" == "1" ]] && MODE="EXECUTE"

log() { printf '%s\n' "$*" | tee -a "$LOG"; }
run() {
  local label="$1"; shift
  log ""
  log "==> [$MODE] $label"
  log "+ $*"
  if [[ "$EXECUTE" == "1" ]]; then
    "$@" 2>&1 | tee -a "$LOG"
  fi
}

need_download_ok() {
  local what="$1"
  if [[ "$ALLOW_DOWNLOAD" != "1" ]]; then
    echo "DOWNLOAD_NOT_APPROVED: $what requires network download; add --allow-download after approval" >&2
    exit 23
  fi
}

find_conda_base() {
  local preferred="${CONDA_BASE:-$TOOLS_ROOT/miniforge3}"
  if [[ -x "$preferred/bin/conda" ]]; then echo "$preferred"; return 0; fi
  return 1
}

set_local_conda_pkg_cache() {
  run "set local conda package cache" bash -lc "mkdir -p '$TOOLS_ROOT/conda_pkgs' && '$CONDA_BASE/bin/conda' config --remove-key pkgs_dirs >/dev/null 2>&1 || true; '$CONDA_BASE/bin/conda' config --add pkgs_dirs '$TOOLS_ROOT/conda_pkgs' >/dev/null 2>&1"
}



ensure_reference_indexes_after_asset_sync() {
  local fasta="$REFERENCE_ROOT/data/ref/hg38/Homo_sapiens_assembly38.fasta"
  if [[ -s "$fasta" && ! -s "$fasta.fai" ]]; then
    if command -v samtools >/dev/null 2>&1; then
      run "index reference FASTA" samtools faidx "$fasta"
    else
      log "WARN: samtools not found; cannot create FASTA index: $fasta.fai"
    fi
  fi
}

sync_assets_if_requested() {
  [[ "$SYNC_ASSETS" == "1" ]] || return 0
  args=(--project-root "$PROJECT_ROOT" --asset-manifest "$ASSET_MANIFEST" --outdir "$OUTDIR/assets")
  [[ -n "$ASSET_SOURCE_HOST" ]] && args+=(--asset-source-host "$ASSET_SOURCE_HOST")
  [[ "$EXECUTE" == "1" ]] && args+=(--execute)
  run "sync large assets from manifest" bash .agents/skills/neoag-remote-deploy/scripts/15_sync_asset_manifest.sh "${args[@]}"
  ensure_reference_indexes_after_asset_sync
  if [[ -f "$REFERENCE_MANIFEST" ]]; then
    verify_ref_args=("$REFERENCE_MANIFEST" --vep-version "$VEP_VERSION")
    [[ "$STRICT_VERIFY" == "1" ]] && verify_ref_args+=(--strict)
    run "verify reference manifest" python3 scripts/verify_reference_manifest.py "${verify_ref_args[@]}"
  else
    log "WARN: reference manifest not found: $REFERENCE_MANIFEST"
  fi
}

stage_bigmhc_models_if_requested() {
  [[ -n "$BIGMHC_MODELS_DIR" ]] || return 0
  local dst="$TOOLS_ROOT/tools/bigmhc/models"
  if [[ -n "$BIGMHC_MODELS_HOST" ]]; then
    run "copy BigMHC models from source host" bash -lc "mkdir -p '$dst' && rsync -a '$BIGMHC_MODELS_HOST:$BIGMHC_MODELS_DIR/' '$dst/'"
  else
    [[ -d "$BIGMHC_MODELS_DIR" ]] || { echo "BIGMHC_MODELS_SOURCE_MISSING: $BIGMHC_MODELS_DIR" >&2; exit 45; }
    run "copy BigMHC models from local directory" bash -lc "mkdir -p '$dst' && cp -a '$BIGMHC_MODELS_DIR/.' '$dst/'"
  fi
}

install_miniforge_if_needed() {
  if CONDA_BASE_FOUND="$(find_conda_base 2>/dev/null)"; then
    CONDA_BASE="$CONDA_BASE_FOUND"
    log "Conda found: $CONDA_BASE"
    set_local_conda_pkg_cache
    return 0
  fi
  CONDA_BASE="${CONDA_BASE:-$TOOLS_ROOT/miniforge3}"
  [[ "$INSTALL_MINIFORGE" == "1" ]] || { echo "CONDA_MISSING: set --conda-base or allow default Miniforge install" >&2; exit 31; }
  need_download_ok "Miniforge3 installer"
  local installer="$OUTDIR/Miniforge3-Linux-x86_64.sh"
  run "download Miniforge3" bash -lc "mkdir -p '$OUTDIR' '$TOOLS_ROOT' && curl -fL --retry 3 -o '$installer' '$MINIFORGE_URL'"
  run "install Miniforge3" bash -lc "bash '$installer' -b -p '$CONDA_BASE'"
  set_local_conda_pkg_cache
}

cd "$PROJECT_ROOT"
[[ -f "pyproject.toml" || -f "setup.py" ]] || { echo "PROJECT_ROOT_INVALID: $PROJECT_ROOT" >&2; exit 30; }

if [[ "$INSTALL_CORE_ENV$INSTALL_VEP$INSTALL_GATK$INSTALL_IMMUNOGENICITY$INSTALL_DEEPIMMUNO$INSTALL_NETMHCSTABPAN$INSTALL_LOHHLA$INSTALL_POLYSOLVER$INSTALL_OPTITYPE$INSTALL_FACETS$INSTALL_ASCAT_PYCLONE$INSTALL_FUSION" =~ 1 ]]; then
  install_miniforge_if_needed
  export NEOAG_CONDA_BASE="$CONDA_BASE"
  export PATH="$CONDA_BASE/bin:$PATH"
fi

export NEOAG_TOOLS_ROOT="$TOOLS_ROOT"
export NEOAG_REF_BUNDLE="$REFERENCE_ROOT"
export NETMHCPAN_HOME="$LICENSED_ROOT/netMHCpan"
export NETMHCpan="$LICENSED_ROOT/netMHCpan"
export NETMHCSTABPAN_HOME="$LICENSED_ROOT/netMHCstabpan"
export PRIME_HOME="$TOOLS_ROOT/tools/prime"
export MIXMHCPRED_HOME="$LICENSED_ROOT/mixMHCpred_install"
export BIGMHC_DIR="$TOOLS_ROOT/tools/bigmhc"
export DEEPIMMUNO_DIR="$TOOLS_ROOT/tools/DeepImmuno"

if [[ -n "$NETMHCPAN_TAR$NETMHCPAN_DIR$NETMHCPAN_URL$MIXMHCPRED_DIR$MIXMHCPRED_ARCHIVE$MIXMHCPRED_URL$NETMHCSTABPAN_DIR$NETMHCSTABPAN_ARCHIVE$NETMHCSTABPAN_URL" ]]; then
  args=(--licensed-root "$LICENSED_ROOT" --outdir "$OUTDIR")
  [[ -n "$NETMHCPAN_TAR" ]] && args+=(--netmhcpan-tar "$NETMHCPAN_TAR")
  [[ -n "$NETMHCPAN_DIR" ]] && args+=(--netmhcpan-dir "$NETMHCPAN_DIR")
  [[ -n "$NETMHCPAN_URL" ]] && args+=(--netmhcpan-url "$NETMHCPAN_URL")
  [[ -n "$MIXMHCPRED_DIR" ]] && args+=(--mixmhcpred-dir "$MIXMHCPRED_DIR")
  [[ -n "$MIXMHCPRED_ARCHIVE" ]] && args+=(--mixmhcpred-archive "$MIXMHCPRED_ARCHIVE")
  [[ -n "$MIXMHCPRED_URL" ]] && args+=(--mixmhcpred-url "$MIXMHCPRED_URL")
  [[ -n "$NETMHCSTABPAN_DIR" ]] && args+=(--netmhcstabpan-dir "$NETMHCSTABPAN_DIR")
  [[ -n "$NETMHCSTABPAN_ARCHIVE" ]] && args+=(--netmhcstabpan-archive "$NETMHCSTABPAN_ARCHIVE")
  [[ -n "$NETMHCSTABPAN_URL" ]] && args+=(--netmhcstabpan-url "$NETMHCSTABPAN_URL")
  [[ "$ALLOW_DOWNLOAD" == "1" ]] && args+=(--allow-download)
  [[ "$EXECUTE" == "1" ]] && args+=(--execute)
  run "install local/downloaded licensed tools" bash .agents/skills/neoag-remote-deploy/scripts/12_install_local_licensed_tools.sh "${args[@]}"
fi

[[ "$INSTALL_CORE_ENV" == "1" ]] && run "install core pVACtools/MHCflurry env" env NEOAG_TOOLS_LITE="$CORE_ENV_LITE" bash scripts/setup_tools_env.sh
[[ "$INSTALL_VEP" == "1" ]] && run "install VEP env" env NEOAG_VEP_VERSION="$VEP_VERSION" bash scripts/install_vep.sh
[[ "$INSTALL_VEP_CACHE" == "1" ]] && { need_download_ok "VEP cache"; run "install VEP cache" env NEOAG_VEP_CACHE_VERSION="$VEP_VERSION" bash scripts/install_vep_cache.sh; }
[[ "$INSTALL_GATK" == "1" ]] && run "install GATK4" bash scripts/install_gatk.sh
IMMUNO_PYTHON="${CONDA_BASE}/envs/neoag-tools/bin/python"
[[ -x "$IMMUNO_PYTHON" ]] || IMMUNO_PYTHON="${CONDA_BASE}/envs/neoag-core/bin/python"
[[ -x "$IMMUNO_PYTHON" ]] || IMMUNO_PYTHON="${CONDA_BASE}/bin/python"
[[ -x "$IMMUNO_PYTHON" ]] || IMMUNO_PYTHON="$(command -v python3)"
[[ "$INSTALL_IMMUNOGENICITY" == "1" ]] && run "install PRIME/MixMHCpred/BigMHC" env NEOAG_SKIP_TORCH_INSTALL="$SKIP_TORCH_INSTALL" NEOAG_IMMUNO_PYTHON="$IMMUNO_PYTHON" bash scripts/install_immunogenicity_tools.sh
if [[ "$INSTALL_NETMHCSTABPAN" == "1" ]]; then
  if [[ -n "$NETMHCSTABPAN_ARCHIVE" && -f "$NETMHCSTABPAN_ARCHIVE" ]]; then
    run "install NetMHCstabpan from archive" bash scripts/install_netmhcstabpan.sh "$NETMHCSTABPAN_ARCHIVE"
  else
    run "install NetMHCstabpan IEDB shim" bash scripts/install_netmhcstabpan.sh --iedb
  fi
fi
if [[ "$INSTALL_DEEPIMMUNO" == "1" ]]; then
  [[ -z "$DEEPIMMUNO_SOURCE" ]] && need_download_ok "DeepImmuno git clone"
  if [[ -n "$DEEPIMMUNO_SOURCE" ]]; then
    run "install DeepImmuno from local source" bash scripts/install_deepimmuno.sh "$DEEPIMMUNO_SOURCE"
  else
    run "install DeepImmuno from official repo" bash scripts/install_deepimmuno.sh "$TOOLS_ROOT/tools/DeepImmuno"
  fi
fi
[[ "$INSTALL_LOHHLA" == "1" ]] && { need_download_ok "LOHHLA git clone"; run "install LOHHLA" bash scripts/install_lohhla.sh; }
if [[ "$INSTALL_POLYSOLVER" == "1" ]]; then
  [[ -n "$POLYSOLVER_HOME_ARG" ]] || { echo "POLYSOLVER_HOME_REQUIRED: pass --polysolver-home" >&2; exit 32; }
  env_cmd="POLYSOLVER_HOME='$POLYSOLVER_HOME_ARG'"
  [[ -n "$NOVOALIGN_LICENSE_FILE_ARG" ]] && env_cmd="$env_cmd NOVOALIGN_LICENSE_FILE='$NOVOALIGN_LICENSE_FILE_ARG'"
  run "configure Polysolver" bash -lc "$env_cmd bash scripts/install_polysolver.sh"
fi
[[ "$INSTALL_OPTITYPE" == "1" ]] && run "install OptiType" bash scripts/install_optitype.sh
[[ "$INSTALL_FACETS" == "1" ]] && run "install FACETS" bash scripts/install_facets.sh
[[ "$INSTALL_ASCAT_PYCLONE" == "1" ]] && run "install ASCAT/PyClone-VI" bash scripts/install_ascat_pyclone.sh
[[ "$INSTALL_FUSION" == "1" ]] && { need_download_ok "fusion tool git clones/conda packages"; run "install fusion tools" bash scripts/install_fusion_tools.sh; }
sync_assets_if_requested
stage_bigmhc_models_if_requested

if [[ "$RUN_VERIFY" == "1" ]]; then
  verify_args=(--smoke "$REFERENCE_ROOT")
  [[ "$STRICT_VERIFY" == "1" ]] && verify_args=(--smoke --strict "$REFERENCE_ROOT")
  run "verify README tools and references" bash scripts/verify_all_tools_and_refs.sh "${verify_args[@]}"
fi

if [[ "$RUN_REAL_VCF_SMOKE" == "1" ]]; then
  real_vcf_outdir="$OUTDIR/real_vcf_smoke"
  real_vcf_args=(
    --project-root "$PROJECT_ROOT"
    --tools-root "$TOOLS_ROOT"
    --licensed-root "$LICENSED_ROOT"
    --conda-base "${CONDA_BASE:-$TOOLS_ROOT/miniforge3}"
    --outdir "$real_vcf_outdir"
    --top-n "$REAL_VCF_SMOKE_TOP_N"
  )
  [[ -n "$REAL_VCF_RAW" ]] && real_vcf_args+=(--raw-vcf "$REAL_VCF_RAW")
  [[ -n "$REAL_VCF_ANNOTATED" ]] && real_vcf_args+=(--annotated-vcf "$REAL_VCF_ANNOTATED")
  [[ -n "$REAL_VCF_HLA_ALLELES" ]] && real_vcf_args+=(--hla-alleles "$REAL_VCF_HLA_ALLELES")
  [[ -n "$REAL_VCF_HLA_FILE" ]] && real_vcf_args+=(--hla-file "$REAL_VCF_HLA_FILE")
  [[ "$REAL_VCF_SMOKE_SKIP_MHCFLURRY" == "1" ]] && real_vcf_args+=(--skip-mhcflurry)
  run "run default real VCF smoke test" bash .agents/skills/neoag-remote-deploy/scripts/14_run_real_vcf_smoke.sh "${real_vcf_args[@]}"
fi

{
  echo "# README tool install report"
  echo
  echo "Mode: \`$MODE\`"
  echo "Project root: \`$PROJECT_ROOT\`"
  echo "Tools root: \`$TOOLS_ROOT\`"
  echo "Licensed root: \`$LICENSED_ROOT\`"
  echo "Reference root: \`$REFERENCE_ROOT\`"
  echo "Conda base: \`${CONDA_BASE:-$TOOLS_ROOT/miniforge3}\`"
  echo "Log: \`$LOG\`"
  echo
  echo "Selected groups:"
  for item in \
    "core-env:$INSTALL_CORE_ENV" "core-env-lite:$CORE_ENV_LITE" "skip-torch-install:$SKIP_TORCH_INSTALL" "vep:$INSTALL_VEP" "vep-cache:$INSTALL_VEP_CACHE" "vep-version:$VEP_VERSION" \
    "gatk:$INSTALL_GATK" "immunogenicity:$INSTALL_IMMUNOGENICITY" \
    "netmhcstabpan:$INSTALL_NETMHCSTABPAN" "deepimmuno:$INSTALL_DEEPIMMUNO" \
    "lohhla:$INSTALL_LOHHLA" "polysolver:$INSTALL_POLYSOLVER" "optitype:$INSTALL_OPTITYPE" \
    "facets:$INSTALL_FACETS" "ascat-pyclone:$INSTALL_ASCAT_PYCLONE" "fusion:$INSTALL_FUSION" \
    "verify:$RUN_VERIFY" "real-vcf-smoke:$RUN_REAL_VCF_SMOKE" "sync-assets:$SYNC_ASSETS" "reference-manifest:${REFERENCE_MANIFEST:+1}" "bigmhc-models:${BIGMHC_MODELS_DIR:+1}"; do
    name="${item%%:*}"; enabled="${item##*:}"
    [[ "$enabled" == "1" ]] && echo "- $name"
  done
  if [[ "$RUN_REAL_VCF_SMOKE" == "1" ]]; then
    echo "- real-vcf-smoke-mhcflurry-default-on"
    [[ "$REAL_VCF_SMOKE_SKIP_MHCFLURRY" == "1" ]] && echo "- real-vcf-smoke-mhcflurry-skipped"
  fi
  echo
  if [[ "$EXECUTE" != "1" ]]; then
    echo "Dry run only. Re-run with \`--execute\`; add \`--allow-download\` for network downloads after approval."
  else
    echo "Install commands completed. Review the log and run production validation before real data."
  fi
} > "$REPORT"

log ""
log "readme_tools_install_report=$REPORT"
