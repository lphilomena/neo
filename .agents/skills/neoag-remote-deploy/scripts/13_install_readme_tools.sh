#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(pwd)"
DEPLOY_ROOT="${NEOAG_DEPLOY_ROOT:-/opt/neoag}"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-$DEPLOY_ROOT/env_tool}"
LICENSED_ROOT="${NEOAG_LICENSED_ROOT:-$DEPLOY_ROOT/licensed_tools}"
REFERENCE_ROOT="${NEOAG_REFERENCE_ROOT:-$DEPLOY_ROOT/refs}"
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
INSTALL_RNA_EXPRESSION=0
INSTALL_IMMUNOGENICITY=0
INSTALL_DEEPIMMUNO=0
INSTALL_SHERPA=0
INSTALL_NETMHCSTABPAN=0
INSTALL_NETCHOP=0
INSTALL_LOHHLA=0
INSTALL_POLYSOLVER=0
INSTALL_OPTITYPE=0
INSTALL_BAM_MATCHER=0
INSTALL_FACETS=0
INSTALL_ASCAT_PYCLONE=0
INSTALL_FUSION=0
INSTALL_SPLICE=0
INSTALL_SNAF=1
INSTALL_SPLICEMUTR=1
INSTALL_SPECHLA=0
INSTALL_HLALA=0
HLALA_VERSION="1.0.4"
INSTALL_SEQUENZA=0
INSTALL_HMF_PURPLE=0
INSTALL_CLAUDE_CODE=0
CLAUDE_CODE_CHANNEL="stable"
CLAUDE_CODE_INSTALLER_URL="https://claude.ai/install.sh"
RUN_VERIFY=0
STRICT_VERIFY=0
RUN_REAL_VCF_SMOKE=0
REAL_VCF_SMOKE_TOP_N=50
REAL_VCF_SMOKE_SKIP_MHCFLURRY=0
REAL_VCF_SMOKE_SKIP_BIGMHC=0
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
ASSET_SSH_KEY="${NEOAG_ASSET_SSH_KEY:-}"
SHARED_ASSET_ROOT="${NEOAG_SHARED_ASSET_ROOT:-}"
CORE_ENV_LITE=1
SKIP_TORCH_INSTALL=1
TORCH_WHEEL_DIR="${TORCH_WHEEL_DIR:-}"
TORCH_INDEX_URL="${TORCH_INDEX_URL:-https://download.pytorch.org/whl/cpu}"

NETMHCPAN_TAR=""
NETMHCPAN_DIR=""
NETMHCPAN_URL=""
MIXMHCPRED_DIR=""
MIXMHCPRED_ARCHIVE=""
MIXMHCPRED_URL=""
NETMHCSTABPAN_DIR=""
NETMHCSTABPAN_ARCHIVE=""
NETMHCSTABPAN_URL=""
NETCHOP_ARCHIVE="${NEOAG_NETCHOP_ARCHIVE:-}"
POLYSOLVER_HOME_ARG=""
NOVOALIGN_LICENSE_FILE_ARG=""
DEEPIMMUNO_SOURCE=""
SHERPA_SOURCE=""
SHERPA_ARCHIVE=""
SHERPA_CONTAINER_IMAGE=""
SHERPA_SMOKE_COMMAND=""
SPECHLA_SOURCE=""

usage() {
  cat <<'USAGE'
Usage: 13_install_readme_tools.sh [options]

Install the external tools listed in README.md by orchestrating the repository's
existing installation scripts. The default mode is dry-run; add --execute to
actually install. Licensed tools still require local archives/directories or an
explicit approved URL plus --allow-download.

Common options:
  --project-root DIR          Project checkout (default: current directory)
  --tools-root DIR            Tool/env root (default: NEOAG_TOOLS_ROOT or /opt/neoag/env_tool)
  --licensed-root DIR         Licensed tool root (default: NEOAG_LICENSED_ROOT or /opt/neoag/licensed_tools)
  --reference-root DIR        Reference root (default: NEOAG_REFERENCE_ROOT or /opt/neoag/refs)
  --conda-base DIR            Miniforge/conda base (default: tools-root/miniforge3)
  --install-miniforge         Install Miniforge3 if conda is missing (default enabled)
  --no-install-miniforge      Do not install Miniforge automatically if conda is missing
  --miniforge-url URL         Miniforge installer URL
  --outdir DIR                Log/report directory (default: work/remote_deploy)
  --allow-download            Permit downloads from official/user-approved URLs
  --execute                   Actually run installation commands
  --full-core-env             Use full core env instead of default lite env
  --install-torch             Let immunogenicity installer install torch if missing
  --skip-torch-install        Skip torch even when installing BigMHC (BigMHC smoke will be partial)
  --torch-wheel-dir DIR       Optional local torch/nvidia wheel directory for offline BigMHC runtime repair
  --torch-index-url URL       PyTorch package index when torch download is approved (default: CPU wheel index)

Tool groups:
  --core-env                  pVACtools/MHCflurry core conda env via scripts/setup_tools_env.sh
  --vep                      VEP conda env via scripts/install_vep.sh
  --vep-cache                VEP cache via scripts/install_vep_cache.sh (large download)
  --vep-version VERSION      Ensembl VEP/cache release to install/use (default: 105)
  --gatk                     GATK4 / Mutect2 via scripts/install_gatk.sh
  --rna-expression           Install Salmon/RSEM in neoag-tools for RNA FASTQ to gene TPM scripts
  --immunogenicity           PRIME + MixMHCpred + BigMHC via scripts/install_immunogenicity_tools.sh
  --deepimmuno               DeepImmuno via scripts/install_deepimmuno.sh
  --sherpa                   Register/install SHERPA-Presentation from a local authorized source/archive/container
  --netmhcstabpan            NetMHCstabpan or IEDB shim via scripts/install_netmhcstabpan.sh (included by --all/--all-open)
  --netchop                  NetChop 3.1d from an authorized local/synchronized archive (included by --all/--all-open)
  --lohhla                   LOHHLA source wrapper via scripts/install_lohhla.sh
  --polysolver               Configure existing Polysolver; requires --polysolver-home
  --optitype                 OptiType via scripts/install_optitype.sh
  --bam-matcher              BAM-matcher plus Python 2 compatibility env via scripts/install_bam_matcher.sh
  --facets                   FACETS via scripts/install_facets.sh
  --ascat-pyclone            ASCAT + PyClone-VI via scripts/install_ascat_pyclone.sh
  --fusion                   Arriba/Nextflow fusion env plus STAR-Fusion/FusionCatcher clones
  --splice                   RegTools + pVACsplice + SNAF + SpliceMutr (defaults) and optional splice tools
  --skip-snaf                Skip SNAF when installing the splice group
  --skip-splicemutr          Skip SpliceMutr when installing the splice group
  --spechla                  Register/load SpecHLA container assets and database if present
  --spechla-source DIR       Existing complete SpecHLA source tree; required when the staged image has runtime only
  --hla-la                   Register/load HLA-LA container assets and PRG graph if present
  --hla-la-version VERSION   Bioconda HLA-LA version (default: 1.0.4)
  --sequenza                 Install Sequenza conda env and reference hooks
  --hmf-purple               Register/load HMF PURPLE/AMBER/COBALT container assets and references
  --claude-code              Install Claude Code with Anthropic's official native installer
  --claude-code-channel V    stable, latest, or exact X.Y.Z version (default: stable)
  --claude-code-installer-url URL
                              Override only with an explicitly approved official URL
  --all-open                 Install open/conda/git tools except very large VEP cache; registers supplied NetMHCstabpan/NetChop assets when present
  --all                      Install supported default groups, including VEP cache and NetMHCstabpan
  --verify                   Run scripts/verify_all_tools_and_refs.sh after installs
  --strict-verify            Treat optional missing tools as failure during verify
  --run-real-vcf-smoke       Run an explicitly configured real VCF top-N smoke test
  --real-vcf-smoke-top-n N   Number of unique peptides for real VCF smoke (default: 50)
  --skip-real-vcf-mhcflurry Skip MHCflurry only for the real VCF smoke fallback
  --skip-real-vcf-bigmhc    Skip BigMHC_IM only for the real VCF smoke fallback
  --real-vcf FILE           Raw VCF path required for real VCF smoke
  --real-annotated-vcf FILE VEP-annotated VCF required for real VCF smoke
  --real-vcf-hla-alleles L  Comma-separated HLA alleles for real VCF smoke
  --real-vcf-hla-file FILE  HLA file for real VCF smoke
  --bigmhc-models-dir DIR   Copy BigMHC models from local/source directory into tools-root
  --bigmhc-models-host HOST Optional source host for --bigmhc-models-dir, e.g. user@source-host
  --asset-manifest FILE    TSV manifest for large assets (default: configs/assets/production_assets.tsv)
  --reference-manifest FILE
                          YAML reference manifest verified after asset sync
  --sync-assets            Sync large assets from manifest (dry-run unless --execute)
  --asset-source-host HOST Default source host for manifest source_path values
  --asset-ssh-key FILE    SSH private key used by rsync for remote assets
  --shared-asset-root DIR Link assets from a locally mounted shared root

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
  --netchop-archive FILE     Authorized netchop-3.1d.Linux.tar.gz; if omitted, search licensed-root/netchop and shared-asset-root
  --polysolver-home DIR      Existing Polysolver distribution
  --novoalign-license-file FILE
  --deepimmuno-source DIR    Existing DeepImmuno checkout; otherwise script clones official repo
  --sherpa-source DIR        Existing SHERPA-Presentation checkout/directory to copy/register
  --sherpa-archive FILE      Existing SHERPA-Presentation tar/zip/wheel archive
  --sherpa-container-image FILE
                              Existing SHERPA-Presentation container image tarball
  --sherpa-smoke-command CMD Optional command to validate the registered SHERPA-Presentation install

Examples:
  bash .agents/skills/neoag-remote-deploy/scripts/13_install_readme_tools.sh \
    --project-root <project-root> \
    --tools-root /opt/neoag/env_tool \
    --conda-base /opt/neoag/env_tool/miniforge3 \
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
    --skip-torch-install) SKIP_TORCH_INSTALL=1; shift ;;
    --torch-wheel-dir) TORCH_WHEEL_DIR="$2"; shift 2 ;;
    --torch-index-url) TORCH_INDEX_URL="$2"; shift 2 ;;
    --core-env) INSTALL_CORE_ENV=1; shift ;;
    --vep) INSTALL_VEP=1; shift ;;
    --vep-cache) INSTALL_VEP=1; INSTALL_VEP_CACHE=1; shift ;;
    --vep-version) VEP_VERSION="$2"; shift 2 ;;
    --gatk) INSTALL_GATK=1; shift ;;
    --rna-expression) INSTALL_RNA_EXPRESSION=1; INSTALL_CORE_ENV=1; shift ;;
    --immunogenicity) INSTALL_IMMUNOGENICITY=1; shift ;;
    --deepimmuno) INSTALL_DEEPIMMUNO=1; shift ;;
    --sherpa) INSTALL_SHERPA=1; shift ;;
    --netmhcstabpan) INSTALL_NETMHCSTABPAN=1; shift ;;
    --netchop) INSTALL_NETCHOP=1; shift ;;
    --lohhla) INSTALL_LOHHLA=1; shift ;;
    --polysolver) INSTALL_POLYSOLVER=1; shift ;;
    --optitype) INSTALL_OPTITYPE=1; shift ;;
    --bam-matcher) INSTALL_BAM_MATCHER=1; shift ;;
    --facets) INSTALL_FACETS=1; shift ;;
    --ascat-pyclone) INSTALL_ASCAT_PYCLONE=1; shift ;;
    --fusion) INSTALL_FUSION=1; shift ;;
    --splice) INSTALL_SPLICE=1; shift ;;
    --skip-snaf) INSTALL_SNAF=0; shift ;;
    --skip-splicemutr) INSTALL_SPLICEMUTR=0; shift ;;
    --spechla) INSTALL_SPECHLA=1; shift ;;
    --spechla-source) SPECHLA_SOURCE="$2"; INSTALL_SPECHLA=1; shift 2 ;;
    --hla-la) INSTALL_HLALA=1; shift ;;
    --hla-la-version) HLALA_VERSION="$2"; INSTALL_HLALA=1; shift 2 ;;
    --sequenza) INSTALL_SEQUENZA=1; shift ;;
    --hmf-purple) INSTALL_HMF_PURPLE=1; shift ;;
    --claude-code) INSTALL_CLAUDE_CODE=1; shift ;;
    --claude-code-channel) CLAUDE_CODE_CHANNEL="$2"; INSTALL_CLAUDE_CODE=1; shift 2 ;;
    --claude-code-installer-url) CLAUDE_CODE_INSTALLER_URL="$2"; INSTALL_CLAUDE_CODE=1; shift 2 ;;
    --all-open)
      INSTALL_CORE_ENV=1; INSTALL_VEP=1; INSTALL_GATK=1; INSTALL_RNA_EXPRESSION=1; INSTALL_IMMUNOGENICITY=1
      INSTALL_DEEPIMMUNO=1; INSTALL_NETMHCSTABPAN=1; INSTALL_NETCHOP=1; INSTALL_LOHHLA=1
      INSTALL_OPTITYPE=1; INSTALL_BAM_MATCHER=1; INSTALL_FACETS=1; INSTALL_ASCAT_PYCLONE=1; INSTALL_FUSION=1; INSTALL_SPLICE=1
      INSTALL_SPECHLA=1; INSTALL_HLALA=1; INSTALL_SEQUENZA=1; INSTALL_HMF_PURPLE=1
      SKIP_TORCH_INSTALL=0
      shift ;;
    --all)
      INSTALL_CORE_ENV=1; INSTALL_VEP=1; INSTALL_VEP_CACHE=1; INSTALL_GATK=1; INSTALL_RNA_EXPRESSION=1; INSTALL_IMMUNOGENICITY=1
      INSTALL_DEEPIMMUNO=1; INSTALL_NETMHCSTABPAN=1; INSTALL_NETCHOP=1; INSTALL_LOHHLA=1
      INSTALL_OPTITYPE=1; INSTALL_BAM_MATCHER=1; INSTALL_FACETS=1; INSTALL_ASCAT_PYCLONE=1; INSTALL_FUSION=1; INSTALL_SPLICE=1
      INSTALL_SPECHLA=1; INSTALL_HLALA=1; INSTALL_SEQUENZA=1; INSTALL_HMF_PURPLE=1
      SKIP_TORCH_INSTALL=0
      shift ;;
    --verify) RUN_VERIFY=1; shift ;;
    --strict-verify) RUN_VERIFY=1; STRICT_VERIFY=1; shift ;;
    --run-real-vcf-smoke) RUN_REAL_VCF_SMOKE=1; shift ;;
    --real-vcf-smoke-top-n) REAL_VCF_SMOKE_TOP_N="$2"; shift 2 ;;
    --skip-real-vcf-mhcflurry) REAL_VCF_SMOKE_SKIP_MHCFLURRY=1; shift ;;
    --skip-real-vcf-bigmhc) REAL_VCF_SMOKE_SKIP_BIGMHC=1; shift ;;
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
    --asset-ssh-key) ASSET_SSH_KEY="$2"; shift 2 ;;
    --shared-asset-root) SHARED_ASSET_ROOT="$2"; shift 2 ;;
    --netmhcpan-tar) NETMHCPAN_TAR="$2"; shift 2 ;;
    --netmhcpan-dir) NETMHCPAN_DIR="$2"; shift 2 ;;
    --netmhcpan-url) NETMHCPAN_URL="$2"; shift 2 ;;
    --mixmhcpred-dir) MIXMHCPRED_DIR="$2"; shift 2 ;;
    --mixmhcpred-archive) MIXMHCPRED_ARCHIVE="$2"; shift 2 ;;
    --mixmhcpred-url) MIXMHCPRED_URL="$2"; shift 2 ;;
    --netmhcstabpan-dir) NETMHCSTABPAN_DIR="$2"; shift 2 ;;
    --netmhcstabpan-archive) NETMHCSTABPAN_ARCHIVE="$2"; shift 2 ;;
    --netmhcstabpan-url) NETMHCSTABPAN_URL="$2"; shift 2 ;;
    --netchop-archive) NETCHOP_ARCHIVE="$2"; INSTALL_NETCHOP=1; shift 2 ;;
    --polysolver-home) POLYSOLVER_HOME_ARG="$2"; shift 2 ;;
    --novoalign-license-file) NOVOALIGN_LICENSE_FILE_ARG="$2"; shift 2 ;;
    --deepimmuno-source) DEEPIMMUNO_SOURCE="$2"; shift 2 ;;
    --sherpa-source) SHERPA_SOURCE="$2"; shift 2 ;;
    --sherpa-archive) SHERPA_ARCHIVE="$2"; shift 2 ;;
    --sherpa-container-image) SHERPA_CONTAINER_IMAGE="$2"; shift 2 ;;
    --sherpa-smoke-command) SHERPA_SMOKE_COMMAND="$2"; shift 2 ;;
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
  if [[ "$EXECUTE" != "1" ]]; then
    log "[DRY_RUN] download approval required during execution: $what"
    return 0
  fi
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
    local samtools_bin=""
    samtools_bin="$(command -v samtools 2>/dev/null || true)"
    if [[ -z "$samtools_bin" && -n "$CONDA_BASE" && -x "$CONDA_BASE/envs/neoag-tools/bin/samtools" ]]; then
      samtools_bin="$CONDA_BASE/envs/neoag-tools/bin/samtools"
    fi
    if [[ -z "$samtools_bin" && -x "$TOOLS_ROOT/miniforge3/envs/neoag-tools/bin/samtools" ]]; then
      samtools_bin="$TOOLS_ROOT/miniforge3/envs/neoag-tools/bin/samtools"
    fi
    if [[ -n "$samtools_bin" ]]; then
      run "index reference FASTA" "$samtools_bin" faidx "$fasta"
    else
      log "WARN: samtools not found; cannot create FASTA index: $fasta.fai"
    fi
  fi
}

register_synced_tool_assets() {
  [[ "$EXECUTE" == "1" ]] || return 0

  local netmhcpan_image="neoag-netmhcpan:4.2c-ubuntu22.04"
  local netmhcpan_tar="$TOOLS_ROOT/container_images/neoag-netmhcpan_4.2c-ubuntu22.04.tar"
  local netmhcpan_wrapper="$TOOLS_ROOT/bin/netMHCpan"

  if [[ -f "$TOOLS_ROOT/tools/DeepImmuno/deepimmuno-cnn.py" ]]; then
    run "register synced DeepImmuno asset" bash scripts/install_deepimmuno.sh "$TOOLS_ROOT/tools/DeepImmuno"
  fi

  if [[ -x "$LICENSED_ROOT/netMHCpan/netMHCpan" ]]; then
    if [[ -s "$netmhcpan_tar" ]]; then
      mkdir -p "$TOOLS_ROOT/bin"
      load_container_image_if_present "NetMHCpan" "$netmhcpan_image" "$netmhcpan_tar"
      cat > "$netmhcpan_wrapper" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export NETMHCPAN_HOME="\${NETMHCPAN_HOME:-$LICENSED_ROOT/netMHCpan}"
export NEOAG_NETMHCPAN_IMAGE="\${NEOAG_NETMHCPAN_IMAGE:-$netmhcpan_image}"
exec "$PROJECT_ROOT/scripts/run_netmhcpan_container.sh" "\$@"
EOF
      chmod +x "$netmhcpan_wrapper"
      log "NetMHCpan container wrapper registered: $netmhcpan_wrapper"
    fi
    if [[ -L "$LICENSED_ROOT/netMHCpan" ]]; then
      log "NetMHCpan licensed asset is a shared symlink; skipping in-place native repair"
    elif [[ -x "${CONDA_BASE:-$TOOLS_ROOT/miniforge3}/envs/neoag-tools/bin/patchelf" || -x "$LICENSED_ROOT/netMHCpan/Linux_x86_64/bin/netMHCpan-4.2" ]]; then
      if ! run "repair/register synced NetMHCpan asset" env NETMHCPAN_HOME="$LICENSED_ROOT/netMHCpan" NEOAG_CONDA_BASE="${CONDA_BASE:-$TOOLS_ROOT/miniforge3}" bash scripts/install_netmhcpan.sh --repair; then
        log "WARN: NetMHCpan asset is present but repair/smoke failed; license asset was left in place for manual validation."
      fi
    fi
  fi

  if [[ -x "$LICENSED_ROOT/netMHCstabpan/netMHCstabpan" ]]; then
    log "NetMHCstabpan asset present: $LICENSED_ROOT/netMHCstabpan/netMHCstabpan"
    write_presentation_tool_env_overrides
  fi

  if [[ -d "$LICENSED_ROOT/polysolver" && -z "$POLYSOLVER_HOME_ARG" ]]; then
    POLYSOLVER_HOME_ARG="$LICENSED_ROOT/polysolver"
  fi
  if [[ -f "$LICENSED_ROOT/novoalign/novoalign.lic" && -z "$NOVOALIGN_LICENSE_FILE_ARG" ]]; then
    NOVOALIGN_LICENSE_FILE_ARG="$LICENSED_ROOT/novoalign/novoalign.lic"
  fi
}

sync_assets_if_requested() {
  [[ "$SYNC_ASSETS" == "1" ]] || return 0
  args=(--project-root "$PROJECT_ROOT" --asset-manifest "$ASSET_MANIFEST" --outdir "$OUTDIR/assets"
    --tools-root "$TOOLS_ROOT" --reference-root "$REFERENCE_ROOT" --licensed-root "$LICENSED_ROOT")
  [[ -n "$ASSET_SOURCE_HOST" ]] && args+=(--asset-source-host "$ASSET_SOURCE_HOST")
  [[ -n "$ASSET_SSH_KEY" ]] && args+=(--asset-ssh-key "$ASSET_SSH_KEY")
  [[ -n "$SHARED_ASSET_ROOT" ]] && args+=(--shared-asset-root "$SHARED_ASSET_ROOT")
  [[ "$EXECUTE" == "1" ]] && args+=(--execute)
  run "sync large assets from manifest" bash .agents/skills/neoag-remote-deploy/scripts/15_sync_asset_manifest.sh "${args[@]}"
  register_synced_tool_assets
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

load_container_image_if_present() {
  local label="$1" image="$2" tarball="$3"
  [[ -s "$tarball" ]] || { log "WARN: $label image tar missing: $tarball"; return 0; }
  if ! command -v docker >/dev/null 2>&1; then
    log "WARN: docker not found; cannot load $label image from $tarball"
    return 0
  fi
  if docker image inspect "$image" >/dev/null 2>&1; then
    log "$label image already loaded: $image"
    return 0
  fi
  run "load $label container image" docker load -i "$tarball"
}

write_presentation_tool_env_overrides() {
  local stab_home="$LICENSED_ROOT/netMHCstabpan"
  local netchop_home="$LICENSED_ROOT/netchop/netchop-3.1"
  local netchop_bin="$TOOLS_ROOT/bin/netChop"
  [[ -x "$stab_home/netMHCstabpan" || -x "$netchop_bin" ]] || return 0
  local local_env="$PROJECT_ROOT/conf/tools.env.local.sh"
  local start="# BEGIN NEOAG LICENSED PRESENTATION TOOLS"
  local end="# END NEOAG LICENSED PRESENTATION TOOLS"
  mkdir -p "$PROJECT_ROOT/conf"
  if [[ -f "$local_env" ]]; then
    sed "/$start/,/$end/d" "$local_env" > "${local_env}.tmp"
    mv "${local_env}.tmp" "$local_env"
  fi
  cat >> "$local_env" <<EOF
$start
export NETMHCSTABPAN_HOME="$stab_home"
export NEOAG_NETMHCSTABPAN_BIN="\${NETMHCSTABPAN_HOME}/netMHCstabpan"
if [[ -x "\${NEOAG_NETMHCSTABPAN_BIN}" ]]; then
  export PATH="\${NETMHCSTABPAN_HOME}:\${PATH}"
fi
export NETCHOP_HOME="$netchop_home"
export NETCHOP="\${NETCHOP_HOME}/Linux_x86_64"
export NEOAG_NETCHOP_BIN="$netchop_bin"
export NETCHOP_BIN="\${NEOAG_NETCHOP_BIN}"
if [[ -x "\${NEOAG_NETCHOP_BIN}" ]]; then
  export PATH="$(dirname "$netchop_bin"):\${PATH}"
fi
$end
EOF
  log "Presentation tool runtime overrides written: $local_env"
}

discover_netchop_archive() {
  local candidate=""
  for candidate in \
    "$LICENSED_ROOT/netchop/netchop-3.1d.Linux.tar.gz" \
    "$TOOLS_ROOT/netchop-3.1d.Linux.tar.gz" \
    "$TOOLS_ROOT/tools/netchop-3.1d.Linux.tar.gz"; do
    [[ -f "$candidate" ]] && { printf '%s\n' "$candidate"; return 0; }
  done
  if [[ -n "$SHARED_ASSET_ROOT" && -d "$SHARED_ASSET_ROOT" ]]; then
    find "$SHARED_ASSET_ROOT" -maxdepth 6 -type f -iname 'netchop-3.1d.Linux.tar.gz' -print -quit 2>/dev/null
  fi
}

register_netmhcstabpan_if_requested() {
  [[ "$INSTALL_NETMHCSTABPAN" == "1" ]] || return 0
  local home="$LICENSED_ROOT/netMHCstabpan"
  local image="neoag-netmhcstabpan:1.0-ubuntu22.04"
  local image_tar="$TOOLS_ROOT/container_images/neoag-netmhcstabpan_1.0-ubuntu22.04.tar"
  local wrapper="$TOOLS_ROOT/bin/netMHCstabpan"

  if [[ -x "$home/netMHCstabpan" && -s "$image_tar" ]]; then
    if [[ "$EXECUTE" != "1" ]]; then
      log ""
      log "==> [DRY_RUN] register licensed NetMHCstabpan container wrapper"
      log "+ load $image_tar and create $wrapper"
      return 0
    fi
    mkdir -p "$TOOLS_ROOT/bin"
    load_container_image_if_present "NetMHCstabpan" "$image" "$image_tar"
    cat > "$wrapper" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export NETMHCSTABPAN_HOME="\${NETMHCSTABPAN_HOME:-$home}"
export NEOAG_NETMHCSTABPAN_IMAGE="\${NEOAG_NETMHCSTABPAN_IMAGE:-$image}"
exec "$PROJECT_ROOT/scripts/run_netmhcstabpan_container.sh" "\$@"
EOF
    chmod +x "$wrapper"
    log "NetMHCstabpan licensed asset and container wrapper registered: $wrapper"
    return 0
  fi

  if [[ -n "$NETMHCSTABPAN_ARCHIVE" && -f "$NETMHCSTABPAN_ARCHIVE" ]]; then
    run "install NetMHCstabpan from archive" bash scripts/install_netmhcstabpan.sh "$NETMHCSTABPAN_ARCHIVE"
  else
    run "install NetMHCstabpan IEDB shim" bash scripts/install_netmhcstabpan.sh --iedb
  fi
}

register_spechla_if_requested() {
  [[ "$INSTALL_SPECHLA" == "1" ]] || return 0
  local home="$TOOLS_ROOT/tools/SpecHLA"
  local image_tar="$TOOLS_ROOT/container_images/neoag-spechla_ubuntu22.04.tar"
  local db="$REFERENCE_ROOT/data/hla/spechla/db"
  [[ -d "$db" ]] || db="$REFERENCE_ROOT/data/hla/spechla_db"
  if [[ "$EXECUTE" != "1" ]]; then
    log ""
    log "==> [DRY_RUN] register SpecHLA container wrappers and DB link"
    log "+ install complete SpecHLA source, repair $home/db -> $db, load/build the dependency-complete image, and create missing novoalign indexes"
    return 0
  fi
  mkdir -p "$home" "$TOOLS_ROOT/bin"
  if [[ -n "$SPECHLA_SOURCE" ]]; then
    [[ -f "$SPECHLA_SOURCE/script/whole/SpecHLA.sh" && -f "$SPECHLA_SOURCE/script/cal.hla.copy.pl" ]] || {
      echo "SPECHLA_SOURCE_INVALID: expected script/whole/SpecHLA.sh and script/cal.hla.copy.pl below $SPECHLA_SOURCE" >&2
      return 1
    }
    rsync -a --exclude '.git/' --exclude 'db/' --exclude 'spechla_env/' "$SPECHLA_SOURCE/" "$home/"
  fi
  [[ -d "$db" ]] || { echo "SPECHLA_DB_MISSING: expected $db" >&2; return 1; }
  if [[ -L "$home/db" ]]; then
    ln -sfn "$db" "$home/db"
  elif [[ ! -e "$home/db" ]]; then
    ln -s "$db" "$home/db"
  elif [[ "$(cd "$home/db" && pwd -P)" != "$(cd "$db" && pwd -P)" ]]; then
    echo "SPECHLA_DB_CONFLICT: $home/db is a real directory and does not resolve to $db" >&2
    return 1
  fi
  [[ -f "$home/script/whole/SpecHLA.sh" && -f "$home/script/cal.hla.copy.pl" ]] || {
    echo "SPECHLA_SOURCE_REQUIRED: provide --spechla-source with a complete official SpecHLA checkout; a database and runtime-only image are insufficient" >&2
    return 1
  }
  if [[ ! -x "$TOOLS_ROOT/bin/SpecHLA" ]]; then
    cat > "$TOOLS_ROOT/bin/SpecHLA" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export SPECHLA_HOME="\${SPECHLA_HOME:-$home}"
exec "$PROJECT_ROOT/scripts/run_spechla_container.sh" "\$@"
EOF
    chmod +x "$TOOLS_ROOT/bin/SpecHLA"
  fi
  load_container_image_if_present "SpecHLA" "neoag-spechla:ubuntu22.04" "$image_tar"
  if command -v docker >/dev/null 2>&1 && docker image inspect neoag-spechla:ubuntu22.04 >/dev/null 2>&1; then
    if ! SPECHLA_HOME="$home" SPECHLA_DB="$db" SPECHLA_ENV="$CONDA_BASE/envs/neoag-tools" \
      SPECHLA_CMD=python3 scripts/run_spechla_container.sh -c \
      'import Bio, networkx, numpy, pandas, pulp, pyfaidx, pysam, scipy, vcf' >/dev/null 2>&1 || \
      ! SPECHLA_HOME="$home" SPECHLA_DB="$db" SPECHLA_ENV="$CONDA_BASE/envs/neoag-tools" \
      SPECHLA_CMD=freebayes scripts/run_spechla_container.sh --version >/dev/null 2>&1; then
      need_download_ok "SpecHLA container runtime dependencies"
      run "rebuild SpecHLA dependency-complete container" bash scripts/build_priority_tool_containers.sh spechla
    fi
  fi

  local resolved_db
  resolved_db=$(cd "$home/db" && pwd -P)
  if [[ -x "$home/bin/novoalign" && -f "$home/bin/novoalign.lic" ]]; then
    [[ -x "$home/bin/novoindex" ]] || { echo "SPECHLA_NOVOINDEX_MISSING: $home/bin/novoindex" >&2; return 1; }
    local prefix fasta index tmp
    for prefix in hla_gen.format.filter.extend.DRB.no26789 hla_gen.format.filter.extend.DRB.no26789.v2; do
      fasta="$resolved_db/ref/$prefix.fasta"
      index="$resolved_db/ref/$prefix.ndx"
      if [[ -s "$fasta" && ! -s "$index" ]]; then
        tmp="$index.tmp.$$"
        "$home/bin/novoindex" "$tmp" "$fasta"
        mv "$tmp" "$index"
      fi
    done
  fi
}

register_hlala_if_requested() {
  [[ "$INSTALL_HLALA" == "1" ]] || return 0
  local home="$TOOLS_ROOT/tools/HLA-LA"
  local env_prefix="$home/.conda"
  local image_tar="$TOOLS_ROOT/container_images/neoag-hla-la_ubuntu22.04.tar"
  if [[ "$EXECUTE" != "1" ]]; then
    log ""
    log "==> [DRY_RUN] install/register HLA-LA and graph path"
    log "+ install hla-la=$HLALA_VERSION at $env_prefix, link real HLA-LA.pl, load $image_tar if present"
    return 0
  fi
  mkdir -p "$home/bin" "$TOOLS_ROOT/bin"
  if [[ ! -x "$env_prefix/bin/HLA-LA.pl" ]]; then
    need_download_ok "HLA-LA Bioconda package"
    "$CONDA_BASE/bin/mamba" create -y -p "$env_prefix" --override-channels -c conda-forge -c bioconda "hla-la=$HLALA_VERSION"
  fi
  if [[ -f "$home/bin/HLA-LA.pl" && ! -L "$home/bin/HLA-LA.pl" ]] && grep -q 'run_hla_la_container.sh' "$home/bin/HLA-LA.pl"; then
    mv "$home/bin/HLA-LA.pl" "$home/bin/HLA-LA.pl.legacy_wrapper"
  fi
  ln -sfn "$env_prefix/bin/HLA-LA.pl" "$home/bin/HLA-LA.pl"
  cat > "$TOOLS_ROOT/bin/HLA-LA.pl" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export HLALA_HOME="\${HLALA_HOME:-$home}"
export HLALA_ENV_PREFIX="\${HLALA_ENV_PREFIX:-$env_prefix}"
export HLALA_GRAPH="\${HLALA_GRAPH:-$REFERENCE_ROOT/data/hla/PRG_MHC_GRCh38_withIMGT}"
exec "$PROJECT_ROOT/scripts/run_hla_la_container.sh" "\$@"
EOF
  chmod +x "$TOOLS_ROOT/bin/HLA-LA.pl"
  load_container_image_if_present "HLA-LA" "neoag-hla-la:ubuntu22.04" "$image_tar"
}

install_sequenza_if_requested() {
  [[ "$INSTALL_SEQUENZA" == "1" ]] || return 0
  local env_path="$CONDA_BASE/envs/neoag-sequenza"
  local helper="$PROJECT_ROOT/.agents/skills/neoag-remote-deploy/scripts/with_bioc_data_cache.sh"
  local -a conda_cmd=("$helper" --conda-base "$CONDA_BASE" \
    --cache-root "$TOOLS_ROOT/install_cache" --package-key genomeinfodbdata-1.2.9 -- \
    "$CONDA_BASE/bin/conda")
  if [[ ! -x "$env_path/bin/sequenza-utils" ]]; then
    [[ -x "$helper" ]] || { echo "BIOC_CACHE_HELPER_MISSING: $helper" >&2; exit 47; }
    if [[ -d "$env_path" ]]; then
      run "repair Sequenza conda env" "${conda_cmd[@]}" env update -n neoag-sequenza \
        -f "$PROJECT_ROOT/conda/env.neoag-sequenza.yml" --prune \
        --override-channels -c conda-forge -c bioconda
    else
      run "install Sequenza conda env" "${conda_cmd[@]}" env create -n neoag-sequenza \
        -f "$PROJECT_ROOT/conda/env.neoag-sequenza.yml" -y \
        --override-channels -c conda-forge -c bioconda
    fi
  else
    log "Sequenza env already present: $env_path"
  fi
}

register_hmf_purple_if_requested() {
  [[ "$INSTALL_HMF_PURPLE" == "1" ]] || return 0
  local image_tar="$TOOLS_ROOT/container_images/neoag-purple-suite_ubuntu22.04.tar"
  if [[ "$EXECUTE" != "1" ]]; then
    log ""
    log "==> [DRY_RUN] register HMF PURPLE/AMBER/COBALT container image"
    log "+ load $image_tar if present"
    return 0
  fi
  load_container_image_if_present "HMF PURPLE suite" "neoag-purple-suite:ubuntu22.04" "$image_tar"
}

ensure_tf_keras_runtime() {
  [[ "$EXECUTE" == "1" ]] || return 0
  local py="$CONDA_BASE/envs/neoag-tools/bin/python"
  [[ -x "$py" ]] || return 0
  if "$py" - <<'PY' >/dev/null 2>&1
import os
os.environ["TF_USE_LEGACY_KERAS"] = "1"
import tf_keras
PY
  then
    log "tf-keras legacy shim already available in neoag-tools"
    return 0
  fi
  run "install tf-keras legacy shim for MHCflurry" bash -lc "source '$CONDA_BASE/etc/profile.d/conda.sh' && conda activate neoag-tools && spec=\$(python - <<'PY'
import tensorflow as tf
major, minor, *_ = tf.__version__.split('.')
print(f'tf-keras>={major}.{minor},<{major}.{int(minor) + 1}')
PY
) && pip install -q \"\$spec\""
}

repair_netmhcpan_frontend() {
  [[ "$EXECUTE" == "1" ]] || return 0
  local nm="$LICENSED_ROOT/netMHCpan/netMHCpan"
  if [[ -L "$LICENSED_ROOT/netMHCpan" ]]; then
    log "NetMHCpan uses an external symlink; leaving the licensed asset unchanged"
    return 0
  fi
  [[ -f "$nm" ]] || return 0
  if grep -qE '/(home|root)/[^/]+/(mini(conda|forge)|mambaforge)' "$nm" || grep -q 'CONDA_BASE=.*mini.*forge' "$nm"; then
    run "repair NetMHCpan frontend conda sysroot path" bash -lc "cp '$nm' '$nm.bak_\$(date +%Y%m%d_%H%M%S)' && perl -0pi -e 's#CONDA_BASE=\"\\\$\\{CONDA_BASE:-[^}]+\\}\"#CONDA_BASE=\"\\\${CONDA_BASE:-$CONDA_BASE}\"#' '$nm'"
  fi
  run "validate NetMHCpan frontend" bash -lc "CONDA_BASE='$CONDA_BASE' '$nm' -h >/dev/null"
}

install_sherpa_if_requested() {
  [[ "$INSTALL_SHERPA" == "1" ]] || return 0
  local home="$TOOLS_ROOT/tools/SHERPA-Presentation"
  local bin="$TOOLS_ROOT/bin/sherpa-presentation"
  if [[ "$EXECUTE" != "1" ]]; then
    log ""
    log "==> [DRY_RUN] register/install SHERPA-Presentation"
    log "+ require one of --sherpa-source, --sherpa-archive, or --sherpa-container-image; install/register under $home"
    return 0
  fi
  if [[ -z "$SHERPA_SOURCE$SHERPA_ARCHIVE$SHERPA_CONTAINER_IMAGE" ]]; then
    echo "SHERPA_PRESENTATION_SOURCE_REQUIRED: provide authorized --sherpa-source, --sherpa-archive, or --sherpa-container-image. SHERPA-Presentation is not publicly auto-downloadable by this installer." >&2
    exit 49
  fi
  mkdir -p "$home" "$TOOLS_ROOT/bin"
  if [[ -n "$SHERPA_SOURCE" ]]; then
    [[ -e "$SHERPA_SOURCE" ]] || { echo "SHERPA_SOURCE_MISSING: $SHERPA_SOURCE" >&2; exit 48; }
    run "copy/register SHERPA-Presentation source" bash -lc "mkdir -p '$home' && rsync -a --delete '$SHERPA_SOURCE/' '$home/'"
    if [[ -f "$home/pyproject.toml" || -f "$home/setup.py" ]]; then
      run "install SHERPA-Presentation Python package from source" bash -lc "source '$CONDA_BASE/etc/profile.d/conda.sh' && conda activate neoag-tools && python -m pip install '$home'"
    fi
  fi
  if [[ -n "$SHERPA_ARCHIVE" ]]; then
    [[ -f "$SHERPA_ARCHIVE" ]] || { echo "SHERPA_ARCHIVE_MISSING: $SHERPA_ARCHIVE" >&2; exit 50; }
    case "$SHERPA_ARCHIVE" in
      *.whl) run "install SHERPA-Presentation wheel" bash -lc "source '$CONDA_BASE/etc/profile.d/conda.sh' && conda activate neoag-tools && python -m pip install '$SHERPA_ARCHIVE'" ;;
      *.tar|*.tar.gz|*.tgz|*.zip) run "extract SHERPA-Presentation archive" bash -lc "rm -rf '$home' && mkdir -p '$home' && tar -xf '$SHERPA_ARCHIVE' -C '$home' --strip-components=1 2>/dev/null || unzip -q '$SHERPA_ARCHIVE' -d '$home'" ;;
      *) echo "SHERPA_ARCHIVE_UNSUPPORTED: $SHERPA_ARCHIVE" >&2; exit 51 ;;
    esac
  fi
  if [[ -n "$SHERPA_CONTAINER_IMAGE" ]]; then
    [[ -f "$SHERPA_CONTAINER_IMAGE" ]] || { echo "SHERPA_CONTAINER_IMAGE_MISSING: $SHERPA_CONTAINER_IMAGE" >&2; exit 52; }
    load_container_image_if_present "SHERPA-Presentation" "neoag-sherpa-presentation:latest" "$SHERPA_CONTAINER_IMAGE"
  fi
  cat > "$bin" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export SHERPA_PRESENTATION_HOME="\${SHERPA_PRESENTATION_HOME:-$home}"
if [[ -x "\$SHERPA_PRESENTATION_HOME/sherpa-presentation" ]]; then
  exec "\$SHERPA_PRESENTATION_HOME/sherpa-presentation" "\$@"
elif [[ -x "\$SHERPA_PRESENTATION_HOME/bin/sherpa-presentation" ]]; then
  exec "\$SHERPA_PRESENTATION_HOME/bin/sherpa-presentation" "\$@"
elif command -v sherpa-presentation >/dev/null 2>&1 && [[ "\$(command -v sherpa-presentation)" != "$bin" ]]; then
  exec sherpa-presentation "\$@"
else
  echo "SHERPA-Presentation is registered at \$SHERPA_PRESENTATION_HOME, but no executable named sherpa-presentation was found." >&2
  exit 127
fi
EOF
  chmod +x "$bin"
  if [[ -n "$SHERPA_SMOKE_COMMAND" ]]; then
    run "validate SHERPA-Presentation smoke command" bash -lc "$SHERPA_SMOKE_COMMAND"
  elif [[ -e "$home" || -f "$SHERPA_CONTAINER_IMAGE" ]]; then
    log "SHERPA-Presentation registered: $home"
  fi
}

ensure_bigmhc_torch_runtime() {
  [[ "$EXECUTE" == "1" ]] || return 0
  [[ "$SKIP_TORCH_INSTALL" == "0" ]] || return 0
  local py="$CONDA_BASE/envs/neoag-tools/bin/python"
  [[ -x "$py" ]] || return 0
  if "$py" - <<'PY' >/dev/null 2>&1
import torch
PY
  then
    log "torch already available in neoag-tools for BigMHC"
    return 0
  fi
  if [[ -n "$TORCH_WHEEL_DIR" && -d "$TORCH_WHEEL_DIR" ]] && compgen -G "$TORCH_WHEEL_DIR/torch-*.whl" >/dev/null; then
    run "install torch from local wheel cache for BigMHC" bash -lc "source '$CONDA_BASE/etc/profile.d/conda.sh' && conda activate neoag-tools && pip install --no-index --no-deps '$TORCH_WHEEL_DIR'/torch-*.whl && if compgen -G '$TORCH_WHEEL_DIR/nvidia_*.whl' >/dev/null || compgen -G '$TORCH_WHEEL_DIR/triton-*.whl' >/dev/null; then pip install --no-index --no-deps '$TORCH_WHEEL_DIR'/nvidia_*.whl '$TORCH_WHEEL_DIR'/triton-*.whl 2>/dev/null || true; fi"
  else
    need_download_ok "PyTorch for BigMHC"
    run "install CPU torch from approved package index for BigMHC" bash -lc "source '$CONDA_BASE/etc/profile.d/conda.sh' && conda activate neoag-tools && pip install --index-url '$TORCH_INDEX_URL' torch"
  fi
  run "install/repair common torch dependencies" bash -lc "source '$CONDA_BASE/etc/profile.d/conda.sh' && conda activate neoag-tools && pip install filelock 'sympy==1.13.1'"
  if ! "$py" - <<'PY' >/dev/null 2>&1
import torch
PY
  then
    if [[ "$ALLOW_DOWNLOAD" == "1" ]]; then
      run "install missing CUDA nvJitLink runtime for torch" bash -lc "source '$CONDA_BASE/etc/profile.d/conda.sh' && conda activate neoag-tools && pip install nvidia-nvjitlink-cu12"
    fi
  fi
  run "validate torch import for BigMHC" "$py" -c "import torch; print(torch.__version__)"
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

if [[ "$INSTALL_CORE_ENV$INSTALL_VEP$INSTALL_GATK$INSTALL_RNA_EXPRESSION$INSTALL_IMMUNOGENICITY$INSTALL_DEEPIMMUNO$INSTALL_SHERPA$INSTALL_NETMHCSTABPAN$INSTALL_NETCHOP$INSTALL_LOHHLA$INSTALL_POLYSOLVER$INSTALL_OPTITYPE$INSTALL_BAM_MATCHER$INSTALL_FACETS$INSTALL_ASCAT_PYCLONE$INSTALL_FUSION$INSTALL_SPLICE$INSTALL_SPECHLA$INSTALL_HLALA$INSTALL_SEQUENZA$INSTALL_HMF_PURPLE" =~ 1 ]]; then
  install_miniforge_if_needed
  export NEOAG_CONDA_BASE="$CONDA_BASE"
  export PATH="$CONDA_BASE/bin:$PATH"
fi

export NEOAG_TOOLS_ROOT="$TOOLS_ROOT"
export NEOAG_REFERENCE_ROOT="$REFERENCE_ROOT"
export NEOAG_REF_BUNDLE="$REFERENCE_ROOT"
if [[ -f "$REFERENCE_ROOT/data/normal/junctions/normal_junctions.GRCh38.tsv.gz" ]]; then
  export NEOAG_NORMAL_JUNCTIONS="$REFERENCE_ROOT/data/normal/junctions/normal_junctions.GRCh38.tsv.gz"
elif [[ -f "$REFERENCE_ROOT/data/normal/junctions/gtex_v8_liver.GRCh38.tsv.gz" ]]; then
  export NEOAG_NORMAL_JUNCTIONS="$REFERENCE_ROOT/data/normal/junctions/gtex_v8_liver.GRCh38.tsv.gz"
fi
export NETMHCPAN_HOME="$LICENSED_ROOT/netMHCpan"
export NETMHCpan="$LICENSED_ROOT/netMHCpan"
export NETMHCSTABPAN_HOME="$LICENSED_ROOT/netMHCstabpan"
export PRIME_HOME="$TOOLS_ROOT/tools/prime"
export MIXMHCPRED_HOME="$LICENSED_ROOT/mixMHCpred_install"
export BIGMHC_DIR="$TOOLS_ROOT/tools/bigmhc"
export BIGMHC_PYTHON="$CONDA_BASE/envs/neoag-tools/bin/python"
export DEEPIMMUNO_DIR="$TOOLS_ROOT/tools/DeepImmuno"
export SHERPA_PRESENTATION_HOME="$TOOLS_ROOT/tools/SHERPA-Presentation"
export SHERPA_PRESENTATION_BIN="$TOOLS_ROOT/bin/sherpa-presentation"
export SPECHLA_HOME="$TOOLS_ROOT/tools/SpecHLA"
export SPECHLA_DB="$REFERENCE_ROOT/data/hla/spechla/db"
export SPECHLA_ENV="$CONDA_BASE/envs/neoag-tools"
export NEOAG_BAM_MATCHER_ENV_PREFIX="$TOOLS_ROOT/conda_envs/neoag-bam-matcher"
export BAM_MATCHER_HOME="$TOOLS_ROOT/tools/bam-matcher"
# Prefer the portable GRCh38 identity panel (full-tier sample_identity_reference).
if [[ -f "$REFERENCE_ROOT/data/sample_identity/bam_matcher.common_snps.hg38.vcf" ]]; then
  export BAM_MATCHER_LOCI="$REFERENCE_ROOT/data/sample_identity/bam_matcher.common_snps.hg38.vcf"
elif [[ -f "$REFERENCE_ROOT/data/facets/reference/bam_matcher.identity.hg38.vcf.gz" ]]; then
  export BAM_MATCHER_LOCI="$REFERENCE_ROOT/data/facets/reference/bam_matcher.identity.hg38.vcf.gz"
else
  export BAM_MATCHER_LOCI="$REFERENCE_ROOT/data/facets/reference/1000G_omni2.5.hg38.biallelic.vcf.gz"
fi
export HLALA_HOME="$TOOLS_ROOT/tools/HLA-LA"
export HLA_LA_HOME="$HLALA_HOME"
export HLALA_ENV_PREFIX="$TOOLS_ROOT/tools/HLA-LA/.conda"
export HLA_LA_ENV_PREFIX="$HLALA_ENV_PREFIX"
export HLALA_BIN="$HLALA_ENV_PREFIX/bin/HLA-LA.pl"
export HLA_LA_BIN="$HLALA_BIN"
export HLALA_GRAPH="$REFERENCE_ROOT/data/hla/PRG_MHC_GRCh38_withIMGT"
export HLA_LA_GRAPH="$HLALA_GRAPH"
export NEOAG_HLALA_BACKEND="auto"
export HMFTOOLS_HOME="$TOOLS_ROOT/tools/HMFTOOLS"
export NEOAG_HMFTOOLS_HOME="$HMFTOOLS_HOME"
export HMF_ENV="$HMFTOOLS_HOME/.conda"
export HMFTOOLS_REFERENCE_ROOT="$REFERENCE_ROOT/data/hmf/purple_reference"
export HMFTOOLS_REFERENCE_FASTA="$REFERENCE_ROOT/data/sequenza/reference/GRCh38.primary_assembly.chr.fa"
export HMFTOOLS_AMBER_LOCI="$REFERENCE_ROOT/data/hmf/purple_reference/amber/GermlineHetPon.38.vcf.gz"
export HMFTOOLS_GC_PROFILE="$REFERENCE_ROOT/data/hmf/purple_reference/cobalt/GC_profile.1000bp.38.cnp"
export HMFTOOLS_ENSEMBL_DATA_DIR="$REFERENCE_ROOT/data/hmf/purple_reference/ensembl_data_cache_38"
export SEQUENZA_FASTA="$REFERENCE_ROOT/data/sequenza/reference/GRCh38.primary_assembly.chr.fa"
export SEQUENZA_GC_WIG="$REFERENCE_ROOT/data/sequenza/reference/Homo_sapiens.GRCh38.dna.primary_assembly.chr.gc50.wig.gz"
export FACETS_R_ENV_PREFIX="$CONDA_BASE/envs/neoag-fusion"
export SNP_PILEUP_BIN="$CONDA_BASE/envs/neoag-tools/bin/snp-pileup"
export SAMTOOLS="$CONDA_BASE/envs/neoag-tools/bin/samtools"
export PATH="$PROJECT_ROOT/bin:$CONDA_BASE/envs/neoag-tools/bin:$PATH"

if [[ -x "$TOOLS_ROOT/conda_pkgs/bedtools-2.31.1-h13024bc_3/bin/bedtools" ]]; then
  mkdir -p "$PROJECT_ROOT/bin"
  cat > "$PROJECT_ROOT/bin/bedtools" <<EOF
#!/usr/bin/env bash
set -euo pipefail
BEDTOOLS_BIN="\${BEDTOOLS_BIN:-$TOOLS_ROOT/conda_pkgs/bedtools-2.31.1-h13024bc_3/bin/bedtools}"
LIB_ROOT="\${BEDTOOLS_LIB_ROOT:-\${NEOAG_CONDA_BASE:-$CONDA_BASE}/envs/neoag-tools/lib}"
export LD_LIBRARY_PATH="\${LIB_ROOT}:\${LD_LIBRARY_PATH:-}"
exec "\${BEDTOOLS_BIN}" "\$@"
EOF
  chmod +x "$PROJECT_ROOT/bin/bedtools"
fi

sync_assets_if_requested

if [[ "$INSTALL_CLAUDE_CODE" == "1" ]]; then
  claude_args=(--outdir "$OUTDIR/claude_code" --channel "$CLAUDE_CODE_CHANNEL" --installer-url "$CLAUDE_CODE_INSTALLER_URL")
  [[ "$ALLOW_DOWNLOAD" == "1" ]] && claude_args+=(--allow-download)
  [[ "$EXECUTE" == "1" ]] && claude_args+=(--execute)
  run "install Claude Code" bash .agents/skills/neoag-remote-deploy/scripts/17_install_claude_code.sh "${claude_args[@]}"
fi

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
ensure_tf_keras_runtime
repair_netmhcpan_frontend
[[ "$INSTALL_VEP" == "1" ]] && run "install VEP env" env NEOAG_VEP_VERSION="$VEP_VERSION" bash scripts/install_vep.sh
[[ "$INSTALL_VEP_CACHE" == "1" ]] && { need_download_ok "VEP cache"; run "install VEP cache" env NEOAG_VEP_CACHE_VERSION="$VEP_VERSION" bash scripts/install_vep_cache.sh; }
[[ "$INSTALL_GATK" == "1" ]] && run "install GATK4" bash scripts/install_gatk.sh
IMMUNO_PYTHON="${CONDA_BASE}/envs/neoag-tools/bin/python"
[[ -x "$IMMUNO_PYTHON" ]] || IMMUNO_PYTHON="${CONDA_BASE}/envs/neoag-core/bin/python"
[[ -x "$IMMUNO_PYTHON" ]] || IMMUNO_PYTHON="${CONDA_BASE}/bin/python"
[[ -x "$IMMUNO_PYTHON" ]] || IMMUNO_PYTHON="$(command -v python3)"
[[ "$INSTALL_IMMUNOGENICITY" == "1" ]] && run "install PRIME/MixMHCpred/BigMHC" env NEOAG_SKIP_TORCH_INSTALL="$SKIP_TORCH_INSTALL" NEOAG_IMMUNO_PYTHON="$IMMUNO_PYTHON" bash scripts/install_immunogenicity_tools.sh
ensure_bigmhc_torch_runtime
register_netmhcstabpan_if_requested
if [[ "$INSTALL_NETCHOP" == "1" ]]; then
  [[ -n "$NETCHOP_ARCHIVE" ]] || NETCHOP_ARCHIVE="$(discover_netchop_archive)"
  [[ -n "$NETCHOP_ARCHIVE" && -f "$NETCHOP_ARCHIVE" ]] || { echo "NETCHOP_ARCHIVE_MISSING: provide --netchop-archive or place netchop-3.1d.Linux.tar.gz under $LICENSED_ROOT/netchop or --shared-asset-root" >&2; exit 46; }
  if [[ "$NETCHOP_ARCHIVE" != "$LICENSED_ROOT/netchop/netchop-3.1d.Linux.tar.gz" ]]; then
    run "stage NetChop archive into licensed root" bash -lc "mkdir -p '$LICENSED_ROOT/netchop' && cp -f '$NETCHOP_ARCHIVE' '$LICENSED_ROOT/netchop/netchop-3.1d.Linux.tar.gz'"
    NETCHOP_ARCHIVE="$LICENSED_ROOT/netchop/netchop-3.1d.Linux.tar.gz"
  fi
  run "install NetChop 3.1d" env NEOAG_LICENSED_ROOT="$LICENSED_ROOT" NEOAG_TOOLS_ROOT="$TOOLS_ROOT" bash scripts/install_netchop.sh "$NETCHOP_ARCHIVE"
  write_presentation_tool_env_overrides
fi
install_sherpa_if_requested
if [[ "$INSTALL_DEEPIMMUNO" == "1" ]]; then
  if [[ -z "$DEEPIMMUNO_SOURCE" && -f "$TOOLS_ROOT/tools/DeepImmuno/deepimmuno-cnn.py" ]]; then
    DEEPIMMUNO_SOURCE="$TOOLS_ROOT/tools/DeepImmuno"
  fi
  [[ -z "$DEEPIMMUNO_SOURCE" ]] && need_download_ok "DeepImmuno git clone"
  if [[ -n "$DEEPIMMUNO_SOURCE" ]]; then
    run "install DeepImmuno from local/source asset" bash scripts/install_deepimmuno.sh "$DEEPIMMUNO_SOURCE"
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
# Doctor / open-neo-run expect lohhla_reference under REFERENCE_ROOT/data/lohhla/polysolver.
# Licensed sync lands Polysolver under LICENSED_ROOT; wire a stable symlink when needed.
wire_lohhla_polysolver_reference() {
  local licensed="${POLYSOLVER_HOME_ARG:-$LICENSED_ROOT/polysolver}"
  local dest="$REFERENCE_ROOT/data/lohhla/polysolver"
  if [[ ! -d "$licensed" ]]; then
    return 0
  fi
  if [[ -e "$dest" || -L "$dest" ]]; then
    return 0
  fi
  mkdir -p "$(dirname "$dest")"
  ln -sfn "$licensed" "$dest"
  log "Linked lohhla_reference -> $dest -> $licensed"
}
if [[ "$INSTALL_LOHHLA" == "1" || -d "${POLYSOLVER_HOME_ARG:-$LICENSED_ROOT/polysolver}" ]]; then
  wire_lohhla_polysolver_reference
fi
[[ "$INSTALL_OPTITYPE" == "1" ]] && run "install OptiType" bash scripts/install_optitype.sh
[[ "$INSTALL_BAM_MATCHER" == "1" ]] && { need_download_ok "BAM-matcher pinned source and compatibility environment"; run "install BAM-matcher" bash scripts/install_bam_matcher.sh; }
[[ "$INSTALL_FACETS" == "1" ]] && run "install FACETS" bash scripts/install_facets.sh
[[ "$INSTALL_ASCAT_PYCLONE" == "1" ]] && run "install ASCAT/PyClone-VI" bash scripts/install_ascat_pyclone.sh
[[ "$INSTALL_FUSION" == "1" ]] && { need_download_ok "fusion tool git clones/conda packages"; run "install fusion tools" bash scripts/install_fusion_tools.sh; }
if [[ "$INSTALL_SPLICE" == "1" ]]; then
  if [[ "$INSTALL_SNAF" == "1" && "$EXECUTE" == "1" ]]; then
    need_download_ok "SNAF pinned Git source"
  fi
  if [[ "$INSTALL_SPLICEMUTR" == "1" && "$EXECUTE" == "1" ]]; then
    need_download_ok "SpliceMutr pinned Git source"
  fi
  run "install splice tools" env NEOAG_INSTALL_SNAF="$INSTALL_SNAF" NEOAG_INSTALL_SPLICEMUTR="$INSTALL_SPLICEMUTR" bash scripts/install_splice_tools.sh
fi
register_spechla_if_requested
register_hlala_if_requested
install_sequenza_if_requested
register_hmf_purple_if_requested
stage_bigmhc_models_if_requested

if [[ "$RUN_VERIFY" == "1" ]]; then
  # NEOAG_REF_BUNDLE is already exported to REFERENCE_ROOT above; pass flags only.
  # Passing the path again used to fail when the env var was already set.
  verify_args=(--smoke)
  [[ "$STRICT_VERIFY" == "1" ]] && verify_args=(--smoke --strict)
  run "verify README tools and references" \
    env NEOAG_REF_BUNDLE="$REFERENCE_ROOT" \
    bash scripts/verify_all_tools_and_refs.sh "${verify_args[@]}"
fi

if [[ "$RUN_REAL_VCF_SMOKE" == "1" ]]; then
  [[ -n "$REAL_VCF_RAW" ]] || { echo "REAL_VCF_REQUIRED: use --real-vcf" >&2; exit 46; }
  [[ -n "$REAL_VCF_ANNOTATED" ]] || { echo "REAL_ANNOTATED_VCF_REQUIRED: use --real-annotated-vcf" >&2; exit 46; }
  [[ -n "$REAL_VCF_HLA_ALLELES" || -n "$REAL_VCF_HLA_FILE" ]] || {
    echo "REAL_VCF_HLA_REQUIRED: use --real-vcf-hla-alleles or --real-vcf-hla-file" >&2
    exit 46
  }
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
  [[ "$REAL_VCF_SMOKE_SKIP_BIGMHC" == "1" ]] && real_vcf_args+=(--skip-bigmhc-im)
  run "run configured real VCF smoke test" bash .agents/skills/neoag-remote-deploy/scripts/14_run_real_vcf_smoke.sh "${real_vcf_args[@]}"
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
    "gatk:$INSTALL_GATK" "rna-expression:$INSTALL_RNA_EXPRESSION" "immunogenicity:$INSTALL_IMMUNOGENICITY" \
    "netmhcstabpan:$INSTALL_NETMHCSTABPAN" "netchop:$INSTALL_NETCHOP" "deepimmuno:$INSTALL_DEEPIMMUNO" "sherpa:$INSTALL_SHERPA" \
    "lohhla:$INSTALL_LOHHLA" "polysolver:$INSTALL_POLYSOLVER" "optitype:$INSTALL_OPTITYPE" "bam-matcher:$INSTALL_BAM_MATCHER" \
    "facets:$INSTALL_FACETS" "ascat-pyclone:$INSTALL_ASCAT_PYCLONE" "fusion:$INSTALL_FUSION" "splice:$INSTALL_SPLICE" \
    "spechla:$INSTALL_SPECHLA" "hla-la:$INSTALL_HLALA" "sequenza:$INSTALL_SEQUENZA" "hmf-purple:$INSTALL_HMF_PURPLE" \
    "claude-code:$INSTALL_CLAUDE_CODE" \
    "verify:$RUN_VERIFY" "real-vcf-smoke:$RUN_REAL_VCF_SMOKE" "sync-assets:$SYNC_ASSETS" "reference-manifest:${REFERENCE_MANIFEST:+1}" "bigmhc-models:${BIGMHC_MODELS_DIR:+1}"; do
    name="${item%%:*}"; enabled="${item##*:}"
    [[ "$enabled" == "1" ]] && echo "- $name"
  done
  if [[ "$RUN_REAL_VCF_SMOKE" == "1" ]]; then
    echo "- real-vcf-smoke-mhcflurry-default-on"
    [[ "$REAL_VCF_SMOKE_SKIP_MHCFLURRY" == "1" ]] && echo "- real-vcf-smoke-mhcflurry-skipped"
    [[ "$REAL_VCF_SMOKE_SKIP_BIGMHC" == "1" ]] && echo "- real-vcf-smoke-bigmhc-skipped"
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
