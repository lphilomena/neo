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
INSTALL_RNA_EXPRESSION=0
INSTALL_IMMUNOGENICITY=0
INSTALL_DEEPIMMUNO=0
INSTALL_SHERPA=0
INSTALL_NETMHCSTABPAN=0
INSTALL_LOHHLA=0
INSTALL_POLYSOLVER=0
INSTALL_OPTITYPE=0
INSTALL_FACETS=0
INSTALL_ASCAT_PYCLONE=0
INSTALL_FUSION=0
INSTALL_SPECHLA=0
INSTALL_HLALA=0
INSTALL_SEQUENZA=0
INSTALL_HMF_PURPLE=0
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
POLYSOLVER_HOME_ARG=""
NOVOALIGN_LICENSE_FILE_ARG=""
DEEPIMMUNO_SOURCE=""
SHERPA_PACKAGE="parameter-sherpa GPy"
SHERPA_SOURCE=""

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
  --sherpa                   SHERPA/parameter-sherpa Python package in existing neoag-tools
  --netmhcstabpan            NetMHCstabpan or IEDB shim via scripts/install_netmhcstabpan.sh
  --lohhla                   LOHHLA source wrapper via scripts/install_lohhla.sh
  --polysolver               Configure existing Polysolver; requires --polysolver-home
  --optitype                 OptiType via scripts/install_optitype.sh
  --facets                   FACETS via scripts/install_facets.sh
  --ascat-pyclone            ASCAT + PyClone-VI via scripts/install_ascat_pyclone.sh
  --fusion                   Arriba/Nextflow fusion env plus STAR-Fusion/FusionCatcher clones
  --spechla                  Register/load SpecHLA container assets and database if present
  --hla-la                   Register/load HLA-LA container assets and PRG graph if present
  --sequenza                 Install Sequenza conda env and reference hooks
  --hmf-purple               Register/load HMF PURPLE/AMBER/COBALT container assets and references
  --all-open                 Install open/conda/git tools except very large VEP cache and licensed packages
  --all                      Install all supported groups, including VEP cache; licensed packages still need sources/URLs
  --verify                   Run scripts/verify_all_tools_and_refs.sh after installs
  --strict-verify            Treat optional missing tools as failure during verify
  --run-real-vcf-smoke       Run default M1ML150017383 real VCF top-N smoke test
  --real-vcf-smoke-top-n N   Number of unique peptides for real VCF smoke (default: 50)
  --skip-real-vcf-mhcflurry Skip MHCflurry only for the real VCF smoke fallback
  --skip-real-vcf-bigmhc    Skip BigMHC_IM only for the real VCF smoke fallback
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
  --sherpa-package SPEC      Pip package spec(s) or approved Git URL for SHERPA (default: parameter-sherpa GPy)
  --sherpa-source DIR        Existing SHERPA source checkout/directory to pip install

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
    --lohhla) INSTALL_LOHHLA=1; shift ;;
    --polysolver) INSTALL_POLYSOLVER=1; shift ;;
    --optitype) INSTALL_OPTITYPE=1; shift ;;
    --facets) INSTALL_FACETS=1; shift ;;
    --ascat-pyclone) INSTALL_ASCAT_PYCLONE=1; shift ;;
    --fusion) INSTALL_FUSION=1; shift ;;
    --spechla) INSTALL_SPECHLA=1; shift ;;
    --hla-la) INSTALL_HLALA=1; shift ;;
    --sequenza) INSTALL_SEQUENZA=1; shift ;;
    --hmf-purple) INSTALL_HMF_PURPLE=1; shift ;;
    --all-open)
      INSTALL_CORE_ENV=1; INSTALL_VEP=1; INSTALL_GATK=1; INSTALL_RNA_EXPRESSION=1; INSTALL_IMMUNOGENICITY=1
      INSTALL_DEEPIMMUNO=1; INSTALL_SHERPA=1; INSTALL_NETMHCSTABPAN=1; INSTALL_LOHHLA=1
      INSTALL_OPTITYPE=1; INSTALL_FACETS=1; INSTALL_ASCAT_PYCLONE=1; INSTALL_FUSION=1
      INSTALL_SPECHLA=1; INSTALL_HLALA=1; INSTALL_SEQUENZA=1; INSTALL_HMF_PURPLE=1
      SKIP_TORCH_INSTALL=0
      shift ;;
    --all)
      INSTALL_CORE_ENV=1; INSTALL_VEP=1; INSTALL_VEP_CACHE=1; INSTALL_GATK=1; INSTALL_RNA_EXPRESSION=1; INSTALL_IMMUNOGENICITY=1
      INSTALL_DEEPIMMUNO=1; INSTALL_SHERPA=1; INSTALL_NETMHCSTABPAN=1; INSTALL_LOHHLA=1
      INSTALL_OPTITYPE=1; INSTALL_FACETS=1; INSTALL_ASCAT_PYCLONE=1; INSTALL_FUSION=1
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
    --sherpa-package) SHERPA_PACKAGE="$2"; shift 2 ;;
    --sherpa-source) SHERPA_SOURCE="$2"; shift 2 ;;
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

register_synced_tool_assets() {
  [[ "$EXECUTE" == "1" ]] || return 0

  if [[ -f "$TOOLS_ROOT/tools/DeepImmuno/deepimmuno-cnn.py" ]]; then
    run "register synced DeepImmuno asset" bash scripts/install_deepimmuno.sh "$TOOLS_ROOT/tools/DeepImmuno"
  fi

  if [[ -x "$LICENSED_ROOT/netMHCpan/netMHCpan" ]]; then
    if [[ -x "${CONDA_BASE:-$TOOLS_ROOT/miniforge3}/envs/neoag-tools/bin/patchelf" || -x "$LICENSED_ROOT/netMHCpan/Linux_x86_64/bin/netMHCpan-4.2" ]]; then
      if ! run "repair/register synced NetMHCpan asset" env NETMHCPAN_HOME="$LICENSED_ROOT/netMHCpan" NEOAG_CONDA_BASE="${CONDA_BASE:-$TOOLS_ROOT/miniforge3}" bash scripts/install_netmhcpan.sh --repair; then
        log "WARN: NetMHCpan asset is present but repair/smoke failed; license asset was left in place for manual validation."
      fi
    fi
  fi

  if [[ -x "$LICENSED_ROOT/netMHCstabpan/netMHCstabpan" ]]; then
    TOOLS_ENV="$PROJECT_ROOT/conf/tools.env.sh"
    run "register synced NetMHCstabpan asset" bash -lc "mkdir -p '$PROJECT_ROOT/conf'; if ! grep -q 'NETMHCSTABPAN_HOME' '$TOOLS_ENV' 2>/dev/null; then printf '\n# NetMHCstabpan (licensed or shim)\nexport NETMHCSTABPAN_HOME=\"$LICENSED_ROOT/netMHCstabpan\"\nexport PATH=\"$LICENSED_ROOT/netMHCstabpan:\$PATH\"\n' >> '$TOOLS_ENV'; fi"
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
  args=(--project-root "$PROJECT_ROOT" --asset-manifest "$ASSET_MANIFEST" --outdir "$OUTDIR/assets")
  [[ -n "$ASSET_SOURCE_HOST" ]] && args+=(--asset-source-host "$ASSET_SOURCE_HOST")
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

register_spechla_if_requested() {
  [[ "$INSTALL_SPECHLA" == "1" ]] || return 0
  local home="$TOOLS_ROOT/tools/SpecHLA"
  local image_tar="$TOOLS_ROOT/container_images/neoag-spechla_ubuntu22.04.tar"
  if [[ "$EXECUTE" != "1" ]]; then
    log ""
    log "==> [DRY_RUN] register SpecHLA container wrappers and DB link"
    log "+ create $home, link $home/db to $REFERENCE_ROOT/data/hla/spechla_db, load $image_tar if present"
    return 0
  fi
  mkdir -p "$home" "$TOOLS_ROOT/bin"
  [[ -d "$REFERENCE_ROOT/data/hla/spechla_db" && ! -e "$home/db" && ! -L "$home/db" ]] && ln -s "$REFERENCE_ROOT/data/hla/spechla_db" "$home/db"
  mkdir -p "$home/script/whole"
  if [[ ! -x "$home/script/whole/SpecHLA.sh" ]]; then
    cat > "$home/script/whole/SpecHLA.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec "$PROJECT_ROOT/scripts/run_spechla_container.sh" "\$@"
EOF
    chmod +x "$home/script/whole/SpecHLA.sh"
  fi
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
}

register_hlala_if_requested() {
  [[ "$INSTALL_HLALA" == "1" ]] || return 0
  local home="$TOOLS_ROOT/tools/HLA-LA"
  local image_tar="$TOOLS_ROOT/container_images/neoag-hla-la_ubuntu22.04.tar"
  if [[ "$EXECUTE" != "1" ]]; then
    log ""
    log "==> [DRY_RUN] register HLA-LA container wrapper and graph path"
    log "+ create $home/bin/HLA-LA.pl, link $TOOLS_ROOT/bin/HLA-LA.pl, load $image_tar if present"
    return 0
  fi
  mkdir -p "$home/bin" "$TOOLS_ROOT/bin"
  cat > "$home/bin/HLA-LA.pl" <<EOF
#!/usr/bin/env bash
set -euo pipefail
export HLALA_HOME="\${HLALA_HOME:-$home}"
export HLALA_GRAPH="\${HLALA_GRAPH:-$REFERENCE_ROOT/data/hla/PRG_MHC_GRCh38_withIMGT}"
exec "$PROJECT_ROOT/scripts/run_hla_la_container.sh" "\$@"
EOF
  chmod +x "$home/bin/HLA-LA.pl"
  ln -sf "$home/bin/HLA-LA.pl" "$TOOLS_ROOT/bin/HLA-LA.pl"
  load_container_image_if_present "HLA-LA" "neoag-hla-la:ubuntu22.04" "$image_tar"
}

install_sequenza_if_requested() {
  [[ "$INSTALL_SEQUENZA" == "1" ]] || return 0
  local env_path="$CONDA_BASE/envs/neoag-sequenza"
  if [[ ! -x "$env_path/bin/sequenza-utils" ]]; then
    if [[ -d "$env_path" ]]; then
      run "repair Sequenza conda env" bash -lc "source '$CONDA_BASE/etc/profile.d/conda.sh' && conda env update -n neoag-sequenza -f '$PROJECT_ROOT/conda/env.neoag-sequenza.yml' --prune"
    else
      run "install Sequenza conda env" bash -lc "source '$CONDA_BASE/etc/profile.d/conda.sh' && conda env create -f '$PROJECT_ROOT/conda/env.neoag-sequenza.yml'"
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
  [[ -f "$nm" ]] || return 0
  if grep -q '/home/na/miniforge3' "$nm" || grep -q 'CONDA_BASE=.*miniforge3' "$nm"; then
    run "repair NetMHCpan frontend conda sysroot path" bash -lc "cp '$nm' '$nm.bak_\$(date +%Y%m%d_%H%M%S)' && perl -0pi -e 's#CONDA_BASE=\"\\\$\\{CONDA_BASE:-[^}]+\\}\"#CONDA_BASE=\"\\\${CONDA_BASE:-$CONDA_BASE}\"#' '$nm'"
  fi
  run "validate NetMHCpan frontend" bash -lc "CONDA_BASE='$CONDA_BASE' '$nm' -h >/dev/null"
}

install_sherpa_if_requested() {
  [[ "$INSTALL_SHERPA" == "1" ]] || return 0
  if [[ "$EXECUTE" != "1" ]]; then
    log ""
    log "==> [DRY_RUN] install SHERPA Python package"
    log "+ install ${SHERPA_SOURCE:-$SHERPA_PACKAGE} into neoag-tools and validate import sherpa"
    return 0
  fi
  local py="$CONDA_BASE/envs/neoag-tools/bin/python"
  [[ -x "$py" ]] || py="$CONDA_BASE/bin/python"
  [[ -x "$py" ]] || { echo "SHERPA_PYTHON_MISSING: install --core-env first" >&2; exit 47; }
  if "$py" - <<'PY' >/dev/null 2>&1
import sherpa
PY
  then
    log "SHERPA already importable in neoag-tools"
    return 0
  fi
  if [[ -n "$SHERPA_SOURCE" ]]; then
    [[ -e "$SHERPA_SOURCE" ]] || { echo "SHERPA_SOURCE_MISSING: $SHERPA_SOURCE" >&2; exit 48; }
    run "install SHERPA from local/source directory" bash -lc "source '$CONDA_BASE/etc/profile.d/conda.sh' && conda activate neoag-tools && python -m pip install '$SHERPA_SOURCE'"
  else
    need_download_ok "SHERPA Python package ($SHERPA_PACKAGE)"
    run "install SHERPA Python package" bash -lc "source '$CONDA_BASE/etc/profile.d/conda.sh' && conda activate neoag-tools && python -m pip install $SHERPA_PACKAGE"
  fi
  run "validate SHERPA import" "$py" -c "import sherpa; print(getattr(sherpa, '__version__', 'unknown'))"
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

if [[ "$INSTALL_CORE_ENV$INSTALL_VEP$INSTALL_GATK$INSTALL_RNA_EXPRESSION$INSTALL_IMMUNOGENICITY$INSTALL_DEEPIMMUNO$INSTALL_SHERPA$INSTALL_NETMHCSTABPAN$INSTALL_LOHHLA$INSTALL_POLYSOLVER$INSTALL_OPTITYPE$INSTALL_FACETS$INSTALL_ASCAT_PYCLONE$INSTALL_FUSION$INSTALL_SPECHLA$INSTALL_HLALA$INSTALL_SEQUENZA$INSTALL_HMF_PURPLE" =~ 1 ]]; then
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
export BIGMHC_PYTHON="$CONDA_BASE/envs/neoag-tools/bin/python"
export DEEPIMMUNO_DIR="$TOOLS_ROOT/tools/DeepImmuno"
export SHERPA_PACKAGE="$SHERPA_PACKAGE"
export SPECHLA_HOME="$TOOLS_ROOT/tools/SpecHLA"
export HLALA_HOME="$TOOLS_ROOT/tools/HLA-LA"
export HLA_LA_HOME="$HLALA_HOME"
export HLALA_GRAPH="$REFERENCE_ROOT/data/hla/PRG_MHC_GRCh38_withIMGT"
export HLA_LA_GRAPH="$HLALA_GRAPH"
export HMFTOOLS_HOME="$TOOLS_ROOT/tools/HMFTOOLS"
export HMFTOOLS_AMBER_LOCI="$REFERENCE_ROOT/data/hmf/purple_reference/amber/GermlineHetPon.38.vcf.gz"
export HMFTOOLS_GC_PROFILE="$REFERENCE_ROOT/data/hmf/purple_reference/cobalt/GC_profile.1000bp.38.cnp"
export HMFTOOLS_ENSEMBL_DATA_DIR="$REFERENCE_ROOT/data/hmf/purple_reference/ensembl_data_cache_38"
export SEQUENZA_FASTA="$REFERENCE_ROOT/data/sequenza/reference/GRCh38.primary_assembly.chr.fa"
export SEQUENZA_GC_WIG="$REFERENCE_ROOT/data/sequenza/reference/gc.wig.gz"

sync_assets_if_requested

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
if [[ "$INSTALL_NETMHCSTABPAN" == "1" ]]; then
  if [[ -n "$NETMHCSTABPAN_ARCHIVE" && -f "$NETMHCSTABPAN_ARCHIVE" ]]; then
    run "install NetMHCstabpan from archive" bash scripts/install_netmhcstabpan.sh "$NETMHCSTABPAN_ARCHIVE"
  else
    run "install NetMHCstabpan IEDB shim" bash scripts/install_netmhcstabpan.sh --iedb
  fi
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
[[ "$INSTALL_OPTITYPE" == "1" ]] && run "install OptiType" bash scripts/install_optitype.sh
[[ "$INSTALL_FACETS" == "1" ]] && run "install FACETS" bash scripts/install_facets.sh
[[ "$INSTALL_ASCAT_PYCLONE" == "1" ]] && run "install ASCAT/PyClone-VI" bash scripts/install_ascat_pyclone.sh
[[ "$INSTALL_FUSION" == "1" ]] && { need_download_ok "fusion tool git clones/conda packages"; run "install fusion tools" bash scripts/install_fusion_tools.sh; }
register_spechla_if_requested
register_hlala_if_requested
install_sequenza_if_requested
register_hmf_purple_if_requested
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
  [[ "$REAL_VCF_SMOKE_SKIP_BIGMHC" == "1" ]] && real_vcf_args+=(--skip-bigmhc-im)
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
    "gatk:$INSTALL_GATK" "rna-expression:$INSTALL_RNA_EXPRESSION" "immunogenicity:$INSTALL_IMMUNOGENICITY" \
    "netmhcstabpan:$INSTALL_NETMHCSTABPAN" "deepimmuno:$INSTALL_DEEPIMMUNO" "sherpa:$INSTALL_SHERPA" \
    "lohhla:$INSTALL_LOHHLA" "polysolver:$INSTALL_POLYSOLVER" "optitype:$INSTALL_OPTITYPE" \
    "facets:$INSTALL_FACETS" "ascat-pyclone:$INSTALL_ASCAT_PYCLONE" "fusion:$INSTALL_FUSION" \
    "spechla:$INSTALL_SPECHLA" "hla-la:$INSTALL_HLALA" "sequenza:$INSTALL_SEQUENZA" "hmf-purple:$INSTALL_HMF_PURPLE" \
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
