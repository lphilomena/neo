#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(pwd)"
DEPLOY_ROOT="${NEOAG_DEPLOY_ROOT:-/opt/neoag}"
TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-$DEPLOY_ROOT/env_tool}"
REFERENCE_ROOT="${NEOAG_REFERENCE_ROOT:-$DEPLOY_ROOT/refs}"
LICENSED_ROOT="${NEOAG_LICENSED_ROOT:-$DEPLOY_ROOT/licensed_tools}"
CONDA_BASE=""
OUTDIR="work/agent_deploy/new_machine_install"
ASSET_MANIFEST="configs/assets/production_assets.tsv"
REFERENCE_MANIFEST="configs/references/reference_manifest.yaml"
ASSET_SOURCE_HOST="${NEOAG_ASSET_SOURCE_HOST:-}"
ASSET_SSH_KEY="${NEOAG_ASSET_SSH_KEY:-}"
SHARED_ASSET_ROOT="${NEOAG_SHARED_ASSET_ROOT:-}"
VEP_VERSION="105"
EXECUTE=0
ALLOW_DOWNLOAD=0
INSTALL_CLAUDE_CODE=0
CLAUDE_CODE_CHANNEL="stable"
CLAUDE_CODE_INSTALLER_URL="https://claude.ai/install.sh"

INSTALL_TOOL_GROUPS=(--core-env --immunogenicity)
SYNC_ASSETS=1
RUN_VERIFY=1
STRICT_VERIFY=0
RUN_RUNTIME_VALIDATE=1
MINI_PRIME=1
RUN_REAL_VCF_SMOKE=0
REAL_VCF_SMOKE_TOP_N=1
SKIP_REAL_VCF_MHCFLURRY=0
REAL_VCF_RAW=""
REAL_VCF_ANNOTATED=""
REAL_VCF_HLA_ALLELES=""
REAL_VCF_HLA_FILE=""

EXTRA_INSTALL_ARGS=()

usage() {
  cat <<'USAGE'
Usage: 16_install_new_machine.sh [options]

One-entry installer for a new NeoAg machine. It orchestrates:
  1) large asset sync from configs/assets/production_assets.tsv,
  2) README-listed tool installation,
  3) activation/wrapper rewrite,
  4) production runtime validation,
  5) optional real VCF smoke test.

Default mode is dry-run. Add --execute to make changes.

Common options:
  --project-root DIR          Project checkout (default: current directory)
  --tools-root DIR            Tool/env root (default: NEOAG_TOOLS_ROOT or /opt/neoag/env_tool)
  --reference-root DIR        Reference root (default: NEOAG_REFERENCE_ROOT or /opt/neoag/refs)
  --licensed-root DIR         Licensed tool root (default: NEOAG_LICENSED_ROOT or /opt/neoag/licensed_tools)
  --conda-base DIR            Miniforge/conda base (default: tools-root/miniforge3)
  --outdir DIR                Work/report directory
  --asset-manifest FILE       Large asset manifest (default: configs/assets/production_assets.tsv)
  --reference-manifest FILE   YAML reference manifest verified after asset sync
  --asset-source-host HOST    Optional source host for remote asset paths (no default)
  --asset-ssh-key FILE        SSH private key used by rsync for remote assets
  --shared-asset-root DIR     Link assets from a locally mounted shared root to save disk space
  --allow-download            Permit official/user-approved network downloads
  --vep-version VERSION       Ensembl VEP/cache release to install/use (default: 105)
  --execute                   Actually run installation/sync/rewrite

Agent tooling:
  --claude-code               Install Claude Code using Anthropic's native installer
  --claude-code-channel V     stable, latest, or exact X.Y.Z version (default: stable)
  --claude-code-installer-url URL
                              Override only with an explicitly approved official URL

Tool group shortcuts:
  --minimal                   Install core env + immunogenicity only (default)
  --standard                  open-neo-run ready: core-env, VEP, GATK, immunogenicity,
                              OptiType, FACETS, splice, LOHHLA, fusion (Nextflow/STAR/
                              EasyFuse family), BAM-matcher, and CPU torch (BigMHC +
                              runtime validate). ASCAT/PyClone remains optional
                              (--add-tool-group --ascat-pyclone).
  --all-open                  Pass --all-open to 13_install_readme_tools.sh; NetMHCstabpan stays skipped
  --add-tool-group FLAG       Add any 13_install_readme_tools.sh group flag, e.g. --vep

Asset / validation toggles:
  --no-sync-assets            Do not sync asset manifest
  --no-verify                 Do not run verify_all_tools_and_refs.sh
  --strict-verify             Treat optional missing tools as verify failure
  --no-runtime-validate       Do not run 11_validate_production_runtime.sh
  --no-mini-prime             Skip PRIME mini smoke inside runtime validation

Real VCF smoke:
  --run-real-vcf-smoke        Run an explicitly configured real VCF smoke after install
  --real-vcf FILE             Raw somatic VCF used as the smoke-test anchor
  --real-annotated-vcf FILE   VEP-annotated VCF used for peptide extraction
  --real-vcf-hla-alleles L    Comma-separated HLA alleles
  --real-vcf-hla-file FILE    File containing HLA alleles
  --real-vcf-smoke-top-n N    Unique peptides for smoke test (default: 1)
  --skip-real-vcf-mhcflurry   Temporary fallback if MHCflurry is broken

Pass-through:
  --                          Remaining args are passed to 13_install_readme_tools.sh.

Examples:
  bash .agents/skills/neoag-remote-deploy/scripts/16_install_new_machine.sh \
    --asset-source-host <user@source-host> \
    --allow-download \
    --execute

  bash .agents/skills/neoag-remote-deploy/scripts/16_install_new_machine.sh \
    --standard \
    --run-real-vcf-smoke \
    --real-vcf <sample.somatic.vcf.gz> \
    --real-annotated-vcf <sample.vep.vcf> \
    --real-vcf-hla-file <sample.hla.txt> \
    --real-vcf-smoke-top-n 1 \
    --allow-download \
    --execute
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --project-root) PROJECT_ROOT="$2"; shift 2 ;;
    --tools-root) TOOLS_ROOT="$2"; shift 2 ;;
    --reference-root) REFERENCE_ROOT="$2"; shift 2 ;;
    --licensed-root) LICENSED_ROOT="$2"; shift 2 ;;
    --conda-base) CONDA_BASE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --asset-manifest) ASSET_MANIFEST="$2"; shift 2 ;;
    --reference-manifest) REFERENCE_MANIFEST="$2"; shift 2 ;;
    --asset-source-host) ASSET_SOURCE_HOST="$2"; shift 2 ;;
    --asset-ssh-key) ASSET_SSH_KEY="$2"; shift 2 ;;
    --shared-asset-root) SHARED_ASSET_ROOT="$2"; shift 2 ;;
    --allow-download) ALLOW_DOWNLOAD=1; shift ;;
    --vep-version) VEP_VERSION="$2"; shift 2 ;;
    --execute) EXECUTE=1; shift ;;
    --claude-code) INSTALL_CLAUDE_CODE=1; shift ;;
    --claude-code-channel) CLAUDE_CODE_CHANNEL="$2"; INSTALL_CLAUDE_CODE=1; shift 2 ;;
    --claude-code-installer-url) CLAUDE_CODE_INSTALLER_URL="$2"; INSTALL_CLAUDE_CODE=1; shift 2 ;;
    --minimal) INSTALL_TOOL_GROUPS=(--core-env --immunogenicity); shift ;;
    # --install-torch: standard includes BigMHC; 11_validate_production_runtime requires torch.
    # Default 13_install skips torch; without this flag Skill1 full/standard fails exit 21.
    # --fusion / --bam-matcher: required by Skill1 full tier groups and by open-neo-run
    # DNA identity + RNA FASTQ paths; previously only available via side-path installs.
    --standard) INSTALL_TOOL_GROUPS=(--core-env --vep --gatk --immunogenicity --optitype --facets --splice --lohhla --fusion --bam-matcher --install-torch); shift ;;
    --all-open) INSTALL_TOOL_GROUPS=(--all-open); shift ;;
    --all) echo "ERROR: --all has been retired for Skill1; use --all-open" >&2; exit 2 ;;
    --add-tool-group) EXTRA_INSTALL_ARGS+=("$2"); shift 2 ;;
    --no-sync-assets) SYNC_ASSETS=0; shift ;;
    --no-verify) RUN_VERIFY=0; shift ;;
    --strict-verify) RUN_VERIFY=1; STRICT_VERIFY=1; shift ;;
    --no-runtime-validate) RUN_RUNTIME_VALIDATE=0; shift ;;
    --no-mini-prime) MINI_PRIME=0; shift ;;
    --run-real-vcf-smoke) RUN_REAL_VCF_SMOKE=1; shift ;;
    --real-vcf) REAL_VCF_RAW="$2"; shift 2 ;;
    --real-annotated-vcf) REAL_VCF_ANNOTATED="$2"; shift 2 ;;
    --real-vcf-hla-alleles) REAL_VCF_HLA_ALLELES="$2"; shift 2 ;;
    --real-vcf-hla-file) REAL_VCF_HLA_FILE="$2"; shift 2 ;;
    --real-vcf-smoke-top-n) REAL_VCF_SMOKE_TOP_N="$2"; shift 2 ;;
    --skip-real-vcf-mhcflurry) SKIP_REAL_VCF_MHCFLURRY=1; shift ;;
    --) shift; EXTRA_INSTALL_ARGS+=("$@"); break ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

mkdir -p "$OUTDIR"
LOG="$OUTDIR/new_machine_install.log"
REPORT="$OUTDIR/new_machine_install_report.md"
: > "$LOG"
MODE="DRY_RUN"
[[ "$EXECUTE" == "1" ]] && MODE="EXECUTE"

log() { printf '%s\n' "$*" | tee -a "$LOG"; }
run_step() {
  local label="$1"; shift
  log ""
  log "==> $label"
  log "+ $*"
  "$@" 2>&1 | tee -a "$LOG"
}

cd "$PROJECT_ROOT"
[[ -f "pyproject.toml" || -f "setup.py" ]] || { echo "PROJECT_ROOT_INVALID: $PROJECT_ROOT" >&2; exit 30; }

CONDA_ARG=()
[[ -n "$CONDA_BASE" ]] && CONDA_ARG=(--conda-base "$CONDA_BASE")

install_args=(
  --project-root "$PROJECT_ROOT"
  --tools-root "$TOOLS_ROOT"
  --licensed-root "$LICENSED_ROOT"
  --reference-root "$REFERENCE_ROOT"
  "${CONDA_ARG[@]}"
  --outdir "$OUTDIR/readme_tools"
  "${INSTALL_TOOL_GROUPS[@]}"
  --vep-version "$VEP_VERSION"
)
[[ "$ALLOW_DOWNLOAD" == "1" ]] && install_args+=(--allow-download)
[[ "$EXECUTE" == "1" ]] && install_args+=(--execute)
if [[ "$INSTALL_CLAUDE_CODE" == "1" ]]; then
  install_args+=(--claude-code --claude-code-channel "$CLAUDE_CODE_CHANNEL" --claude-code-installer-url "$CLAUDE_CODE_INSTALLER_URL")
fi
if [[ "$SYNC_ASSETS" == "1" ]]; then
  install_args+=(--asset-manifest "$ASSET_MANIFEST" --reference-manifest "$REFERENCE_MANIFEST" --sync-assets)
  [[ -n "$ASSET_SOURCE_HOST" ]] && install_args+=(--asset-source-host "$ASSET_SOURCE_HOST")
  [[ -n "$ASSET_SSH_KEY" ]] && install_args+=(--asset-ssh-key "$ASSET_SSH_KEY")
  [[ -n "$SHARED_ASSET_ROOT" ]] && install_args+=(--shared-asset-root "$SHARED_ASSET_ROOT")
fi
if [[ "$RUN_VERIFY" == "1" ]]; then
  if [[ "$STRICT_VERIFY" == "1" ]]; then install_args+=(--strict-verify); else install_args+=(--verify); fi
fi
if [[ "$RUN_REAL_VCF_SMOKE" == "1" ]]; then
  [[ -n "$REAL_VCF_RAW" ]] || { echo "REAL_VCF_REQUIRED: use --real-vcf" >&2; exit 46; }
  [[ -n "$REAL_VCF_ANNOTATED" ]] || { echo "REAL_ANNOTATED_VCF_REQUIRED: use --real-annotated-vcf" >&2; exit 46; }
  [[ -n "$REAL_VCF_HLA_ALLELES" || -n "$REAL_VCF_HLA_FILE" ]] || {
    echo "REAL_VCF_HLA_REQUIRED: use --real-vcf-hla-alleles or --real-vcf-hla-file" >&2
    exit 46
  }
  install_args+=(--run-real-vcf-smoke --real-vcf-smoke-top-n "$REAL_VCF_SMOKE_TOP_N")
  install_args+=(--real-vcf "$REAL_VCF_RAW" --real-annotated-vcf "$REAL_VCF_ANNOTATED")
  [[ -n "$REAL_VCF_HLA_ALLELES" ]] && install_args+=(--real-vcf-hla-alleles "$REAL_VCF_HLA_ALLELES")
  [[ -n "$REAL_VCF_HLA_FILE" ]] && install_args+=(--real-vcf-hla-file "$REAL_VCF_HLA_FILE")
  [[ "$SKIP_REAL_VCF_MHCFLURRY" == "1" ]] && install_args+=(--skip-real-vcf-mhcflurry)
fi
install_args+=("${EXTRA_INSTALL_ARGS[@]}")

run_step "install tools and sync assets" bash .agents/skills/neoag-remote-deploy/scripts/13_install_readme_tools.sh "${install_args[@]}"

rewrite_args=(
  --project-root "$PROJECT_ROOT"
  --tools-root "$TOOLS_ROOT"
  --reference-root "$REFERENCE_ROOT"
  --licensed-root "$LICENSED_ROOT"
)
[[ "$EXECUTE" == "1" ]] && rewrite_args+=(--write)
run_step "rewrite activation and wrappers" bash .agents/skills/neoag-remote-deploy/scripts/10_rewrite_production_activation.sh "${rewrite_args[@]}"

if [[ "$RUN_RUNTIME_VALIDATE" == "1" ]]; then
  validate_args=(--project-root "$PROJECT_ROOT" --tools-root "$TOOLS_ROOT" --outdir "$OUTDIR/production_runtime")
  [[ "$MINI_PRIME" == "1" ]] && validate_args+=(--mini-prime)
  if [[ "$EXECUTE" == "1" ]]; then
    run_step "validate production runtime" bash .agents/skills/neoag-remote-deploy/scripts/11_validate_production_runtime.sh "${validate_args[@]}"
  else
    log ""
    log "==> [DRY_RUN] validate production runtime after installation"
    log "+ bash .agents/skills/neoag-remote-deploy/scripts/11_validate_production_runtime.sh ${validate_args[*]}"
  fi
fi

{
  echo "# New machine install report"
  echo
  echo "Mode: \`$MODE\`"
  echo "Project root: \`$PROJECT_ROOT\`"
  echo "Tools root: \`$TOOLS_ROOT\`"
  echo "Reference root: \`$REFERENCE_ROOT\`"
  echo "Licensed root: \`$LICENSED_ROOT\`"
  echo "Asset manifest: \`$ASSET_MANIFEST\`"
  echo "Reference manifest: \`$REFERENCE_MANIFEST\`"
  echo "Asset source host: \`${ASSET_SOURCE_HOST:-none}\`"
  echo "Shared asset root: \`${SHARED_ASSET_ROOT:-none}\`"
  echo "Claude Code: \`$INSTALL_CLAUDE_CODE\`"
  [[ "$INSTALL_CLAUDE_CODE" == "1" ]] && echo "Claude Code channel/version: \`$CLAUDE_CODE_CHANNEL\`"
  echo "Log: \`$LOG\`"
  echo
  echo "Next step: review logs under \`$OUTDIR\`."
} > "$REPORT"

log ""
log "new_machine_install_report=$REPORT"
