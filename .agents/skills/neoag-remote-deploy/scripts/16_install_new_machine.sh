#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(pwd)"
TOOLS_ROOT="/root/neo/env_tool"
REFERENCE_ROOT="/root/neo/neodata4git"
LICENSED_ROOT="/root/neo/licensed_tools"
CONDA_BASE=""
OUTDIR="work/agent_deploy/new_machine_install"
ASSET_MANIFEST="configs/assets/production_assets.tsv"
ASSET_SOURCE_HOST="na@10.200.50.134"
VEP_VERSION="105"
EXECUTE=0
ALLOW_DOWNLOAD=0

INSTALL_TOOL_GROUPS=(--core-env --immunogenicity)
SYNC_ASSETS=1
RUN_VERIFY=1
STRICT_VERIFY=0
RUN_RUNTIME_VALIDATE=1
MINI_PRIME=1
RUN_REAL_VCF_SMOKE=0
REAL_VCF_SMOKE_TOP_N=1
SKIP_REAL_VCF_MHCFLURRY=0

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
  --tools-root DIR            Tool/env root (default: /root/neo/env_tool)
  --reference-root DIR        Reference root (default: /root/neo/neodata4git)
  --licensed-root DIR         Licensed tool root (default: /root/neo/licensed_tools)
  --conda-base DIR            Miniforge/conda base (default: tools-root/miniforge3)
  --outdir DIR                Work/report directory
  --asset-manifest FILE       Large asset manifest (default: configs/assets/production_assets.tsv)
  --asset-source-host HOST    Source host for asset paths (default: na@10.200.50.134)
  --allow-download            Permit official/user-approved network downloads
  --vep-version VERSION       Ensembl VEP/cache release to install/use (default: 105)
  --execute                   Actually run installation/sync/rewrite

Tool group shortcuts:
  --minimal                   Install core env + immunogenicity only (default)
  --standard                  Add VEP, GATK, OptiType, FACETS, ASCAT/PyClone
  --all-open                  Pass --all-open to 13_install_readme_tools.sh
  --all                      Pass --all to 13_install_readme_tools.sh
  --add-tool-group FLAG       Add any 13_install_readme_tools.sh group flag, e.g. --vep

Asset / validation toggles:
  --no-sync-assets            Do not sync asset manifest
  --no-verify                 Do not run verify_all_tools_and_refs.sh
  --strict-verify             Treat optional missing tools as verify failure
  --no-runtime-validate       Do not run 11_validate_production_runtime.sh
  --no-mini-prime             Skip PRIME mini smoke inside runtime validation

Real VCF smoke:
  --run-real-vcf-smoke        Run default real VCF smoke after install
  --real-vcf-smoke-top-n N    Unique peptides for smoke test (default: 1)
  --skip-real-vcf-mhcflurry   Temporary fallback if MHCflurry is broken

Pass-through:
  --                          Remaining args are passed to 13_install_readme_tools.sh.

Examples:
  bash .agents/skills/neoag-remote-deploy/scripts/16_install_new_machine.sh \
    --asset-source-host na@10.200.50.134 \
    --allow-download \
    --execute

  bash .agents/skills/neoag-remote-deploy/scripts/16_install_new_machine.sh \
    --standard \
    --run-real-vcf-smoke \
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
    --asset-source-host) ASSET_SOURCE_HOST="$2"; shift 2 ;;
    --allow-download) ALLOW_DOWNLOAD=1; shift ;;
    --vep-version) VEP_VERSION="$2"; shift 2 ;;
    --execute) EXECUTE=1; shift ;;
    --minimal) INSTALL_TOOL_GROUPS=(--core-env --immunogenicity); shift ;;
    --standard) INSTALL_TOOL_GROUPS=(--core-env --vep --gatk --immunogenicity --optitype --facets --ascat-pyclone); shift ;;
    --all-open) INSTALL_TOOL_GROUPS=(--all-open); shift ;;
    --all) INSTALL_TOOL_GROUPS=(--all); shift ;;
    --add-tool-group) EXTRA_INSTALL_ARGS+=("$2"); shift 2 ;;
    --no-sync-assets) SYNC_ASSETS=0; shift ;;
    --no-verify) RUN_VERIFY=0; shift ;;
    --strict-verify) RUN_VERIFY=1; STRICT_VERIFY=1; shift ;;
    --no-runtime-validate) RUN_RUNTIME_VALIDATE=0; shift ;;
    --no-mini-prime) MINI_PRIME=0; shift ;;
    --run-real-vcf-smoke) RUN_REAL_VCF_SMOKE=1; shift ;;
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
if [[ "$SYNC_ASSETS" == "1" ]]; then
  install_args+=(--asset-manifest "$ASSET_MANIFEST" --sync-assets)
  [[ -n "$ASSET_SOURCE_HOST" ]] && install_args+=(--asset-source-host "$ASSET_SOURCE_HOST")
fi
if [[ "$RUN_VERIFY" == "1" ]]; then
  if [[ "$STRICT_VERIFY" == "1" ]]; then install_args+=(--strict-verify); else install_args+=(--verify); fi
fi
if [[ "$RUN_REAL_VCF_SMOKE" == "1" ]]; then
  install_args+=(--run-real-vcf-smoke --real-vcf-smoke-top-n "$REAL_VCF_SMOKE_TOP_N")
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
  run_step "validate production runtime" bash .agents/skills/neoag-remote-deploy/scripts/11_validate_production_runtime.sh "${validate_args[@]}"
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
  echo "Asset source host: \`${ASSET_SOURCE_HOST:-none}\`"
  echo "Log: \`$LOG\`"
  echo
  echo "Next step: review logs under \`$OUTDIR\`."
} > "$REPORT"

log ""
log "new_machine_install_report=$REPORT"
