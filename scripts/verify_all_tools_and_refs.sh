#!/usr/bin/env bash
# Unified acceptance checks for neoag-v03 external tools and reference data.
# This is a lightweight readiness check: it verifies executables, environment
# variables, and key reference files. It does not run patient-scale workflows.
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SMOKE=0
STRICT=0
SKIP_REFERENCE=0
REF_BUNDLE="${NEOAG_REF_BUNDLE:-}"

usage() {
  cat <<USAGE
Usage: bash scripts/verify_all_tools_and_refs.sh [--smoke] [--strict] [--skip-reference] [/path/to/neodata4git]

Options:
  --smoke           Also pass --smoke to scripts/verify_external_tools.sh.
  --strict          Treat missing optional/specialized tools as failures.
  --skip-reference  Skip scripts/verify_reference_bundle.sh.

Environment:
  NEOAG_REF_BUNDLE=/path/to/neodata4git  Portable reference bundle root.
  NEOAG_STRICT_VERIFY=1                  Same as --strict.

This script checks the core deployment tools plus VEP, GATK, NetMHCpan,
NetMHCstabpan, PRIME/MixMHCpred/BigMHC, EasyFuse, SpecHLA, HLA-LA, PURPLE/AMBER/COBALT, and Sequenza.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --smoke) SMOKE=1 ;;
    --strict) STRICT=1 ;;
    --skip-reference) SKIP_REFERENCE=1 ;;
    -h|--help) usage; exit 0 ;;
    *)
      if [[ -z "$REF_BUNDLE" && -d "$1" ]]; then
        REF_BUNDLE="$1"
      else
        echo "ERROR: unknown argument or missing directory: $1" >&2
        usage >&2
        exit 2
      fi
      ;;
  esac
  shift
done
[[ "${NEOAG_STRICT_VERIFY:-0}" == "1" ]] && STRICT=1

FAILED=0
WARNED=0
pass() { echo "[OK] $*"; }
warn() { echo "[WARN] $*" >&2; WARNED=1; }
fail() { echo "[FAIL] $*" >&2; FAILED=1; }
soft_fail() { [[ "$STRICT" == "1" ]] && fail "$*" || warn "$*"; }

have_cmd() { command -v "$1" >/dev/null 2>&1; }
first_existing_file() {
  local p
  for p in "$@"; do
    [[ -n "$p" && -f "$p" ]] && { echo "$p"; return 0; }
  done
  return 1
}
check_file() { [[ -f "$1" || -L "$1" ]] && pass "$2: $1" || soft_fail "$2 missing: $1"; }
check_dir() { [[ -d "$1" || -L "$1" ]] && pass "$2: $1" || soft_fail "$2 missing: $1"; }

source_if_exists() {
  local f="$1"
  if [[ -f "$f" ]]; then
    # shellcheck source=/dev/null
    source "$f" || warn "source returned non-zero: $f"
    pass "sourced $f"
  fi
}

export NEOAG_PROJECT_ROOT="$ROOT"
if [[ -n "$REF_BUNDLE" ]]; then
  export NEOAG_REF_BUNDLE="$REF_BUNDLE"
  source_if_exists "$REF_BUNDLE/neodata4git.env.sh"
fi
source_if_exists "$ROOT/conf/tools.env.sh"
if [[ -n "${NEOAG_REF_BUNDLE:-}" ]]; then
  source_if_exists "$NEOAG_REF_BUNDLE/neodata4git.env.sh"
fi
# Local overrides must win over portable bundle defaults on a configured server.
source_if_exists "$ROOT/conf/tools.env.local.sh"
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

run_check_script() {
  local label="$1"; shift
  echo
  echo "==> $label"
  "$@"
  local rc=$?
  if [[ "$rc" == "0" ]]; then
    pass "$label"
  else
    soft_fail "$label failed with exit code $rc"
  fi
}

echo "==> Unified tool/reference acceptance"
echo "ROOT=$ROOT"
echo "STRICT=$STRICT SMOKE=$SMOKE"
echo "NEOAG_REF_BUNDLE=${NEOAG_REF_BUNDLE:-unset}"

external_args=()
[[ "$SMOKE" == "1" ]] && external_args+=(--smoke)
run_check_script "Core external tools" bash "$ROOT/scripts/verify_external_tools.sh" "${external_args[@]}"

if [[ "$SKIP_REFERENCE" != "1" ]]; then
  if [[ -n "${NEOAG_REF_BUNDLE:-}" && -d "${NEOAG_REF_BUNDLE}" ]]; then
    run_check_script "Reference bundle" bash "$ROOT/scripts/verify_reference_bundle.sh" "${NEOAG_REF_BUNDLE}"
  else
    soft_fail "NEOAG_REF_BUNDLE is unset or missing; skip reference bundle verification"
  fi
fi

echo

echo "==> VEP"
vep_bin="${NEOAG_VEP_BIN:-}"
[[ -z "$vep_bin" ]] && vep_bin="$(command -v vep 2>/dev/null || true)"
if [[ -n "$vep_bin" && -x "$vep_bin" ]]; then
  pass "VEP executable: $vep_bin"
  "$vep_bin" --help >/dev/null 2>&1 || warn "VEP executable found but --help returned non-zero"
else
  soft_fail "VEP executable missing; set NEOAG_VEP_BIN or install with scripts/install_vep.sh"
fi
if [[ -n "${NEOAG_VEP_CACHE:-}" ]]; then
  check_dir "$NEOAG_VEP_CACHE/homo_sapiens/${NEOAG_VEP_CACHE_VERSION:-105}_GRCh38" "VEP GRCh38 cache"
  check_file "$NEOAG_VEP_CACHE/homo_sapiens/${NEOAG_VEP_CACHE_VERSION:-105}_GRCh38/info.txt" "VEP cache info"
else
  soft_fail "NEOAG_VEP_CACHE unset"
fi
if [[ -n "${NEOAG_VEP_PLUGINS:-}" ]]; then
  check_file "$NEOAG_VEP_PLUGINS/Wildtype.pm" "VEP Wildtype plugin"
  check_file "$NEOAG_VEP_PLUGINS/Frameshift.pm" "VEP Frameshift plugin"
else
  soft_fail "NEOAG_VEP_PLUGINS unset"
fi

echo

echo "==> GATK / Mutect2"
gatk_bin="$(command -v gatk 2>/dev/null || true)"
if [[ -n "$gatk_bin" ]]; then
  pass "GATK executable: $gatk_bin"
  gatk --version >/dev/null 2>&1 || warn "gatk --version returned non-zero"
else
  soft_fail "gatk missing from PATH; install with scripts/install_gatk.sh or source the GATK env"
fi
if [[ -n "${NEOAG_REFERENCE_FASTA:-}" ]]; then
  check_file "$NEOAG_REFERENCE_FASTA" "GRCh38 FASTA"
  check_file "$NEOAG_REFERENCE_FASTA.fai" "GRCh38 FASTA index"
else
  soft_fail "NEOAG_REFERENCE_FASTA unset"
fi

echo

echo "==> RNA expression quantification"
salmon_bin="${SALMON_BIN:-}"
[[ -z "$salmon_bin" ]] && salmon_bin="$(command -v salmon 2>/dev/null || true)"
if [[ -n "$salmon_bin" && -x "$salmon_bin" ]]; then
  pass "Salmon executable: $salmon_bin"
  "$salmon_bin" --version >/dev/null 2>&1 || warn "salmon --version returned non-zero"
else
  warn "Salmon executable missing; RNA FASTQ to gene TPM via scripts/run_salmon_fastq_to_tpm.sh will not run"
fi
if [[ -n "${SALMON_INDEX:-}" ]]; then
  check_dir "$SALMON_INDEX" "Salmon index"
else
  [[ -d "${NEOAG_REF_BUNDLE:-}/data/rna/salmon_index" ]] && pass "Salmon index: ${NEOAG_REF_BUNDLE}/data/rna/salmon_index" || warn "SALMON_INDEX unset and default data/rna/salmon_index missing"
fi
if [[ -n "${SALMON_TX2GENE:-}" ]]; then
  check_file "$SALMON_TX2GENE" "Salmon tx2gene"
else
  [[ -f "${NEOAG_REF_BUNDLE:-}/data/rna/tx2gene.tsv" ]] && pass "Salmon tx2gene: ${NEOAG_REF_BUNDLE}/data/rna/tx2gene.tsv" || warn "SALMON_TX2GENE unset and default data/rna/tx2gene.tsv missing"
fi
rsem_bin="${RSEM_BIN:-}"
[[ -z "$rsem_bin" ]] && rsem_bin="$(command -v rsem-calculate-expression 2>/dev/null || true)"
if [[ -n "$rsem_bin" && -x "$rsem_bin" ]]; then
  pass "RSEM executable: $rsem_bin"
else
  warn "RSEM executable missing; scripts/run_rsem_fastq_to_tpm.sh will not run"
fi

echo

echo "==> NetMHCpan"
netmhcpan_bin="${NEOAG_NETMHCPAN_BIN:-${NETMHCPAN_HOME:-}/netMHCpan}"
[[ ! -x "$netmhcpan_bin" ]] && netmhcpan_bin="$(command -v netMHCpan 2>/dev/null || true)"
if [[ -n "$netmhcpan_bin" && -x "$netmhcpan_bin" ]]; then
  pass "NetMHCpan executable: $netmhcpan_bin"
  "$netmhcpan_bin" -h >/dev/null 2>&1 || warn "NetMHCpan -h returned non-zero; check license/data paths"
else
  soft_fail "NetMHCpan executable missing; set NETMHCPAN_HOME or NEOAG_NETMHCPAN_BIN"
fi
[[ -n "${NETMHCPAN_HOME:-}" ]] && check_dir "$NETMHCPAN_HOME/data" "NetMHCpan data directory" || soft_fail "NETMHCPAN_HOME unset"

echo

echo "==> NetMHCstabpan"
netmhcstab_bin="${NETMHCSTABPAN_HOME:-}/netMHCstabpan"
[[ ! -x "$netmhcstab_bin" ]] && netmhcstab_bin="$(command -v netMHCstabpan 2>/dev/null || true)"
if [[ -n "$netmhcstab_bin" && -x "$netmhcstab_bin" ]]; then
  pass "NetMHCstabpan executable: $netmhcstab_bin"
  "$netmhcstab_bin" -h >/dev/null 2>&1 || warn "NetMHCstabpan -h returned non-zero; this may be normal for some wrappers"
else
  soft_fail "NetMHCstabpan executable missing; install with scripts/install_netmhcstabpan.sh"
fi
[[ -n "${NETMHCSTABPAN_HOME:-}" && -d "${NETMHCSTABPAN_HOME}/data" ]] && pass "NetMHCstabpan data directory: ${NETMHCSTABPAN_HOME}/data" || warn "NetMHCstabpan data directory not found; OK for IEDB shim, not OK for DTU binary"

echo

echo "==> PRIME / MixMHCpred / BigMHC"
prime_bin="${NEOAG_PRIME_BIN:-${PRIME_HOME:-}/PRIME}"
[[ ! -x "$prime_bin" ]] && prime_bin="$(command -v PRIME 2>/dev/null || true)"
if [[ -n "$prime_bin" && -x "$prime_bin" ]]; then
  pass "PRIME executable: $prime_bin"
else
  soft_fail "PRIME executable missing; set PRIME_HOME or NEOAG_PRIME_BIN"
fi
mix_bin="${MIXMHCPRED_BIN:-${MIXMHCPRED_HOME:-}/MixMHCpred}"
[[ ! -x "$mix_bin" ]] && mix_bin="$(command -v MixMHCpred 2>/dev/null || true)"
if [[ -n "$mix_bin" && -x "$mix_bin" ]]; then
  pass "MixMHCpred executable: $mix_bin"
else
  soft_fail "MixMHCpred executable missing; set MIXMHCPRED_HOME or MIXMHCPRED_BIN"
fi
bigmhc_dir="${BIGMHC_DIR:-${NEOAG_TOOLS_ROOT:-$ROOT}/tools/bigmhc}"
check_file "$bigmhc_dir/src/predict.py" "BigMHC predict.py"
bigmhc_bin="${BIGMHC_BIN:-}"
[[ -z "$bigmhc_bin" ]] && bigmhc_bin="$(command -v bigmhc_predict 2>/dev/null || true)"
if [[ -n "$bigmhc_bin" && -x "$bigmhc_bin" ]]; then
  pass "BigMHC wrapper: $bigmhc_bin"
else
  soft_fail "BigMHC wrapper missing; set BIGMHC_BIN or put bigmhc_predict on PATH"
fi
[[ -d "$bigmhc_dir/models" ]] && pass "BigMHC models: $bigmhc_dir/models" || warn "BigMHC models directory missing; predictions may fail"

sherpa_home="${SHERPA_PRESENTATION_HOME:-${NEOAG_TOOLS_ROOT:-$ROOT}/tools/SHERPA-Presentation}"
sherpa_bin="${SHERPA_PRESENTATION_BIN:-}"
[[ -z "$sherpa_bin" ]] && sherpa_bin="$(command -v sherpa-presentation 2>/dev/null || true)"
if [[ -n "$sherpa_bin" && -x "$sherpa_bin" ]]; then
  pass "SHERPA-Presentation wrapper: $sherpa_bin"
elif [[ -d "$sherpa_home" ]]; then
  warn "SHERPA-Presentation directory registered without executable wrapper: $sherpa_home"
else
  warn "SHERPA-Presentation not installed; provide authorized source/archive/container and install with --sherpa"
fi
echo

echo "==> Splice tools"
pvacsplice_bin="${NEOAG_PVACSPLICE_BIN:-}"
[[ -z "$pvacsplice_bin" ]] && pvacsplice_bin="$(command -v pvacsplice-neoag 2>/dev/null || command -v pvacsplice 2>/dev/null || true)"
if [[ -n "$pvacsplice_bin" && -x "$pvacsplice_bin" ]]; then
  pass "pVACsplice executable: $pvacsplice_bin"
else
  soft_fail "pVACsplice missing; install core pVACtools env and/or scripts/install_splice_tools.sh"
fi
regtools_bin="${NEOAG_REGTOOLS_BIN:-}"
[[ -z "$regtools_bin" ]] && regtools_bin="$(command -v regtools-neoag 2>/dev/null || command -v regtools 2>/dev/null || true)"
if [[ -n "$regtools_bin" && -x "$regtools_bin" ]]; then
  pass "RegTools executable: $regtools_bin"
else
  warn "RegTools missing; RNA BAM to junction extraction will not run"
fi
if command -v snaf-neoag >/dev/null 2>&1; then
  pass "SNAF wrapper: $(command -v snaf-neoag)"
else
  soft_fail "SNAF wrapper missing; splice installation includes SNAF by default (use --skip-snaf only intentionally)"
fi

echo

echo "==> EasyFuse / Nextflow"
if [[ -n "${NEOAG_EASYFUSE_HOME:-}" ]]; then
  check_dir "$NEOAG_EASYFUSE_HOME" "EasyFuse home"
  check_file "$NEOAG_EASYFUSE_HOME/main.nf" "EasyFuse main.nf"
else
  soft_fail "NEOAG_EASYFUSE_HOME unset"
fi
if [[ -n "${NEOAG_EASYFUSE_REF:-}" ]]; then
  check_dir "$NEOAG_EASYFUSE_REF" "EasyFuse reference"
  [[ -f "$NEOAG_EASYFUSE_REF/BEFORE_EXECUTING_EASYFUSE" ]] && pass "EasyFuse reference marker" || warn "EasyFuse reference marker BEFORE_EXECUTING_EASYFUSE missing"
else
  soft_fail "NEOAG_EASYFUSE_REF unset"
fi
nextflow_bin="${NEOAG_NEXTFLOW:-}"
[[ -z "$nextflow_bin" ]] && nextflow_bin="$(command -v nextflow 2>/dev/null || true)"
if [[ -n "$nextflow_bin" && -x "$nextflow_bin" ]]; then
  pass "Nextflow executable: $nextflow_bin"
  nxf_home="${NXF_HOME:-${HOME}/.nextflow}"
  nxf_ver="${NXF_VER:-}"
  if [[ -n "$nxf_ver" ]]; then
    [[ -s "$nxf_home/framework/$nxf_ver/nextflow-$nxf_ver-one.jar" ]] || warn "Nextflow executable exists; framework jar may still download on first run unless cache is pre-staged"
  fi
else
  soft_fail "Nextflow missing; required for EasyFuse"
fi

echo

echo "==> SpecHLA"
spechla_home="${SPECHLA_HOME:-${NEOAG_SPECHLA_HOME:-${NEOAG_TOOLS_ROOT:-$ROOT}/tools/SpecHLA}}"
spechla_env="${SPECHLA_ENV:-${spechla_home}/spechla_env}"
if [[ -d "$spechla_home" ]]; then
  pass "SpecHLA home: $spechla_home"
  check_dir "$spechla_home/db" "SpecHLA database"
  check_dir "$spechla_home/script" "SpecHLA scripts"
else
  soft_fail "SpecHLA home missing: $spechla_home"
fi
[[ -d "$spechla_env/bin" ]] && pass "SpecHLA env: $spechla_env" || warn "SpecHLA env missing: $spechla_env"

echo

echo "==> HLA-LA"
hlala_graph="${HLALA_GRAPH:-${HLA_LA_GRAPH:-${NEOAG_REF_BUNDLE:-}/data/hla/PRG_MHC_GRCh38_withIMGT}}"
hlala_bin="${HLALA_BIN:-${HLA_LA_BIN:-}}"
[[ -z "$hlala_bin" ]] && hlala_bin="$(command -v HLA-LA.pl 2>/dev/null || command -v hla-la 2>/dev/null || true)"
if [[ -n "$hlala_bin" && -x "$hlala_bin" ]]; then
  pass "HLA-LA executable: $hlala_bin"
else
  soft_fail "HLA-LA executable missing; set HLALA_BIN or HLA_LA_BIN"
fi
[[ -n "$hlala_graph" ]] && check_dir "$hlala_graph" "HLA-LA PRG graph" || soft_fail "HLA-LA graph path unset"

echo

echo "==> PURPLE / AMBER / COBALT"
hmf_home="${HMFTOOLS_HOME:-${NEOAG_HMFTOOLS_HOME:-}}"
if [[ -n "$hmf_home" ]]; then
  check_dir "$hmf_home" "HMFTOOLS home"
fi
for tool in amber cobalt purple; do
  upper_tool="${tool^^}"
  if have_cmd "$tool"; then
    tool_path="$(command -v "$tool")"
    pass "${upper_tool} wrapper: $tool_path"
    "$tool" -version >/dev/null 2>&1 || warn "${upper_tool} wrapper found but -version returned non-zero"
    continue
  fi
  jar_var="HMFTOOLS_${upper_tool}_JAR"
  jar="${!jar_var:-}"
  if [[ -n "$jar" ]]; then
    check_file "$jar" "${upper_tool} jar"
  else
    found="$(first_existing_file \
      "${hmf_home:+$hmf_home/${tool}.jar}" \
      "${hmf_home:+$hmf_home/${tool}-latest.jar}" \
      "${hmf_home:+$hmf_home/${tool}/${tool}.jar}" \
      "${hmf_home:+$hmf_home/${tool}/${tool}-latest.jar}" 2>/dev/null || true)"
    [[ -n "$found" ]] && pass "${upper_tool} jar: $found" || soft_fail "${upper_tool} jar or wrapper missing; set $jar_var, HMFTOOLS_HOME, or put $tool on PATH"
  fi
done
[[ -n "${HMFTOOLS_AMBER_LOCI:-}" ]] && check_file "$HMFTOOLS_AMBER_LOCI" "AMBER loci VCF" || warn "HMFTOOLS_AMBER_LOCI unset"
[[ -n "${HMFTOOLS_GC_PROFILE:-}" ]] && check_file "$HMFTOOLS_GC_PROFILE" "COBALT GC profile" || warn "HMFTOOLS_GC_PROFILE unset"
[[ -n "${HMFTOOLS_ENSEMBL_DATA_DIR:-}" ]] && check_dir "$HMFTOOLS_ENSEMBL_DATA_DIR" "PURPLE Ensembl data" || warn "HMFTOOLS_ENSEMBL_DATA_DIR unset"

echo

echo "==> Sequenza"
if have_cmd sequenza-utils; then
  pass "sequenza-utils: $(command -v sequenza-utils)"
  sequenza-utils --help >/dev/null 2>&1 || warn "sequenza-utils --help returned non-zero"
else
  soft_fail "sequenza-utils missing from PATH"
fi
if have_cmd Rscript; then
  if Rscript -e 'suppressPackageStartupMessages(library(sequenza)); cat(as.character(packageVersion("sequenza")), "\n")' >/tmp/neoag_sequenza_verify.txt 2>&1; then
    pass "R sequenza package: $(tail -1 /tmp/neoag_sequenza_verify.txt)"
  else
    soft_fail "R package sequenza not loadable"
  fi
else
  soft_fail "Rscript missing; required for Sequenza fit"
fi
seq_ref="${SEQUENZA_FASTA:-${NEOAG_REF_BUNDLE:-}/data/sequenza/reference/GRCh38.primary_assembly.chr.fa}"
[[ -n "$seq_ref" ]] && check_file "$seq_ref" "Sequenza FASTA" || warn "SEQUENZA_FASTA unset"
seq_gc="${SEQUENZA_GC_WIG:-${NEOAG_REF_BUNDLE:-}/data/sequenza/reference/gc.wig.gz}"
[[ -n "$seq_gc" ]] && check_file "$seq_gc" "Sequenza GC wiggle" || warn "SEQUENZA_GC_WIG unset"

echo
if [[ "$FAILED" != "0" ]]; then
  echo "==> Unified acceptance failed" >&2
  exit 1
fi
if [[ "$WARNED" != "0" ]]; then
  echo "==> Unified acceptance passed with warnings"
  exit 0
fi
echo "==> Unified acceptance passed"
