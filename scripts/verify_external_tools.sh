#!/usr/bin/env bash
# Verify closed-loop deployment for LOHHLA, FACETS, ASCAT, Arriba, PRIME, and OptiType.
#
# Usage:
#   source conf/tools.env.sh
#   bash scripts/verify_external_tools.sh
#   bash scripts/verify_external_tools.sh --smoke
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SMOKE=0
if [[ "${1:-}" == "--smoke" ]]; then
  SMOKE=1
fi

# shellcheck source=/dev/null
source "${ROOT}/conf/tools.env.sh"
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"

pass() { echo "[OK] $*"; }
warn() { echo "[WARN] $*" >&2; }
fail() { echo "[FAIL] $*" >&2; FAILED=1; }

FAILED=0

echo "==> Registry check-tools"
python -m neoag.cli check-tools | grep -E 'lohhla|facets|ascat|arriba|prime|optitype' || warn "check-tools does not report every external helper; continuing with direct checks"

if command -v LOHHLA >/dev/null 2>&1; then
  pass "LOHHLA executable: $(command -v LOHHLA)"
  [[ -f "${LOHHLA_HOME:-}/LOHHLAscript.R" ]] || warn "LOHHLAscript.R missing under LOHHLA_HOME=${LOHHLA_HOME:-unset}"
  [[ -n "${POLYSOLVER_HOME:-}" ]] || warn "POLYSOLVER_HOME unset; production LOHHLA runs need Polysolver"
else
  fail "LOHHLA not on PATH"
fi

if command -v runFACETS.R >/dev/null 2>&1; then
  pass "FACETS wrapper: $(runFACETS.R --version 2>/dev/null | head -1 || echo runFACETS.R)"
  [[ -x "${ROOT}/bin/snp-pileup" || -n "${NEOAG_DBSNP_VCF:-}" ]] || warn "bin/snp-pileup and NEOAG_DBSNP_VCF both missing; real FACETS pileup runs need one SNP reference path"
else
  fail "runFACETS.R not on PATH"
fi

if command -v ascat.R >/dev/null 2>&1; then
  pass "ASCAT wrapper: $(ascat.R --version 2>/dev/null | head -1 || echo ascat.R)"
else
  fail "ascat.R not on PATH"
fi

if command -v arriba >/dev/null 2>&1; then
  pass "Arriba: $(arriba --version 2>/dev/null | head -1 || echo arriba)"
else
  fail "arriba not on PATH; run bash scripts/install_fusion_tools.sh"
fi

if [[ -n "${NEOAG_PRIME_BIN:-}" && -x "${NEOAG_PRIME_BIN}" ]]; then
  pass "PRIME: ${NEOAG_PRIME_BIN}"
  [[ -x "${MIXMHCPRED_BIN:-}" ]] || warn "MIXMHCPRED_BIN missing; PRIME runs require MixMHCpred"
else
  fail "PRIME not executable; run bash scripts/install_immunogenicity_tools.sh"
fi

if command -v optitype >/dev/null 2>&1; then
  optitype check-deps >/dev/null 2>&1 \
    && pass "OptiType: $(command -v optitype)" \
    || fail "OptiType is on PATH but optitype check-deps failed"
elif [[ "${SKIP_OPTITYPE:-0}" == "1" ]]; then
  warn "OptiType skipped by SKIP_OPTITYPE=1"
else
  fail "OptiType not on PATH; run bash scripts/install_optitype.sh"
fi

if [[ "${SMOKE}" == "1" ]]; then
  echo "==> Optional smoke tests"
  if [[ -f "${LOHHLA_HOME:-}/example-file/bam/example_tumor_sorted.bam" ]]; then
    pass "LOHHLA bundled example BAM present"
  else
    warn "LOHHLA example BAM missing; use bash scripts/run_lohhla_example.sh only after staging example-file data"
  fi
  if [[ -f "${PRIME_HOME:-}/test/test.txt" && -x "${NEOAG_PRIME_BIN:-}" && -x "${MIXMHCPRED_BIN:-}" ]]; then
    "${NEOAG_PRIME_BIN}" -i "${PRIME_HOME}/test/test.txt" -o /tmp/neoag_prime_verify.tsv \
      -a A0101 -mix "${MIXMHCPRED_BIN}" >/tmp/neoag_prime_verify.log 2>&1 \
      && pass "PRIME smoke test passed" \
      || warn "PRIME smoke test failed; see /tmp/neoag_prime_verify.log"
  fi
fi

if [[ "${FAILED}" != "0" ]]; then
  echo "==> Verification failed" >&2
  exit 1
fi
echo "==> External tool verification passed"
