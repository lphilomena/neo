#!/usr/bin/env bash
# Inspect external tool availability without claiming biological validation.
set -euo pipefail
printf 'tool\tstatus\tversion_or_path\n'
VERIFY_TIMEOUT="${NEOAG_TOOL_VERIFY_TIMEOUT:-120}"
check() {
 local label="$1" cmd="$2"; shift 2
 if command -v "$cmd" >/dev/null 2>&1; then
  local value rc tmp
  tmp="$(mktemp)"
  set +e
  timeout "$VERIFY_TIMEOUT" "$cmd" "$@" >"$tmp" 2>&1
  rc=$?
  set -e
  value="$(head -1 "$tmp" || true)"
  rm -f "$tmp"
  if [[ "$rc" == 0 ]]; then
   [[ -n "$value" ]] || value="$(command -v "$cmd")"
   printf '%s\tPASS\t%s\n' "$label" "$value"
  elif [[ "$rc" == 124 ]]; then
   printf '%s\tTIMEOUT\t%s seconds\n' "$label" "$VERIFY_TIMEOUT"
  else
   printf '%s\tFAIL\texit %s: %s\n' "$label" "$rc" "$value"
  fi
 else printf '%s\tMISSING\t\n' "$label"; fi
}
MOPEPGEN_BIN="${NEOAG_MOPEPGEN_BIN:-moPepGen}"
BP_QUANT_BIN="${NEOAG_BP_QUANT_BIN:-bp_quant}"
PVACSPLICE_BIN="${NEOAG_PVACSPLICE_BIN:-pvacsplice}"
PVACBIND_BIN="${NEOAG_PVACBIND_BIN:-pvacbind}"
RSCRIPT_BIN="${NEOAG_SPLICE2NEO_RSCRIPT:-Rscript}"

if command -v "$MOPEPGEN_BIN" >/dev/null 2>&1; then
 py="${NEOAG_PYTHON:-$(dirname "$(command -v "$MOPEPGEN_BIN")")/python}"
 version="$($py -c 'import moPepGen; print(moPepGen.__version__)' 2>/dev/null || true)"
 timeout "$VERIFY_TIMEOUT" "$MOPEPGEN_BIN" callAltTranslation -h >/dev/null 2>&1 || version=""
 [[ -n "$version" ]] && printf 'moPepGen\tPASS\t%s\n' "$version" || printf 'moPepGen\tFAIL\tcommand smoke failed\n'
else printf 'moPepGen\tMISSING\t\n'; fi
check EasyQuant "$BP_QUANT_BIN" --help
check pVACsplice "$PVACSPLICE_BIN" --help
check pVACbind "$PVACBIND_BIN" --help
K4NEO_BIN="${NEOAG_K4NEO_BIN:-k4neo-annotator}"
check k4neo "$K4NEO_BIN" --version
if [[ -x "$RSCRIPT_BIN" ]] || command -v "$RSCRIPT_BIN" >/dev/null 2>&1; then
 value="$("$RSCRIPT_BIN" -e "cat(if (requireNamespace('splice2neo', quietly=TRUE)) as.character(packageVersion('splice2neo')) else 'MISSING')" 2>/dev/null || true)"
 [[ "$value" == "MISSING" || -z "$value" ]] && printf 'splice2neo\tMISSING\t\n' || printf 'splice2neo\tPASS\t%s\n' "$value"
else printf 'splice2neo\tMISSING\tconfigured Rscript unavailable\n'; fi
