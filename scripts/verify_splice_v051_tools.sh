#!/usr/bin/env bash
# Inspect external tool availability without claiming biological validation.
set -euo pipefail
printf 'tool\tstatus\tversion_or_path\n'
check() {
 local label="$1" cmd="$2"; shift 2
 if command -v "$cmd" >/dev/null 2>&1; then
  local value
  value="$("$cmd" "$@" 2>&1 | head -1 || true)"
  [[ -n "$value" ]] || value="$(command -v "$cmd")"
  printf '%s\tPASS\t%s\n' "$label" "$value"
 else printf '%s\tMISSING\t\n' "$label"; fi
}
check moPepGen moPepGen --version
check EasyQuant bp_quant --version
check pVACsplice pvacsplice --version
check pVACbind pvacbind --version
check k4neo k4neo-annotator --version
if command -v Rscript >/dev/null 2>&1; then
 value="$(Rscript -e "cat(if (requireNamespace('splice2neo', quietly=TRUE)) as.character(packageVersion('splice2neo')) else 'MISSING')" 2>/dev/null || true)"
 [[ "$value" == "MISSING" || -z "$value" ]] && printf 'splice2neo\tMISSING\t\n' || printf 'splice2neo\tPASS\t%s\n' "$value"
else printf 'splice2neo\tMISSING\tRscript unavailable\n'; fi
