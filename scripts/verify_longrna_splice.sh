#!/usr/bin/env bash
set -euo pipefail
WORKDIR=${1:?Usage: verify_longrna_splice.sh WORKDIR}
STATUS="$WORKDIR/longrna_splice_status.tsv"
fail=0
check() {
  if [[ -s "$1" ]]; then printf 'PASS\t%s\n' "$1"; else printf 'MISSING\t%s\n' "$1"; fail=1; fi
}
check "$STATUS"
check "$WORKDIR/README.workflow.txt"
for f in "$WORKDIR"/isoquant/*.transcript_models.gtf "$WORKDIR"/sqanti3/*_classification.txt "$WORKDIR"/translation/*.pep "$WORKDIR"/comprehensive_evidence.tsv; do
  [[ -e "$f" ]] && check "$f"
done
grep -q $'\tFAILED\t' "$STATUS" && fail=1 || true
if (( fail == 0 )); then echo workflow_status=PASS; else echo workflow_status=INCOMPLETE; exit 1; fi
