#!/usr/bin/env bash
set -euo pipefail

status=0
check_file() { if [[ -s "$2" ]]; then printf "PASS\t%s\t%s\n" "$1" "$2"; else printf "MISSING\t%s\t%s\n" "$1" "$2"; status=1; fi; }

SNAF_DB="${NEOAG_SNAF_DB:-}"
[[ -n "$SNAF_DB" ]] || { echo -e "MISSING\tsnaf_db\tNEOAG_SNAF_DB"; exit 1; }
check_file snaf_exon_table "$SNAF_DB/Alt91_db/Hs_Ensembl_exon_add_col.txt"
check_file snaf_transcript_db "$SNAF_DB/Alt91_db/mRNA-ExonIDs.txt"
check_file snaf_gene_fasta "$SNAF_DB/Alt91_db/Hs_gene-seq-2000_flank.fa"
check_file snaf_gtex "$SNAF_DB/controls/GTEx_junction_counts.h5ad"
if command -v docker >/dev/null && docker image inspect "${NEOAG_ALTANALYZE_IMAGE:-neoag-altanalyze:snaf}" >/dev/null 2>&1; then
  echo -e "PASS\taltanalyze_image\t${NEOAG_ALTANALYZE_IMAGE:-neoag-altanalyze:snaf}"
else
  echo -e "MISSING\taltanalyze_image\t${NEOAG_ALTANALYZE_IMAGE:-neoag-altanalyze:snaf}"
  status=1
fi
if [[ -n "${NEOAG_SPLICEMUTR_HOME:-}" && -d "${NEOAG_SPLICEMUTR_HOME}" ]]; then
  echo -e "PASS\tsplicemutr_home\t${NEOAG_SPLICEMUTR_HOME}"
else
  echo -e "UNASSESSED\tsplicemutr_home\tNEOAG_SPLICEMUTR_HOME"
fi
exit "$status"
