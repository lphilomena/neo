#!/usr/bin/env bash
set -euo pipefail

SAMPLE=260F400881011
WORK=/mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data/work/longrna/${SAMPLE}_20260722
REF=/root/neo/neodata4git/data/ref/hg38/Homo_sapiens_assembly38.fasta
GTF=/root/neo/neodata4git/data/ref/hg38/gencode.gtf
SQANTI_HOME=/root/neo/env_tool/tools/SQANTI3
SQANTI_ENV=/root/neo/env_tool/miniforge3/envs/neoag-sqanti3
STAR_DIR=/mnt/zjl-bgi-zzb/peixunban/gl/data/chenxiaoliang_data/work/snaf/FP500004780_L01_203_20260721/star
STAR_SJ="$STAR_DIR/FP500004780_L01_203.SJ.out.tab"
THREADS="${LONGRNA_DOWNSTREAM_THREADS:-4}"

mkdir -p "$WORK"/{logs,sqanti3,orf,fusion}
exec > >(tee -a "$WORK/logs/downstream.driver.log") 2>&1

stage() {
  printf '[%s]\t%s\n' "$(date '+%F %T')" "$1"
  printf 'stage\t%s\nupdated_at\t%s\n' "$1" "$(date -Iseconds)" > "$WORK/downstream_status.tsv"
}

fail() {
  code=$?
  printf 'status\tFAILED\nexit_code\t%s\nfailed_at\t%s\n' "$code" "$(date -Iseconds)" >> "$WORK/downstream_status.tsv"
  exit "$code"
}
trap fail ERR

stage WAITING_FOR_ISOQUANT
while [[ ! -e "$WORK/.isoquant_complete" ]]; do
  if [[ -e "$WORK/run_status.tsv" ]] && grep -q $'status\tFAILED' "$WORK/run_status.tsv"; then
    echo 'IsoQuant driver failed; downstream analysis will not start.' >&2
    exit 1
  fi
  sleep 300
done

ISOFORM_GTF="$(find "$WORK/isoquant" -type f -name '*.transcript_models.gtf' -size +0c -print | sort | head -1)"
test -n "$ISOFORM_GTF"
test -s "$REF"
test -s "$GTF"

stage SQANTI3_QC_ORF
if [[ ! -e "$WORK/sqanti3/.complete" ]]; then
  export PATH="$SQANTI_ENV/bin:$PATH"
  args=(
    --isoforms "$ISOFORM_GTF"
    --refGTF "$GTF"
    --refFasta "$REF"
    --include_ORF
    --isoform_hits
    --report html
    --cpus "$THREADS"
    --chunks 1
    --output "$SAMPLE"
    --dir "$WORK/sqanti3"
  )
  if [[ -s "$STAR_SJ" ]]; then
    args+=(--coverage "$STAR_SJ")
  else
    printf 'WARN\tSTAR junction file missing; SQANTI3 will run without short-read junction coverage: %s\n' "$STAR_SJ"
  fi
  "$SQANTI_ENV/bin/python" "$SQANTI_HOME/sqanti3_qc.py" "${args[@]}"
  test -s "$WORK/sqanti3/${SAMPLE}_classification.txt"
  test -s "$WORK/sqanti3/${SAMPLE}_corrected.gtf"
  test -s "$WORK/sqanti3/${SAMPLE}_corrected.fasta"
  touch "$WORK/sqanti3/.complete"
fi

stage SQANTI3_OUTPUT_VALIDATION
CLASSIFICATION="$WORK/sqanti3/${SAMPLE}_classification.txt"
awk -F '\t' 'NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
  {n++; c[$h["structural_category"]]++}
  END {print "metric\tvalue"; print "classified_isoforms\t" n;
       for(k in c) print "structural_category:" k "\t" c[k]}' \
  "$CLASSIFICATION" | sort > "$WORK/sqanti3/${SAMPLE}.classification_summary.tsv"

printf 'status\tSQANTI3_ORF_COMPLETE\ncompleted_at\t%s\n' "$(date -Iseconds)" >> "$WORK/downstream_status.tsv"
touch "$WORK/.sqanti3_complete"
stage READY_FOR_JAFFAL_AND_NEOANTIGEN_TRANSLATION
