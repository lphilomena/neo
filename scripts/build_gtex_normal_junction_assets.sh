#!/usr/bin/env bash
set -euo pipefail

BASE_URL="${RECOUNT3_BASE_URL:-https://recount-opendata.s3.amazonaws.com/recount3/release/human/data_sources/gtex}"
ASSET_ROOT="${NEOAG_NORMAL_JUNCTION_ROOT:-/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git/data/normal/junctions}"
PROJECT_ROOT="${NEOAG_PROJECT_ROOT:-/root/neo/src/na0707_upload_release}"
PYTHON="${NEOAG_PYTHON:-/root/neo/env_tool/miniforge3/envs/neoag-tools/bin/python}"
RAW_ROOT="$ASSET_ROOT/recount3_gtex_v8/raw"
TISSUE_ROOT="$ASSET_ROOT/recount3_gtex_v8/per_tissue"
METADATA="$RAW_ROOT/gtex.recount_project.MD.gz"
STATUS="$ASSET_ROOT/recount3_gtex_v8/build_status.tsv"

mkdir -p "$RAW_ROOT" "$TISSUE_ROOT"
download_gzip() {
  local url="$1" destination="$2"
  if [[ -s "$destination" ]] && gzip -t "$destination" 2>/dev/null; then
    return 0
  fi
  curl -fsSL --retry 30 --retry-all-errors --retry-delay 5 --connect-timeout 30 -C - -o "$destination" "$url"
  gzip -t "$destination"
}

download_gzip "$BASE_URL/metadata/gtex.recount_project.MD.gz" "$METADATA"

if (( $# )); then
  tissues=("$@")
else
  mapfile -t tissues < <(gzip -cd "$METADATA" | awk -F '\t' 'NR > 1 && $4 != "STUDY_NA" {print $4}' | sort -u)
fi

printf 'tissue\tstage\tupdated_at\n' > "$STATUS"
for tissue in "${tissues[@]}"; do
  suffix="${tissue: -2}"
  raw_dir="$RAW_ROOT/$tissue"
  out="$TISSUE_ROOT/${tissue}.GRCh38.tsv.gz"
  mkdir -p "$raw_dir"
  if [[ -s "$out" && -s "$out.meta.json" ]] && gzip -t "$out" 2>/dev/null; then
    printf '%s\tDONE\t%s\n' "$tissue" "$(date -Is)" >> "$STATUS"
    continue
  fi
  printf '%s\tDOWNLOAD\t%s\n' "$tissue" "$(date -Is)" >> "$STATUS"
  for ext in RR MM ID; do
    file="gtex.junctions.${tissue}.ALL.${ext}.gz"
    download_gzip "$BASE_URL/junctions/$suffix/$tissue/$file" "$raw_dir/$file"
  done
  sha256sum "$raw_dir"/*.gz > "$raw_dir/SHA256SUMS"
  printf '%s\tCONVERT\t%s\n' "$tissue" "$(date -Is)" >> "$STATUS"
  "$PYTHON" "$PROJECT_ROOT/scripts/build_recount3_normal_junctions.py" \
    --rr "$raw_dir/gtex.junctions.${tissue}.ALL.RR.gz" \
    --mm "$raw_dir/gtex.junctions.${tissue}.ALL.MM.gz" \
    --tissue "$tissue" \
    --output "$out"
  printf '%s\tDONE\t%s\n' "$tissue" "$(date -Is)" >> "$STATUS"
done

if [[ -f "$TISSUE_ROOT/LIVER.GRCh38.tsv.gz" ]]; then
  cp -f "$TISSUE_ROOT/LIVER.GRCh38.tsv.gz" "$ASSET_ROOT/gtex_v8_liver.GRCh38.tsv.gz"
  cp -f "$TISSUE_ROOT/LIVER.GRCh38.tsv.gz.meta.json" "$ASSET_ROOT/gtex_v8_liver.GRCh38.tsv.gz.meta.json"
fi

if (( ${#tissues[@]} > 1 )); then
  mapfile -t tables < <(find "$TISSUE_ROOT" -maxdepth 1 -type f -name '*.GRCh38.tsv.gz' | sort)
  "$PYTHON" "$PROJECT_ROOT/scripts/merge_normal_junction_tissues.py" \
    --inputs "${tables[@]}" \
    --output "$ASSET_ROOT/normal_junctions.GRCh38.tsv.gz"
fi

sha256sum "$ASSET_ROOT"/*.GRCh38.tsv.gz > "$ASSET_ROOT/SHA256SUMS"
