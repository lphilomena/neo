#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: build_bam_matcher_grch38_loci.sh \
  --source-panel BAM_MATCHER_GRCH37_PANEL.vcf \
  --grch38-common-vcf COMMON_SNP.hg38.vcf.gz \
  --reference-fasta GRCh38.fa \
  --output BAM_MATCHER.hg38.vcf \
  [--min-matched 1000] [--min-fraction-percent 80]

Builds a GRCh38 BAM-matcher identity panel by exact dbSNP ID lookup. It does
not lift coordinates heuristically. Only biallelic SNVs whose REF agrees with
the supplied GRCh38 FASTA are retained.
EOF
}

SOURCE_PANEL=""
COMMON_VCF=""
REFERENCE_FASTA=""
OUTPUT=""
MIN_MATCHED=1000
MIN_FRACTION_PERCENT=80
SAMTOOLS="${SAMTOOLS:-$(command -v samtools 2>/dev/null || true)}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --source-panel) SOURCE_PANEL="$2"; shift 2 ;;
    --grch38-common-vcf) COMMON_VCF="$2"; shift 2 ;;
    --reference-fasta) REFERENCE_FASTA="$2"; shift 2 ;;
    --output) OUTPUT="$2"; shift 2 ;;
    --min-matched) MIN_MATCHED="$2"; shift 2 ;;
    --min-fraction-percent) MIN_FRACTION_PERCENT="$2"; shift 2 ;;
    --samtools) SAMTOOLS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

for path in "$SOURCE_PANEL" "$COMMON_VCF" "$REFERENCE_FASTA"; do
  [[ -s "$path" ]] || { echo "ERROR: missing input: $path" >&2; exit 2; }
done
[[ -n "$OUTPUT" ]] || { usage; exit 2; }
[[ -x "$SAMTOOLS" ]] || { echo "ERROR: samtools is required" >&2; exit 127; }
[[ -s "$REFERENCE_FASTA.fai" ]] || { echo "ERROR: FASTA index missing: $REFERENCE_FASTA.fai" >&2; exit 2; }

mkdir -p "$(dirname "$OUTPUT")"
work="$(mktemp -d "${TMPDIR:-/tmp}/bam_matcher_loci.XXXXXX")"
trap 'rm -rf "$work"' EXIT
ids="$work/source_ids.txt"
matched="$work/matched_ids.txt"
raw="$work/raw.hg38.vcf"

awk '!/^#/ && $3 != "." {print $3}' "$SOURCE_PANEL" | sort -u > "$ids"
source_count="$(wc -l < "$ids" | tr -d ' ')"

gzip -cd "$COMMON_VCF" | awk -v ids="$ids" -v matched="$matched" '
  BEGIN { while ((getline id < ids) > 0) wanted[id] = 1; close(ids) }
  /^##fileformat=/ { print; next }
  /^#CHROM/ {
    print "##reference=GRCh38"
    print "##neoag_build_method=exact_dbSNP_ID_from_BAM_matcher_1kg_panel"
    print; next
  }
  /^#/ { print; next }
  ($3 in wanted) && length($4) == 1 && length($5) == 1 {
    print
    if (!seen[$3]++) print $3 > matched
  }
' > "$raw"

matched_count="$(sort -u "$matched" | wc -l | tr -d ' ')"
if (( matched_count < MIN_MATCHED || matched_count * 100 < source_count * MIN_FRACTION_PERCENT )); then
  echo "ERROR: only $matched_count/$source_count source dbSNP IDs mapped to GRCh38; require >=$MIN_MATCHED and >=${MIN_FRACTION_PERCENT}%" >&2
  exit 4
fi

grep '^#' "$raw" > "$OUTPUT.tmp"
while IFS=$'\t' read -r chrom pos id ref alt qual filter info; do
  output_chrom="$chrom"
  if ! grep -q -F "$chrom"$'\t' "$REFERENCE_FASTA.fai"; then
    output_chrom="${chrom#chr}"
  fi
  observed="$($SAMTOOLS faidx "$REFERENCE_FASTA" "${output_chrom}:${pos}-${pos}" | tail -n 1 | tr '[:lower:]' '[:upper:]')"
  if [[ "$observed" == "${ref^^}" ]]; then
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$output_chrom" "$pos" "$id" "$ref" "$alt" "$qual" "$filter" "$info" >> "$OUTPUT.tmp"
  else
    printf '%s\t%s\t%s\t%s\t%s\n' "$chrom" "$pos" "$id" "$ref" "$observed" >> "$OUTPUT.ref_mismatch.tsv"
  fi
done < <(grep -v '^#' "$raw")
mv "$OUTPUT.tmp" "$OUTPUT"
comm -23 "$ids" <(sort -u "$matched") > "$OUTPUT.missing_ids.txt"
record_count="$(grep -vc '^#' "$OUTPUT")"

source_sha="$(sha256sum "$SOURCE_PANEL" | awk '{print $1}')"
common_sha="$(sha256sum "$COMMON_VCF" | awk '{print $1}')"
fasta_sha="$(sha256sum "$REFERENCE_FASTA" | awk '{print $1}')"
output_sha="$(sha256sum "$OUTPUT" | awk '{print $1}')"
cat > "$OUTPUT.meta.json" <<EOF
{
  "schema_version": "neoag-bam-matcher-loci-v1",
  "genome_build": "GRCh38",
  "selection": "exact dbSNP ID match; biallelic SNV; REF validated against GRCh38 FASTA",
  "source_panel": "$SOURCE_PANEL",
  "source_panel_sha256": "$source_sha",
  "grch38_common_vcf": "$COMMON_VCF",
  "grch38_common_vcf_sha256": "$common_sha",
  "reference_fasta": "$REFERENCE_FASTA",
  "reference_fasta_sha256": "$fasta_sha",
  "source_id_count": $source_count,
  "matched_id_count": $matched_count,
  "output_record_count": $record_count,
  "output": "$OUTPUT",
  "output_sha256": "$output_sha"
}
EOF
sha256sum "$OUTPUT" "$OUTPUT.meta.json" "$OUTPUT.missing_ids.txt" > "$OUTPUT.SHA256SUMS"
echo "Built $OUTPUT with $record_count records ($matched_count/$source_count source IDs mapped)."
