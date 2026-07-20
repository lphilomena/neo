#!/usr/bin/env bash
# check_peptides_similarity.sh — Compare two variant_peptides.annotated.tsv files
#
# Usage:
#   bash tests/check_peptides_similarity.sh <new_file> <reference_file>
#
# Example:
#   bash tests/check_peptides_similarity.sh \
#     ~/working/neoantigen/upstream/tools/variant_peptides.annotated.tsv \
#     ~/working/neoantigen_saved/sbb/upstream/tools/variant_peptides.annotated.tsv

set -euo pipefail

NEW="${1:?Missing new file path}"
REF="${2:?Missing reference file path}"

if [[ ! -f "$NEW" ]]; then echo "ERROR: new file not found: $NEW"; exit 1; fi
if [[ ! -f "$REF" ]]; then echo "ERROR: reference file not found: $REF"; exit 1; fi

SEP="============================================================"

# --- 1. Basic stats -----------------------------------------------------------
echo "$SEP"
echo "1. BASIC STATS"
echo "$SEP"
echo "  New file:      $(du -h "$NEW" | cut -f1), $(wc -l < "$NEW") lines"
echo "  Reference:     $(du -h "$REF" | cut -f1), $(wc -l < "$REF") lines"

# --- 2. MD5 ------------------------------------------------------------------
echo
echo "$SEP"
echo "2. MD5 CHECK"
echo "$SEP"
new_md5=$(md5sum "$NEW" | awk '{print $1}')
ref_md5=$(md5sum "$REF" | awk '{print $1}')
echo "  New:       $new_md5"
echo "  Reference: $ref_md5"
if [[ "$new_md5" == "$ref_md5" ]]; then
  echo "  VERDICT:   EXACT MATCH (binary identical)"
  exit 0
fi
echo "  VERDICT:   DIFFER (checking details below)"

# --- 3. Header comparison ----------------------------------------------------
echo
echo "$SEP"
echo "3. HEADER COMPARISON"
echo "$SEP"
if diff <(head -1 "$NEW" | tr '\t' '\n') <(head -1 "$REF" | tr '\t' '\n') >/dev/null 2>&1; then
  ncols=$(head -1 "$NEW" | tr '\t' '\n' | wc -l)
  echo "  Headers: IDENTICAL ($ncols columns)"
else
  echo "  Headers: DIFFER"
  diff <(head -1 "$NEW" | tr '\t' '\n') <(head -1 "$REF" | tr '\t' '\n') || true
fi

# --- 4. Peptide ID overlap ----------------------------------------------------
echo
echo "$SEP"
echo "4. PEPTIDE ID OVERLAP"
echo "$SEP"

awk -F'\t' 'NR>1{print $1}' "$NEW" | sort > /tmp/_cps_new_ids
awk -F'\t' 'NR>1{print $1}' "$REF" | sort > /tmp/_cps_ref_ids

new_uniq=$(wc -l < /tmp/_cps_new_ids)
ref_uniq=$(wc -l < /tmp/_cps_ref_ids)
inter=$(comm -12 /tmp/_cps_new_ids /tmp/_cps_ref_ids | wc -l)
only_new=$(comm -23 /tmp/_cps_new_ids /tmp/_cps_ref_ids | wc -l)
only_ref=$(comm -13 /tmp/_cps_new_ids /tmp/_cps_ref_ids | wc -l)

echo "  New file unique IDs:     $new_uniq"
echo "  Reference unique IDs:    $ref_uniq"
echo "  Intersection:            $inter"
echo "  Only in new:             $only_new"
echo "  Only in reference:       $only_ref"
pct=$(awk "BEGIN {printf \"%.2f\", ($inter / ($new_uniq + 0.0001)) * 100}")
echo "  Overlap:                 ${pct}%"

# --- 5. Column-by-column comparison (Python) -----------------------------------
echo
echo "$SEP"
echo "5. COLUMN-BY-COLUMN COMPARISON"
echo "$SEP"

python3 << PYEOF
import sys

# Read headers
with open("$NEW") as f:
    header_new = f.readline().rstrip('\n').split('\t')
with open("$REF") as f:
    header_ref = f.readline().rstrip('\n').split('\t')

assert header_new == header_ref, "Headers differ!"
ncols = len(header_new)

# Load data rows into dicts keyed by peptide_id (column 0)
def load_dict(path):
    d = {}
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.rstrip('\n').split('\t')
            d[parts[0]] = parts
    return d

new_dict = load_dict("$NEW")
ref_dict = load_dict("$REF")

diff_count = [0] * ncols
first_new = [""] * ncols
first_ref = [""] * ncols

common_ids = set(new_dict) & set(ref_dict)

for pid in common_ids:
    nr = new_dict[pid]
    rr = ref_dict[pid]
    for i in range(ncols):
        vn = nr[i] if i < len(nr) else ""
        vr = rr[i] if i < len(rr) else ""
        if vn != vr:
            if diff_count[i] == 0:
                first_new[i] = vn
                first_ref[i] = vr
            diff_count[i] += 1

same = 0
differ = 0
for i in range(ncols):
    if diff_count[i] > 0:
        print(f"  [{header_new[i]}] {diff_count[i]} rows differ")
        print(f"      new:  {first_new[i][:80]}")
        print(f"      ref:  {first_ref[i][:80]}")
        differ += 1
    else:
        same += 1

print()
print(f"  Identical columns:  {same}")
print(f"  Differing columns:  {differ}")
PYEOF

# --- 6. Summary ----------------------------------------------------------------
echo
echo "$SEP"
echo "6. SUMMARY"
echo "$SEP"
echo "  File size:      $(du -h "$NEW" | cut -f1) vs $(du -h "$REF" | cut -f1)"
echo "  Rows:           $(wc -l < "$NEW") vs $(wc -l < "$REF")"
echo "  MD5 match:      $([[ "$new_md5" == "$ref_md5" ]] && echo YES || echo NO)"
echo "  Peptide overlap: ${pct}% ($inter / $new_uniq)"
echo

# Cleanup
rm -f /tmp/_cps_new_ids /tmp/_cps_ref_ids
