#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: run_wes_qc.sh --bam FILE --outdir DIR [options]

Options:
  --sample-id ID       Sample label (default: BAM basename)
  --capture-bed FILE   Assay-specific capture BED; enables formal capture rate
  --gencode-gtf FILE   GENCODE GTF used only to build a CDS coverage proxy
  --samtools FILE      samtools executable
  --threads N          Worker threads (default: 16)
  --reuse-existing     Reuse completed raw metrics in OUTDIR
  --skip-stats         Skip optional insert-size/error-rate samtools stats pass
EOF
}

BAM=""; OUTDIR=""; SAMPLE_ID=""; CAPTURE_BED=""; GTF=""
SAMTOOLS="${SAMTOOLS:-samtools}"; THREADS="${THREADS:-16}"
REUSE_EXISTING=false; SKIP_STATS=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM=$2; shift 2 ;;
    --outdir) OUTDIR=$2; shift 2 ;;
    --sample-id) SAMPLE_ID=$2; shift 2 ;;
    --capture-bed) CAPTURE_BED=$2; shift 2 ;;
    --gencode-gtf) GTF=$2; shift 2 ;;
    --samtools) SAMTOOLS=$2; shift 2 ;;
    --threads) THREADS=$2; shift 2 ;;
    --reuse-existing) REUSE_EXISTING=true; shift ;;
    --skip-stats) SKIP_STATS=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -f "$BAM" ]] || { echo "Missing BAM: $BAM" >&2; exit 2; }
[[ -n "$OUTDIR" ]] || { echo "--outdir is required" >&2; exit 2; }
command -v "$SAMTOOLS" >/dev/null 2>&1 || [[ -x "$SAMTOOLS" ]] || {
  echo "samtools is not executable: $SAMTOOLS" >&2; exit 2;
}
SAMPLE_ID=${SAMPLE_ID:-$(basename "$BAM" .bam)}
mkdir -p "$OUTDIR/raw"

TARGET_BED=""
TARGET_DEFINITION="UNASSESSED"
CAPTURE_STATUS="UNASSESSED_CAPTURE_BED_MISSING"
if [[ -n "$CAPTURE_BED" ]]; then
  [[ -f "$CAPTURE_BED" ]] || { echo "Missing capture BED: $CAPTURE_BED" >&2; exit 2; }
  TARGET_BED="$CAPTURE_BED"
  TARGET_DEFINITION="ASSAY_CAPTURE_BED"
  CAPTURE_STATUS="ASSESSED"
elif [[ -n "$GTF" ]]; then
  [[ -f "$GTF" ]] || { echo "Missing GENCODE GTF: $GTF" >&2; exit 2; }
  TARGET_BED="$OUTDIR/gencode_cds_proxy.merged.bed"
  if [[ ! -s "$TARGET_BED" ]]; then
    awk 'BEGIN{OFS="\t"} $0 !~ /^#/ && $3=="CDS" {
      c=$1; if (c=="MT") c="M";
      if (c ~ /^([0-9]+|X|Y|M)$/) print "chr" c,$4-1,$5
    }' "$GTF" | LC_ALL=C sort -k1,1V -k2,2n -k3,3n | \
      awk 'BEGIN{OFS="\t"}
        NR==1 {c=$1;s=$2;e=$3;next}
        $1==c && $2<=e {if($3>e)e=$3;next}
        {print c,s,e;c=$1;s=$2;e=$3}
        END{if(NR)print c,s,e}' > "$TARGET_BED"
  fi
  TARGET_DEFINITION="GENCODE_CDS_PROXY_NOT_ASSAY_CAPTURE_BED"
fi

"$SAMTOOLS" quickcheck -v "$BAM"
if ! $REUSE_EXISTING || [[ ! -s "$OUTDIR/raw/samtools.flagstat.txt" ]]; then
  "$SAMTOOLS" flagstat -@ "$THREADS" "$BAM" > "$OUTDIR/raw/samtools.flagstat.txt"
fi
if $SKIP_STATS; then
  : > "$OUTDIR/raw/samtools.stats.txt"
elif ! $REUSE_EXISTING || [[ ! -s "$OUTDIR/raw/samtools.stats.txt" ]]; then
  "$SAMTOOLS" stats -@ "$THREADS" "$BAM" > "$OUTDIR/raw/samtools.stats.txt"
fi

if [[ -n "$TARGET_BED" ]]; then
  if ! $REUSE_EXISTING || [[ ! -s "$OUTDIR/raw/on_target_primary_mapped_reads.txt" ]]; then
    "$SAMTOOLS" view -@ "$THREADS" -c -F 0x904 -L "$TARGET_BED" "$BAM" \
      > "$OUTDIR/raw/on_target_primary_mapped_reads.txt"
  fi
  if ! $REUSE_EXISTING || [[ ! -s "$OUTDIR/raw/target_depth_summary.tsv" ]]; then
    "$SAMTOOLS" depth -a -Q 20 -q 20 -G 0xF04 -b "$TARGET_BED" "$BAM" | \
      awk -v out="$OUTDIR/raw/target_depth_summary.tsv" 'BEGIN{OFS="\t"; print "target_bases","mean_coverage","pct_1x","pct_10x","pct_20x","pct_30x","pct_50x","pct_100x" > out}
        {n++;s+=$3;if($3>=1)a++;if($3>=10)b++;if($3>=20)c++;if($3>=30)d++;if($3>=50)e++;if($3>=100)f++}
        END{if(n==0)n=1;printf "%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",n,s/n,100*a/n,100*b/n,100*c/n,100*d/n,100*e/n,100*f/n >> out}'
  fi
fi

python3 - "$SAMPLE_ID" "$BAM" "$OUTDIR" "$TARGET_BED" "$TARGET_DEFINITION" "$CAPTURE_STATUS" <<'PY'
import csv, json, re, sys
from datetime import datetime, timezone
from pathlib import Path

sample, bam, outdir, target_bed, target_definition, capture_status = sys.argv[1:]
out = Path(outdir)
flagstat = (out / "raw/samtools.flagstat.txt").read_text()
stats_lines = (out / "raw/samtools.stats.txt").read_text().splitlines()

def flag_value(label):
    match = re.search(rf"^(\d+) \+ (\d+) {re.escape(label)}(?:\s|$)", flagstat, re.M)
    return sum(map(int, match.groups())) if match else 0

sn = {}
for line in stats_lines:
    if line.startswith("SN\t"):
        parts = line.split("\t", 2)
        if len(parts) == 3:
            sn[parts[1].rstrip(":")] = parts[2]

total = flag_value("in total (QC-passed reads + QC-failed reads)")
primary = flag_value("primary") or total
primary_mapped = flag_value("primary mapped") or flag_value("mapped")
primary_duplicates = flag_value("primary duplicates") or flag_value("duplicates")
properly_paired = flag_value("properly paired")
on_target = 0
on_target_path = out / "raw/on_target_primary_mapped_reads.txt"
if on_target_path.exists():
    on_target = int((on_target_path.read_text().strip() or "0"))

depth = {}
depth_path = out / "raw/target_depth_summary.tsv"
if depth_path.exists():
    with depth_path.open(newline="") as handle:
        depth = next(csv.DictReader(handle, delimiter="\t"), {})

def pct(num, den):
    return round(100.0 * num / den, 4) if den else None

row = {
    "sample_id": sample,
    "data_type": "tumor_WES",
    "bam": bam,
    "bam_size_bytes": Path(bam).stat().st_size,
    "qc_generated_at": datetime.now(timezone.utc).isoformat(),
    "qc_status": "PASS_WITH_CAPTURE_RATE_UNASSESSED" if capture_status != "ASSESSED" else "ASSESSED",
    "total_reads": total,
    "primary_reads": primary,
    "primary_mapped_reads": primary_mapped,
    "primary_mapping_rate_pct": pct(primary_mapped, primary),
    "properly_paired_reads": properly_paired,
    "properly_paired_rate_pct": pct(properly_paired, primary),
    "primary_duplicate_reads": primary_duplicates,
    "duplicate_rate_pct": pct(primary_duplicates, primary_mapped),
    "target_definition": target_definition,
    "target_bed": target_bed,
    "on_target_primary_mapped_reads": on_target if target_bed else None,
    "on_target_rate_pct": pct(on_target, primary_mapped) if target_bed else None,
    "capture_rate_status": capture_status,
    "formal_capture_rate_pct": pct(on_target, primary_mapped) if capture_status == "ASSESSED" else None,
    "coverage_filter": "MAPQ>=20;BaseQ>=20;exclude_unmapped,secondary,qcfail,duplicate,supplementary",
    "target_bases": depth.get("target_bases"),
    "mean_target_coverage": depth.get("mean_coverage"),
    "pct_target_bases_1x": depth.get("pct_1x"),
    "pct_target_bases_10x": depth.get("pct_10x"),
    "pct_target_bases_20x": depth.get("pct_20x"),
    "pct_target_bases_30x": depth.get("pct_30x"),
    "pct_target_bases_50x": depth.get("pct_50x"),
    "pct_target_bases_100x": depth.get("pct_100x"),
    "insert_size_average": sn.get("insert size average"),
    "insert_size_standard_deviation": sn.get("insert size standard deviation"),
    "error_rate": sn.get("error rate"),
}

with (out / "wes_qc.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerow({k: "" if v is None else v for k, v in row.items()})
(out / "wes_qc.json").write_text(json.dumps(row, indent=2) + "\n")

def show(value, suffix=""):
    return "UNASSESSED" if value in (None, "") else f"{value}{suffix}"

summary = f"""# Independent tumor WES QC

- Sample: `{sample}`
- QC status: **{row['qc_status']}**
- Primary mapping rate: **{show(row['primary_mapping_rate_pct'], '%')}**
- Duplicate rate: **{show(row['duplicate_rate_pct'], '%')}**
- Mean target coverage: **{show(row['mean_target_coverage'], 'x')}**
- Target bases >=20x: **{show(row['pct_target_bases_20x'], '%')}**
- Target bases >=30x: **{show(row['pct_target_bases_30x'], '%')}**
- Proxy on-target rate: **{show(row['on_target_rate_pct'], '%')}**
- Formal assay capture rate: **{show(row['formal_capture_rate_pct'], '%')}**

The current coverage and on-target values use `{target_definition}`. They are
coverage proxies, not formal assay capture metrics. Supply the original capture
kit BED and rerun this script to populate `formal_capture_rate_pct`.
"""
(out / "wes_qc_summary.md").write_text(summary)
PY

echo "WES QC completed: $OUTDIR/wes_qc.tsv"
