#!/usr/bin/env bash
# SNAF: AltAnalyze junction counts -> snaf.initialize / prediction
#
# Single-sample rules:
#   - Mount the real STAR BAM once. Never create sample_replicate from the
#     same file (that is not independent RNA support).
#   - After docker exits, replace bind-mount leftovers with a symlink to the
#     real BAM so the work dir is not 0-byte placeholders.
#   - Prefer counts.original.full (pruned is often empty without two groups).
#   - If STAR SJ.out.tab is available, keep only junctions whose AltAnalyze
#     coords match STAR splice junctions (off-by-one aware). GTEx comparison
#     still uses AltAnalyze counts; STAR is the RNA-support gate.
#   - snaf_candidates.tsv genomic start<=end (donor-acceptor order is swapped).
set -euo pipefail

usage() {
  echo "Usage: $0 --bam BAM --hla-file FILE --sample-id ID --db-dir DIR --outdir DIR [--threads N] [--star-sj SJ.out.tab]" >&2
}

BAM=""; HLA=""; SAMPLE=""; DB=""; OUT=""; THREADS=8
IMAGE="${NEOAG_ALTANALYZE_IMAGE:-neoag-altanalyze:snaf}"
SNAF_PY="${SNAF_PYTHON:-python}"
SKIP_ALTANALYZE="${SKIP_ALTANALYZE:-0}"
STAR_SJ="${STAR_SJ:-}"
KEEP_UNVERIFIED="${SNAF_KEEP_UNVERIFIED_JUNCTIONS:-0}"
RUN_AS_HOST_USER="${NEOAG_ALTANALYZE_RUN_AS_HOST_USER:-1}"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2;;
    --hla-file) HLA="$2"; shift 2;;
    --sample-id) SAMPLE="$2"; shift 2;;
    --db-dir) DB="$2"; shift 2;;
    --outdir) OUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --altanalyze-image) IMAGE="$2"; shift 2;;
    --star-sj) STAR_SJ="$2"; shift 2;;
    --skip-altanalyze) SKIP_ALTANALYZE=1; shift;;
    --keep-unverified-junctions) KEEP_UNVERIFIED=1; shift;;
    *) usage; exit 2;;
  esac
done
for value in "$BAM" "$HLA" "$SAMPLE" "$DB" "$OUT"; do [[ -n "$value" ]] || { usage; exit 2; }; done
[[ -s "$BAM" && -s "$HLA" ]] || { echo "ERROR: missing BAM or HLA file" >&2; exit 2; }
[[ -s "$DB/controls/GTEx_junction_counts.h5ad" ]] || { echo "ERROR: incomplete SNAF database: $DB" >&2; exit 2; }

if [[ -z "$STAR_SJ" ]]; then
  _sj_guess="$(dirname "$BAM")/SJ.out.tab"
  if [[ -s "$_sj_guess" ]]; then
    STAR_SJ="$_sj_guess"
  fi
fi
if [[ -n "$STAR_SJ" && ! -s "$STAR_SJ" ]]; then
  echo "ERROR: --star-sj is not a readable file: $STAR_SJ" >&2
  exit 2
fi

WORK="$OUT/altanalyze_work"
FULL="$WORK/altanalyze_output/ExpressionInput/counts.original.full.txt"
PRUNED="$WORK/altanalyze_output/ExpressionInput/counts.original.pruned.txt"
ORIG="$WORK/altanalyze_output/ExpressionInput/counts.original.txt"
SNAF_INPUT="$OUT/altanalyze_counts.snaf_input.tsv"
STAR_REPORT="$OUT/star_sj_overlap.tsv"

mkdir -p "$WORK/bam" "$OUT/assets"

normalize_work_ownership() {
  # Older runs may have root-owned files because Docker used its default user.
  # Repair those files from inside the same image before the host resumes.
  if find "$WORK" -maxdepth 4 ! -user "$(id -u)" -print -quit 2>/dev/null | grep -q .; then
    docker run --rm -v "$WORK:/mnt" --entrypoint /bin/sh "$IMAGE" -c \
      "chown -R $(id -u):$(id -g) /mnt 2>/dev/null || chmod -R a+rwX /mnt"
  fi
}

provenance_symlink_bam() {
  mkdir -p "$WORK/bam"
  rm -f \
    "$WORK/bam/${SAMPLE}.bam" \
    "$WORK/bam/${SAMPLE}.bam.bai" \
    "$WORK/bam/${SAMPLE}_replicate.bam" \
    "$WORK/bam/${SAMPLE}_replicate.bam.bai"
  # AltAnalyze runs in a container and may recreate this file as root.
  # Remove it through the caller-owned work directory before rewriting it.
  rm -f "$WORK/samples.txt"
  ln -sfn "$BAM" "$WORK/bam/${SAMPLE}.bam"
  if [[ -s "${BAM}.bai" ]]; then
    ln -sfn "${BAM}.bai" "$WORK/bam/${SAMPLE}.bam.bai"
  fi
  printf '%s.bam\n' "$SAMPLE" > "$WORK/samples.txt"
}

if [[ "$SKIP_ALTANALYZE" != "1" ]]; then
  command -v docker >/dev/null || { echo "ERROR: docker is required for AltAnalyze" >&2; exit 2; }
  docker image inspect "$IMAGE" >/dev/null 2>&1 || { echo "ERROR: missing AltAnalyze image: $IMAGE" >&2; exit 2; }
  normalize_work_ownership
  # Drop leftover replicate / empty bind-mount files before identify bam.
  rm -f \
    "$WORK/bam/${SAMPLE}.bam" \
    "$WORK/bam/${SAMPLE}.bam.bai" \
    "$WORK/bam/${SAMPLE}_replicate.bam" \
    "$WORK/bam/${SAMPLE}_replicate.bam.bai"
  rm -f "$WORK/samples.txt"
  mounts=(
    -v "$WORK:/mnt"
    -v "$BAM:/mnt/bam/${SAMPLE}.bam:ro"
  )
  if [[ -s "${BAM}.bai" ]]; then
    mounts+=(
      -v "${BAM}.bai:/mnt/bam/${SAMPLE}.bam.bai:ro"
    )
  fi
  docker_user=()
  if [[ "$RUN_AS_HOST_USER" == "1" ]]; then
    docker_user=(--user "$(id -u):$(id -g)" -e HOME=/tmp)
  fi
  echo "==> AltAnalyze identify bam sample=${SAMPLE} (single BAM, no fake replicate)"
  set +e
  docker run --rm "${docker_user[@]}" "${mounts[@]}" "$IMAGE" identify bam "$THREADS"
  aa_rc=$?
  set -e
  normalize_work_ownership
  provenance_symlink_bam
  if [[ ! -s "$FULL" ]] || [[ "$(wc -l < "$FULL")" -le 2 ]]; then
    echo "ERROR: AltAnalyze docker rc=${aa_rc} and no usable counts.original.full" >&2
    exit 1
  fi
  if [[ "$aa_rc" != "0" ]]; then
    echo "==> WARN: AltAnalyze docker rc=${aa_rc} (single-sample PSI/EventAnnotation often fails); continuing with counts.original.full"
  fi
else
  # Resume path: replace 0-byte leftovers with a real BAM symlink.
  command -v docker >/dev/null || { echo "ERROR: docker is required to repair AltAnalyze work ownership" >&2; exit 2; }
  docker image inspect "$IMAGE" >/dev/null 2>&1 || { echo "ERROR: missing AltAnalyze image: $IMAGE" >&2; exit 2; }
  normalize_work_ownership
  provenance_symlink_bam
fi

# Prefer full matrix. pruned is often header-only for a single sample (dPSI=0).
MATRIX=""
if [[ -s "$FULL" ]] && [[ "$(wc -l < "$FULL")" -gt 2 ]]; then
  MATRIX="$FULL"
elif [[ -s "$PRUNED" ]] && [[ "$(wc -l < "$PRUNED")" -gt 2 ]]; then
  MATRIX="$PRUNED"
elif [[ -s "$ORIG" ]] && [[ "$(wc -l < "$ORIG")" -gt 2 ]]; then
  MATRIX="$ORIG"
else
  echo "ERROR: no usable AltAnalyze counts matrix under $WORK/altanalyze_output/ExpressionInput" >&2
  ls -la "$WORK/altanalyze_output/ExpressionInput" 2>/dev/null || true
  exit 1
fi
echo "==> SNAF matrix source: $MATRIX ($(wc -l < "$MATRIX") lines)"
if [[ -n "$STAR_SJ" ]]; then
  echo "==> STAR SJ gate: $STAR_SJ keep_unverified=$KEEP_UNVERIFIED"
else
  echo "==> STAR SJ gate: skipped (no SJ.out.tab next to BAM and --star-sj not set)"
fi

# Build single-sample snaf_input.tsv and optionally gate on STAR SJ.out.tab.
"$SNAF_PY" - "$MATRIX" "$SNAF_INPUT" "$SAMPLE" "$ORIG" "$STAR_SJ" "$STAR_REPORT" "$KEEP_UNVERIFIED" <<'PY'
import sys
from pathlib import Path
import pandas as pd

src, dst, sample, orig_path, star_path, report_path, keep_unverified = sys.argv[1:8]
keep_unverified = keep_unverified == "1"


def parse_coord(text):
    text = str(text)
    if ":" not in text or "-" not in text:
        return None
    chrom, rest = text.split(":", 1)
    rest = rest.split("(", 1)[0]
    start_s, end_s = (rest.split("-", 1) + [""])[:2]
    try:
        start, end = int(start_s), int(end_s)
    except ValueError:
        return None
    return chrom, min(start, end), max(start, end)


def load_coord_map(path):
    mapping = {}
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        return mapping
    df = pd.read_csv(p, sep="\t", index_col=0)
    for raw in df.index.astype(str):
        uid, _, coord = raw.partition("=")
        parsed = parse_coord(coord)
        if parsed:
            mapping[uid.split("=", 1)[0]] = parsed
    return mapping


def load_star(path):
    keys = {}
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        return keys
    with p.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 7:
                continue
            try:
                start, end, unique = int(parts[1]), int(parts[2]), int(parts[6])
            except ValueError:
                continue
            chrom = parts[0]
            a, b = min(start, end), max(start, end)
            keys[(chrom, a, b)] = keys.get((chrom, a, b), 0) + unique
    return keys


def star_hit(star, chrom, start, end):
    a, b = min(start, end), max(start, end)
    for ds in (-1, 0, 1):
        for de in (-1, 0, 1):
            key = (chrom, a + ds, b + de)
            if key in star:
                return key, star[key]
    return None, 0


df = pd.read_csv(src, sep="\t", index_col=0)
df.index = [str(i).split("=", 1)[0] for i in df.index]
df = df.loc[~pd.Index(df.index).duplicated(keep="first")]
cols = [c for c in df.columns if "replicate" not in str(c).lower()]
if not cols:
    cols = list(df.columns)
if not cols:
    raise SystemExit(f"ERROR: matrix has no sample columns: {src}")
df = df.loc[:, [cols[0]]].copy()
df.columns = [sample]
df = df.loc[(df.fillna(0) > 0).any(axis=1)]
if df.empty:
    raise SystemExit(f"ERROR: matrix empty after filtering: {src}")

coord_map = load_coord_map(orig_path)
n_in = len(df)
n_star = n_no_coord = n_drop = 0
star = load_star(star_path) if star_path else {}
report_rows = []

if star:
    keep = []
    for uid in df.index.astype(str):
        parsed = coord_map.get(uid)
        if parsed is None:
            n_no_coord += 1
            report_rows.append((uid, "", "", "", "no_coord", 0))
            if keep_unverified:
                keep.append(uid)
            else:
                n_drop += 1
            continue
        chrom, start, end = parsed
        key, unique = star_hit(star, chrom, start, end)
        if key is None:
            n_drop += 1
            report_rows.append((uid, chrom, start, end, "no_star_sj", 0))
            if keep_unverified:
                keep.append(uid)
            continue
        n_star += 1
        report_rows.append((uid, chrom, start, end, "star_sj", unique))
        keep.append(uid)
    df = df.loc[keep]
    if df.empty:
        raise SystemExit(
            "ERROR: STAR SJ gate removed every junction. "
            "Check AltAnalyze coords vs SJ.out.tab, or pass --keep-unverified-junctions."
        )

dst = Path(dst)
dst.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(dst, sep="\t")
report = Path(report_path)
if report_rows:
    pd.DataFrame(
        report_rows,
        columns=["uid", "chrom", "start", "end", "star_status", "star_unique_reads"],
    ).to_csv(report, sep="\t", index=False)
print(
    f"wrote {dst} shape={df.shape} input_nonzero={n_in} "
    f"star_kept={n_star} dropped={n_drop} no_coord={n_no_coord}"
)
PY

# Matrix changed: do not resume an old after_prediction.p from the fake-replicate run.
rm -f "$OUT/after_prediction.p"

export NEOAG_SNAF_OUTDIR="$OUT"
export NEOAG_SNAF_MATRIX="$SNAF_INPUT"
export NEOAG_SNAF_DB="$DB"
export NEOAG_SNAF_HLA_FILE="$HLA"
export NEOAG_SNAF_SAMPLE_ID="$SAMPLE"
export NEOAG_SNAF_CORES="$THREADS"
SNAF_WORKFLOW="$(cd "$(dirname "$0")" && pwd -P)/snaf_sample_workflow.py"
# SNAF's dashboard cleanup expects an assets directory under the current
# working directory. Run from the caller-owned output directory so cleanup
# never targets the repository or fails after an otherwise successful run.
(
  cd "$OUT"
  "$SNAF_PY" "$SNAF_WORKFLOW"
)
# Header-only candidates still count as success (stage-3 may be empty).
[[ -f "$OUT/snaf_candidates.tsv" ]] || { echo "ERROR: missing snaf_candidates.tsv" >&2; exit 1; }
date -Is > "$OUT/.snaf_complete"
date -Is > "$OUT/.snaf.done"
