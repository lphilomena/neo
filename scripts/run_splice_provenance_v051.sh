#!/usr/bin/env bash
# NeoAg v0.5.1 production driver: RNA-driven, DNA-causal and normal-background chains.
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/run_splice_provenance_v051.sh \
    --sample-id SAMPLE001 --outdir results/SAMPLE001/splice_v051 \
    --junctions regtools_junctions.tsv \
    --immunopepper-meta immunopepper_meta.tsv \
    --mopepgen-fasta mopepgen_peptides.fasta \
    --mopepgen-provenance-map mopepgen_map.tsv \
    --splice2neo splice2neo_export.tsv [options]

Three-pass contract:
  1. phase1_registry: construct exact events, peptide origins and query IDs.
  2. external evidence: EasyQuant, k4neo and pVACsplice use phase-1 maps.
  3. pre_pvacbind/final: rebuild all three chains, then optionally import pVACbind.

Core repeatable inputs:
  --spladder-gff3 PATH        --spladder-txt PATH
  --irfinder PATH             --immunopepper-meta PATH
  --immunopepper-kmers PATH   --mopepgen-fasta PATH
  --mopepgen-gvf PATH         --mopepgen-provenance-map PATH
  --splice2neo PATH
  --normal-junctions PATH     --normal-coverage PATH

Precomputed external outputs:
  --easyquant PATH
  --k4neo-healthy-sample-rate PATH
  --k4neo-annotated PATH
  --k4neo-uniqueness PATH
  --pvacsplice PATH
  --pvacbind PATH

Run EasyQuant from phase-1 query table:
  --run-easyquant --easyquant-bam RNA.bam
    or --run-easyquant --easyquant-fq1 R1.fq.gz --easyquant-fq2 R2.fq.gz

Run k4neo from phase-1 query table:
  --run-k4neo --k4neo-database DB --k4neo-index MANIFEST --k4neo-license-accepted

Run pVACsplice:
  --run-pvacsplice --pvacsplice-junctions REGTOOLS_CIS.tsv \
  --annotated-vcf VEP.vcf.gz --ref-fasta GRCh38.fa --gtf gencode.gtf \
  --hla HLA-A*02:01,HLA-B*07:02

Run pVACbind on all resolved ORFs:
  --hla LIST | --hla-file PATH
  --pvacbind-algorithms MHCflurry,MHCflurryEL

General:
  --genome-build GRCh38       --disease-profile default
  --tool-version TOOL=VERSION (repeatable)
  --critical-tissue NAME      (repeatable)
  --strict                    --overwrite
  --skip-pvacbind
USAGE
}

SAMPLE=""; OUTDIR=""; BUILD="GRCh38"; PROFILE="default"; STRICT=0; OVERWRITE=0
JUNCTIONS=""; JUNCTION_COORD="auto"; STAR_JUNCTIONS=""; IR_COORD="intron_1based_closed"; NORMAL_COORD="auto"
SPL_GFF3=(); SPL_TXT=(); IR=(); IMM_META=(); IMM_KMER=(); MO_FASTA=(); MO_GVF=(); MO_MAP=(); S2N=()
NORMAL_JUNC=(); NORMAL_COV=(); CRITICAL=(); VERSIONS=()
EASY=(); K4_HEALTHY=(); K4_ANNOTATED=(); K4_UNIQUE=(); PVS=(); PVB=()
RUN_EASY=0; EASY_BAM=""; EASY_FQ1=""; EASY_FQ2=""; EASY_THREADS="4"
RUN_K4=0; K4_DB=""; K4_INDEX=""; K4_ACCEPT=0; K4_PREFIX="neoag"
RUN_PVS=0; PVS_JUNCTIONS=""; ANNOTATED_VCF=""; REF_FASTA=""; GTF=""; PVS_ALG="MHCflurry"; PVS_THREADS="4"
HLA=""; HLA_FILE=""; PVB_ALG="MHCflurry"; PVB_THREADS="4"; REF_PROTEOME=""; SKIP_PVB=0
while [[ $# -gt 0 ]]; do
 case "$1" in
  --sample-id) SAMPLE="$2"; shift 2;; --outdir) OUTDIR="$2"; shift 2;;
  --genome-build) BUILD="$2"; shift 2;; --disease-profile) PROFILE="$2"; shift 2;;
  --junctions) JUNCTIONS="$2"; shift 2;; --junction-coordinate-system) JUNCTION_COORD="$2"; shift 2;;
  --star-junctions) STAR_JUNCTIONS="$2"; shift 2;;
  --spladder-gff3) SPL_GFF3+=("$2"); shift 2;; --spladder-txt) SPL_TXT+=("$2"); shift 2;;
  --irfinder) IR+=("$2"); shift 2;; --irfinder-coordinate-system) IR_COORD="$2"; shift 2;;
  --immunopepper-meta) IMM_META+=("$2"); shift 2;; --immunopepper-kmers) IMM_KMER+=("$2"); shift 2;;
  --mopepgen-fasta) MO_FASTA+=("$2"); shift 2;; --mopepgen-gvf) MO_GVF+=("$2"); shift 2;;
  --mopepgen-provenance-map) MO_MAP+=("$2"); shift 2;; --splice2neo) S2N+=("$2"); shift 2;;
  --normal-junctions) NORMAL_JUNC+=("$2"); shift 2;; --normal-coverage) NORMAL_COV+=("$2"); shift 2;;
  --normal-coordinate-system) NORMAL_COORD="$2"; shift 2;; --critical-tissue) CRITICAL+=("$2"); shift 2;;
  --tool-version) VERSIONS+=("$2"); shift 2;;
  --easyquant) EASY+=("$2"); shift 2;; --run-easyquant) RUN_EASY=1; shift;;
  --easyquant-bam) EASY_BAM="$2"; shift 2;; --easyquant-fq1) EASY_FQ1="$2"; shift 2;; --easyquant-fq2) EASY_FQ2="$2"; shift 2;;
  --easyquant-threads) EASY_THREADS="$2"; shift 2;;
  --k4neo-healthy-sample-rate) K4_HEALTHY+=("$2"); shift 2;; --k4neo-annotated) K4_ANNOTATED+=("$2"); shift 2;;
  --k4neo-uniqueness) K4_UNIQUE+=("$2"); shift 2;; --run-k4neo) RUN_K4=1; shift;;
  --k4neo-database) K4_DB="$2"; shift 2;; --k4neo-index) K4_INDEX="$2"; shift 2;;
  --k4neo-prefix) K4_PREFIX="$2"; shift 2;; --k4neo-license-accepted) K4_ACCEPT=1; shift;;
  --pvacsplice) PVS+=("$2"); shift 2;; --run-pvacsplice) RUN_PVS=1; shift;;
  --pvacsplice-junctions) PVS_JUNCTIONS="$2"; shift 2;; --annotated-vcf) ANNOTATED_VCF="$2"; shift 2;;
  --ref-fasta) REF_FASTA="$2"; shift 2;; --gtf) GTF="$2"; shift 2;;
  --pvacsplice-algorithms) PVS_ALG="$2"; shift 2;; --pvacsplice-threads) PVS_THREADS="$2"; shift 2;;
  --pvacbind) PVB+=("$2"); shift 2;; --hla) HLA="$2"; shift 2;; --hla-file) HLA_FILE="$2"; shift 2;;
  --pvacbind-algorithms) PVB_ALG="$2"; shift 2;; --pvacbind-threads) PVB_THREADS="$2"; shift 2;;
  --reference-proteome) REF_PROTEOME="$2"; shift 2;; --skip-pvacbind) SKIP_PVB=1; shift;;
  --strict) STRICT=1; shift;; --overwrite) OVERWRITE=1; shift;;
  -h|--help) usage; exit 0;; *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2;;
 esac
done
[[ -n "$SAMPLE" && -n "$OUTDIR" ]] || { echo "ERROR: --sample-id and --outdir required" >&2; exit 2; }
if [[ "$RUN_K4" == 1 || ${#K4_HEALTHY[@]} -gt 0 || ${#K4_ANNOTATED[@]} -gt 0 || ${#K4_UNIQUE[@]} -gt 0 ]]; then
 [[ "$K4_ACCEPT" == 1 ]] || { echo "ERROR: k4neo use requires --k4neo-license-accepted" >&2; exit 2; }
fi
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${NEOAG_PYTHON:-$(command -v python3 || command -v python || true)}"
[[ -n "$PY" ]] || { echo "ERROR: python3 or python is required" >&2; exit 3; }
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
PHASE1="$OUTDIR/phase1_registry"; EXT="$OUTDIR/external"; PRE="$OUTDIR/pre_pvacbind"; FINAL="$OUTDIR/splice_provenance"
if [[ -e "$FINAL" && "$OVERWRITE" != 1 ]]; then echo "ERROR: final output exists; use --overwrite: $FINAL" >&2; exit 2; fi
if [[ "$OVERWRITE" == 1 ]]; then rm -rf "$PHASE1" "$PRE" "$FINAL"; fi
mkdir -p "$OUTDIR" "$EXT"

make_build_args() {
 local target="$1"; BUILD_ARGS=(-m neoag.splice.cli build --sample-id "$SAMPLE" --outdir "$target" --genome-build "$BUILD" --disease-profile "$PROFILE" --junction-coordinate-system "$JUNCTION_COORD" --irfinder-coordinate-system "$IR_COORD" --normal-coordinate-system "$NORMAL_COORD")
 [[ -n "$JUNCTIONS" ]] && BUILD_ARGS+=(--junctions "$JUNCTIONS")
 [[ -n "$STAR_JUNCTIONS" ]] && BUILD_ARGS+=(--star-junctions "$STAR_JUNCTIONS")
 local x
 for x in "${SPL_GFF3[@]}"; do BUILD_ARGS+=(--spladder-gff3 "$x"); done
 for x in "${SPL_TXT[@]}"; do BUILD_ARGS+=(--spladder-txt "$x"); done
 for x in "${IR[@]}"; do BUILD_ARGS+=(--irfinder "$x"); done
 for x in "${IMM_META[@]}"; do BUILD_ARGS+=(--immunopepper-meta "$x"); done
 for x in "${IMM_KMER[@]}"; do BUILD_ARGS+=(--immunopepper-kmers "$x"); done
 for x in "${MO_FASTA[@]}"; do BUILD_ARGS+=(--mopepgen-fasta "$x"); done
 for x in "${MO_GVF[@]}"; do BUILD_ARGS+=(--mopepgen-gvf "$x"); done
 for x in "${MO_MAP[@]}"; do BUILD_ARGS+=(--mopepgen-provenance-map "$x"); done
 for x in "${S2N[@]}"; do BUILD_ARGS+=(--splice2neo "$x"); done
 for x in "${NORMAL_JUNC[@]}"; do BUILD_ARGS+=(--normal-junctions "$x"); done
 for x in "${NORMAL_COV[@]}"; do BUILD_ARGS+=(--normal-coverage "$x"); done
 for x in "${CRITICAL[@]}"; do BUILD_ARGS+=(--critical-tissue "$x"); done
 for x in "${VERSIONS[@]}"; do BUILD_ARGS+=(--tool-version "$x"); done
 [[ "$STRICT" == 1 ]] && BUILD_ARGS+=(--strict)
}

# Phase 1 establishes exact IDs and query maps. No externally generated result may precede this step.
make_build_args "$PHASE1"
"$PY" "${BUILD_ARGS[@]}" > "$OUTDIR/phase1.outputs.json"

if [[ "$RUN_EASY" == 1 ]]; then
 eq=(bash "$ROOT/scripts/run_easyquant_sample.sh" --query-table "$PHASE1/splice_easyquant_input.tsv" --outdir "$EXT/easyquant" --threads "$EASY_THREADS")
 if [[ -n "$EASY_BAM" ]]; then eq+=(--bam "$EASY_BAM"); else eq+=(--fq1 "$EASY_FQ1" --fq2 "$EASY_FQ2"); fi
 eq_list="$("${eq[@]}")"
 while IFS= read -r x; do [[ -n "$x" ]] && EASY+=("$x"); done < "$eq_list"
fi
if [[ "$RUN_K4" == 1 ]]; then
 k4=(bash "$ROOT/scripts/run_k4neo_sample.sh" --query-table "$PHASE1/splice_k4neo_input.tsv" --database "$K4_DB" --index "$K4_INDEX" --outdir "$EXT/k4neo" --prefix "$K4_PREFIX" --license-accepted)
 mapfile -t k4_lists < <("${k4[@]}")
 [[ ${#k4_lists[@]} -ge 1 && -s "${k4_lists[0]}" ]] && while IFS= read -r x; do [[ -n "$x" ]] && K4_HEALTHY+=("$x"); done < "${k4_lists[0]}"
 [[ ${#k4_lists[@]} -ge 2 && -s "${k4_lists[1]}" ]] && while IFS= read -r x; do [[ -n "$x" ]] && K4_ANNOTATED+=("$x"); done < "${k4_lists[1]}"
 [[ ${#k4_lists[@]} -ge 3 && -s "${k4_lists[2]}" ]] && while IFS= read -r x; do [[ -n "$x" ]] && K4_UNIQUE+=("$x"); done < "${k4_lists[2]}"
fi
if [[ "$RUN_PVS" == 1 ]]; then
 pvs=(bash "$ROOT/scripts/run_pvacsplice_sample.sh" --junctions "$PVS_JUNCTIONS" --annotated-vcf "$ANNOTATED_VCF" --sample-id "$SAMPLE" --algorithms "$PVS_ALG" --outdir "$EXT/pvacsplice" --ref-fasta "$REF_FASTA" --gtf "$GTF" --threads "$PVS_THREADS")
 [[ -n "$HLA" ]] && pvs+=(--hla "$HLA")
 [[ -n "$HLA_FILE" ]] && pvs+=(--hla-file "$HLA_FILE")
 pvs_list="$("${pvs[@]}")"
 while IFS= read -r x; do [[ -n "$x" ]] && PVS+=("$x"); done < "$pvs_list"
fi

# Phase 2 imports external evidence strictly against the phase-1 maps.
make_build_args "$PRE"
for x in "${EASY[@]}"; do BUILD_ARGS+=(--easyquant "$x"); done
[[ ${#EASY[@]} -gt 0 ]] && BUILD_ARGS+=(--easyquant-query-map "$PHASE1/splice_easyquant_query_map.tsv")
for x in "${PVS[@]}"; do BUILD_ARGS+=(--pvacsplice "$x"); done
for x in "${K4_HEALTHY[@]}"; do BUILD_ARGS+=(--k4neo-healthy-sample-rate "$x"); done
for x in "${K4_ANNOTATED[@]}"; do BUILD_ARGS+=(--k4neo-annotated "$x"); done
for x in "${K4_UNIQUE[@]}"; do BUILD_ARGS+=(--k4neo-uniqueness "$x"); done
if [[ ${#K4_HEALTHY[@]} -gt 0 || ${#K4_ANNOTATED[@]} -gt 0 || ${#K4_UNIQUE[@]} -gt 0 ]]; then
 BUILD_ARGS+=(--k4neo-query-map "$PHASE1/splice_k4neo_query_map.tsv" --k4neo-license-accepted)
fi
"$PY" "${BUILD_ARGS[@]}" > "$OUTDIR/pre_pvacbind.outputs.json"

if [[ "$SKIP_PVB" != 1 && ${#PVB[@]} -eq 0 && ( -n "$HLA" || -n "$HLA_FILE" ) && -s "$PRE/splice_pvacbind_input.fasta" ]]; then
 pvb=(bash "$ROOT/scripts/run_pvacbind_sample.sh" --fasta "$PRE/splice_pvacbind_input.fasta" --sample-id "$SAMPLE" --algorithms "$PVB_ALG" --outdir "$EXT/pvacbind" --threads "$PVB_THREADS")
 [[ -n "$HLA" ]] && pvb+=(--hla "$HLA")
 [[ -n "$HLA_FILE" ]] && pvb+=(--hla-file "$HLA_FILE")
 [[ -n "$REF_PROTEOME" ]] && pvb+=(--reference-proteome "$REF_PROTEOME")
 pvb_list="$("${pvb[@]}")"
 while IFS= read -r x; do [[ -n "$x" ]] && PVB+=("$x"); done < "$pvb_list"
fi

if [[ ${#PVB[@]} -gt 0 ]]; then
 make_build_args "$FINAL"
 for x in "${EASY[@]}"; do BUILD_ARGS+=(--easyquant "$x"); done
 [[ ${#EASY[@]} -gt 0 ]] && BUILD_ARGS+=(--easyquant-query-map "$PHASE1/splice_easyquant_query_map.tsv")
 for x in "${PVS[@]}"; do BUILD_ARGS+=(--pvacsplice "$x"); done
 for x in "${K4_HEALTHY[@]}"; do BUILD_ARGS+=(--k4neo-healthy-sample-rate "$x"); done
 for x in "${K4_ANNOTATED[@]}"; do BUILD_ARGS+=(--k4neo-annotated "$x"); done
 for x in "${K4_UNIQUE[@]}"; do BUILD_ARGS+=(--k4neo-uniqueness "$x"); done
 if [[ ${#K4_HEALTHY[@]} -gt 0 || ${#K4_ANNOTATED[@]} -gt 0 || ${#K4_UNIQUE[@]} -gt 0 ]]; then BUILD_ARGS+=(--k4neo-query-map "$PHASE1/splice_k4neo_query_map.tsv" --k4neo-license-accepted); fi
 for x in "${PVB[@]}"; do BUILD_ARGS+=(--pvacbind "$x"); done
 BUILD_ARGS+=(--pvacbind-fasta-map "$PRE/splice_pvacbind_fasta_map.tsv")
 "$PY" "${BUILD_ARGS[@]}" > "$OUTDIR/final.outputs.json"
else
 rm -rf "$FINAL"; mkdir -p "$FINAL"; cp -a "$PRE/." "$FINAL/"; cp "$OUTDIR/pre_pvacbind.outputs.json" "$OUTDIR/final.outputs.json"
fi

val=(-m neoag.splice.cli validate --layer-dir "$FINAL" --report "$OUTDIR/validation.json")
[[ "$STRICT" == 1 ]] && val+=(--strict)
"$PY" "${val[@]}"
"$PY" - "$OUTDIR" "$SAMPLE" "$PHASE1" "$PRE" "$FINAL" <<'PY'
from __future__ import annotations
import hashlib,json,sys
from datetime import datetime,timezone
from pathlib import Path
out,sample,phase1,pre,final=map(Path,sys.argv[1:])
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for c in iter(lambda:f.read(1024*1024),b''): h.update(c)
 return h.hexdigest()
def item(p): return {'path':str(p.resolve()),'sha256':sha(p)} if p.is_file() else None
payload={'schema_version':'neoag-v051-driver-v1','created_at':datetime.now(timezone.utc).isoformat(),'sample_id':sample.name,
 'phase1_manifest':item(phase1/'provenance_manifest.json'),'pre_pvacbind_manifest':item(pre/'provenance_manifest.json'),
 'final_manifest':item(final/'provenance_manifest.json'),'validation':item(out/'validation.json')}
(out/'v051_driver_manifest.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')
PY
printf '%s\n' "$FINAL"
