# SpliceMutr normal-RNA first-layer P0 fixes

## Scope

This release hardens the single-tumor SpliceMutr/SNAF path when no compatible
normal RNA cohort is available.  It does not manufacture a cohort comparison,
delta PSI, q-value, outlier statistic, local normal coverage, or tumor-specific
claim.  The schema version is `0.5.3-splicemutr-normal-p0`.

## Evidence contract

The pipeline keeps three statements separate:

1. `ALTERED_JUNCTION_SPANNING_SEQUENCE`: structural sequence novelty was
   reconstructed across the exact canonical junction.
2. `NOT_LISTED_IN_NORMAL_CATALOG`: the exact junction was not listed in the
   selected presence-only catalog.  This is not an adequate-coverage negative.
3. `UNASSESSED_NO_COMPATIBLE_NORMAL_RNA_COHORT`: tumor specificity was not
   statistically assessed against a compatible normal RNA cohort.

When statement 3 applies, normal safety is at most `N1` and the final evidence
tier is at most `R3`, regardless of HLA prediction strength.

## Implemented fixes

### Normal catalog metrics

`src/neoag/splice/normal_background.py` now parses and preserves:

- `normal_reads`: maximum observed normal-junction count;
- `normal_total_reads`: cohort-wide junction-read sum;
- `normal_samples`: positive sample count;
- `total_samples`: cohort denominator when supplied;
- `sample_prevalence`;
- normal tissue names and tissue count;
- dataset and release provenance.

A catalog row with any positive read/sample evidence is `DETECTED`, even when a
generic caller adapter would otherwise report zero.  A catalog negative remains
unavailable unless an explicit coverage table reports both `NOT_DETECTED` and
adequate locus coverage.

### Catalog non-membership

Historical `ABSENT_GTEX_V11` and `SEEN_GTEX_V11` labels are replaced by:

```text
NOT_LISTED_IN_NORMAL_CATALOG
SEEN_IN_NORMAL_CATALOG
```

Catalog non-membership is still allowed as a screening filter, but contributes
no positive composite-score term.

### Structural novelty versus tumor specificity

SpliceMutr altered paths now emit:

```text
structural_novelty_status=ALTERED_JUNCTION_SPANNING_SEQUENCE
tumor_specificity_status=UNASSESSED_NO_COMPATIBLE_NORMAL_RNA_COHORT
cohort_analysis_status=UNASSESSED_NO_COMPATIBLE_NORMAL_RNA_COHORT
```

They no longer receive `mutant_specificity_status=MT_SPECIFIC` solely because
the path was labelled modified.  The legacy mutant-specificity gate becomes
`UNASSESSED/REVIEW_REQUIRED` and is capped at `R3`.

### Reference provenance

The normal-junction build scripts identify the resource as:

```text
source=recount3_GTEx_v8
dataset=GTEx_v8
reference_release=recount3_GTEx_v8_GRCh38
```

This is distinct from any separately maintained GTEx v11 expression resource.
For a merged pan-tissue table, `total_samples` counts sample-tissue records;
donors can contribute more than one tissue.  This denominator definition is
written to the sidecar metadata.

## Rebuilding the normal catalog

Per tissue:

```bash
python scripts/build_recount3_normal_junctions.py \
  --rr references/liver.junctions.tsv.gz \
  --mm references/liver.counts.mtx.gz \
  --tissue Liver \
  --output references/normal_junctions.recount3_gtex_v8_liver_grch38.tsv.gz
```

Merge tissues:

```bash
python scripts/merge_normal_junction_tissues.py \
  --inputs references/normal_junctions.recount3_gtex_v8_*_grch38.tsv.gz \
  --output references/normal_junctions.recount3_gtex_v8_grch38.tsv.gz
```

The same GRCh38 coordinate contract and exact strand-aware canonical junction
identity must be used by tumor and normal resources.

The main normalization CLI always emits the missing-cohort state.  A future
LeafCutter or LeafCutterMD adapter must import the actual delta PSI, q-value,
outlier, prevalence, and coverage records; changing a label is not an allowed
substitute for those data.

## Migrating historical output tables

The migration is downgrade-only and retains the input file unchanged:

```bash
python scripts/migrate_splicemutr_normal_p0.py \
  --input old_raw_peptides.tsv \
  --output migrated_raw_peptides.tsv \
  --kind peptides
```

Supported kinds are `candidates`, `peptides`, and `events`.  Only splice rows
are downgraded; SNV/indel `MT_SPECIFIC` calls are not modified.

## Running SpliceMutr origin reconstruction

```bash
PYTHONPATH=src python scripts/rebuild_splice_origins_from_splicemutr.py \
  --sample-id CASE001 \
  --genome-build GRCh38 \
  --candidates results/snaf_candidates.tsv \
  --formal-events results/splice_events.tsv \
  --splicemutr-glob 'results/splicemutr/formed_transcripts/**/*_data_splicemutr_cp_corrected.txt' \
  --outdir results/splicemutr_normal_p0
```

## Production interpretation

- Suitable for discovery and manual research review.
- Catalog detection is negative safety evidence and can hard-fail a candidate.
- Catalog non-membership is not tumor-specific evidence.
- Without a compatible normal RNA cohort, candidates cannot exceed `N1/R3`.
- Peptide synthesis or clinical release still requires a compatible normal
  comparison and/or targeted experimental validation.

## Tests

The standard-library regression suite is runnable without pytest:

```bash
PYTHONPATH=src python tests/test_splicemutr_normal_p0_unittest.py
```

It covers normal-read parsing, sample denominators, no catalog-absence scoring,
N1/R3 capping, SpliceMutr specificity downgrading, historical migration, and
recount3 release provenance.
