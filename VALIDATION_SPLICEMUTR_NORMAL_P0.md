# Validation: SpliceMutr normal-RNA first-layer P0

## Release

- Version: `0.5.3-splicemutr-normal-p0`
- Base: `0.5.2-p0`
- Validation date: 2026-08-27
- Scope: normal-catalog parsing and evidence semantics when no compatible
  normal RNA cohort exists

## Results

| Validation | Result |
|---|---:|
| Dedicated standard-library regression tests | 8 passed |
| Existing Splice P0 direct tests | 8 passed |
| Existing v0.5.0/v0.5.1 direct splice tests | 33 passed |
| Total directly executed tests | 49 passed |
| End-to-end normal-catalog pipeline smoke | passed |
| Python compileall | passed |
| v0.5.3 changed-table schema parity | passed |
| JSON syntax checks | passed |
| Shell syntax check | passed |

The environment does not contain pytest, so the complete historical pytest
suite was not executed and is not claimed as passing.  The 33 legacy tests were
run directly with standard temporary-directory fixtures; parameterized pytest
cases outside this execution route were not counted.

## Regression assertions

1. `normal_reads=19`, `normal_samples=7`, and `normal_total_reads=83` produce a
   real normal detection rather than a zero-read record.
2. Normal-catalog positive evidence becomes `DETECTED_BROAD_NORMAL`.
3. Catalog non-membership does not increase the SNAF composite score.
4. Catalog non-membership without local coverage remains `N1` and caps final
   evidence at `R3`.
5. The main normalization path emits neutral numeric tumor specificity (`0.5`)
   and `UNASSESSED_NO_COMPATIBLE_NORMAL_RNA_COHORT`.
6. A structurally altered SpliceMutr path is no longer labelled `MT_SPECIFIC`.
7. Historical splice tables are migrated by downgrade only; SNV/indel rows are
   outside the migration rule.
8. Per-tissue recount3 tables retain sample denominators and the merged table
   documents sample-tissue-record denominator semantics.
9. Dataset provenance remains `recount3_GTEx_v8`, `GTEx_v8`, and
   `recount3_GTEx_v8_GRCh38` through formal normal-background output.
10. Existing exact-junction, read-QC, no-cross-junction borrowing, pVACsplice
    residue-boundary, normal-detection, fallback-event, and partial-ORF P0 tests
    remain passing.

## End-to-end smoke

The production splice provenance builder was run with:

- one exact GRCh38 tumor junction;
- complete RNA read-QC metrics;
- one matching recount3/GTEx catalog row containing `normal_reads=9`,
  `normal_samples=3`, `normal_total_reads=17`, and `total_samples=251`.

The emitted `splice_normal_background.tsv` retained the counts and release and
reported `DETECTED_BROAD_NORMAL`.

## Remaining limitations

- This release does not perform LeafCutter/LeafCutterMD cohort analysis.
- It does not infer local normal coverage from a presence-only catalog.
- It does not calculate delta PSI, q-values, outlier scores, or protocol/batch
  corrections.
- Therefore it supports discovery/manual review, not automated synthesis or
  clinical release.
