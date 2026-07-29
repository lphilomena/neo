# v0.4.4 Splice Junction Normalization and Provenance Contract

## 1. Purpose

This contract prevents junction evidence from crossing between biologically different splice events and makes every merge auditable. It applies to RegTools, STAR junctions, pVACsplice, SNAF, SpliceMutr, and generic normalized splice tables.

## 2. Canonical coordinate model

The authoritative identity is:

```text
SJ|BUILD|CHR|INTRON_START_1BASED|INTRON_END_1BASED|STRAND
```

Example:

```text
SJ|GRCh38|chr1|154589929|154590184|-
```

For the minus strand, donor is the higher genomic coordinate and acceptor is the lower coordinate. Caller coordinates remain in source fields and must never replace the canonical fields during matching.

## 3. Required output semantics

### `canonical_junction_id`

Resolved exact identity. Empty means the source row is unresolved or conflicted.

### `source_junction_id`

Caller-local identifier such as `JUNC00000017`. It is not globally unique. An alias that maps to more than one canonical junction is `AMBIGUOUS` and transfers zero support.

### `provided_rna_junction_reads`

The count stated by the source row. It is provenance only until exact resolution.

### `rna_junction_reads`

Verified support from the exact primary RNA-junction entity. Unresolved, ambiguous, source-only, or unstranded-primary records receive zero.

### `junction_support_status`

Common states:

- `SUPPORTED_EXACT_JUNCTION`
- `MATCHED_ZERO_READS`
- `RESOLVED_WITHOUT_PRIMARY_SUPPORT`
- `UNSTRANDED_PRIMARY_UNVERIFIED`
- `PROVIDED_UNVERIFIED`
- `AMBIGUOUS`
- `UNRESOLVED`

## 4. Matching order

1. canonical junction ID;
2. exact build/chromosome/intron interval/strand;
3. exact unique source alias;
4. explicit unique variant-to-junction relation;
5. unique unstranded coordinate with caution, without verified read transfer when the primary evidence itself is unstranded.

The implementation must not use gene, nearby position, maximum count, or approximate interval matching.

## 5. Merge contract

Event merge key:

```text
event_id
```

For canonical splice events, `event_id` equals `canonical_junction_id`.

Peptide merge key:

```text
event_id + peptide + hla_allele
```

Caller-specific peptide IDs are provenance and do not create duplicate biological entities.

Every input row generates a provenance record containing stage, tool, file, row number, source record ID, and SHA-256. Conflicts are materialized instead of discarded.

## 6. Running the normalizer

```bash
PYTHONPATH=src python scripts/normalize_rna_fusion_splice.py \
  --sample-id SAMPLE001 \
  --profile rna_fusion_splice_v1 \
  --genome-build GRCh38 \
  --junctions regtools_junctions.tsv \
  --snaf snaf_candidates.tsv \
  --splicemutr splicemutr_candidates.tsv \
  --normal-junctions normal_junctions.tsv \
  --outdir results/SAMPLE001/splice/intermediates
```

`--snaf`, `--splicemutr`, and `--normal-junctions` are optional. The primary `--junctions` table is required.

## 7. Safety invariants

- An unresolved source count never increases a candidate's RNA support score.
- Two rows sharing only a gene never merge.
- Two rows sharing an ambiguous alias never exchange reads.
- Normal non-observation without locus coverage is not a strong negative.
- Cross-domain confirmation requires the same canonical junction.
- All source records remain recoverable after entity deduplication.
