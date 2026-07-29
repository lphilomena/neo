# NeoAg v0.4.4 — Exact Junction, Evidence Non-Leakage, and Provenance Repair

Release date: 2026-07-28

## Problem corrected

v0.4.3 could transfer RNA junction support through a gene-level or nearby-locus fallback. In a gene containing multiple junctions, a highly expressed canonical junction could therefore supply reads to a weak or unrelated neo-junction. The production runner also retained only the first duplicate event/peptide row, discarding later caller provenance. Finally, splice consensus used loose string keys, so different events in the same gene could be incorrectly treated as cross-domain confirmation.

## Canonical identity

v0.4.4 defines one internal junction identity:

```text
SJ|genome_build|chromosome|intron_start_1based|intron_end_1based|strand
```

The canonical interval is the intron, represented as 1-based closed coordinates. Genome build and strand are mandatory parts of the identity. Source coordinates, source coordinate system, conversion method, file, row, tool, and row hash are retained separately.

Supported conversions include:

- RegTools annotated tables: zero-based start/end to 1-based closed intron (`start + 1`, `end`).
- RegTools BED12 extraction output: intron derived from block sizes and block starts, not from the outer BED interval.
- STAR `SJ.out.tab`: direct 1-based intron interval.
- SNAF-style outer splice-boundary UIDs: explicit boundary-to-intron conversion.
- Already canonical `SJ|...` identifiers.

## Evidence non-leakage policy

Verified `rna_junction_reads` may be transferred only through:

1. an exact canonical junction ID;
2. exact build/chromosome/intron-start/intron-end/strand;
3. an exact unique source-junction alias already registered to one canonical junction;
4. an explicit unique variant-to-junction relation recorded by the source;
5. a unique unstranded interval only with an explicit caution state; unstranded primary counts are not transferable as verified support.

Removed fallbacks:

- same-gene maximum junction reads;
- nearest genomic locus;
- approximate coordinate windows;
- gene-only caller agreement;
- ambiguous source aliases.

Caller-provided but unresolved values are retained in `provided_rna_junction_reads` and contribute zero to `rna_junction_reads`.

## Provenance-preserving merge

The production runner now groups biologically equivalent entities while writing one provenance row per input record. It emits:

- `merged/event_provenance.tsv`;
- `merged/peptide_provenance.tsv`;
- `merged/evidence_conflicts.tsv`.

For the splice normalizer it additionally emits:

- `splice_junctions.tsv`;
- `splice_tool_evidence.long.tsv`;
- `splice_peptide_provenance.tsv`;
- `splice_event_merge_provenance.tsv`;
- `splice_peptide_merge_provenance.tsv`;
- `splice_merge_conflicts.tsv`;
- `junction_aliases.tsv`;
- `splice_consensus.tsv`;
- `splice_consensus_provenance.tsv`;
- `splice_consensus_conflicts.tsv`;
- `evidence_conflicts.tsv`;
- `splice_qc.tsv`;
- `provenance_manifest.json`.

Scalar output tables remain backward compatible. `source_tools`, `source_records`, and `provenance_record_count` preserve the full source set. Conflicting scalar observations are retained in the conflict table rather than silently overwritten.

## Exact consensus

`CROSS_DOMAIN_CONFIRMED_EXACT_JUNCTION` now requires RNA-junction evidence and neoantigen evidence to resolve to the same canonical junction. Unresolved or ambiguous rows cannot produce cross-domain confirmation. Same-gene but different-junction records remain separate.

## Normal background semantics

A normal-panel absence is reported as `NOT_DETECTED_COVERAGE_UNASSESSED` unless locus-level coverage is available. The pipeline no longer converts panel-level non-observation into a strong negative assertion.

## Compatibility

- `raw_events.tsv` and `raw_peptides.tsv` remain the standard intermediate tables.
- New schema fields are additive.
- `rna_junction_reads` now consistently means verified exact-junction support.
- `provided_rna_junction_reads` means upstream caller-provided support not yet independently resolved.
- Existing Mode C and pVACsplice interfaces remain available.
