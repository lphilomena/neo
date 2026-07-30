# NeoAg v0.5.0 Formal Splice Provenance Layer — Validation Report

- Release date: 2026-07-30
- Baseline: NeoAg v0.4.4 exact-junction/provenance repair
- Scope: source-level, schema, referential-integrity, compatibility, packaging, and synthetic two-pass presentation validation

## 1. Implemented production contract

The release materializes the following authoritative chain:

```text
canonical junction
  → biological splice event
  → transcript hypothesis
  → ORF / translated segment
  → peptide origin
  → peptide-HLA presentation
```

It preserves the v0.4.4 exact-junction non-leakage rules and adds stable foreign-key relations, normal-background states, independent evidence-group consensus, pVACbind exact FASTA-index mapping, conflict materialization, and compatibility projections.

## 2. Automated validation results

| Check | Result |
| --- | --- |
| Python `compileall` over `src/` and `scripts/` | PASS |
| Shell syntax for v0.5.0 production drivers | PASS |
| CLI help for source installation | PASS |
| v0.4.4 + v0.5.0 splice regression suite | 33 passed, 0 failed |
| Full default repository test suite | 437 passed, 105 skipped, 0 failed |
| Unstranded-junction downgrade and strict rejection | PASS |
| pVACbind event/transcript/ORF/hash-chain rejection tests | PASS |
| Cryptic-exon/exitron event normalization | PASS |
| Strict external-tool version locks | PASS |
| Seven-state normal-background fixture | PASS |
| SplAdder explicit reference/alternative path roles | PASS |
| IRFinder-S missing coordinate declaration rejection | PASS |
| High-order evidence exact-entity linking and E3/O3 fixture | PASS |
| Same-assay RNA source collapse / independent-assay separation | PASS |
| Consensus hard-fail and priority-cap conflict materialization | PASS |
| Install Skill cross-machine inaccessible-path fallback | PASS |
| Doctor inaccessible tool/reference reporting without crash | PASS |
| Full-tier install-check plan and verify execution | PASS (`BLOCKED` readiness result as expected) |
| Project wrapper, conda environment, and portable asset auto-discovery | PASS |
| Prediction/full tier default installer profile | PASS (`standard`; explicit override preserved) |
| Full-tier required tool/reference readiness after portable installation | READY (0 required gaps) |
| BAM-matcher GRCh38 loci provenance and FASTA REF validation | PASS (1,303 exact rsID-mapped loci) |
| Headered/headerless SplAdder fixture parsing | PASS |
| RegTools exact-junction non-leakage fixture | PASS |
| IRFinder-S retained/spliced hypothesis fixture | PASS |
| ImmunoPepper semicolon-coordinate and partial-ORF fixture | PASS |
| pVACbind unresolved Index rejection | PASS |
| pVACbind epitope/ORF mismatch rejection | PASS |
| Strict no-pVAC two-pass driver smoke | PASS |
| Synthetic pVACbind two-pass driver smoke | PASS |
| Final manifest/hash and referential-integrity validation | PASS |
| JSON example config and machine-readable schema parsing | PASS |
| Wheel content, isolated import, and console CLI smoke | PASS |

The 105 skipped tests remain conditional tests in the existing project configuration. They are not counted as successful executions. Integration validation also confirmed that the dependency-light manifest parser preserves nested scalar lists when PyYAML is unavailable.

## 3. Safety properties explicitly tested

1. A high-read junction in the same gene cannot transfer reads to another junction.
2. Event identity is deterministic and independent of junction input order.
3. Event-junction, transcript-event, ORF-transcript, and peptide-origin foreign keys remain valid.
4. An ImmunoPepper local translated segment is not relabelled as a confirmed full-length transcript.
5. pVACbind results contribute presentation evidence only through a unique generated FASTA Index.
6. A pVACbind epitope absent from the mapped ORF is rejected and recorded as a conflict.
7. Unknown normal coverage remains incomplete; it is not promoted to a strong normal-negative result.
8. Single-generator ORFs are capped below the highest evidence tier.
9. Compatibility `raw_events.tsv`, `raw_peptides.tsv`, and `rna_junction_evidence.tsv` retain formal provenance identifiers.
10. SplAdder path roles remain unresolved unless the source explicitly identifies reference or alternative paths.
11. STAR and RegTools records sharing one `source_assay_id` remain one RNA source; different assay IDs remain distinct.
12. High-order evidence with an unknown formal entity ID cannot upgrade E3/O3.
13. Every consensus hard-fail and cap reason is present as an auditable conflict record.
14. Tool/reference paths copied from another machine but inaccessible to the current account are reported as missing and do not crash auto-configuration, Doctor, or tier assessment.

## 4. Synthetic presentation smoke

The production driver was run with a synthetic pVACbind-compatible executable. The smoke test verified:

- pass 1 generated an ORF FASTA and exact FASTA map;
- the synthetic predictor returned one exact epitope for the generated Index;
- pass 2 mapped the result to one `ORF|…`, one `STH|…`, one `SEV|…`, and one `POR|…`;
- the accepted record received a stable `PRE|…` presentation identifier;
- final validation returned `PASS`.

This validates orchestration and provenance. It is not a biological benchmark of pVACbind or an HLA predictor.

## 5. Validation boundaries

The following were not claimed or inferred from the automated tests:

- clinical validity or treatment utility;
- experimental translation of a predicted ORF;
- endogenous HLA presentation;
- T-cell recognition;
- calibration on a clinical truth set;
- successful execution of licensed predictors or site-local reference assets;
- complete compatibility with every historical output variant of each external tool.

Production deployment must lock external tool versions, references, genome build, coordinate contracts, HLA inputs, and normal-background resources, then run site-specific fixtures before patient-data analysis.
