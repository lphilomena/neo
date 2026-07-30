# NeoAg v0.5.0 — Formal Splice Provenance Layer

Release date: 2026-07-30

## Added

- Formal five-level splice model: junction → event → transcript hypothesis → ORF → peptide origin.
- Stable deterministic IDs and exact foreign-key link tables.
- SplAdder GFF3 and TXT adapters, including headerless build-mode event files.
- IRFinder-S retained-intron adapter with explicit coordinate-system handling.
- ImmunoPepper metadata and junction-kmer adapters with support for official semicolon coordinate output.
- pVACbind adapter using exact generated FASTA indices and exact epitope-to-ORF validation.
- Coverage-aware normal-background state model.
- Independent evidence-group consensus with E/O/N/P grades and R1–R4 caps.
- Formal compatibility projection into existing `raw_events.tsv`, `raw_peptides.tsv`, and `rna_junction_evidence.tsv`.
- `neoag-splice-layer` CLI and two-pass production shell driver.
- Machine-readable table contract at `resources/splice_provenance_v050.schema.json`.
- Seven-state, coverage-aware normal-background assessment contract.
- Exact high-order evidence adapter for long-read, DNA-causal, protein-validation, and ligandome evidence.
- RNA `source_assay_id` tracking so STAR and RegTools from the same BAM are not counted as independent assays.

## Safety changes

- No gene or nearest-locus fallback is used anywhere in the formal layer.
- Caller-provided read counts do not transfer without exact canonical junction identity.
- Unstranded junctions are retained only as unresolved review evidence; strict mode rejects them and they cannot contribute verified exact-junction support.
- pVACbind presentation requires an intact event→transcript→ORF chain and matching FASTA-map, stored-ORF, and computed sequence SHA-256 values.
- pVACbind positions inconsistent with the mapped ORF are corrected only through a unique exact sequence match; otherwise the row is rejected to the conflict table.
- Normal non-detection without explicit adequate coverage remains incomplete.
- Unknown peptide novelty receives a cap rather than being asserted as true or hard-failed as false.
- ImmunoPepper `isJunctionList` is not used as a crossing flag.
- Cryptic-exon, exitron, novel-junction, and complex-splice event labels are normalized explicitly.
- Strict production builds require a locked version for every external tool whose input is consumed.
- IRFinder-S input now requires an explicit coordinate-system declaration at the production entrypoint.
- SplAdder reference/alternative path classification is populated only from explicit source roles.
- High-order evidence with an unknown entity ID, unsupported group/status, or missing strict-mode version is rejected to `splice_conflicts.tsv`.
- Every consensus cap and hard-fail reason is materialized as an auditable conflict row.

## Compatibility

- v0.4.4 exact-junction outputs remain available.
- Existing ranked/report layers continue to consume the compatibility tables.
- New formal fields are additive.
- `rna_fusion_splice_v1` remains unchanged; the formal layer is invoked through `neoag-splice-layer` or `scripts/run_splice_provenance_v050.sh`.
