# NeoAg v0.5.1 change log

## Added

- moPepGen exact peptide/full-ORF provenance adapter.
- splice2neo mutation-junction-transcript/peptide adapter.
- EasyQuant exact query-ID targeted re-quantification adapter.
- pVACsplice exact variant + strand-aware junction branch.
- k4neo exact cts_id normal-background adapter with explicit license gate.
- Canonical variant IDs and exact event index.
- Three independent evidence-chain table.
- External query registry and generated EasyQuant/k4neo inputs/maps.
- Multi-pass production driver and auditable external-tool wrappers.

## Changed

- Schema version 0.5.0 → 0.5.1.
- Consensus now distinguishes exact peptide consensus from exact full-ORF consensus.
- DNA-causal support can raise event evidence to E3.
- k4neo-only negative evidence cannot masquerade as coverage-aware normal negativity.
- Manifest records the matching and evidence-chain policies.

## Preserved

- v0.4.4 exact junction normalization and no-evidence-leakage rules.
- v0.5.0 junction→event→transcript hypothesis→ORF→peptide origin→HLA foreign-key chain.
- Compatibility projections for `raw_events.tsv`, `raw_peptides.tsv`, and `rna_junction_evidence.tsv`.
