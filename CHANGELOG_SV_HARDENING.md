# SV hardening patch

- Added canonical, build-aware, order-independent SV adjacency identity.
- Removed silent tumor/normal sample-order inference for multi-sample VCFs.
- Added default exact-PASS VCF filtering with explicit research override.
- Replaced gene-pair RNA lookup with exact oriented-breakpoint lookup.
- Added RNA junction QC fields and thresholds.
- Added confirmed expressed-product input and peptide-generation gate.
- Demoted DNA/CDS heuristic reconstructions to unresolved hypotheses.
- Removed junction-read-to-TPM conversion.
- Added stable SHA-256 peptide identifiers.
- Made WES capture BED mandatory and capture-unassessed events uninterpretable.
- Wired capture BED, sample names, genome build and expressed products through
  the Nextflow SV build process.
- Removed fixture defaults from the hardened WES workflow entry point.
- Added a DNA_SV-specific C1-C4 source-chain requirement builder.
- Added focused hardening tests and exact-evidence fixtures.

