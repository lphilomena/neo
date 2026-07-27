# Deterministic routing

Input detection precedence is: sample manifest, explicit CLI fill-ins, directory scan, then extension/header inference.

Execution routes are selected in this order:

1. Explicit production manifest.
2. Existing result/evidence inputs.
3. BAM/FASTQ production inputs; execution requires a reviewed production manifest.
4. Standard raw events + raw peptides.
5. Somatic VCF, fusion table, splice table, peptide table, and/or SV inputs.

Raw SV is first normalized by the WGS/WES SV adapter. In a multi-entry case its standard raw tables are merged with the VCF/fusion/splice catalogue before the shared evidence and ranking layers.

Multiple non-SV entry branches may be merged into a single standard event/peptide catalogue before the shared evidence layers and rankings are run.
