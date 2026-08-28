# Deterministic routing

Input detection precedence is: sample manifest, explicit CLI fill-ins, directory scan, then extension/header inference.

Execution routes are selected in this order:

1. Explicit production manifest.
2. Completed `case_root` plus `somatic_vcf`: use `scripts/run_production_case.sh` to generate and execute a production-results manifest from existing upstream outputs.
3. Existing result/evidence inputs.
4. BAM/FASTQ production inputs; execution requires a reviewed production manifest.
5. Standard raw events + raw peptides.
6. Somatic VCF, fusion table, splice table, peptide table, and/or SV inputs.

Raw SV is first normalized by the WGS/WES SV adapter. In a multi-entry case its standard raw tables are merged with the VCF/fusion/splice catalogue before the shared evidence and ranking layers.

Multiple non-SV entry branches may be merged into a single standard event/peptide catalogue before the shared evidence layers and rankings are run.

For completed case roots, fusion inputs are cumulative rather than mutually exclusive: EasyFuse pass/unfiltered evidence, STAR-Fusion, Arriba, FusionCatcher, JAFFAL and configured caller-root directories are normalized into one fusion union while preserving caller provenance and single-caller labels.

Normal background references should be routed from explicit inputs first, then from the fixed asset root: normal junctions with a reusable sqlite index, normal expression, HLA ligandome, and reference proteome. Missing background is reported as safety partial/UNASSESSED and must not be interpreted as negative evidence.

RNA allele evidence has exactly one active mode per run: RNA FASTQ-derived STAR pileup, supplied RNA BAM pileup, or existing RNA VAF table.
