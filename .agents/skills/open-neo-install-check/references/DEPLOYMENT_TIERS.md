# Deployment tiers

- **review**: read event/peptide results, compare rankings, and generate reports.
- **core**: run the Python core, fixtures, and evidence-consensus ranking on precomputed inputs.
- **prediction**: additionally requires VEP, pVACseq, NetMHCpan, MHCflurry,
  GRCh38 FASTA plus indexes/dictionary, GENCODE GTF, VEP cache and normal
  proteome. At least one immunogenicity-like implementation must be available.
- **full**: additionally requires Java/Nextflow, BWA, STAR and GATK; RNA
  alignment/quantification and safety references; and at least one working
  implementation in each HLA typing, HLA LOH, purity/CNV, fusion, splice and
  CCF capability group. Tool-specific optional alternatives are reported but
  do not all have to be installed simultaneously.

`tier_requirements.tsv` is the authoritative readiness matrix. `READY` is
allowed only when every required item and every any-of capability group is
satisfied. A missing higher-tier dependency remains `PARTIAL`, while a broken
Python/Open-Neo core is `BLOCKED`.

A higher-tier missing tool does not make a lower tier unusable. Status is evaluated against the tier requested by the user.
