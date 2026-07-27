# Deployment tiers

- **review**: read event/peptide results, compare rankings, and generate reports.
- **core**: run the Python core, fixtures, and evidence-consensus ranking on precomputed inputs.
- **prediction**: additionally run VEP, pVACtools, NetMHCpan, MHCflurry, and prediction assets.
- **full**: additionally process BAM/FASTQ and production HLA, fusion, splice, CNV/purity, and licensed-tool workflows.

A higher-tier missing tool does not make a lower tier unusable. Status is evaluated against the tier requested by the user.
