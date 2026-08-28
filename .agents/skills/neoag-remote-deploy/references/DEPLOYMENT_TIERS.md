# Deployment tiers

## Tier 0: result review

Use when the target machine only needs to inspect existing outputs such as
`ranked_peptides`, `evidence_report`, APPM, CCF, or patient/technical reports.

Acceptance:

```bash
python -m pip install -e .
neoag-skill validate --root . --outdir work/skill_validate
neoag run-demo --outdir work/demo_v043 --sample-id DEMO001
```

## Tier 1: core scoring

Use when standardized intermediate tables already exist and the target machine
only needs recombination scoring and reports.

Acceptance:

```bash
neoag run-demo --outdir work/demo_v043 --sample-id DEMO001
pytest -q tests/test_skills_taxonomy_abcd.py
```

## Tier 2: containerized prediction subset

Use when the target machine must run selected external tools such as VEP,
NetMHCpan, MHCflurry, HLA tools, fusion tools, or CNV tools.

Acceptance:

```bash
neoag-doctor --project-root . \
  --tools-manifest configs/local/tools_manifest.yaml \
  --reference-manifest configs/local/reference_manifest.yaml \
  --outdir work/doctor --dry-run --mini-smoke
```

## Tier 3: full production/HPC

Use only when the target machine must run from BAM/FASTQ/VCF to final reports.
Full execution requires explicit human approval.

Pre-execution acceptance:

```bash
neoag pipeline-full \
  --sample-manifest configs/local/sample_manifest.yaml \
  --tools-manifest configs/local/tools_manifest.yaml \
  --reference-manifest configs/local/reference_manifest.yaml \
  --outdir results/SAMPLE001 \
  --profile slurm
```
