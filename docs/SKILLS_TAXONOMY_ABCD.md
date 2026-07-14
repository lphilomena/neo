# Project B Skills Taxonomy A/B/C/D

This document describes the upgraded Skills system for Project B. It follows the principle that **Agent decides what to do**, **Skill defines a stable SOP**, and **Tool/CLI/Nextflow performs computation**.

## A. Entry adapter Skills

These skills normalize heterogeneous user inputs into Project B standard `raw_events.tsv`, `raw_peptides.tsv`, or evidence tables.

- `neoag-vcf`: somatic VCF / VEP VCF → raw events.
- `neoag-fusion`: fusion caller output → fusion events and junction peptide inputs.
- `neoag-splice`: splice junction table → splice/exon-junction events.
- `neoag-sv-wgs`: WGS SV VCF → SV reconstruction tasks.
- `neoag-sv-wes`: WES/capture-limited SV VCF → conservative SV tasks with confidence caps.
- `neoag-peptide-csv`: existing peptide-HLA table → normalized raw peptides and presentation evidence.

## B. Public evidence analysis Skills

These skills are shared by all entry routes.

- `neoag-hla-typing-loh`
- `neoag-presentation`
- `neoag-expression`
- `neoag-rna-evidence`
- `neoag-ccf`
- `neoag-appm-escape`
- `neoag-safety`
- `neoag-ranking`

## C. Review/report/experiment-design Skills

- `neoag-ranking-compare`
- `neoag-experiment-design`
- `neoag-patient-report`
- `neoag-technical-report`
- `neoag-concept-explainer`

## D. Governance/execution-control Skills

- `neoag-input-qc`
- `neoag-doctor`
- `neoag-tool-reference-qc`
- `neoag-run-demo-and-smoke`
- `neoag-pipeline-full`
- `neoag-release-qc`
- `neoag-gateway-submit`
- `neoag-hpc-runner`
- `neoag-remote-deploy`

## Safety boundaries

- Skills do not make clinical decisions.
- Missing evidence is marked as `missing` or `unassessed`; it is not interpreted as negative evidence.
- HPC submission, installation, deletion, overwrite, and large reference downloads require human approval.
- Skills do not contain patient BAM/FASTQ/VCF, VEP cache, GRCh38, NetMHCpan license, LOHHLA reference, or large conda environments.

## CLI

For deployment on a new machine, do not ask an agent to infer the workflow from
README examples. Use the skill-first migration guide and bootstrap:

```bash
bash scripts/bootstrap_agent_deploy.sh
```

See `docs/SKILL_FIRST_MIGRATION.md` for the recommended agent prompt, manifest
rules, risk policy and success criteria.

List skills:

```bash
neoag-skill list
```

Describe a skill:

```bash
neoag-skill describe neoag-vcf
```

Validate repo-scoped skill directories:

```bash
neoag-skill validate --root . --outdir work/skill_validate
```

Run a dry-run:

```bash
neoag-skill run neoag-vcf --outdir work/neoag-vcf --dry-run --arg vcf=sample.vcf.gz
```

Run a normalized peptide table:

```bash
neoag-skill run neoag-peptide-csv --outdir work/peptides --arg peptide_csv=peptides.tsv
```
