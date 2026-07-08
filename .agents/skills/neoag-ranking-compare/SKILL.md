---
name: neoag-ranking-compare
description: Compare NetMHCpan42 ranking with recommendation ranking and identify concordant/conflicting candidates.
---

# neoag-ranking-compare

## Use when

- The user request matches this skill description.
- The required inputs are available or can be discovered from a result directory.

## Do not use when

- The user asks for a clinical diagnosis or guaranteed treatment effect.
- Required inputs are absent and cannot be requested from the user.
- A high-impact operation would overwrite data, submit HPC jobs, or install tools without explicit confirmation.

## Required inputs

See the command wrapper and `references/INPUTS.md`. Missing inputs must be reported as missing evidence, not interpreted as negative biological evidence.

## Primary command

```bash
python -m neoag_v03.agent_skills.ranking_compare --netmhcpan42 ranked_peptides.netmhcpan42.tsv --recommendation ranked_peptides.recommendation.tsv --outdir <OUTDIR>
```

## Outputs

Each skill writes a Markdown report plus TSV/JSON sidecars under the requested output directory.

## Safety and interpretation boundaries

- Computational triage only. Candidate neoantigens are not confirmed neoantigens.
- Missing RNA/HLA LOH/APPM/CCF evidence must be labelled as missing or partial, not negative.
- Patient-facing outputs must avoid clinical efficacy promises.
