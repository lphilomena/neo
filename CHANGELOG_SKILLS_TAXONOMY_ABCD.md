# Changelog: Skills Taxonomy A/B/C/D

## Added

- Added `neoag.skill_taxonomy` package.
- Added A/B/C/D skill taxonomy registry with 27 skills.
- Added `neoag-skill` CLI.
- Added repo-scoped `.agents/skills/*/SKILL.md` definitions for all skills.
- Added entry adapter implementations for VCF, fusion, splice, WGS SV, WES SV and peptide CSV.
- Added public evidence skill implementations for HLA typing/LOH, presentation, expression, RNA evidence, CCF, APPM/escape, safety and ranking.
- Added review/report/design skill implementations for ranking compare, experiment design, patient report, technical report and concept explainer.
- Added governance wrappers for input QC, doctor, tool/reference QC, run-demo smoke, pipeline-full, release QC, gateway submit and HPC runner.

## Boundaries

- Skills are SOP wrappers and do not make clinical decisions.
- Heavy operations remain dry-run or require explicit approval.
