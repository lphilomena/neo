# Project B Agent + Skills P0/P1 Implementation

This release adds a lightweight, repo-scoped Skills Pack and a Coordinator Agent for Project B.

## P0 core skills

1. `neoag-input-qc` — input/result directory inspection and workflow recommendation.
2. `neoag-tool-and-reference-qc` — external tool and reference readiness checks.
3. `neoag-run-demo-and-smoke` — pytest, run-demo, and optional Nextflow smoke planning/execution.
4. `neoag-evidence-scoring` — evidence scoring wrapper / command plan.
5. `neoag-hla-loh-appm-review` — HLA LOH, APPM, and immune escape interpretation.
6. `neoag-ccf-clonality-review` — purity, CCF, clonality, and CCF modifier interpretation.
7. `neoag-ranking-compare` — NetMHCpan42 vs recommendation ranking comparison.
8. `neoag-patient-report` — patient-facing report draft generation.

Each skill is stored under `.agents/skills/<skill-name>/` and contains a `SKILL.md`, wrapper scripts, and reference notes. Skills are SOPs and should not store patient data, BAM/FASTQ/VCF files, large references, VEP cache, or licensed NetMHCpan bundles.

## P1 Coordinator Agent

The Coordinator Agent is available as:

```bash
neoag-agent --message "比较两个排序文件" --result-dir results/sample --outdir work/agent_plan
```

Default mode is dry-run planning. Use `--execute` to run skills that do not require high-impact approval.

## Examples

### Compare rankings

```bash
neoag-agent \
  --message "比较 recommendation 和 NetMHCpan42 排序差异" \
  --file ranked_peptides.recommendation.tsv \
  --file ranked_peptides.netmhcpan42.tsv \
  --outdir work/agent_ranking \
  --execute
```

### Update patient report

```bash
neoag-agent \
  --message "根据最新结果更新患者沟通版报告" \
  --result-dir results/sample \
  --outdir work/agent_patient_report \
  --execute
```

## Interpretation boundary

All outputs are computational triage. Candidate neoantigens are not confirmed neoantigens and are not treatment prescriptions.
