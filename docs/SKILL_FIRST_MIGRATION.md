# Skill-first migration package

This document is the hand-off guide for moving Project B to a new machine with
a programming agent. The agent should read the repo-scoped skills first, then
use manifests, Doctor and dry-run planning. It should not infer patient paths,
tool locations or reference locations from old local runs.

## Migration principle

Use this order on every new machine:

1. Read the skill registry and selected `SKILL.md` files.
2. Create local manifests from templates.
3. Install the Python package entry points.
4. Run skill validation.
5. Run read-only Doctor.
6. Run `pipeline-full` dry-run.
7. Execute only after the user approves high-risk or heavy tasks.

The source package must not contain large references, licensed binaries, patient
data, conda packs, `work/`, `results/`, Nextflow cache or machine-specific
absolute paths.

## Files the new-machine agent should read first

Start with these files from the project root:

```text
.agents/config/skills_registry.abcd.json
docs/SKILLS_TAXONOMY_ABCD.md
docs/CONTROLLED_EXECUTION_PHASE0_2.md
docs/INSTALL_AND_DATA.md
docs/TOOL_INVENTORY.md
```

Then open only the relevant skill SOPs under `.agents/skills/<skill>/SKILL.md`.
For example, VCF-to-ranking work should usually read `neoag-vcf`,
`neoag-presentation`, `neoag-ranking`, `neoag-doctor` and
`neoag-pipeline-full`.

## One-command bootstrap

For a target-machine programming agent, first read the dedicated deployment skill:

```text
.agents/skills/neoag-remote-deploy/SKILL.md
```

This skill provides the fixed deployment SOP and bundled scripts for checksum, preflight, core install, smoke tests, manifest generation, Doctor and deployment reporting.

For a fresh checkout, run:

```bash
bash scripts/bootstrap_agent_deploy.sh
```

The script is intentionally conservative. It performs local setup and checks,
but does not download references, install licensed tools, submit HPC jobs,
delete files or run heavy production workflows.

Useful options:

```bash
bash scripts/bootstrap_agent_deploy.sh --prefix conf --outdir work/agent_deploy
bash scripts/bootstrap_agent_deploy.sh --skip-install
bash scripts/bootstrap_agent_deploy.sh --mini-smoke
```

Expected outputs:

```text
conf/tools_manifest.yaml
conf/reference_manifest.yaml
conf/sample_manifest.yaml
work/agent_deploy/skill_validate/
work/agent_deploy/doctor/
work/agent_deploy/pipeline_full/
work/agent_deploy/agent_deploy_summary.md
```

## Manual bootstrap steps

If the agent must run the steps manually:

```bash
python -m pip install -e '.[test]'

mkdir -p conf
cp -n configs/controlled_execution/tools_manifest.example.yaml conf/tools_manifest.yaml
cp -n configs/controlled_execution/reference_manifest.example.yaml conf/reference_manifest.yaml
cp -n configs/controlled_execution/sample_manifest.example.yaml conf/sample_manifest.yaml

neoag-skill validate --root . --outdir work/skill_validate

neoag-doctor \
  --project-root . \
  --tools-manifest conf/tools_manifest.yaml \
  --reference-manifest conf/reference_manifest.yaml \
  --sample-manifest conf/sample_manifest.yaml \
  --outdir work/doctor \
  --dry-run

neoag-pipeline-full \
  --sample-manifest conf/sample_manifest.yaml \
  --tools-manifest conf/tools_manifest.yaml \
  --reference-manifest conf/reference_manifest.yaml \
  --outdir work/pipeline_full \
  --profile local
```

Only after reviewing `doctor/recommended_fixes.md`,
`doctor/blocking_issues.tsv` and `pipeline_full/pipeline_plan.md` should the
agent suggest tool installation, reference staging or real execution.

## Manifest rules

All machine-specific paths belong in manifests or local env overrides, not in
skills or source code:

```text
conf/tools_manifest.yaml
conf/reference_manifest.yaml
conf/sample_manifest.yaml
conf/tools.env.local.sh
```

Do not hard-code paths such as:

```text
/home/na/...
/mnt/zjl-bgi-zzb/...
/root/miniconda3/...
patient-specific BAM/FASTQ/VCF paths
```

## Risk policy for agents

Low-risk actions can be run by default:

- read docs and skills;
- create local manifest copies;
- install editable Python entry points into the active environment;
- run `neoag-skill validate`;
- run `neoag-doctor --dry-run`;
- run `neoag-pipeline-full` without `--execute`.

High-risk actions require explicit human approval:

- installing external tools or system packages;
- downloading large references;
- using licensed tools such as NetMHCpan, NetMHCstabpan, LOHHLA or Novoalign;
- submitting HPC jobs;
- running `pipeline-full --execute` on real data;
- overwriting existing results;
- deleting files.

## New-machine agent prompt

Give this prompt to another programming agent:

```text
You are deploying Project B on a new machine. Start by reading
.agents/config/skills_registry.abcd.json, docs/SKILLS_TAXONOMY_ABCD.md,
docs/CONTROLLED_EXECUTION_PHASE0_2.md and docs/INSTALL_AND_DATA.md.

Do not guess old server paths. Do not run heavy workflows yet. First run
bash scripts/bootstrap_agent_deploy.sh. Then inspect the generated Doctor and
pipeline-full dry-run outputs. Update conf/tools_manifest.yaml,
conf/reference_manifest.yaml and conf/sample_manifest.yaml for this machine.
Only propose high-risk actions after showing the missing tools/references and
asking for approval.
```

## Success criteria

A machine is migration-ready when:

- `neoag-skill validate` returns `PASS`;
- Doctor is `READY` or an explicitly accepted `PARTIAL`;
- `pipeline-full` dry-run writes a complete DAG and run manifest;
- no source file or skill SOP contains machine-specific patient paths;
- required external tools and references are declared in manifests.

