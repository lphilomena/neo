---
name: neoag-sliding-run
description: Run NeoAg sliding-window SNV/InDel variant peptide extraction and ranking from a somatic VCF. Use when the user asks to run sliding-window variant peptides, run-full snv_indel, VEP annotation, short peptide generation, NetMHCpan/MHCflurry scoring, or ranked event/peptide reports from a VCF.
---

# NeoAg Sliding-window SNV/InDel Run

## Use This Skill When

Use this skill for NeoAg SNV/InDel workflows that start from a somatic VCF and need sliding-window mutant peptides, binding prediction, ranking, validation plans, and reports.

Prefer the staged workflow when the user wants debugging, production traceability, or partial reruns. Use the one-command `run-full` path for smoke tests and straightforward repeat runs.

## Required Inputs

Before running, infer or ask for:

- `sample_id`
- `variants_vcf`
- `tumor_sample_name`
- `hla_alleles`
- `outdir`
- `profile` (`default` unless specified)
- whether the VCF already has VEP `CSQ`
- whether to use real tools or stub mode
- optional `normal_expression`, `normal_hla_ligands`, and `normal_proteome_fasta`

Never invent real patient paths or HLA alleles. Ask if they are missing.

## Environment Bootstrap

Use this section only when the environment is not installed, `neoag` is missing, imports fail, or the user explicitly asks to prepare the machine. Do not reinstall tools on every run.

Before installing licensed or large tools, confirm the user's intended mode:

- **Fixture/stub mode**: install only the Python package and lightweight env.
- **Production local predictor mode**: also install/configure VEP, VEP cache/plugins, NetMHCpan, and MHCflurry model data.
- **Precomputed predictor mode**: skip local NetMHCpan/MHCflurry execution if the user already has predictor outputs.

### Minimal bootstrap

```bash
source conf/tools.env.sh
python -m pip install -e '.[test]'
neoag --help
```

### Lightweight conda bootstrap

Use this for development or smoke tests:

```bash
NEOAG_TOOLS_LITE=1 bash scripts/setup_tools_env.sh
source conf/tools.env.sh
python -m pip install -e '.[test]'
neoag check-tools
```

### Production bootstrap

Use this when the run needs real local VEP and binding prediction:

```bash
bash scripts/setup_tools_env.sh
bash scripts/install_vep.sh
# VEP cache may be large; use an existing cache if available.
bash scripts/install_vep_cache.sh
# After install, point NEOAG_VEP_CACHE at the cache root (~/.vep by default):
# export NEOAG_VEP_CACHE="${HOME}/.vep"
# export NEOAG_VEP_CACHE_VERSION=105
# NetMHCpan requires a licensed DTU tarball.
bash scripts/install_netmhcpan.sh /path/to/netMHCpan-4.2c.Linux.tar.gz
source conf/tools.env.sh
python -m pip install -e '.[test]'
neoag check-tools
```

If the site already has references or licensed tools, prefer setting local paths in `conf/tools.env.local.sh` over reinstalling.

## Preflight

Run from the repository root:

```bash
source conf/tools.env.sh
python -m pip install -e .
neoag check-tools
```

For production runs with auto VEP annotation, verify:

- `NEOAG_VEP_BIN`
- `NEOAG_VEP_CACHE` (cache root; release dir is `$NEOAG_VEP_CACHE/homo_sapiens/105_GRCh38/`)
- `NEOAG_VEP_CACHE_VERSION` (default `105`)
- `NEOAG_VEP_PLUGINS`
- `NEOAG_REFERENCE_FASTA`
- `NETMHCPAN_HOME` or `NEOAG_NETMHCPAN_BIN`

For VEP plugins, `NEOAG_VEP_PLUGINS` should contain `Wildtype.pm` and `Frameshift.pm`.

## Workflow Choice

### One-command path

Use `neoag run-full` when the user wants a single command and accepts automatic VEP annotation if `CSQ` is missing.

Create a private config:

```bash
cat > conf/run.<sample_id>.sliding.private.toml <<'TOML'
[sample]
id = "<sample_id>"
profile = "default"

[tools]
stub = false
enabled = ["netmhcpan", "mhcflurry"]
immunogenicity_stub = false

[inputs]
entry_mode = "snv_indel"
variant_peptide_extraction = true
variants_vcf = "<variants_vcf>"
tumor_sample_name = "<tumor_sample_name>"
hla_alleles = ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02"]
extract_appm_from_vcf = false
normal_expression = "resources/normal_expression.example.tsv"
normal_hla_ligands = "resources/normal_hla_ligands.example.tsv"
TOML
```

Then run:

```bash
neoag run-full \
  --config conf/run.<sample_id>.sliding.private.toml \
  --outdir results/<sample_id>_sliding
```

### Staged production/debug path

Use the staged path when the user wants the workflow split into:

1. VEP annotation
2. short peptide generation
3. tool scoring, ranking, validation plan, and reports

Read `staged-workflow.md` for the exact command templates before executing staged runs.

## Key Outputs

Report these paths at the end:

- `<outdir>/upstream/tools/variant_peptides.tsv`
- `<outdir>/upstream/tools/variant_peptides.annotated.tsv` (created during peptide extraction, then refreshed after tool scoring)
- `<outdir>/upstream/parsed/raw_events.tsv`
- `<outdir>/upstream/parsed/raw_peptides.tsv`
- `<outdir>/scoring/ranked_events.tsv`
- `<outdir>/scoring/ranked_peptides.tsv`
- `<outdir>/scoring/validation_plan.tsv`
- `<outdir>/reports/evidence_report.html`
- `<outdir>/reports/evidence_report.patient.html`
- `<outdir>/reports/evidence_report.technical.html`

## Validation

After running:

```bash
test -s <outdir>/upstream/parsed/raw_events.tsv
test -s <outdir>/upstream/parsed/raw_peptides.tsv
test -s <outdir>/scoring/ranked_peptides.tsv
test -s <outdir>/reports/evidence_report.technical.html
```

If VEP fails, debug only the VEP stage first. If NetMHCpan/MHCflurry fails, do not rerun VEP or peptide extraction unless the inputs changed. After successful scoring in staged mode, refresh `variant_peptides.annotated.tsv` with `.agents/skills/neoag-sliding-run/scripts/refresh_variant_peptides_annotated.py`.

## References

- For staged commands, read `staged-workflow.md`.
- For the full skill decomposition diagram, read `workflow-breakdown.md`.
- For refreshing the final annotated peptide catalog, use `scripts/refresh_variant_peptides_annotated.py`.
- For project setup and tool paths, see `README.md` and `docs/TOOLS_SETUP.md`.
