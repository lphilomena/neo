# NeoAg Event Pipeline v0.4.3 Online Release

NeoAg Event Pipeline is a research-oriented neoantigen prioritization pipeline. It converts SNV/InDel, fusion, splice, structural-variant, and peptide-only candidates into standardized event and peptide-HLA tables, then layers presentation, APPM, CCF, safety, immune-escape, validation-plan, and report evidence.

This package is a lightweight online release. It includes source code, CLI entry points, Nextflow workflows, tests, fixtures, profiles, setup scripts, and documentation. It does not bundle large references, licensed tools, conda environments, cached work directories, real patient data, or production results.

Important boundary: the pipeline produces computational triage and validation-planning outputs. It does not make clinical diagnoses, clinical resistance calls, or validated treatment recommendations.

## What It Does

The pipeline can:

- Parse pVACtools-like SNV/fusion/splice outputs into `raw_events.tsv` and `raw_peptides.tsv`.
- Generate sliding-window variant peptides from VEP-annotated VCFs, with optional automatic VEP annotation when CSQ annotations are missing.
- Score MHC presentation evidence from NetMHCpan, MHCflurry, and optional stability/immunogenicity tools.
- Build APPM 2.0 evidence, including input completeness, conflicts, peptide modifiers, and immune-context annotations.
- Estimate CCF/clonality from purity, CNV, and VAF context.
- Build peptide safety evidence from normal expression, normal ligandome, normal junction, matched-normal, and reference-proteome context.
- Build immune-escape evidence from HLA LOH, APPM, CCF, B2M/JAK/APM context, and related evidence tables.
- Generate long-peptide and minigene validation designs for frameshift, splice, exon-junction, fusion, and SV candidates.
- Produce both patient-facing and technical HTML reports.
- Run fixture workflows through the CLI or the included Nextflow wrappers.

The `.v03.tsv` suffix in ranked outputs is a schema-compatibility label. It is not the software version. The current release is v0.4.3 and writes v03-compatible tables so older downstream scripts can keep reading the same filenames.

## Agent Skills And Coordinator

This release includes a repo-scoped agent skills pack under `.agents/skills/` and a lightweight coordinator CLI:

```bash
neoag-agent --message "compare recommendation and NetMHCpan42 rankings" --result-dir results/sample --outdir work/agent_plan
```

Default mode is dry-run planning. Add `--execute` for supported low-impact skills. See `docs/AGENT_SKILLS_P0_P1.md` for the skill list, expected inputs, outputs, and interpretation boundaries.

For moving this package to a new machine with another programming agent, use the skill-first migration flow: read `.agents/config/skills_registry.abcd.json`, create local manifests, run Doctor, then run `pipeline-full` dry-run. See `docs/SKILL_FIRST_MIGRATION.md` and the helper:

```bash
bash scripts/bootstrap_agent_deploy.sh
```

For stricter hand-off to a target-machine programming agent, use the dedicated deployment skill at `.agents/skills/neoag-remote-deploy/SKILL.md`. For a fresh production-like machine, the preferred consolidated entrypoint is:

```bash
bash .agents/skills/neoag-remote-deploy/scripts/16_install_new_machine.sh \
  --asset-source-host na@10.200.50.134 \
  --allow-download \
  --execute
```

Add `--standard` for a broader common tool set, and add `--run-real-vcf-smoke --real-vcf-smoke-top-n 1` when a post-install real VCF smoke test is approved. The default installer pins VEP to Ensembl release 105 (`--vep-version 105`) to match the `homo_sapiens/105_GRCh38` cache.

## New Machine Migration

Use the machine-readable manifests as the source of truth. README is only the
human navigation layer:

- `configs/assets/production_assets.tsv`: large assets to synchronize, including
  references, VEP cache, BigMHC/DeepImmuno models, EasyFuse references, LOHHLA,
  and licensed-tool assets.
- `configs/references/reference_manifest.yaml`: reference paths, genome build,
  required markers, and VEP cache version checks.
- `configs/tools/tools_manifest.yaml`: tools, environments, containers,
  license/distribution flags, and smoke commands.

For a fresh target machine, clone this branch and run the consolidated installer
from the checkout:

```bash
mkdir -p /root/neo/src

git clone --branch na0707_upload_release \
  https://github.com/lphilomena/neo.git \
  /root/neo/src/na0707_upload_release

cd /root/neo/src/na0707_upload_release
```

Recommended production-like install:

```bash
bash .agents/skills/neoag-remote-deploy/scripts/13_install_readme_tools.sh \
  --project-root /root/neo/src/na0707_upload_release \
  --tools-root /root/neo/env_tool \
  --licensed-root /root/neo/licensed_tools \
  --reference-root /root/neo/neodata4git \
  --conda-base /root/neo/env_tool/miniforge3 \
  --asset-manifest configs/assets/production_assets.tsv \
  --reference-manifest configs/references/reference_manifest.yaml \
  --sync-assets \
  --asset-source-host na@10.200.50.134 \
  --all-open \
  --verify \
  --allow-download \
  --execute
```

To run the default real VCF smoke test after installation, add
`--run-real-vcf-smoke`:

```bash
bash .agents/skills/neoag-remote-deploy/scripts/13_install_readme_tools.sh \
  --project-root /root/neo/src/na0707_upload_release \
  --tools-root /root/neo/env_tool \
  --licensed-root /root/neo/licensed_tools \
  --reference-root /root/neo/neodata4git \
  --conda-base /root/neo/env_tool/miniforge3 \
  --asset-manifest configs/assets/production_assets.tsv \
  --reference-manifest configs/references/reference_manifest.yaml \
  --sync-assets \
  --asset-source-host na@10.200.50.134 \
  --all-open \
  --verify \
  --run-real-vcf-smoke \
  --allow-download \
  --execute
```

The default asset source is:

```bash
/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git
```

The installer defaults to Miniforge3 and pins VEP to Ensembl release 105 to
match `homo_sapiens/105_GRCh38` in the VEP cache. Licensed tools such as
NetMHCpan, NetMHCstabpan, Polysolver, and Novoalign are represented in the asset
and tools manifests, but the operator must confirm the license permits use on
the new machine.

After installation, verify the machine explicitly:

```bash
source /root/neo/src/na0707_upload_release/conf/tools.env.sh

neoag-v03 check-tools

python3 scripts/verify_reference_manifest.py \
  configs/references/reference_manifest.yaml \
  --vep-version 105
```

The older wrapper remains available for short-form installs:

```bash
bash .agents/skills/neoag-remote-deploy/scripts/16_install_new_machine.sh \
  --standard \
  --asset-source-host na@10.200.50.134 \
  --allow-download \
  --execute
```

Prefer the explicit `13_install_readme_tools.sh` command above when recording a
production migration, because it names the tool root, licensed-tool root,
reference root, and manifest files directly.

## Quick Start

Run these commands from the project root:

```bash
python -m pip install -e .
neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001
```

Important demo outputs include:

- `work/demo_v043/scoring/ranked_peptides.v03.tsv`
- `work/demo_v043/scoring/ranked_events.v03.tsv`
- `work/demo_v043/scoring/validation_plan.v03.tsv`
- `work/demo_v043/reports/evidence_report.v03.html`
- `work/demo_v043/reports/evidence_report.patient.html`
- `work/demo_v043/reports/evidence_report.technical.html`
- `work/demo_v043/appm/appm_summary.tsv`
- `work/demo_v043/appm/appm_peptide_modifiers.tsv`
- `work/demo_v043/clonality/ccf_lite.tsv`
- `work/demo_v043/safety/peptide_safety.tsv`
- `work/demo_v043/immune_escape/peptide_escape_flags.tsv`

For tests:

```bash
python -m pip install -e '.[test]'
pytest -q
```

The default test command intentionally skips integration, benchmark, and external-tool tests.

## Common Run Commands

### Prepare The Environment

For fixture-only development:

```bash
python -m pip install -e '.[test]'
pytest -q
neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001
```

For runs that need external tools:

```bash
bash scripts/setup_tools_env.sh
source conf/tools.env.sh
python -m pip install -e '.[test]'
neoag-v03 check-tools
```

For a smaller development/test environment:

```bash
NEOAG_TOOLS_LITE=1 bash scripts/setup_tools_env.sh
source conf/tools.env.sh
python -m pip install -e '.[test]'
pytest -q
```

### Run From Existing pVAC-like Tables

Use this when you already have pVACseq/pVACfuse/pVACsplice-like aggregated tables:

```bash
neoag-v03 run-v03 \
  --outdir results/sample \
  --sample-id SAMPLE001 \
  --profile default \
  --pvac data/fixtures/pvacseq_aggregated.tsv \
  --immunogenicity-stub
```

### Run From Pre-built Raw Intermediates

Use this when `parsed/raw_events.tsv` and `parsed/raw_peptides.tsv` already exist:

```bash
neoag-v03 run-v03 \
  --outdir results/sample \
  --sample-id SAMPLE001 \
  --profile default \
  --raw-events results/sample/parsed/raw_events.tsv \
  --raw-peptides results/sample/parsed/raw_peptides.tsv \
  --netmhcpan results/sample/presentation/netmhcpan.xls \
  --mhcflurry results/sample/presentation/mhcflurry.csv \
  --expression results/sample/parsed/expression.tsv \
  --hla-loh results/sample/tools/hla_loh.tsv \
  --purity results/sample/tools/purity.tsv \
  --cnv results/sample/tools/cnv_segments.tsv
```

### Sliding-window Variant Peptides To Ranking

Use this path when you have a somatic SNV/InDel VCF and want to generate mutant peptides by sliding window, predict peptide-HLA presentation, and rank event/peptide candidates.

If the VCF already contains VEP `CSQ` annotations, the pipeline uses it directly. If `CSQ` is missing, `run-full` will run VEP annotation first when VEP, cache, reference FASTA, and plugins are configured.

```bash
cat > conf/run.sliding.private.toml <<'TOML'
[sample]
id = "SAMPLE001"
profile = "default"

[tools]
stub = false
enabled = ["netmhcpan", "mhcflurry"]
immunogenicity_stub = false

[inputs]
entry_mode = "snv_indel"
variant_peptide_extraction = true
variants_vcf = "/path/to/sample.somatic.pass.vcf.gz"
tumor_sample_name = "TUMOR"
hla_alleles = ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:02"]
extract_appm_from_vcf = false
normal_expression = "resources/normal_expression.example.tsv"
normal_hla_ligands = "resources/normal_hla_ligands.example.tsv"
TOML

neoag-v03 run-full \
  --config conf/run.sliding.private.toml \
  --outdir results/SAMPLE001_sliding
```

Key outputs:

- `results/SAMPLE001_sliding/upstream/tools/variant_peptides.tsv`
- `results/SAMPLE001_sliding/upstream/tools/variant_peptides.annotated.tsv`
- `results/SAMPLE001_sliding/upstream/parsed/raw_events.tsv`
- `results/SAMPLE001_sliding/upstream/parsed/raw_peptides.tsv`
- `results/SAMPLE001_sliding/scoring/ranked_events.v03.tsv`
- `results/SAMPLE001_sliding/scoring/ranked_peptides.v03.tsv`
- `results/SAMPLE001_sliding/scoring/validation_plan.v03.tsv`
- `results/SAMPLE001_sliding/reports/evidence_report.v03.html`
- `results/SAMPLE001_sliding/reports/evidence_report.patient.html`
- `results/SAMPLE001_sliding/reports/evidence_report.technical.html`

Manual debug path for variant peptide extraction:

```bash
neoag-v03 extract-variant-peptides \
  --input-vcf /path/to/sample.vep.annotated.vcf.gz \
  --output results/SAMPLE001_sliding/upstream/tools/variant_peptides.tsv \
  --sample-id SAMPLE001 \
  --lengths 8,9,10,11 \
  --mini-len 27 \
  --hla-alleles HLA-A*02:01,HLA-B*07:02,HLA-C*07:02 \
  --tumor-sample-name TUMOR \
  --normal-proteome-fasta /path/to/Homo_sapiens.GRCh38.pep.all.fa \
  --filter-normal-proteome
```

For smoke tests without licensed predictors, set `stub = true` in the TOML or add `--immunogenicity-stub` to direct `run-v03` calls. For production ranking, use real NetMHCpan/MHCflurry outputs and real normal-expression/normal-ligand evidence instead of fixture resources.

### Build Standard Evidence Sidecars

```bash
neoag-v03 build-evidence-layer \
  --outdir results/sample \
  --profile default \
  --sample-id SAMPLE001 \
  --raw-events results/sample/parsed/raw_events.tsv \
  --raw-peptides results/sample/parsed/raw_peptides.tsv \
  --expression results/sample/parsed/gene_expression.tsv \
  --rna-vaf results/sample/parsed/rna_vaf.tsv \
  --rna-junction results/sample/parsed/rna_junctions.tsv \
  --fusion-evidence results/sample/parsed/fusion_evidence.tsv \
  --normal-expression resources/normal_expression.example.tsv \
  --normal-hla-ligands resources/normal_hla_ligands.example.tsv
```

### HLA LOH Conversion And Cross-check

```bash
neoag-v03 convert-lohhla \
  -i results/sample/tools/LOHHLA.HLAlossPrediction_CI.xls \
  -o results/sample/tools/lohhla.hla_loh.tsv

neoag-v03 convert-spechla \
  -i results/sample/tools/merge.hla.copy.txt \
  -o results/sample/tools/spechla.hla_loh.tsv

neoag-v03 crosscheck-hla-loh \
  --lohhla-hla-loh results/sample/tools/lohhla.hla_loh.tsv \
  --spechla-hla-loh results/sample/tools/spechla.hla_loh.tsv \
  --out results/sample/tools/hla_loh.crosscheck.tsv \
  --consensus-out results/sample/tools/hla_loh.consensus.tsv
```

### Generate Reports

Generate the default combined report plus patient and technical audience-specific reports:

```bash
neoag-v03 report-v03 \
  --profile default \
  --ranked-events results/sample/scoring/ranked_events.v03.tsv \
  --ranked-peptides results/sample/scoring/ranked_peptides.v03.tsv \
  --appm-summary results/sample/appm/appm_summary.tsv \
  --validation-plan results/sample/scoring/validation_plan.v03.tsv \
  --outdir results/sample \
  --audience both \
  --out results/sample/reports/evidence_report.v03.html
```

### Nextflow Fixture Run

Use the project wrapper rather than calling `nextflow` directly. The wrapper prioritizes the current checkout's `bin/neoag-v03`, sets project paths, and avoids writing Nextflow metadata into a root-owned location.

```bash
export NXF_HOME=/path/to/writable/nextflow_cache
bin/neoag-nextflow -version
bin/neoag-nextflow run workflows/main.nf \
  -w /tmp/neoag_nf_work \
  --pvac_files data/fixtures/pvacseq_aggregated.tsv \
  --outdir results/demo_nf \
  --sample_id NF_DEMO
```

For command-specific options:

```bash
neoag-v03 <command> --help
```

## Configuration Files

Real deployment paths should be kept in local/private files. Start with the templates below and copy them before editing site-specific values.

| File | Purpose | Edit for real data? | Commit/package? |
| --- | --- | --- | --- |
| `conf/tools.env.sh` | Main environment entry point. Sets project paths, conda env names, tool roots, VEP cache fallback, and wrapper `PATH`. | Usually no; override locally instead. | Yes |
| `conf/tools.env.local.example.sh` | Template for private site paths such as patient data roots, shared references, licensed tool installs, and cache directories. | Copy to `conf/tools.env.local.sh` and edit. | Example yes; copied local file no |
| `conf/site.config.example` | Site/cluster/Nextflow executor template. | Copy to `conf/site.config` and edit. | Example yes; copied local file no |
| `conf/run.private.example.toml` | Private real-sample run configuration template. | Copy to a private TOML and edit. | Example yes; copied local file no |
| `conf/run.snv_wes.example.toml` | WES SNV workflow example with Mutect2/annotation inputs. | Copy before using with real BAM/VCF paths. | Example yes |
| `conf/run.stub.toml` | Lightweight stub/demo upstream config. | No for production; useful for smoke tests. | Yes |
| `conf/*.example.toml` | Workflow-specific examples for SV, fusion, splice, peptide-only, or site modes. | Copy and edit. | Example yes |

Typical setup pattern:

```bash
cp conf/tools.env.local.example.sh conf/tools.env.local.sh
cp conf/run.private.example.toml conf/run.sample.private.toml
# Edit both files with site paths and sample inputs.
source conf/tools.env.sh
neoag-v03 check-tools
neoag-v03 run-full --config conf/run.sample.private.toml --outdir results/sample
```

Important variables commonly set by `conf/tools.env.sh` or local overrides:

| Variable | Meaning |
| --- | --- |
| `NEOAG_PROJECT_ROOT` | Project checkout root. |
| `NEOAG_TOOLS_ROOT` | Root for external tools, wrappers, and local artifact bundles. |
| `NEOAG_CONDA_BASE` | Miniforge/Mambaforge install path. |
| `NEOAG_CONDA_ENV` | Main Python CLI environment. |
| `NEOAG_VEP_ENV` | VEP conda environment name/path. |
| `NEOAG_VEP_BIN` | VEP executable or wrapper path. |
| `NEOAG_VEP_CACHE` | VEP offline cache **root** directory. Must contain `homo_sapiens/<version>_GRCh38/` (this project uses `105_GRCh38`). Do **not** point this variable at the `105_GRCh38` subdirectory itself. |
| `NEOAG_VEP_CACHE_VERSION` | Ensembl cache release passed to VEP as `--cache_version` (default `105`). |
| `NEOAG_VEP_PLUGINS` | VEP plugin directory containing `Wildtype.pm` and `Frameshift.pm`. |
| `NEOAG_REFERENCE_FASTA` | GRCh38 FASTA used by VEP/GATK/SV peptide workflows. |
| `NEOAG_NORMAL_PROTEOME_FASTA` | Normal/reference proteome FASTA used by peptide safety filtering. |
| `NETMHCPAN_HOME` / `NEOAG_NETMHCPAN_BIN` | NetMHCpan install and executable path. |
| `NEOAG_NETMHCPAN_TMPDIR` | Short temporary directory for NetMHCpan. |
| `NETMHCSTABPAN_HOME` | NetMHCstabpan install path. |
| `LOHHLA_HOME`, `POLYSOLVER_HOME`, `NOVOALIGN_LICENSE_FILE` | LOHHLA and dependency paths/licenses. |
| `FACETS_HOME`, `NEOAG_DBSNP_VCF` | FACETS scripts and dbSNP/common SNP VCF for `snp-pileup`. |
| `NEOAG_ASCAT_ENV`, `ASCAT_HOME` | ASCAT conda env and wrapper path. |
| `NXF_HOME` | Writable Nextflow cache; required for clean online/offline workflow runs. |

Do not commit or package local/private files such as:

- `conf/tools.env.local.sh`
- `conf/site.config`
- `conf/private/*`
- `conf/*.private.toml`
- files containing patient identifiers, absolute clinical data paths, cluster credentials, or licensed tool paths

## Installation, Tools, And Data

### Base System Environment

Recommended baseline:

| Component | Recommended | Notes |
| --- | --- | --- |
| OS | Linux x86_64 | The included wrappers and most upstream bioinformatics tools assume Linux. |
| Shell | Bash | Scripts are written for Bash. |
| CPU/RAM | 8+ CPU, 32+ GB RAM for real data | Fixture demo is small; patient-scale WES/WGS/fusion runs need more. |
| Disk | 200 GB+ for real-data deployment | VEP cache, hg38 references, fusion references, and Nextflow work dirs are large. |
| Network | Required for first install/download | Offline installs should pre-stage tarballs, conda packages, references, and Nextflow cache. |
| Java | Java 11 or newer | Required by Nextflow and some tools. |
| Conda/Mamba | Miniforge/Mambaforge recommended | Tool environments are defined under `conda/`. |

Suggested Ubuntu/Debian packages:

```bash
sudo apt-get update
sudo apt-get install -y \
  bash coreutils curl wget git tar gzip unzip bzip2 xz-utils \
  ca-certificates build-essential openjdk-17-jre-headless rsync file
```

If a migrated archive lost executable permissions:

```bash
find bin -maxdepth 1 -type f -exec chmod +x {} \;
find scripts -maxdepth 1 -type f -name '*.sh' -exec chmod +x {} \;
```

### External Tool Installation Table

Tools are optional for the fixture demo but required by specific real-data modes. Install only the tools needed for your workflow.

| Tool | Needed for | Required? | Install/download command | Key variables | Verify |
| --- | --- | --- | --- | --- | --- |
| pVACtools (`pvacseq`, `pvacfuse`, `pvacsplice`) | Upstream SNV/fusion/splice candidate generation | Optional unless using pVAC upstream modes | `bash scripts/setup_tools_env.sh` | `NEOAG_PVAC_DOCKER`, `NEOAG_PVAC_WORKDIR` | `neoag-v03 check-tools` |
| NetMHCpan 4.2 | Binding/presentation prediction | Required for real local NetMHCpan runs unless using fallback/stub | `bash scripts/install_netmhcpan.sh /path/to/netMHCpan-4.2*.tar.gz` | `NETMHCPAN_HOME`, `NETMHCpan`, `NEOAG_NETMHCPAN_BIN`, `NEOAG_NETMHCPAN_BACKEND` | `neoag-v03 check-tools` |
| MHCflurry | Binding/presentation prediction | Optional alternative/complement to NetMHCpan | `bash scripts/setup_tools_env.sh`; then `mhcflurry-downloads fetch` if needed | `NEOAG_CONDA_ENV`, `NEOAG_FORCE_CPU` | `neoag-v03 check-tools` |
| NetMHCstabpan | pMHC stability evidence | Optional | `bash scripts/install_netmhcstabpan.sh --iedb` or licensed tarball install | `NETMHCSTABPAN_HOME` | `neoag-v03 check-tools` |
| PRIME / MixMHCpred / BigMHC | Immunogenicity evidence | Optional | `bash scripts/install_immunogenicity_tools.sh` | `PRIME_HOME`, `MIXMHCPRED_HOME`, `BIGMHC_DIR`, `NEOAG_PRIME_JOBS` | `neoag-v03 check-tools` |
| DeepImmuno | Optional immunogenicity evidence | Optional | `bash scripts/install_deepimmuno.sh` | `DEEPIMMUNO_DIR` | `neoag-v03 check-tools` |
| SHERPA-Presentation | MHC-I presentation ranking; restricted/authorized asset required | Optional | `bash .agents/skills/neoag-remote-deploy/scripts/13_install_readme_tools.sh --sherpa --sherpa-source /path/to/SHERPA-Presentation --execute` | `SHERPA_PRESENTATION_HOME`, `SHERPA_PRESENTATION_BIN` | `scripts/verify_all_tools_and_refs.sh` |
| VEP | VCF annotation and peptide extraction | Required for `vep-annotate` / auto-annotation in `run-full` | `bash scripts/install_vep.sh`; cache with `bash scripts/install_vep_cache.sh` | `NEOAG_VEP_ENV`, `NEOAG_VEP_BIN`, `NEOAG_VEP_CACHE`, `NEOAG_VEP_PLUGINS`, `NEOAG_VEP_ONLINE` | `neoag-v03 check-tools` |
| GATK4 / Mutect2 | WES/WGS SNV calling | Required if starting from BAMs | `bash scripts/install_gatk.sh` | `NEOAG_GATK_ENV` | `neoag-v03 check-tools` |
| LOHHLA | HLA LOH evidence | Optional but recommended for immune-escape evidence | `bash scripts/install_lohhla.sh`; configure Polysolver/Novoalign separately | `LOHHLA_HOME`, `POLYSOLVER_HOME`, `NOVOALIGN_LICENSE_FILE` | `neoag-v03 check-tools` |
| SpecHLA | HLA copy/LOH conversion | Optional | Install externally and provide output to `convert-spechla` | Site-specific | `neoag-v03 convert-spechla --help` |
| OptiType | HLA-A/B/C typing from DNA/RNA FASTQ or BAM | Optional HLA typing cross-check | `bash scripts/install_optitype.sh` or `mamba create -n neoag-optitype -c conda-forge -c bioconda optitype glpk coincbc razers3` | `OPTITYPE_ENV`, `OPTITYPE_BIN`, `OPTITYPE_REFERENCE` | `optitype check-deps` |
| FACETS | Purity/CNV/LOH evidence | Optional but recommended for CCF/escape | `bash scripts/install_facets.sh` | `FACETS_HOME`, `NEOAG_DBSNP_VCF` | `neoag-v03 check-tools` |
| ASCAT | CNV/LOH evidence | Optional | `bash scripts/install_ascat_pyclone.sh` | `NEOAG_ASCAT_ENV`, `ASCAT_HOME` | `neoag-v03 check-tools` |
| PyClone-VI | Clonality context | Optional | `bash scripts/install_ascat_pyclone.sh` | `NEOAG_PYCLONE_ENV`, `NEOAG_PYCLONE_BIN` | `neoag-v03 check-tools` |
| STAR-Fusion / FusionCatcher / Arriba / EasyFuse | Fusion discovery | Optional; required for corresponding fusion workflows | Install/mount externally; seed EasyFuse envs only when a Nextflow conda cache exists | `NEOAG_FUSION_ENV`, `NEOAG_STAR_FUSION_HOME`, `NEOAG_CTAT_LIB_DIR`, `NEOAG_EASYFUSE_HOME`, `NEOAG_EASYFUSE_REF` | `neoag-v03 check-tools` |
| Manta / GRIDSS / SvABA / Sniffles2 | SV discovery | Optional upstream SV callers | Install externally or via site conda/modules | `NEOAG_SV_ENV`, `NEOAG_MANTA_ENV` | `neoag-v03 check-tools` |
| PURPLE / AMBER / COBALT | Purity, ploidy, CNV, LOH evidence | Optional | See `docs/TOOLS_SETUP.md` and local wrappers | `HMFTOOLS_HOME`, site-specific references | Tool-specific wrapper `--help` |
| DASH | HLA LOH / allele-specific deletion evidence | Optional | See `docs/TOOLS_SETUP.md`; model may need to be provided separately | DASH env/model path | Tool-specific wrapper |

Licensed tools such as NetMHCpan, NetMHCstabpan, LOHHLA, and Novoalign/Polysolver components may require academic or institutional approval. Do not redistribute their binaries inside the online release.

### NetMHCpan Notes

Install from the DTU-licensed Linux tarball:

```bash
mkdir -p vendor
cp /path/to/netMHCpan-4.2c.Linux.tar.gz vendor/
export NEOAG_CONDA_BASE="$(conda info --base)"
bash scripts/install_netmhcpan.sh vendor/netMHCpan-4.2c.Linux.tar.gz
source conf/tools.env.sh
neoag-v03 check-tools
netMHCpan -h | head
```

If NetMHCpan is already installed but its wrapper is broken:

```bash
bash scripts/install_netmhcpan.sh --repair
```

For patched NetMHCpan binaries on older host glibc, the `neoag-tools` environment must retain `sysroot_linux-64` and `patchelf`. The lite environment file includes these dependencies.

### VEP Notes

Install VEP and configure cache/plugins:

```bash
bash scripts/install_vep.sh
source conf/tools.env.sh
neoag-v03 check-tools
```

#### VEP cache layout

`NEOAG_VEP_CACHE` is the **cache root**, not the release subdirectory. VEP is invoked with:

```bash
vep --cache --offline \
  --dir_cache "$NEOAG_VEP_CACHE" \
  --cache_version "$NEOAG_VEP_CACHE_VERSION"
```

Expected directory layout:

```text
$NEOAG_VEP_CACHE/
└── homo_sapiens/
    └── 105_GRCh38/
        ├── info.txt
        ├── 1/ 2/ ... 22/          # per-chromosome indexed cache files
        └── ...                     # optional FASTA / auxiliary files from vep_install
```

For this project:

| Item | Value |
| --- | --- |
| Species | `homo_sapiens` |
| Assembly | `GRCh38` |
| Cache release | `105` (`NEOAG_VEP_CACHE_VERSION=105`) |
| Data package | Ensembl indexed cache `homo_sapiens_vep_105_GRCh38.tar.gz` |
| Typical size | about 12–16 GB after extraction |
| Download source | `https://ftp.ensembl.org/pub/release-105/variation/indexed_vep_cache/homo_sapiens_vep_105_GRCh38.tar.gz` |

Install a new cache, or point to an existing one:

```bash
# Fresh install (downloads to ~/.vep by default).
bash scripts/install_vep_cache.sh
export NEOAG_VEP_CACHE="${HOME}/.vep"
export NEOAG_VEP_CACHE_VERSION=105

# Or reuse a site-managed cache root.
export NEOAG_VEP_CACHE=/path/to/data/vep
export NEOAG_VEP_CACHE_VERSION=105
source conf/tools.env.sh
```

Verify:

```bash
test -f "$NEOAG_VEP_CACHE/homo_sapiens/105_GRCh38/info.txt"
```

Recommended VEP cache resolution for a portable release:

1. Set `NEOAG_REF_BUNDLE=/path/to/neodata4git` and source `$NEOAG_REF_BUNDLE/neodata4git.env.sh`.
2. Or set `NEOAG_VEP_CACHE=/path/to/data/vep` directly.
3. Keep server-specific cache paths in `conf/tools.env.local.sh`, not in tracked release files.

```bash
export NEOAG_REF_BUNDLE=/path/to/neodata4git
source "$NEOAG_REF_BUNDLE/neodata4git.env.sh"
# sets NEOAG_VEP_CACHE=$NEOAG_REF_BUNDLE/data/vep when the bundle env file is configured
```

`NEOAG_VEP_PLUGINS` should point to a directory containing `Wildtype.pm` and `Frameshift.pm`. Plain VEP can run without these plugins, but this project uses them for pVACseq-compatible WT/frameshift information and more complete peptide extraction.

### Reference Data Table

Large data should live under `NEOAG_TOOLS_ROOT`, `NEOAG_SHARED_REF_DIR`, or another site-managed reference area, not inside the source checkout.

| Data/reference | Needed for | Download/setup command | Expected variable/path | Verify |
| --- | --- | --- | --- | --- |
| VEP cache, GRCh38 release 105 | Offline VEP annotation | `bash scripts/install_vep_cache.sh` or reuse site cache | `NEOAG_VEP_CACHE=/path/to/data/vep` (cache root); release dir is `$NEOAG_VEP_CACHE/homo_sapiens/105_GRCh38/`; `NEOAG_VEP_CACHE_VERSION=105` | `test -f "$NEOAG_VEP_CACHE/homo_sapiens/105_GRCh38/info.txt"` |
| VEP plugins | WT and frameshift plugin annotations | Installed by VEP/pVAC tooling or copied from a site bundle | `NEOAG_VEP_PLUGINS=/path/to/work/vep_plugins` | `test -f "$NEOAG_VEP_PLUGINS/Wildtype.pm"` |
| GRCh38 FASTA and indices | VEP peptide extraction, GATK, SV peptide building | `bash scripts/download_ref_hg38.sh /path/to/ref/hg38` | `NEOAG_REFERENCE_FASTA=/path/to/Homo_sapiens_assembly38.fasta` | `test -f "$NEOAG_REFERENCE_FASTA"` |
| dbSNP/common SNP VCF | FACETS `snp-pileup`, some CNV workflows | Included in site reference bundle or downloaded with hg38 bundle where available | `NEOAG_DBSNP_VCF=/path/to/dbsnp_chr.vcf.gz` | `test -f "$NEOAG_DBSNP_VCF"` |
| gnomAD AF VCF and PoN | GATK Mutect2 filtering | `bash scripts/download_ref_hg38.sh /path/to/ref/hg38` or site bundle | Paths inside selected run config | `test -f /path/to/af-only-gnomad.hg38.vcf.gz` |
| Ensembl protein FASTA | Peptide safety normal/reference proteome screen | Download Ensembl GRCh38 peptide FASTA manually or from site bundle | `NEOAG_NORMAL_PROTEOME_FASTA=/path/to/Homo_sapiens.GRCh38.pep.all.fa` | `test -f "$NEOAG_NORMAL_PROTEOME_FASTA"` |
| Normal expression table | Peptide safety evidence | Site-generated TSV or fixture example | CLI argument or run config path | Check expected TSV header |
| Normal HLA ligand table | Peptide safety evidence | Site-generated TSV or fixture example | CLI argument or run config path | Check expected TSV header |
| CTAT genome lib | STAR-Fusion | Download per STAR-Fusion/CTAT docs or mount site bundle | `CTAT_GENOME_LIB`, `NEOAG_CTAT_LIB_DIR`, `NEOAG_SHARED_REF_DIR` | `test -d "$CTAT_GENOME_LIB"` |
| EasyFuse reference | EasyFuse workflow | Download per EasyFuse docs or mount site bundle | `NEOAG_EASYFUSE_REF`, `NEOAG_SHARED_REF_DIR` | `test -d "$NEOAG_EASYFUSE_REF"` |
| GTF annotation | SV/fusion peptide generation | Use GENCODE/Ensembl GTF matching reference FASTA | CLI argument `--gencode-gtf` | `test -f /path/to/genes.gtf` |
| Capture BED | WES SV Phase 1.5 | Use panel/exome capture BED | CLI argument `--capture-bed` | `test -f /path/to/capture.bed` |
| HLA allele file | Peptide prediction and SV workflows | Site HLA typing output converted to one allele per line | CLI argument `--hla` | `head /path/to/hla.txt` |


### Reference Files By Tool And Module

The table below lists the highest-priority reference files by tool/module. Put large references in `NEOAG_TOOLS_ROOT`, `NEOAG_SHARED_REF_DIR`, or another site-managed reference bundle. Do not commit licensed tools, public reference bundles, or patient data to Git.

| Tool/module | Required reference files | Optional reference files | Environment variables / recommended paths | Purpose |
| --- | --- | --- | --- | --- |
| Core GRCh38 reference | GRCh38 FASTA, `.fai`, and sequence dictionary, for example `Homo_sapiens_assembly38.fasta`, `Homo_sapiens_assembly38.fasta.fai`, `Homo_sapiens_assembly38.dict` | Capture BED or interval BED for WES/panel workflows | `NEOAG_REFERENCE_FASTA=$NEOAG_TOOLS_ROOT/data/ref/hg38/Homo_sapiens_assembly38.fasta` | Shared genome reference for GATK, VEP peptide extraction, SV/fusion peptide generation, and read-based evidence workflows. |
| VEP | GRCh38 VEP cache, expected as `homo_sapiens/105_GRCh38/info.txt`; matching reference FASTA | VEP plugins such as `Wildtype.pm` and `Frameshift.pm` | `NEOAG_VEP_CACHE=$NEOAG_TOOLS_ROOT/data/vep`; `NEOAG_VEP_CACHE_VERSION=105`; `NEOAG_VEP_PLUGINS=$NEOAG_TOOLS_ROOT/work/vep_plugins` | Offline annotation and pVACseq-compatible WT/frameshift peptide context. |
| GATK / Mutect2 | GRCh38 FASTA bundle, gnomAD AF-only VCF + `.tbi`, Panel of Normals VCF + `.tbi` | Exome/panel intervals, contamination resources | Run-config paths; recommended under `$NEOAG_TOOLS_ROOT/data/ref/hg38/` | Somatic SNV/indel calling and filtering. |
| NetMHCpan | Licensed NetMHCpan executable and its `data/` directory | None | `NETMHCPAN_HOME=$NEOAG_TOOLS_ROOT/tools/netMHCpan`; `NETMHCpan=$NETMHCPAN_HOME`; `NEOAG_NETMHCPAN_BIN=$NETMHCPAN_HOME/netMHCpan` | Primary peptide-HLA binding and presentation prediction. |
| MHCflurry | Installed MHCflurry environment and downloaded models | Custom model cache | `NEOAG_MHCFLURRY_ENV`; model cache managed by MHCflurry | Optional binding/presentation cross-check. |
| PRIME / MixMHCpred / BigMHC / DeepImmuno | Tool-specific executables or model directories | Model/version-specific data files | `PRIME_HOME`, `MIXMHCPRED_HOME`, `BIGMHC_DIR`, `DEEPIMMUNO_DIR` | Immunogenicity and presentation evidence layers. |
| NetMHCstabpan / IEDB shim | Licensed NetMHCstabpan executable or configured IEDB-compatible shim | None | `NETMHCSTABPAN_HOME=$NEOAG_TOOLS_ROOT/tools/netMHCstabpan` | Optional peptide stability evidence. |
| FACETS | Common biallelic SNP/dbSNP VCF + `.tbi` matching the BAM reference build | Downsampled SNP pileups or curated SNP-only panels such as Omni2.5/common SNP resources | `NEOAG_DBSNP_VCF=$NEOAG_TOOLS_ROOT/data/facets/reference/common_snp.hg38.vcf.gz` | Tumor purity, ploidy, CNV, and LOH evidence. |
| ASCAT | hg38 loci and alleles resources; GC and RT correction files for ASCAT v3/prepareHTS modes | Platform-specific SNP panels | Site variables or run config paths such as `ASCAT_LOCI_PREFIX`, `ASCAT_ALLELES_PREFIX`, `ASCAT_GC_FILE`, `ASCAT_RT_FILE` | Independent purity/ploidy/CNV/LOH cross-check. |
| Sequenza | GRCh38 FASTA, GC wiggle/file, matched tumor/normal BAMs | Precomputed per-chromosome `.seqz.gz` blocks | Run script arguments and Sequenza config paths | Independent purity, ploidy, and CNV estimation. |
| PURPLE / AMBER / COBALT | HMF reference bundle matching GRCh38; tumor/normal BAMs; AMBER and COBALT resources | Somatic VCF for richer PURPLE interpretation | HMF/PURPLE config paths; keep the bundle outside Git | Purity, ploidy, CNV, LOH, and QC cross-check. |
| STAR-Fusion / CTAT | CTAT genome library directory | RNA FASTQ QC and caller-specific caches | `CTAT_GENOME_LIB` or `$NEOAG_TOOLS_ROOT/data/ref/ctat/current` | RNA fusion discovery. |
| Arriba | GRCh38 FASTA, matching GTF, STAR index/reference files | Blacklist/known fusion resources | `GTF=/path/to/gencode.gtf`; Arriba/STAR paths in run config | RNA fusion discovery and SV/fusion peptide generation. |
| EasyFuse | EasyFuse reference bundle | Prebuilt conda/Docker environment cache | `NEOAG_EASYFUSE_REF` or `$NEOAG_SHARED_REF_DIR/easyfuse_ref_v4`; `NEOAG_EASYFUSE_HOME` | RNA fusion discovery and fusion evidence cross-check. |
| LOHHLA / Polysolver | Polysolver distribution, Novoalign license, patient HLA FASTA or Polysolver HLA output | FACETS/ASCAT/PURPLE purity/CNV evidence | `POLYSOLVER_HOME=/path/to/polysolver`; `NOVOALIGN_LICENSE_FILE=/path/to/novoalign.lic` | HLA LOH evidence and HLA typing support. |
| SpecHLA | SpecHLA database/reference files and tumor/normal BAM or extracted HLA reads | User-supplied purity/ploidy and HLA typing | SpecHLA install path plus run-script arguments | HLA typing and HLA LOH cross-check. |
| HLA-LA | `PRG_MHC_GRCh38_withIMGT` graph directory | Long-read BAM input for higher-confidence HLA typing | HLA-LA graph argument/path; keep graph outside Git | Independent HLA typing cross-check, especially useful for long-read data. |
| OptiType | HLA-I reference bundled with OptiType: `hla_reference_dna.fasta`, `hla_reference_rna.fasta`, and allele parquet files | RNA-derived FASTQ input | `OPTITYPE_ENV=$NEOAG_CONDA_BASE/envs/neoag-optitype`; `OPTITYPE_BIN=$OPTITYPE_ENV/bin/optitype`; `OPTITYPE_REFERENCE=$OPTITYPE_ENV/share/optitype/data`; optional bundle link: `$NEOAG_REF_BUNDLE/data/hla/optitype_reference` | HLA-A/B/C typing cross-check from DNA or RNA short reads. |
| Normal safety evidence | Normal/reference proteome FASTA; normal expression TSV; normal HLA ligand TSV | Site-specific normal tissue evidence tables | `NEOAG_NORMAL_PROTEOME_FASTA=/path/to/Homo_sapiens.GRCh38.pep.all.fa`; CLI/run-config paths for TSV evidence | Filtering or annotation against normal protein and normal tissue/ligand evidence. |
| RNA evidence | RNA allele-count/RNA VAF TSV; RNA expression matrix or TPM table | Targeted junction/fusion validation tables | CLI arguments such as `--rna-vaf`; generated by `scripts/rna_allele_counts_pysam.py` or RNA workflows | RNA support, RNA VAF, and expression evidence for candidate events. |
| SV/long-read module | GRCh38 FASTA, matching GTF, HLA allele file, SV VCF or long-read caller output | Capture BED for WES-like SV workflows | `--gencode-gtf`, `--hla`, `--capture-bed`; run config paths | SV/fusion event parsing, peptide generation, and ranking. |

### RNA FASTQ to TPM

RNA expression evidence can be reviewed through the `neoag-rna-fastq-to-tpm` agent skill. When paired RNA FASTQ files and Salmon references are configured, the skill can also execute Salmon and generate a normalized gene TPM table.

```bash
export SALMON_BIN=salmon
export SALMON_INDEX="$NEOAG_DATA_ROOT/data/ref/rna/salmon_index"
export SALMON_TX2GENE="$NEOAG_DATA_ROOT/data/ref/rna/tx2gene.tsv"

bin/neoag-llm-agent \
  --message "Generate RNA TPM with Salmon, sample_id=SAMPLE01" \
  --file tumor_R1.fastq.gz \
  --file tumor_R2.fastq.gz \
  --outdir results/llm_agent_rna_tpm \
  --mode execute-with-approval \
  --allow-high-risk
```

Direct script entry points are also available:

```bash
bash scripts/run_salmon_fastq_to_tpm.sh \
  --fastq1 tumor_R1.fastq.gz \
  --fastq2 tumor_R2.fastq.gz \
  --sample-id SAMPLE01 \
  --outdir results/rna_tpm/SAMPLE01

bash scripts/run_rsem_fastq_to_tpm.sh --help
```

Required Salmon references are a transcriptome index and a transcript-to-gene mapping table. RSEM requires `RSEM_REFERENCE` built for the same genome/GTF version and writes both native `*.genes.results` and normalized `gene_tpm.tsv` outputs.

For the full installation and data checklist, including acceptance commands and common fixes, see `docs/INSTALL_AND_DATA.md`.

OptiType local acceptance:

```bash
source conf/tools.env.sh
optitype --version
optitype check-deps
# DNA paired FASTQ example
optitype run -i tumor_R1.fastq.gz -i tumor_R2.fastq.gz --dna -o results/optitype_sample --solver cbc --threads 8
# RNA single/paired FASTQ uses --rna instead of --dna. BAM input is also supported by OptiType.
```

The bundled OptiType references are also linked into the portable reference layout:

```text
$NEOAG_REF_BUNDLE/data/hla/optitype_reference
```

### Portable `neodata4git` Reference Bundle

For a new machine, stage large public references and licensed/local tool references outside Git using this layout. Example configuration:

```bash
export NEOAG_REF_BUNDLE=/path/to/neodata4git
source "$NEOAG_REF_BUNDLE/neodata4git.env.sh"
bash scripts/verify_reference_bundle.sh "$NEOAG_REF_BUNDLE"
NEOAG_REF_BUNDLE="$NEOAG_REF_BUNDLE" bash scripts/deploy_external_tools.sh --smoke
```

Recommended layout:

```text
neodata4git/
  data/
    ref/
      hg38/
        Homo_sapiens_assembly38.fasta
        Homo_sapiens_assembly38.fasta.fai
        Homo_sapiens_assembly38.dict
        gencode.gtf
        capture.bed                  # optional; required for WES/panel SV workflows
      ctat/current/
    facets/reference/
      common_snp.hg38.vcf.gz
      common_snp.hg38.vcf.gz.tbi
    vep/homo_sapiens/105_GRCh38/
    easyfuse/easyfuse_ref_v4/
    ascat/reference/WGS_hg38/
    sequenza/reference/
      GRCh38.primary_assembly.chr.fa
      GRCh38.primary_assembly.chr.fa.fai
      gc.wig.gz                      # optional until Sequenza is used
    hla/
      spechla_db/
      PRG_MHC_GRCh38_withIMGT/
      optitype_reference/
    lohhla/
      polysolver/
      novoalign.lic
    predictors/
      netMHCpan/
      netMHCstabpan/
      prime/
      mixMHCpred_install/
      bigmhc/models/
      DeepImmuno/
    hmf/purple_reference/            # optional until PURPLE is used
    normal/proteome/
  work/
    vep_plugins/
    nextflow_cache/
  neodata4git.env.sh
  tool_reference_manifest.tsv
```

Current known gaps in this bundle are `data/ref/hg38/capture.bed`, `data/sequenza/reference/gc.wig.gz`, and `data/hmf/purple_reference`. They are not needed for every workflow, but must be staged before WES/panel SV, Sequenza, or PURPLE production runs respectively.

For this deployment line, the canonical large-asset source root is:

```bash
/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git
```


Current synchronized assets larger than 1G are:

| Asset | Current size | Source path | Install target | Required by default |
| --- | ---: | --- | --- | --- |
| EasyFuse reference | 104G | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git/data/easyfuse/easyfuse_ref_v4` | `/root/neo/neodata4git/data/easyfuse/easyfuse_ref_v4` | No |
| HLA-LA graph | 29G | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git/data/hla/PRG_MHC_GRCh38_withIMGT` | `/root/neo/neodata4git/data/hla/PRG_MHC_GRCh38_withIMGT` | No |
| VEP cache | 16G | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git/data/vep` | `/root/neo/neodata4git/data/vep` | Yes |
| Sequenza reference | 3.2G | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git/data/sequenza/reference` | `/root/neo/neodata4git/data/sequenza/reference` | No |
| GRCh38 FASTA | 3.0G | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git/data/easyfuse/easyfuse_ref_v4/Homo_sapiens.GRCh38.dna.primary_assembly.fa` | `/root/neo/neodata4git/data/ref/hg38/Homo_sapiens_assembly38.fasta` | Yes |
| BigMHC models | 2.4G | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git/data/predictors/bigmhc/models` | `/root/neo/env_tool/tools/bigmhc/models` | Yes |
| Polysolver distribution | 1.7G | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git/data/lohhla/polysolver` | `/root/neo/licensed_tools/polysolver` | No |
| GENCODE GTF | 1.4G | `/mnt/zjl-bgi-zzb/peixunban/gl/liup/neodata4git/data/easyfuse/easyfuse_ref_v4/Homo_sapiens.GRCh38.110.gtf` | `/root/neo/neodata4git/data/ref/hg38/gencode.gtf` | Yes |

Smaller assets are still kept in `production_assets.tsv` when they are needed
for automated installation or validation, but they are intentionally omitted
from this size-focused checklist.

The installer reads `configs/assets/production_assets.tsv` and can sync large
assets from that root via:

```bash
bash .agents/skills/neoag-remote-deploy/scripts/13_install_readme_tools.sh \
  --asset-manifest configs/assets/production_assets.tsv \
  --sync-assets \
  --asset-source-host na@10.200.50.134 \
  --execute
```

`production_assets.tsv` stores paths and markers only; it does not put BigMHC
models, EasyFuse references, real VCFs, licensed tools, or patient data into Git.


For the full portable reference-bundle layout and acceptance command, use the `neodata4git` section above. Minimal configuration is:

```bash
export NEOAG_REF_BUNDLE=/path/to/neodata4git
source "$NEOAG_REF_BUNDLE/neodata4git.env.sh"
bash scripts/verify_reference_bundle.sh "$NEOAG_REF_BUNDLE"
neoag-v03 check-tools
```

Full release acceptance entry point:

```bash
# Lightweight complete check: reports missing specialized tools as warnings.
NEOAG_REF_BUNDLE=/path/to/neodata4git bash scripts/verify_all_tools_and_refs.sh --smoke

# Release-gate check: missing specialized tools or references fail the command.
NEOAG_REF_BUNDLE=/path/to/neodata4git bash scripts/verify_all_tools_and_refs.sh --strict
```

`verify_all_tools_and_refs.sh` calls the core external-tool check and reference-bundle check, then adds dedicated checks for VEP, GATK/Mutect2, NetMHCpan, NetMHCstabpan, EasyFuse/Nextflow, SpecHLA, HLA-LA, PURPLE/AMBER/COBALT, and Sequenza.

For a detailed local path inventory, see `docs/PROJECT_DATA_PATHS.md`.

## Workflow Dependency Matrix

| Workflow / command | Minimal inputs | Tools | Reference/data |
| --- | --- | --- | --- |
| Fixture demo: `neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001` | Bundled fixtures | None beyond Python package | Bundled fixtures/resources |
| Parsed pVAC results: `neoag-v03 run-v03 --outdir results/sample --sample-id SAMPLE001 --pvac data/fixtures/pvacseq_aggregated.tsv --immunogenicity-stub` | pVAC-like TSVs | None if inputs already exist | Optional normal expression/ligand tables |
| Raw intermediates: `neoag-v03 run-v03 --outdir results/sample --raw-events raw_events.tsv --raw-peptides raw_peptides.tsv` | `raw_events.tsv`, `raw_peptides.tsv` | NetMHCpan/MHCflurry outputs if provided; optional evidence tools | Optional expression, LOH, purity, CNV, normal evidence |
| Full upstream run: `neoag-v03 run-full --config conf/run.sample.private.toml --outdir results/sample` | Run config | Depends on enabled tools | Depends on enabled tools |
| Binding prediction only: `peptide-predict` | Peptide/HLA table | NetMHCpan, MHCflurry, PRIME/BigMHC/DeepImmuno as selected | HLA alleles; predictor model data |
| VEP annotation: `vep-annotate` | VCF | VEP | VEP cache, reference FASTA, plugins |
| Variant peptide extraction: `extract-variant-peptides` | VEP-annotated VCF | Python; optional VEP pre-step | Reference FASTA, optional normal proteome |
| WES SNV calling: `snv-call-wes` | Tumor/normal BAM | GATK4 | GRCh38 FASTA, gnomAD AF VCF, PoN, intervals as needed |
| WES SNV full: `snv-run-full-wes` | Somatic VCF or BAMs | GATK if BAM mode; pVAC/binding tools if enabled | GRCh38 FASTA, HLA, optional normal evidence |
| SV WGS raw build: `sv-build-raw` | SV VCF, FASTA, GTF, HLA | Python | Reference FASTA, GTF, HLA file |
| SV WES raw build: `sv-build-raw-wes` | SV VCF, FASTA, GTF, HLA, capture BED | Python | Reference FASTA, GTF, capture BED, HLA file |
| SV score: `sv-score-v03` | Raw events/peptides | NetMHCpan/MHCflurry unless `--binding-stub` | HLA alleles, optional evidence tables |
| Long-read SV wrapper | FASTQ/BAM or Sniffles2 VCF | minimap2/samtools/Sniffles2 as selected | Reference FASTA, GTF, HLA |
| Fusion discovery | FASTQ/BAM or caller outputs | STAR-Fusion, FusionCatcher, Arriba, EasyFuse as selected | CTAT/EasyFuse/fusion caller references |
| Immune escape evidence: `immune-escape` | Raw peptides, APPM/CCF/LOH evidence | Optional LOHHLA/FACETS upstream | HLA LOH, CNV, VEP/APM/JAK/B2M evidence |
| Nextflow fixture | Bundled pVAC fixture | Java/Nextflow runtime | Bundled fixtures; writable `NXF_HOME` |

## Tests

Default pytest runs fast unit tests only:

```bash
pytest -q
```

Run broader groups explicitly:

```bash
pytest -q --run-integration
pytest -q --run-benchmark
pytest -q --run-external
pytest -q --run-all
```

Marker form is also supported:

```bash
pytest -q -m unit
pytest -q -m integration --run-integration
pytest -q -m benchmark --run-benchmark
pytest -q -m external --run-external
```

This split prevents lightweight release users from accidentally running long Nextflow, benchmark, or external-tool tests with plain `pytest`.

## Installation Acceptance Commands

Run these from the project root after installation.

### Basic Package Acceptance

```bash
source conf/tools.env.sh
python -m pip install -e '.[test]'
pytest -q
neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001
```

### Tool Visibility Acceptance

```bash
source conf/tools.env.sh
neoag-v03 check-tools
bash scripts/check_tools_env.sh
```

`check-tools` may report optional tools as missing if your selected workflow does not need them. For production runs, every tool required by the selected workflow should be `OK`.

### Nextflow Acceptance

```bash
export NXF_HOME=/path/to/writable/nextflow_cache
bin/neoag-nextflow -version
bin/neoag-nextflow run workflows/main.nf \
  -w /tmp/neoag_nf_work \
  --pvac_files data/fixtures/pvacseq_aggregated.tsv \
  --outdir results/demo_nf \
  --sample_id NF_DEMO
```

Expected outputs include:

- `results/demo_nf/scoring/ranked_peptides.v03.tsv`
- `results/demo_nf/scoring/ranked_events.v03.tsv`
- `results/demo_nf/reports/evidence_report.v03.html`
- `results/demo_nf/provenance/workflow_provenance.yml`

### Reference File Acceptance

```bash
test -f "$NEOAG_REFERENCE_FASTA"
test -f "$NEOAG_VEP_CACHE/homo_sapiens/105_GRCh38/info.txt"
test -f "$NEOAG_VEP_PLUGINS/Wildtype.pm"
test -f "$NEOAG_NORMAL_PROTEOME_FASTA"
```

Run only the checks relevant to your selected workflow and configured paths.

## Common Errors And Fixes

| Symptom | Likely cause | Fix |
| --- | --- | --- |
| `neoag-v03: command not found` | Package not installed or project `bin/` not on `PATH` | Run `source conf/tools.env.sh`, then `python -m pip install -e '.[test]'`. |
| `No module named neoag_v03` | `PYTHONPATH` or editable install missing | Run `python -m pip install -e .` or execute with `PYTHONPATH=src python -m neoag_v03.cli ...`. |
| `pytest: command not found` | Test extra not installed | Run `python -m pip install -e '.[test]'`. |
| `bin/neoag-nextflow: Permission denied` | Executable bit lost during archive/migration | `find bin -maxdepth 1 -type f -exec chmod +x {} \;`. |
| `conda not found` | Miniforge/Mambaforge not installed or not initialized | Install Miniforge and open a new shell, or source its `etc/profile.d/conda.sh`. |
| `mamba: unrecognized arguments -n ...` | Pip package `mamba` was installed instead of conda-forge mamba | Do not use `pip install mamba`; use the scripts' default conda mode or install real conda/mamba. |
| `CXXABI_1.3.15 not found` | MHCflurry/scipy loaded an old system `libstdc++` | `conda install -n neoag-tools -c conda-forge 'libstdcxx-ng>=13'` and ensure conda `lib` is first in `LD_LIBRARY_PATH`. |
| `mhcflurry-downloads fetch failed` | Network/model download issue | Activate the env and rerun `mhcflurry-downloads fetch`; for offline deploys, pre-stage model data. |
| `NetMHCpan MISSING` | Licensed tarball not installed or `NETMHCPAN_HOME` wrong | Install with `bash scripts/install_netmhcpan.sh /path/to/tar.gz`, then `source conf/tools.env.sh`. |
| NetMHCpan wrapper looks for a stale conda path | Wrapper was created with an old conda base | Set `NEOAG_CONDA_BASE="$(conda info --base)"` and run `bash scripts/install_netmhcpan.sh --repair`. |
| NetMHCpan reports missing `ld-linux-x86-64.so.2` | Conda sysroot was removed from `neoag-tools` | Install `sysroot_linux-64` and `patchelf` in `neoag-tools`; the lite env file now keeps them. |
| `VEP cache not found` | Offline cache missing, `NEOAG_VEP_CACHE` points at the wrong directory, or release `105_GRCh38` is absent | Set `NEOAG_VEP_CACHE` to the cache root (not `.../105_GRCh38`), run `bash scripts/install_vep_cache.sh`, or verify `test -f "$NEOAG_VEP_CACHE/homo_sapiens/105_GRCh38/info.txt"`. |
| `vep MISSING` but VEP is installed | VEP env/wrapper not on configured path | Run `bash scripts/install_vep.sh`, source `conf/tools.env.sh`, and verify `NEOAG_VEP_BIN`. |
| `Can't locate DBI.pm` during VEP | Perl environment from another conda env polluted VEP | Use `bin/vep-neoag`, which clears conflicting Perl environment variables. |
| `No CSQ annotations` | Input VCF was not VEP-annotated | Use `run-full` with VEP configured for auto-annotation, or run `neoag-v03 vep-annotate` first. |
| `.nextflow/history.lock (Permission denied)` | Root-owned `.nextflow` metadata | Use `export NXF_HOME=/path/to/writable/cache` and run `bin/neoag-nextflow`. |
| `Downloading nextflow dependencies` hangs | First launch without cache or blocked network | Pre-populate `NXF_HOME`, use a shared cache, or allow network until download completes. |
| `Java not found` or unsupported Java | Java missing/old | Install OpenJDK 11+; verify with `java -version`. |
| `Permission denied` under `work/`, `results/`, or `tools/` | Directory owned by another user/root | Use a user-writable output/work directory or ask an administrator to fix ownership. |
| `GATK reference dictionary missing` | FASTA index/dict missing | Run `bash scripts/download_ref_hg38.sh /path/to/ref/hg38` or create `.fai`/`.dict` with samtools/picard. |
| Real-data workflow runs with fixture paths | Private run config not edited | Copy an example config to a private local config and update all paths before production. |
| Optional tool is missing but demo works | Tool not needed for fixture demo | Install only if selected workflow requires it; see the dependency matrix. |

## Release Boundary

Do not commit or package:

- `.git`, `.venv`, `.nextflow`, `.pytest_cache`
- `tools/`, `results/`, `work/`, `dist/`, `conda_packs/`
- `conf/tools.env.local.sh`
- `conf/site.config`
- `conf/private/*`
- `conf/*.private.toml`
- real patient data or sample identifiers
- licensed tool binaries
- large references such as `data/ref` and `data/vep`

Use `scripts/check_release_boundary.sh` before preparing an online release.

## Additional Documentation

- `docs/TOOLS_SETUP.md`: external tool installation details.
- `docs/PROJECT_DATA_PATHS.md`: project reference/data path inventory.
- `docs/INSTALL_AND_DATA.md`: installation and data setup guide.
- `docs/V043_CCF21.md`: CCF 2.1 notes.
- `docs/V042_P1_APPM_EXPLAINABILITY.md`: APPM explainability notes.
- `docs/V04_EVIDENCE_SAFETY_ESCAPE.md`: safety and immune-escape evidence notes.
- `RELEASE.md`: release boundary and test summary.

## Interpretation Boundary

This pipeline is for research triage and validation planning. Ranked candidates should be reviewed with assay validation, disease context, HLA typing, tumor purity, expression/protein support, safety evidence, immune-escape context, and appropriate clinical or research governance.

### NetMHCpan 4.2c container runtime

For servers where the official NetMHCpan 4.2c binary cannot run because of `tcsh` or glibc compatibility, use the Docker/Apptainer runtime documented in [docs/NETMHCPAN_CONTAINER.md](docs/NETMHCPAN_CONTAINER.md). The image contains only OS dependencies; the licensed official `tools/netMHCpan` directory is mounted at runtime.

### Priority tool containers

Docker/Apptainer runtimes for NetMHCpan, NetMHCstabpan, HLA-LA, SpecHLA, PURPLE/AMBER/COBALT, and EasyFuse are documented in [docs/PRIORITY_TOOL_CONTAINERS.md](docs/PRIORITY_TOOL_CONTAINERS.md). These images contain only runtime dependencies; licensed tools and large reference data are mounted from host paths.

## LLM-assisted Coordinator P1

This release adds an optional LLM-assisted Coordinator layer on top of the P0 Skills Pack. The default mode is dependency-free and rule-based; installing the optional `agent-llm` extra enables LiteLLM/LangGraph integration.

Plan only:

```bash
neoag-llm-agent --message "compare recommendation and NetMHCpan42 rankings" \
  --file ranked_peptides.recommendation.tsv \
  --file ranked_peptides.netmhcpan42.tsv \
  --outdir work/llm_plan --mode plan
```

Execute safe Skills:

```bash
neoag-llm-agent --message "compare recommendation and NetMHCpan42 rankings" \
  --file ranked_peptides.recommendation.tsv \
  --file ranked_peptides.netmhcpan42.tsv \
  --outdir work/llm_execute --mode execute-safe
```

Local Qwen/vLLM through LiteLLM/OpenAI-compatible API:

```bash
neoag-llm-agent --message "update patient report" \
  --file evidence_report.v04x_latest.html \
  --file ranked_peptides.recommendation.tsv \
  --file ranked_peptides.netmhcpan42.tsv \
  --outdir work/llm_report --mode execute-safe \
  --llm-provider litellm --model openai/qwen3-32b \
  --api-base http://localhost:8000/v1 --api-key-env LOCAL_VLLM_API_KEY
```

The Coordinator does not replace Project B CLI/Nextflow. It plans and calls registered Skills; high-impact operations such as HPC submission, installation, deletion, and overwrite require explicit approval.

See `docs/LLM_COORDINATOR_P1.md` and `docs/MODEL_API_AND_AGENT_FRAMEWORK_SELECTION.md`.

- [Tool inventory](docs/TOOL_INVENTORY.md): external tools, Docker images, environment variables, references, and licensing boundaries.

## Skills Taxonomy A/B/C/D

This release includes an upgraded NeoAg Skills taxonomy organized into four categories:

- **A Entry adapter skills**: `neoag-vcf`, `neoag-fusion`, `neoag-splice`, `neoag-sv-wgs`, `neoag-sv-wes`, `neoag-peptide-csv`.
- **B Public evidence analysis skills**: `neoag-hla-typing-loh`, `neoag-presentation`, `neoag-expression`, `neoag-rna-evidence`, `neoag-ccf`, `neoag-appm-escape`, `neoag-safety`, `neoag-ranking`.
- **C Review/report/design skills**: `neoag-ranking-compare`, `neoag-experiment-design`, `neoag-patient-report`, `neoag-technical-report`, `neoag-concept-explainer`.
- **D Governance/execution-control skills**: `neoag-input-qc`, `neoag-doctor`, `neoag-tool-reference-qc`, `neoag-run-demo-and-smoke`, `neoag-pipeline-full`, `neoag-release-qc`, `neoag-gateway-submit`, `neoag-hpc-runner`.

Use:

```bash
neoag-skill list
neoag-skill describe neoag-vcf
neoag-skill validate --root . --outdir work/skill_validate
neoag-skill run neoag-peptide-csv --outdir work/peptides --arg peptide_csv=peptides.tsv
```

Skills are SOP wrappers. They do not make clinical decisions, do not include patient BAM/FASTQ/VCF or large references, and high-risk execution paths remain dry-run or human-approval gated.
