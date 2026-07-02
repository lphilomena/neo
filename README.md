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
| `NEOAG_VEP_CACHE` | VEP offline cache root, expected to contain `homo_sapiens/105_GRCh38` or equivalent. |
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
| VEP | VCF annotation and peptide extraction | Required for `vep-annotate` / auto-annotation in `run-full` | `bash scripts/install_vep.sh`; cache with `bash scripts/install_vep_cache.sh` | `NEOAG_VEP_ENV`, `NEOAG_VEP_BIN`, `NEOAG_VEP_CACHE`, `NEOAG_VEP_PLUGINS`, `NEOAG_VEP_ONLINE` | `neoag-v03 check-tools` |
| GATK4 / Mutect2 | WES/WGS SNV calling | Required if starting from BAMs | `bash scripts/install_gatk.sh` | `NEOAG_GATK_ENV` | `neoag-v03 check-tools` |
| LOHHLA | HLA LOH evidence | Optional but recommended for immune-escape evidence | `bash scripts/install_lohhla.sh`; configure Polysolver/Novoalign separately | `LOHHLA_HOME`, `POLYSOLVER_HOME`, `NOVOALIGN_LICENSE_FILE` | `neoag-v03 check-tools` |
| SpecHLA | HLA copy/LOH conversion | Optional | Install externally and provide output to `convert-spechla` | Site-specific | `neoag-v03 convert-spechla --help` |
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

Install or point to a VEP cache:

```bash
# Online cache install. This can be slow and large.
bash scripts/install_vep_cache.sh

# Or use an existing cache.
export NEOAG_VEP_CACHE=/path/to/vep_cache
export NEOAG_VEP_CACHE_VERSION=105
source conf/tools.env.sh
```

`NEOAG_VEP_PLUGINS` should point to a directory containing `Wildtype.pm` and `Frameshift.pm`. Plain VEP can run without these plugins, but this project uses them for pVACseq-compatible WT/frameshift information and more complete peptide extraction.

### Reference Data Table

Large data should live under `NEOAG_TOOLS_ROOT`, `NEOAG_SHARED_REF_DIR`, or another site-managed reference area, not inside the source checkout.

| Data/reference | Needed for | Download/setup command | Expected variable/path | Verify |
| --- | --- | --- | --- | --- |
| VEP cache, GRCh38 release 105 | Offline VEP annotation | `bash scripts/install_vep_cache.sh` | `NEOAG_VEP_CACHE=/path/to/data/vep`, contains `homo_sapiens/105_GRCh38` or equivalent | `test -d "$NEOAG_VEP_CACHE/homo_sapiens"` |
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

Recommended external bundle layout:

```text
/path/to/neoag_artifact_bundle/
  tools/
    netMHCpan/
    netMHCstabpan/
    DeepImmuno/
    prime/
    mixMHCpred_install/
  data/
    ref/hg38/
    vep/
    ref/ctat/
  work/
    vep_plugins/
```

Then configure:

```bash
export NEOAG_TOOLS_ROOT=/path/to/neoag_artifact_bundle
export NEOAG_SHARED_REF_DIR=/path/to/shared_refs
source conf/tools.env.sh
neoag-v03 check-tools
```

For a detailed local path inventory, see `docs/PROJECT_DATA_PATHS.md`.

## Workflow Dependency Matrix

| Workflow / command | Minimal inputs | Tools | Reference/data |
| --- | --- | --- | --- |
| Fixture demo: `neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001` | Bundled fixtures | None beyond Python package | Bundled fixtures/resources |
| Parsed pVAC results: `neoag-v03 run-v03 --outdir results/sample --sample-id SAMPLE001 --pvac data/fixtures/pvacseq_aggregated.tsv --immunogenicity-stub` | pVAC-like TSVs | None if inputs already exist | Optional normal expression/ligand tables |
| Raw intermediates: `neoag-v03 run-v03 --raw-events ... --raw-peptides ...` | `raw_events.tsv`, `raw_peptides.tsv` | NetMHCpan/MHCflurry outputs if provided; optional evidence tools | Optional expression, LOH, purity, CNV, normal evidence |
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
test -d "$NEOAG_VEP_CACHE/homo_sapiens"
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
| `VEP cache not found` | Offline cache missing or wrong `NEOAG_VEP_CACHE` | Run `bash scripts/install_vep_cache.sh` or set `NEOAG_VEP_CACHE` in `conf/tools.env.local.sh`. |
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
