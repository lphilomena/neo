# Installation, Tool, and Data Setup

This guide is the deployment checklist for NeoAg Event Pipeline v0.4.3. It covers the base environment, Python/conda setup, Nextflow cache, external tools, reference data, workflow-specific dependencies, acceptance checks, and common fixes.

The online release is intentionally lightweight. It includes source code, workflows, tests, configuration examples, and small fixtures. It does not include licensed tools, large references, real patient data, results, work directories, conda packs, or Nextflow caches.

## 1. Base System Environment

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

Suggested system packages:

```bash
sudo apt-get update
sudo apt-get install -y bash coreutils curl wget git tar gzip unzip bzip2 xz-utils ca-certificates build-essential openjdk-17-jre-headless
```

If `sudo` is not available, ask the site administrator to provide these tools on `PATH`.

## 2. Python And Conda Environment

### 2.1 Create The Main Python Environment

Install Miniforge or Mambaforge first, then create the main tool environment:

```bash
bash scripts/setup_tools_env.sh
source conf/tools.env.sh
python -m pip install -e '.[test]'
neoag-v03 check-tools
```

For a smaller development/test environment that skips heavier optional stacks:

```bash
NEOAG_TOOLS_LITE=1 bash scripts/setup_tools_env.sh
source conf/tools.env.sh
python -m pip install -e '.[test]'
pytest -q
```

Default `pytest -q` is the release-safe quick check. `pytest -q --run-all` is a maintainer-level check that opts into integration, benchmark, external-tool, and longer workflow tests; run it only after tools, references, and Nextflow cache are installed.

Important environment variables:

| Variable | Purpose | Typical value |
| --- | --- | --- |
| `NEOAG_PROJECT_ROOT` | Source checkout path | Set by `conf/tools.env.sh` |
| `NEOAG_TOOLS_ROOT` | External artifact/tool/reference root | `/path/to/neoag_artifact_bundle` |
| `NEOAG_CONDA_BASE` | Miniforge/Mambaforge root | `/opt/miniforge3` or `/home/user/miniforge3` |
| `NEOAG_CONDA_ENV` | Main environment name | `neoag-tools` |
| `NEOAG_FORCE_CPU` | Disable GPU discovery | `1` on CPU-only nodes |

For site-specific paths, copy the local override template:

```bash
cp conf/tools.env.local.example.sh conf/tools.env.local.sh
```

Then edit `conf/tools.env.local.sh`. This file is intentionally gitignored and should not be included in online releases.

### 2.2 Verify The Lightweight CLI

```bash
neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001
# RNA VAF / junction evidence acceptance on your own raw tables
neoag-v03 build-evidence-layer --outdir results/sample --profile default \
  --raw-events results/sample/parsed/raw_events.tsv \
  --raw-peptides results/sample/parsed/raw_peptides.tsv \
  --rna-vaf results/sample/parsed/rna_vaf.tsv \
  --rna-junction results/sample/parsed/rna_junctions.tsv
# HLA LOH cross-check acceptance after converting LOHHLA and SpecHLA outputs
neoag-v03 crosscheck-hla-loh \
  --lohhla-hla-loh results/sample/tools/lohhla.hla_loh.tsv \
  --spechla-hla-loh results/sample/tools/spechla.hla_loh.tsv \
  --out results/sample/tools/hla_loh.crosscheck.tsv \
  --consensus-out results/sample/tools/hla_loh.consensus.tsv
pytest -q
```

Expected demo outputs include:

- `work/demo_v043/scoring/ranked_peptides.v03.tsv`
- `work/demo_v043/scoring/ranked_events.v03.tsv`
- `work/demo_v043/scoring/validation_plan.v03.tsv`
- `work/demo_v043/reports/evidence_report.v03.html`
- `work/demo_v043/provenance.v03.json`

## 3. Nextflow, Java, And Cache Setup

Use the project wrapper rather than calling `nextflow` directly. The wrapper prioritizes the current checkout's `bin/neoag-v03`, sets `PYTHONPATH=src`, and stores Nextflow metadata outside the repository root.

```bash
export NXF_HOME=/path/to/writable/nextflow_cache
bin/neoag-nextflow -version
```

Run the fixture workflow:

```bash
bin/neoag-nextflow run workflows/main.nf -w /tmp/neoag_nf_work --pvac_files data/fixtures/pvacseq_aggregated.tsv --outdir results/demo_nf --sample_id NF_DEMO
```

Notes:

- Online mode: the package includes the lightweight `bin/nextflow` launcher. First launch may download Nextflow runtime dependencies into `NXF_HOME`.
- Offline mode: pre-stage Java, the Nextflow runtime cache, conda/container assets, and all external references before launching. Set `NXF_HOME` to that pre-populated writable cache.
- On shared clusters, set `NXF_HOME` to a writable shared cache to avoid repeated downloads.
- Avoid launching from a root-owned `.nextflow` directory. If you see `.nextflow/history.lock (Permission denied)`, set `NXF_HOME` to a writable path and use `bin/neoag-nextflow`.

## 4. External Tool Installation Table

Tools are optional for fixture demo but required by specific real-data modes. Install only the tools needed for your workflow.

| Tool | Needed for | Required? | Install/download command | Key variables | Verify |
| --- | --- | --- | --- | --- | --- |
| pVACtools (`pvacseq`, `pvacfuse`, `pvacsplice`) | Upstream SNV/fusion/splice candidate generation | Optional unless using `run-upstream` with pVACtools | `bash scripts/setup_tools_env.sh` or Docker via `scripts/pull_docker_tools.sh` | `NEOAG_PVAC_DOCKER`, `NEOAG_PVAC_WORKDIR` | `neoag-v03 check-tools` |
| NetMHCpan | Binding/presentation prediction | Required for real binding runs unless using stub/API fallback | `bash scripts/install_netmhcpan.sh /path/to/netMHCpan-4.2*.tar.gz` | `NETMHCPAN_HOME`, `NETMHCpan`, `NEOAG_NETMHCPAN_BIN`, `NEOAG_NETMHCPAN_BACKEND` | `neoag-v03 check-tools` |
| MHCflurry | Binding/presentation prediction | Optional alternative/complement to NetMHCpan | `bash scripts/setup_tools_env.sh`; then run `mhcflurry-downloads fetch` if needed | `NEOAG_CONDA_ENV`, `NEOAG_FORCE_CPU` | `neoag-v03 check-tools` |
| NetMHCstabpan | Stability evidence | Optional | `bash scripts/install_netmhcstabpan.sh --iedb` or `bash scripts/install_netmhcstabpan.sh /path/to/netMHCstabpan*.tar.gz` | `NETMHCSTABPAN_HOME` | `neoag-v03 check-tools` |
| PRIME / MixMHCpred / BigMHC | Immunogenicity evidence | Optional | `bash scripts/install_immunogenicity_tools.sh` | `PRIME_HOME`, `MIXMHCPRED_HOME`, `BIGMHC_DIR`, `NEOAG_PRIME_JOBS` | `neoag-v03 check-tools` |
| DeepImmuno | Optional immunogenicity evidence | Optional | `bash scripts/install_deepimmuno.sh` | `DEEPIMMUNO_DIR` | `neoag-v03 check-tools` |
| VEP | Variant annotation and peptide extraction | Required for `vep-annotate` / `extract-variant-peptides` | `bash scripts/install_vep.sh`; cache with `bash scripts/install_vep_cache.sh` | `NEOAG_VEP_ENV`, `NEOAG_VEP_BIN`, `NEOAG_VEP_CACHE`, `NEOAG_VEP_PLUGINS`, `NEOAG_VEP_ONLINE` | `neoag-v03 check-tools` |
| GATK4 / Mutect2 | WES SNV calling | Required for `snv-call-wes` if not providing a somatic VCF | `bash scripts/install_gatk.sh` | `NEOAG_GATK_ENV` | `neoag-v03 check-tools` |
| LOHHLA | HLA LOH evidence | Optional but recommended for immune escape evidence | Install per LOHHLA license/docs; set path in local env | `LOHHLA_HOME`, `POLYSOLVER_HOME`, `NOVOALIGN_LICENSE_FILE` | `neoag-v03 check-tools` |
| SpecHLA | HLA copy/LOH conversion and LOHHLA cross-check | Optional | Install externally; provide output to `convert-spechla` | Site-specific | `neoag-v03 convert-spechla ...`; `neoag-v03 crosscheck-hla-loh ...` |
| FACETS | Purity/CNV/LOH evidence | Optional but recommended for CCF/escape | `bash scripts/install_ascat_pyclone.sh` and/or site R install | `FACETS_HOME`, `NEOAG_DBSNP_VCF` | `neoag-v03 check-tools` |
| ASCAT 2.5.2 | CNV/LOH evidence; legacy baseline | Optional | `bash scripts/install_ascat_pyclone.sh` or `conda env create -f conda/env.neoag-ascat.yml` | `NEOAG_ASCAT_ENV`, `ASCAT_HOME` | `neoag-v03 check-tools` |
| ASCAT 3.2.0 | CNV/LOH evidence cross-check and newer ASCAT runs | Optional | `conda env create -f conda/env.neoag-ascat-v3.yml` | `NEOAG_ASCAT_V3_ENV`, `NEOAG_ASCAT_V3_BIN` | `bin/ascat-v3 --check` |
| PyClone-VI | Clonality context | Optional | `bash scripts/install_ascat_pyclone.sh` | `NEOAG_PYCLONE_ENV`, `NEOAG_PYCLONE_BIN` | `neoag-v03 check-tools` |
| STAR-Fusion / FusionCatcher / Arriba / EasyFuse | Fusion discovery | Optional, required for corresponding fusion workflows | Install/mount externally; seed EasyFuse envs with `bash scripts/seed_easyfuse_conda_envs.sh` when applicable | `NEOAG_FUSION_ENV`, `NEOAG_STAR_FUSION_HOME`, `NEOAG_CTAT_LIB_DIR`, `NEOAG_EASYFUSE_HOME`, `NEOAG_EASYFUSE_REF` | `neoag-v03 check-tools` |
| Manta / GRIDSS / SvABA | SV discovery | Optional upstream SV callers | Install externally or via site conda/modules | `NEOAG_SV_ENV`, `NEOAG_MANTA_ENV` | `neoag-v03 check-tools` |

Licensed tools such as NetMHCpan, NetMHCstabpan, LOHHLA, and Novoalign/Polysolver components may require academic or institutional approval. Do not redistribute their binaries inside the online release.

## 5. Reference Data Download Table

Large data should live under `NEOAG_TOOLS_ROOT` or another site-managed reference area, not in the source checkout.

| Data/reference | Needed for | Download/setup command | Expected variable/path | Verify |
| --- | --- | --- | --- | --- |
| VEP cache, GRCh38 release 105 | Offline VEP annotation | `bash scripts/install_vep_cache.sh` | `NEOAG_VEP_CACHE=/path/to/data/vep`, contains `homo_sapiens/105_GRCh38` or equivalent | `test -d "$NEOAG_VEP_CACHE/homo_sapiens"` |
| GRCh38 FASTA and indices | VEP peptide extraction, GATK, SV peptide building | `bash scripts/download_ref_hg38.sh /path/to/ref/hg38` | `NEOAG_REFERENCE_FASTA=/path/to/Homo_sapiens_assembly38.fasta` | `test -f "$NEOAG_REFERENCE_FASTA"` |
| dbSNP/common SNP VCF | FACETS `snp-pileup`, some CNV workflows | Included in site reference bundle or downloaded with hg38 bundle where available | `NEOAG_DBSNP_VCF=/path/to/dbsnp_chr.vcf.gz` | `test -f "$NEOAG_DBSNP_VCF"` |
| gnomAD AF VCF and PoN | GATK Mutect2 filtering | `bash scripts/download_ref_hg38.sh /path/to/ref/hg38` | Paths inside selected run config | `test -f /path/to/af-only-gnomad.hg38.vcf.gz` |
| Ensembl protein FASTA | Peptide safety normal/reference proteome screen | Download Ensembl GRCh38 peptide FASTA manually or from site bundle | `NEOAG_NORMAL_PROTEOME_FASTA=/path/to/Homo_sapiens.GRCh38.pep.all.fa` | `test -f "$NEOAG_NORMAL_PROTEOME_FASTA"` |
| Normal expression table | Peptide safety evidence | Site-generated TSV or `resources/normal_expression.example.tsv` for fixtures | CLI argument or run config path | Header should match expected TSV schema |
| Normal HLA ligand table | Peptide safety evidence | Site-generated TSV or `resources/normal_hla_ligands.example.tsv` for fixtures | CLI argument or run config path | Header should match expected TSV schema |
| RNA allele-count / RNA VAF TSV | RNA variant support in `rna_junction_evidence.tsv` | Generate with `scripts/rna_allele_counts_pysam.py` or site RNA genotyper | CLI argument `--rna-vaf` to `build-evidence-layer` | Columns may include `event_id`, `gene`, `chrom`, `pos`, `ref`, `alt`, `rna_ref_reads`, `rna_alt_reads`, `rna_depth`, `rna_vaf` |
| CTAT genome lib | STAR-Fusion | Download per STAR-Fusion/CTAT docs or mount site bundle | `CTAT_GENOME_LIB`, `NEOAG_CTAT_LIB_DIR`, `NEOAG_SHARED_REF_DIR` | `test -d "$CTAT_GENOME_LIB"` |
| EasyFuse reference | EasyFuse workflow | Download per EasyFuse docs or mount site bundle | `NEOAG_EASYFUSE_REF`, `NEOAG_SHARED_REF_DIR` | `test -d "$NEOAG_EASYFUSE_REF"` |
| GTF annotation | SV/fusion peptide generation | Use GENCODE/Ensembl GTF matching reference FASTA | CLI argument `--gencode-gtf` | `test -f /path/to/genes.gtf` |
| Capture BED | WES SV Phase 1.5 | Use panel/exome capture BED | CLI argument `--capture-bed` | `test -f /path/to/capture.bed` |
| HLA allele file | peptide prediction and SV workflows | Site HLA typing output converted to one allele per line | CLI argument `--hla` | `head /path/to/hla.txt` |

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

## 6. Workflow Dependency Matrix

| Workflow / command | Minimal inputs | Tools | Reference/data |
| --- | --- | --- | --- |
| Fixture demo: `neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001` | Bundled fixtures | None beyond Python package | Bundled fixtures/resources |
| Parsed pVAC results: `neoag-v03 run-v03 --outdir results/sample --sample-id SAMPLE001 --pvac data/fixtures/pvacseq_aggregated.tsv --immunogenicity-stub` | pVAC-like TSVs | None if inputs already exist | Optional normal expression/ligand tables |
| Full upstream run: `neoag-v03 run-upstream --config conf/run.stub.toml --outdir results/upstream` | Run config | Depends on enabled tools | Depends on enabled tools |
| Binding prediction only: `peptide-predict` | Peptide/HLA table | NetMHCpan, MHCflurry, PRIME/BigMHC/DeepImmuno as selected | HLA alleles; predictor model data |
| VEP annotation: `vep-annotate` | VCF | VEP | VEP cache, reference FASTA, plugins |
| Variant peptide extraction: `extract-variant-peptides` | VEP-annotated VCF | Python; optional VEP pre-step | Reference FASTA, optional normal proteome |
| WES SNV calling: `snv-call-wes` | Tumor/normal BAM | GATK4 | GRCh38 FASTA, gnomAD AF VCF, PoN, intervals as needed |
| WES SNV full: `snv-run-full-wes` | Somatic VCF or BAMs | GATK if BAM mode; pVAC/binding tools if enabled | GRCh38 FASTA, HLA, optional normal evidence |
| SV WGS raw build: `sv-build-raw` | SV VCF, FASTA, GTF, HLA | Python | Reference FASTA, GTF, HLA file |
| SV WES raw build: `sv-build-raw-wes` | SV VCF, FASTA, GTF, HLA, capture BED | Python | Reference FASTA, GTF, capture BED, HLA file |
| SV score: `sv-score-v03` | raw events/peptides | NetMHCpan/MHCflurry unless `--binding-stub` | HLA alleles, optional evidence tables |
| Fusion discovery | FASTQ/BAM or caller outputs | STAR-Fusion, FusionCatcher, Arriba, EasyFuse as selected | CTAT/EasyFuse/fusion caller references |
| Immune escape evidence: `immune-escape` | raw peptides, APPM/CCF/LOH evidence | Optional LOHHLA/SpecHLA/FACETS upstream | HLA LOH consensus, CNV, VEP/APM/JAK/B2M evidence |
| HLA LOH cross-check: `crosscheck-hla-loh` | normalized LOHHLA and/or SpecHLA `hla_loh.tsv` | LOHHLA and SpecHLA outputs already converted | Optional `hla_loh.consensus.tsv` for downstream APPM/immune escape |
| Nextflow fixture: `bin/neoag-nextflow run workflows/main.nf -w /tmp/neoag_nf_work --pvac_files data/fixtures/pvacseq_aggregated.tsv --outdir results/demo_nf --sample_id NF_DEMO` | Bundled pVAC fixture | Java/Nextflow runtime | Bundled fixtures; writable `NXF_HOME` |

Use `neoag-v03 <command> --help` locally for full argument details.

## 7. Installation Acceptance Commands

Run these from the project root after installation.

### 7.1 Basic Package Acceptance

```bash
source conf/tools.env.sh
python -m pip install -e '.[test]'
pytest -q
neoag-v03 run-demo --outdir work/demo_v043 --sample-id DEMO001
# RNA VAF / junction evidence acceptance on your own raw tables
neoag-v03 build-evidence-layer --outdir results/sample --profile default \
  --raw-events results/sample/parsed/raw_events.tsv \
  --raw-peptides results/sample/parsed/raw_peptides.tsv \
  --rna-vaf results/sample/parsed/rna_vaf.tsv \
  --rna-junction results/sample/parsed/rna_junctions.tsv
# HLA LOH cross-check acceptance after converting LOHHLA and SpecHLA outputs
neoag-v03 crosscheck-hla-loh \
  --lohhla-hla-loh results/sample/tools/lohhla.hla_loh.tsv \
  --spechla-hla-loh results/sample/tools/spechla.hla_loh.tsv \
  --out results/sample/tools/hla_loh.crosscheck.tsv \
  --consensus-out results/sample/tools/hla_loh.consensus.tsv
```

### 7.2 Tool Visibility Acceptance

```bash
source conf/tools.env.sh
neoag-v03 check-tools
bash scripts/check_tools_env.sh
```

`check-tools` may report optional tools as missing if your selected workflow does not need them. For production runs, every tool required by the selected workflow should be `OK`.

### 7.3 Nextflow Acceptance

```bash
export NXF_HOME=/path/to/writable/nextflow_cache
bin/neoag-nextflow -version
bin/neoag-nextflow run workflows/main.nf -w /tmp/neoag_nf_work --pvac_files data/fixtures/pvacseq_aggregated.tsv --outdir results/demo_nf --sample_id NF_DEMO
```

Expected outputs include:

- `results/demo_nf/scoring/ranked_peptides.v03.tsv`
- `results/demo_nf/scoring/ranked_events.v03.tsv`
- `results/demo_nf/reports/evidence_report.v041.html`
- `results/demo_nf/provenance/workflow_provenance.yml`

### 7.4 Reference File Acceptance

```bash
test -f "$NEOAG_REFERENCE_FASTA"
test -d "$NEOAG_VEP_CACHE/homo_sapiens"
test -f "$NEOAG_NORMAL_PROTEOME_FASTA"
```

Run only the checks relevant to your selected workflow and configured paths.

## 8. Common Errors And Fixes

| Symptom | Likely cause | Fix |
| --- | --- | --- |
| `neoag-v03: command not found` | Package not installed or project `bin/` not on `PATH` | Run `source conf/tools.env.sh`, then `python -m pip install -e '.[test]'`. |
| `No module named neoag_v03` | `PYTHONPATH` or editable install missing | Run `python -m pip install -e .` or execute with `PYTHONPATH=src python -m neoag_v03.cli ...`. |
| `pytest: command not found` | Test extra not installed | Run `python -m pip install -e '.[test]'`. |
| `conda not found` | Miniforge/Mambaforge not installed or not initialized | Install Miniforge and open a new shell, or source its `etc/profile.d/conda.sh`. |
| `mhcflurry-downloads fetch failed` | Network/model download issue | Activate the env and rerun `mhcflurry-downloads fetch`; for offline deploys, pre-stage model data. |
| `NetMHCpan MISSING` | Licensed tarball not installed or `NETMHCPAN_HOME` wrong | Install with `bash scripts/install_netmhcpan.sh /path/to/tar.gz`, then `source conf/tools.env.sh`. |
| `VEP cache not found` | Offline cache missing or wrong `NEOAG_VEP_CACHE` | Run `bash scripts/install_vep_cache.sh` or set `NEOAG_VEP_CACHE` in `conf/tools.env.local.sh`. |
| `.nextflow/history.lock (Permission denied)` | Root-owned `.nextflow` metadata | Use `export NXF_HOME=/path/to/writable/cache` and run `bin/neoag-nextflow`. |
| `Downloading nextflow dependencies` hangs | First launch without cache or blocked network | Pre-populate `NXF_HOME`, use a shared cache, or allow network until download completes. |
| `Java not found` or unsupported Java | Java missing/old | Install OpenJDK 11+; verify with `java -version`. |
| `Permission denied` under `work/`, `results/`, or `tools/` | Directory owned by another user/root | Use a user-writable output/work directory or ask admin to fix ownership. |
| `GATK reference dictionary missing` | FASTA index/dict missing | Run `bash scripts/download_ref_hg38.sh /path/to/ref/hg38` or create `.fai`/`.dict` with samtools/picard. |
| Real-data workflow runs with fixture paths | Private run config not edited | Copy an example config to a private local config and update all paths before production. |
| Optional tool is missing but demo works | Tool not needed for fixture demo | Install only if selected workflow requires it; see the dependency matrix above. |

## 9. Release Boundary Reminder

Do not commit or package:

- `conf/tools.env.local.sh`
- `conf/site.config`
- `conf/private/*`
- real patient data, sample identifiers, patient-specific scripts, or site-local absolute paths
- licensed tool binaries
- large references
- `tools/`, `results/`, `work/`, `dist/`, `conda_packs/`, `.nextflow*`

Use `scripts/check_release_boundary.sh` before preparing an online release.
