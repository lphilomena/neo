# Tool Inventory

This document is the project-level inventory for external tools used by NeoAg Event Pipeline. It is intended for release review, new-machine deployment, and agent/tool registration checks.

Large reference bundles, licensed binaries, patient data, Conda caches, Nextflow work directories, Docker image archives, and generated workflow outputs must stay outside Git.

## Runtime Modes

| Mode | How to enable | Notes |
| --- | --- | --- |
| Conda/local | default | Uses executables from `PATH`, `conf/tools.env.sh`, or run config `executables.<tool>`. |
| Docker | `NEOAG_RUNNER_MODE=docker` | Uses known Docker images from the tool registry when available; tools without a registered image fall back to local/Conda mode. |
| Apptainer/Singularity | Build from Docker images with `scripts/build_priority_tool_apptainer.sh` | Intended for HPC sites; licensed tools and references are still mounted from host paths. |

Priority container details are documented in `docs/PRIORITY_TOOL_CONTAINERS.md`.

## Tool Summary

| Category | Tool | Project registry key / executable | Main use | Install path or command | Docker image | Key variables | Required references/data |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Variant annotation | VEP | `vep` / `vep` | Variant consequence annotation and peptide-context extraction | `bash scripts/install_vep.sh`; cache with `bash scripts/install_vep_cache.sh` | `ensemblorg/ensembl-vep:release_115.2` | `NEOAG_VEP_BIN`, `NEOAG_VEP_CACHE`, `NEOAG_VEP_CACHE_VERSION`, `NEOAG_VEP_PLUGINS`, `NEOAG_REFERENCE_FASTA` | GRCh38 FASTA, VEP cache `homo_sapiens/105_GRCh38`, optional VEP plugins. |
| Variant calling | GATK4 / Mutect2 | `gatk` / `gatk` | Tumor/normal SNV/InDel calling from BAM | `bash scripts/install_gatk.sh` | `broadinstitute/gatk:4.6.2.0` | `NEOAG_GATK_ENV`, run-config GATK resource paths | GRCh38 FASTA bundle, gnomAD AF VCF, Panel of Normals, optional intervals/BED. |
| Neoantigen upstream | pVACseq | `pvacseq` / `pvacseq` | SNV/InDel upstream neoantigen prediction | `bash scripts/setup_tools_env.sh` | `griffithlab/pvactools:6.1.1` | `NEOAG_PVAC_DOCKER`, `NEOAG_PVAC_WORKDIR` | VCF, HLA alleles, reference annotations required by pVACtools mode. |
| Neoantigen upstream | pVACfuse | `pvacfuse` / `pvacfuse` | Fusion neoantigen prediction | `bash scripts/setup_tools_env.sh` | `griffithlab/pvactools:6.1.1` | `NEOAG_PVAC_DOCKER`, `NEOAG_PVAC_WORKDIR` | Fusion calls and HLA alleles. |
| Neoantigen upstream | pVACsplice | `pvacsplice` / `pvacsplice` | Splice-junction neoantigen prediction | `bash scripts/setup_tools_env.sh` | `griffithlab/pvactools:6.1.1` | `NEOAG_PVAC_DOCKER`, `NEOAG_PVAC_WORKDIR` | VCF, RegTools junctions, HLA alleles. |
| Presentation | NetMHCpan 4.2 | `netmhcpan` / `netMHCpan` | Primary MHC-I binding and EL prediction | `bash scripts/install_netmhcpan.sh /path/to/netMHCpan-4.2*.tar.gz` | `neoag-netmhcpan:4.2c-ubuntu22.04` | `NETMHCPAN_HOME`, `NETMHCpan`, `NEOAG_NETMHCPAN_BIN`, `NEOAG_NETMHCPAN_BACKEND` | Licensed DTU NetMHCpan directory including `data/`. Do not redistribute. |
| Presentation | MHCflurry | `mhcflurry` / `mhcflurry-predict` | Alternative/complementary MHC-I prediction | `bash scripts/setup_tools_env.sh`; then `mhcflurry-downloads fetch` if needed | `neoag-base-bioinfo:ubuntu22.04` | `MHCFLURRY_DOWNLOADS_DIR`, `TF_USE_LEGACY_KERAS`, `NEOAG_FORCE_CPU` | MHCflurry downloaded models. |
| Stability | NetMHCstabpan | `netmhcstabpan` / `netMHCstabpan` | Optional pMHC stability evidence | `bash scripts/install_netmhcstabpan.sh --iedb` or licensed tarball install | `neoag-netmhcstabpan:1.0-ubuntu22.04` | `NETMHCSTABPAN_HOME` | Licensed DTU package or IEDB-compatible shim. Do not redistribute licensed package. |
| Immunogenicity | PRIME / MixMHCpred | `prime` / `PRIME` | MHC-I immunogenicity and presentation evidence | `bash scripts/install_immunogenicity_tools.sh` | `neoag-base-bioinfo:ubuntu22.04` | `PRIME_HOME`, `MIXMHCPRED_HOME`, `NEOAG_PRIME_JOBS` | PRIME/MixMHCpred model and executable directories. |
| Immunogenicity | BigMHC_IM | `bigmhc_im` / `bigmhc_predict` | Neoepitope immunogenicity evidence | `bash scripts/install_immunogenicity_tools.sh` | `neoag-base-bioinfo:ubuntu22.04` | `BIGMHC_DIR`, optional `BIGMHC_PYTHON` | BigMHC repository/model directory. |
| Immunogenicity | DeepImmuno-CNN | `deepimmuno` / `deepimmuno-cnn.py` | Optional 9/10-mer immunogenicity evidence | `bash scripts/install_deepimmuno.sh` | `neoag-base-bioinfo:ubuntu22.04` | `DEEPIMMUNO_DIR` | DeepImmuno data/model files. |
| HLA LOH | LOHHLA | `lohhla` / `LOHHLA` | HLA allele-specific LOH evidence | `bash scripts/install_lohhla.sh`; configure Polysolver/Novoalign separately | `quay.io/biocontainers/lohhla:20171108--hdfd78af_3` | `LOHHLA_HOME`, `POLYSOLVER_HOME`, `NOVOALIGN_LICENSE_FILE` | Polysolver distribution, Novoalign license, HLA/purity inputs. |
| HLA typing/LOH | SpecHLA | `spechla` / site wrapper | HLA typing, HLA-region analysis, HLA LOH cross-check | Site install or container wrapper | `neoag-spechla:ubuntu22.04` | `SPECHLA_HOME`, `NEOAG_SPECHLA_HOME`, `SPECHLA_ENV` | SpecHLA database/reference files; BAM or extracted HLA reads. |
| HLA typing | HLA-LA | `hla-la` / site executable | Graph-based HLA typing, useful for long reads | Site install; graph outside Git | `neoag-hla-la:ubuntu22.04` | `HLALA_HOME`, `HLALA_BIN`, `HLALA_GRAPH`, `HLA_LA_GRAPH` | `PRG_MHC_GRCh38_withIMGT` graph directory and aligned BAM. |
| HLA typing | OptiType | `optitype` / `optitype` | HLA-A/B/C typing from DNA/RNA FASTQ or BAM | `bash scripts/install_optitype.sh` | `fred2/optitype:release-v1.3.1` | `OPTITYPE_ENV`, `OPTITYPE_BIN`, `OPTITYPE_REFERENCE` | OptiType bundled HLA-I reference. |
| Purity/CNV | FACETS | `facets` / `runFACETS.R` | Tumor purity, ploidy, CNV, LOH evidence | `bash scripts/install_facets.sh` | `neoag-purple-suite:ubuntu22.04` | `FACETS_HOME`, `NEOAG_DBSNP_VCF` | Common SNP/dbSNP VCF and index; tumor/normal BAM-derived pileup. |
| Purity/CNV | ASCAT 2.5.2 | `ascat` / `ascat.R` | Allele-specific purity/CNV baseline | `bash scripts/install_ascat_pyclone.sh` | `quay.io/biocontainers/ascat:2.5.2--r40hdfd78af_3` | `NEOAG_ASCAT_ENV`, `ASCAT_HOME`, run-specific ASCAT loci/alleles variables | hg38 loci/alleles resources, GC/RT correction files where applicable. |
| Purity/CNV | ASCAT 3.2.0 | site wrapper / `ascat-v3` | Newer ASCAT cross-check and prepareHTS workflows | `conda env create -f conda/env.neoag-ascat-v3.yml` | site/container dependent | `NEOAG_ASCAT_V3_ENV`, `NEOAG_ASCAT_V3_BIN` | hg38 loci/alleles resources, GC/RT correction files. |
| Purity/CNV | PURPLE / AMBER / COBALT | site wrappers/JARs | Purity, ploidy, CNV, LOH, QC cross-check | See `docs/TOOLS_SETUP.md` | `neoag-purple-suite:ubuntu22.04` | `HMFTOOLS_HOME`, `PURPLE_JAR`, `AMBER_JAR`, `COBALT_JAR` | HMF GRCh38 reference bundle; tumor/normal BAM; optional somatic VCF. |
| Purity/CNV | Sequenza | site R/Python tools | Independent purity/ploidy/CNV estimate | Conda/R environment; see Sequenza scripts | none registered | Sequenza script arguments; `sequenza-utils` on `PATH` | GRCh38 FASTA, GC wiggle/file, tumor/normal BAM or existing seqz blocks. |
| Clonality | PyClone-VI | `pyclone` / `pyclone` | Clonal cluster and CCF context | `bash scripts/install_ascat_pyclone.sh` | `neoag-base-bioinfo:ubuntu22.04` | `NEOAG_PYCLONE_ENV`, `NEOAG_PYCLONE_BIN` | Mutation cellular prevalence / copy-number input. |
| RNA fusion | STAR-Fusion | `star_fusion` / `star-fusion-neoag` | RNA fusion discovery | `bash scripts/install_fusion_tools.sh` or external install | `quay.io/biocontainers/star-fusion:1.15.1--hdfd78af_1` | `NEOAG_STAR_FUSION_HOME`, `CTAT_GENOME_LIB` | CTAT genome library. |
| RNA fusion | FusionCatcher | `fusioncatcher` / `fusioncatcher-neoag` | RNA fusion discovery | `bash scripts/install_fusion_tools.sh` or external install | `quay.io/biocontainers/fusioncatcher:1.33b--hdfd78af_0` | `NEOAG_FUSIONCATCHER_HOME` | FusionCatcher reference data. |
| RNA fusion | Arriba | `arriba` / `arriba` | RNA fusion discovery | Conda/site install or external module | `uhrigs/arriba:2.5.1` | `NEOAG_FUSION_ENV`, `GTF`, STAR/reference paths | GRCh38 FASTA, GTF, STAR index, optional blacklist resources. |
| RNA fusion | EasyFuse | `easyfuse` / `easyfuse-neoag` | Fusion metacaller and evidence cross-check | External install or prebuilt Nextflow/Conda cache | `tronbioinformatics/easyfuse:1.3.7` | `NEOAG_EASYFUSE_HOME`, `NEOAG_EASYFUSE_REF`, `NEOAG_CONDA_BASE` | EasyFuse reference bundle and Nextflow/Conda cache as needed. |
| RNA expression | Salmon | script-level | RNA FASTQ to gene/transcript TPM | `scripts/run_salmon_fastq_to_tpm.sh` | none registered | `SALMON_BIN`, `SALMON_INDEX`, `SALMON_TX2GENE`, `SALMON_THREADS` | Salmon index and tx2gene mapping. |
| RNA expression | RSEM | script-level | RNA FASTQ to TPM with RSEM | `scripts/run_rsem_fastq_to_tpm.sh` | none registered | `RSEM_BIN`, `RSEM_REFERENCE`, `RSEM_THREADS` | RSEM reference prefix. |
| SV discovery | Manta | site tool | Structural variant discovery | Site install/module | `quay.io/biocontainers/manta:1.6.0--h9ee0642_3` | `NEOAG_MANTA_ENV`, run config paths | GRCh38 FASTA and BAM inputs. |
| SV discovery | GRIDSS2 | site tool | Structural variant discovery | Site install/module | `gridss/gridss:2.13.2` | run config paths | GRCh38 FASTA and BAM inputs. |
| SV discovery | SvABA | site tool | Structural variant discovery | Site install/module | `quay.io/biocontainers/svaba:1.2.0--h69ac913_1` | run config paths | GRCh38 FASTA and BAM inputs. |
| SV discovery | DELLY | site tool | Structural variant discovery | Site install/module | `dellytools/delly:v2.3.0` | run config paths | GRCh38 FASTA and BAM inputs. |
| Long-read SV | Sniffles2 / minimap2 / samtools | site tools | Long-read alignment and SV calling | Site install/module | varies | run config paths | Long-read FASTQ/BAM, GRCh38 FASTA, GTF, HLA file. |
| Basic utilities | samtools / tabix | `samtools`, `tabix` | BAM/VCF indexing and utility operations | Conda/site install | `quay.io/biocontainers/samtools:1.23.1--ha83d96e_0`, `quay.io/biocontainers/tabix:1.11--hdfd78af_0` | `PATH` | BAM/VCF inputs. |
| Workflow engine | Nextflow / Java | site executable | Fusion workflows and optional workflow orchestration | Install Java 11+ and Nextflow; configure cache | none registered | `NXF_HOME`, `NXF_CONDA_CACHEDIR`, `JAVA_HOME` | Writable cache/work directories; prebuilt cache for offline use. |

## Agent-Level Skills

The LLM agent uses higher-level skills that may call one or more tools. These skills are registered under `.agents/config/skills_registry.json` and `.agents/config/tools_registry.json`.

| Skill | Role | Typical tools |
| --- | --- | --- |
| `neoag-input-qc` | Validate input files and recommend workflow mode | Python, file checks |
| `neoag-sliding-run` | Run VCF + HLA sliding-window neoantigen workflow | VEP, NetMHCpan/MHCflurry, scoring pipeline |
| `neoag-result-inspector` | Summarize existing workflow results | Python parsers |
| `neoag-purity-cnv-review` | Compare purity/CNV calls | FACETS, PURPLE, Sequenza, ASCAT |
| `neoag-hla-typing-compare` | Compare HLA typing calls | OptiType, SpecHLA, HLA-LA |
| `neoag-fusion-rna-run` | Fusion RNA workflow coordination | STAR-Fusion, FusionCatcher, Arriba, EasyFuse |
| `neoag-rna-fastq-to-tpm` | RNA FASTQ to TPM workflow | Salmon or RSEM |
| `neoag-tool-and-reference-qc` | Tool/reference acceptance checks | `verify_all_tools_and_refs.sh`, tool-specific checks |

## Verification Commands

```bash
# Fast Python/unit-test check
pytest -q

# Core tool availability check
neoag check-tools

# Full tool/reference acceptance, warnings allowed by default
NEOAG_REF_BUNDLE=/path/to/neodata4git bash scripts/verify_all_tools_and_refs.sh

# Strict release-gate check
NEOAG_REF_BUNDLE=/path/to/neodata4git bash scripts/verify_all_tools_and_refs.sh --strict

# Priority container smoke checks
bash scripts/verify_priority_tool_containers.sh
```

## Licensing Boundary

Do not commit or redistribute the following inside the public release branch unless the license explicitly allows it:

- DTU NetMHCpan and NetMHCstabpan official packages.
- Novoalign/Polysolver licensed components.
- HLA-LA graph bundles such as `PRG_MHC_GRCh38_withIMGT`.
- Large VEP/CTAT/EasyFuse/HMF/ASCAT/FACETS reference bundles.
- Patient BAM/FASTQ/VCF files or generated analysis outputs.
- Docker image tarballs or Apptainer `.sif` files unless they contain only redistributable software and are intentionally released outside Git.
