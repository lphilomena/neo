# Project Data Paths

This document records the preferred host-side data layout for real deployments. The repository should remain lightweight: large references, licensed databases, BAM/FASTQ/VCF inputs, Conda/Nextflow caches, and generated outputs stay outside git and are mounted or configured through environment variables.

## Current Staged Bundle

Preferred staged bundle root on the current server:

```bash
${NEOAG_DATA_ROOT}
```

Recommended environment setup:

```bash
source /home/na/project/neoantigen/neoag_event_pipeline_v03_rc/conf/tools.env.sh
```

## Reference Data

| Category | Preferred path | Environment variable | Used for |
| --- | --- | --- | --- |
| GRCh38 FASTA | `${NEOAG_DATA_ROOT}/data/ref/hg38/Homo_sapiens_assembly38.fasta` | `NEOAG_REFERENCE_FASTA` | VEP, GATK, SV peptide building |
| GRCh38 FASTA index | `${NEOAG_DATA_ROOT}/data/ref/hg38/Homo_sapiens_assembly38.fasta.fai` | derived from FASTA | samtools/GATK/VEP checks |
| GRCh38 dictionary | `${NEOAG_DATA_ROOT}/data/ref/hg38/Homo_sapiens_assembly38.dict` | derived from FASTA | GATK |
| GTF | `${NEOAG_DATA_ROOT}/data/ref/hg38/gencode.gtf` | workflow-specific | RNA TPM, SV/fusion annotation |
| Capture BED | `${NEOAG_DATA_ROOT}/data/ref/hg38/capture.bed` | workflow-specific | WES/panel workflows |
| VEP cache | `${NEOAG_DATA_ROOT}/data/vep/homo_sapiens/105_GRCh38` | `NEOAG_VEP_CACHE`, `NEOAG_VEP_CACHE_VERSION` | Offline VEP |
| CTAT library | `${NEOAG_DATA_ROOT}/data/ctat/current` | `CTAT_GENOME_LIB` | STAR-Fusion |
| EasyFuse reference | `${NEOAG_DATA_ROOT}/data/easyfuse/current` | `NEOAG_EASYFUSE_REF` | EasyFuse |
| FACETS common SNP VCF | `${NEOAG_DATA_ROOT}/data/facets/reference/common_snp.hg38.vcf.gz` | workflow-specific | FACETS common SNP mode |
| HLA-LA graph | `${NEOAG_TOOLS_ROOT}/hla-la/env/opt/hla-la/graphs/PRG_MHC_GRCh38_withIMGT` | `HLALA_GRAPH` / `HLA_LA_GRAPH` | HLA-LA |

## Tool Install Roots

| Tool | Preferred host path | Notes |
| --- | --- | --- |
| NetMHCpan | `tools/netMHCpan` | Licensed official package; use container wrapper for runtime compatibility. |
| NetMHCstabpan | `tools/netMHCstabpan` | Licensed package or IEDB shim. |
| SpecHLA | `tools/SpecHLA` | Keep `db/` and `script/` outside release archives if large/licensed. |
| HLA-LA | `${NEOAG_TOOLS_ROOT}/hla-la/env/opt/hla-la` | External installation; graph mounted at runtime. |
| HMFTOOLS | `tools/HMFTOOLS` | PURPLE/AMBER/COBALT jars discovered under `.conda/share`. |
| EasyFuse | `tools/EasyFuse` | Code in tools root; reference bundle mounted from `neodata4git`. |

## Container Mount Boundary

Priority tool containers are documented in [PRIORITY_TOOL_CONTAINERS.md](PRIORITY_TOOL_CONTAINERS.md). They mount host-side paths instead of baking licensed tools or large data into images.

Typical mounted roots:

```bash
/home/na/project/neoantigen/neoag_event_pipeline_v03_rc
${NEOAG_DATA_ROOT}
${NEOAG_TOOLS_ROOT}
```

## Git Boundary

Do not commit:

- patient BAM/FASTQ/VCF files;
- VEP cache, HLA-LA graph, EasyFuse/CTAT references;
- official NetMHCpan/NetMHCstabpan licensed tarballs;
- Nextflow work/cache, Conda packs, Docker/Apptainer SIF images;
- generated reports and full workflow outputs.

Keep only lightweight fixtures, configuration examples, wrappers, verification scripts, and documentation in git.
