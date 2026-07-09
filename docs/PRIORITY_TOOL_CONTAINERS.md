# Priority Tool Docker/Apptainer Runtimes

This project provides container runtimes for the external tools that are hardest to reproduce on a new machine:

- NetMHCpan 4.2c
- NetMHCstabpan
- HLA-LA
- SpecHLA
- PURPLE / AMBER / COBALT
- EasyFuse

The containers hold the operating-system runtime only. Licensed tools, Java jars, Conda caches, and large reference databases remain outside the image and are mounted at runtime.

## Build Docker images

```bash
cd /home/na/project/neoantigen/neoag_event_pipeline_v03_rc
./scripts/build_priority_tool_containers.sh all
```

Build one image only:

```bash
./scripts/build_priority_tool_containers.sh spechla
./scripts/build_priority_tool_containers.sh purple-suite
```

Docker images built by default:

| Tool | Image |
| --- | --- |
| shared bioinfo base | `neoag-base-bioinfo:ubuntu22.04` |
| NetMHCpan | `neoag-netmhcpan:4.2c-ubuntu22.04` |
| NetMHCstabpan | `neoag-netmhcstabpan:1.0-ubuntu22.04` |
| HLA-LA | `neoag-hla-la:ubuntu22.04` |
| SpecHLA | `neoag-spechla:ubuntu22.04` |
| PURPLE / AMBER / COBALT | `neoag-purple-suite:ubuntu22.04` |
| EasyFuse | `neoag-easyfuse:ubuntu22.04` |

## Verify all priority containers

```bash
./scripts/verify_priority_tool_containers.sh
```

This verifies that each image starts and that configured host-side tool/reference paths are visible. It does not run full production workflows.

## Run wrappers

| Tool | Wrapper |
| --- | --- |
| NetMHCpan | `scripts/run_netmhcpan_container.sh` |
| NetMHCstabpan | `scripts/run_netmhcstabpan_container.sh` |
| HLA-LA | `scripts/run_hla_la_container.sh` |
| SpecHLA | `scripts/run_spechla_container.sh` |
| PURPLE / AMBER / COBALT | `scripts/run_purple_suite_container.sh` |
| EasyFuse | `scripts/run_easyfuse_container.sh` |

Examples:

```bash
./scripts/run_netmhcpan_container.sh -p peptides.txt -a HLA-A*02:06
./scripts/run_spechla_container.sh -n SAMPLE -1 hla_R1.fq.gz -2 hla_R2.fq.gz -o outdir
./scripts/run_purple_suite_container.sh amber -sample SAMPLE [amber args]
./scripts/run_purple_suite_container.sh cobalt -sample SAMPLE [cobalt args]
./scripts/run_purple_suite_container.sh purple -sample SAMPLE [purple args]
```

## Host-side paths

The wrappers source `conf/tools.env.sh` when present. Important variables:

| Variable | Meaning |
| --- | --- |
| `NETMHCPAN_HOME` / `NETMHCpan` | Official NetMHCpan installation directory. |
| `NETMHCSTABPAN_HOME` | NetMHCstabpan installation or IEDB shim directory. |
| `HLALA_HOME`, `HLALA_BIN`, `HLALA_GRAPH` | HLA-LA installation, executable, and PRG graph. |
| `SPECHLA_HOME` | SpecHLA source/database directory. |
| `HMFTOOLS_HOME`, `PURPLE_JAR`, `AMBER_JAR`, `COBALT_JAR` | HMF tools jars. |
| `NEOAG_EASYFUSE_HOME`, `NEOAG_EASYFUSE_REF`, `NEOAG_CONDA_BASE` | EasyFuse code, reference bundle, and Conda/Mamba root. |

## Apptainer/Singularity

On HPC systems with Apptainer/Singularity, first build Docker images, then convert them to SIF files:

```bash
./scripts/build_priority_tool_containers.sh all
./scripts/build_priority_tool_apptainer.sh all
```

Default SIF output directory:

```bash
containers/sif/
```

The SIF directory is ignored by git.

## Licensing and data boundary

Do not commit licensed software or large reference data into git. In particular, keep NetMHCpan/NetMHCstabpan official packages, HLA-LA graphs, SpecHLA databases, PURPLE references, EasyFuse references, BAM/FASTQ/VCF files, and Nextflow/Conda caches outside the repository or under ignored runtime directories.
