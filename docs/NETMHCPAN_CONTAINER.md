# NetMHCpan 4.2c Container Runtime

This project keeps the official licensed NetMHCpan package outside the image and mounts it at runtime.
The container image only provides a compatible operating-system runtime: Ubuntu 22.04, glibc 2.35, `tcsh`, and minimal shared libraries.

## Why this is needed

The current server runs Ubuntu 20.04 with glibc 2.31. The official NetMHCpan 4.2c Linux binary requires newer glibc symbols (`GLIBC_2.33` / `GLIBC_2.34`) and the official launcher requires `/bin/tcsh`. Running through the previous conda/sysroot wrapper can trigger buffer overflow errors. The container runtime avoids both issues.

## Required host files

The official NetMHCpan directory must already exist and be configured with the correct `NMHOME` in the official launcher:

```bash
/home/na/project/neoantigen/neoag_event_pipeline/tools/netMHCpan
```

The runtime uses:

```bash
export NETMHCPAN_HOME=/home/na/project/neoantigen/neoag_event_pipeline/tools/netMHCpan
export NETMHCpan=$NETMHCPAN_HOME
```

## Build Docker image

```bash
cd /home/na/project/neoantigen/neoag_event_pipeline
./scripts/build_netmhcpan_container.sh docker
```

Default image tag:

```bash
neoag-netmhcpan:4.2c-ubuntu22.04
```

Override if needed:

```bash
NEOAG_NETMHCPAN_IMAGE=my-netmhcpan:4.2c ./scripts/build_netmhcpan_container.sh docker
```

## Run NetMHCpan through Docker

```bash
cd /home/na/project/neoantigen/neoag_event_pipeline
./scripts/run_netmhcpan_container.sh -p peptides.txt -a HLA-A02:06
```

The wrapper mounts:

- the official NetMHCpan directory as read-only;
- the project directory;
- the current working directory;
- `/mnt` when present;
- a writable temporary directory at `work/netmhcpan_tmp`.

If input or output files are outside these paths, add mounts with:

```bash
NEOAG_NETMHCPAN_EXTRA_MOUNTS=/data:/data,/scratch:/scratch \
  ./scripts/run_netmhcpan_container.sh -p /data/peptides.txt -a HLA-A02:06
```

## Verify runtime

```bash
cd /home/na/project/neoantigen/neoag_event_pipeline
./scripts/verify_netmhcpan_container.sh
```

Expected result:

```text
PASS: NetMHCpan container runtime produced HLA-A02:06 output
```

## Apptainer/Singularity

A definition file is provided for HPC systems:

```bash
./scripts/build_netmhcpan_container.sh apptainer
NEOAG_NETMHCPAN_ENGINE=apptainer ./scripts/run_netmhcpan_container.sh -p peptides.txt -a HLA-A02:06
```

Default SIF path:

```bash
containers/netmhcpan/netmhcpan-4.2c-ubuntu22.04.sif
```

## Important boundary

Do not commit or bake the official NetMHCpan package into this repository or image. It is licensed software and should remain under `tools/netMHCpan` or another site-local path mounted at runtime.

## HLA allele format

The official `tcsh` launcher can treat `*` as a shell wildcard. The wrapper normalizes allele arguments such as `HLA-A*02:06` to the NetMHCpan-safe form `HLA-A02:06` before entering the container. Direct calls to `tools/netMHCpan/netMHCpan` should still use `HLA-A02:06` format.
