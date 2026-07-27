# NetMHCpan runtime-only container

This compatibility directory intentionally does not contain or copy the
licensed NetMHCpan package. The canonical definitions are under
`containers/netmhcpan/`; the official site-local installation is mounted
read-only when `scripts/run_netmhcpan_container.sh` executes.

Build and test:

```bash
bash deploy/build_containers/netmhcpan/build.sh docker
bash deploy/build_containers/netmhcpan/test.sh
```

Do not add the official tarball, installation directory, a prebuilt image tar,
or a SIF containing NetMHCpan to Git, GHCR, Docker Hub, release archives, or
cross-site migration assets. Each user or institution must obtain and use the
software under its own applicable DTU license.
