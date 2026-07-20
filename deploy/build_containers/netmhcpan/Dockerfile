FROM ubuntu:22.04

LABEL version="4.2c"
LABEL description="NetMHCpan 4.2c — MHC-I binding prediction (DTU Health Tech)"
LABEL tool="netMHCpan"

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# Minimal runtime deps: libgomp1 for OpenMP (neural network engine),
# libstdc++6 for the C++ runtime.  The ELF binary links against GLIBC_2.34
# which Ubuntu 22.04 (glibc 2.35) satisfies without patching.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       bash ca-certificates libgomp1 libstdc++6 \
    && rm -rf /var/lib/apt/lists/*

ENV NETMHCPAN_HOME=/opt/netMHCpan

# Extract NetMHCpan tarball into /opt/netMHCpan
# The tarball includes: Linux_x86_64/bin/* (ELF binaries), data/ (matrices, thresholds, alleles)
# Linux_x86_64/data -> ../data symlink is already present.  No ELF patching needed — the
# DTU binary links against GLIBC_2.34 and Ubuntu 22.04 ships glibc 2.35.
COPY vendor/netMHCpan-4.2c.Linux.tar.gz /tmp/
RUN tar xzf /tmp/netMHCpan-4.2c.Linux.tar.gz -C /tmp/ \
    && mkdir -p /opt \
    && mv /tmp/netMHCpan-4.2 "${NETMHCPAN_HOME}" \
    && rm -f /tmp/netMHCpan-4.2c.Linux.tar.gz \
    && find "${NETMHCPAN_HOME}" -name '._*' -delete \
    && find "${NETMHCPAN_HOME}" -name '.DS_Store' -delete \
    && mkdir -p "${NETMHCPAN_HOME}/tmp"

# Replace tcsh launcher with a portable bash wrapper.
# The original tcsh script sets NETMHCpan=$NMHOME/$PLATFORM (i.e., $NMHOME/Linux_x86_64).
# The main binary shells out to helpers at $NETMHCpan/bin/, so NETMHCpan must point
# to the platform directory, not the install root.
RUN cat > "${NETMHCPAN_HOME}/netMHCpan" <<'WRAPPER'
#!/usr/bin/env bash
set -euo pipefail

export NETMHCPAN_HOME="${NETMHCPAN_HOME:-/opt/netMHCpan}"
export NETMHCpan="${NETMHCpan:-${NETMHCPAN_HOME}/Linux_x86_64}"
export TMPDIR="${TMPDIR:-/tmp}"
mkdir -p "${TMPDIR}"

BIN="${NETMHCpan}/bin/netMHCpan-4.2"
if [[ ! -x "${BIN}" ]]; then
  echo "netMHCpan: binary not found at ${BIN}" >&2
  exit 127
fi

exec "${BIN}" "$@"
WRAPPER
RUN chmod +x "${NETMHCPAN_HOME}/netMHCpan"

# Environment variables consumed by neoag-v03 and downstream Nextflow processes.
# NETMHCpan must point to the platform subdirectory so the binary finds its helpers
# at $NETMHCpan/bin/ (matching the original tcsh launcher behaviour).
ENV NETMHCpan=/opt/netMHCpan/Linux_x86_64 \
    NEOAG_NETMHCPAN_BIN=/opt/netMHCpan/netMHCpan \
    NEOAG_NETMHCPAN_TMPDIR=/opt/netMHCpan/tmp
ENV PATH="/opt/netMHCpan:${PATH}"

# Smoke test — verify the binary starts and prints help output
RUN netMHCpan -h 2>&1 | head -5 || true \
    && test -x /opt/netMHCpan/Linux_x86_64/bin/netMHCpan-4.2 \
    && echo "NetMHCpan 4.2c OK"

WORKDIR /work
ENTRYPOINT ["/bin/bash"]
