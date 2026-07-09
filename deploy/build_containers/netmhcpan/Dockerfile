FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive     LANG=C.UTF-8     LC_ALL=C.UTF-8

RUN apt-get update     && apt-get install -y --no-install-recommends        bash ca-certificates coreutils gawk grep        libgomp1 libstdc++6 procps tcsh zlib1g     && rm -rf /var/lib/apt/lists/*

WORKDIR /work
ENTRYPOINT ["/bin/bash", "-lc"]
