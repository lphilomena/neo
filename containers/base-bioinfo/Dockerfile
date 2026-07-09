FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       bash ca-certificates coreutils curl file findutils gawk grep gzip less \
       libgomp1 libstdc++6 locales openssl procps sed tar unzip wget xz-utils \
       bzip2 pigz zip \
       python3 python3-pip python3-setuptools python3-venv \
       perl default-jre-headless \
       samtools bcftools bedtools bwa bowtie2 minimap2 ncbi-blast+ \
       tcsh \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work
ENTRYPOINT ["/bin/bash", "-lc"]
