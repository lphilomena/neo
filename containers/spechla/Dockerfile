FROM neoag-base-bioinfo:ubuntu22.04
RUN apt-get update \
    && apt-get install -y --no-install-recommends freebayes \
    && rm -rf /var/lib/apt/lists/*
RUN python3 -m pip install --no-cache-dir \
    biopython==1.81 \
    networkx==3.1 \
    numpy==1.24.4 \
    pandas==1.5.3 \
    pulp==2.7.0 \
    pyfaidx==0.8.1.1 \
    pysam==0.22.1 \
    PyVCF3==1.0.3 \
    scipy==1.10.1
WORKDIR /work
ENTRYPOINT ["/bin/bash", "-lc"]
