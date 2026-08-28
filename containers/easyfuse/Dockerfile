FROM neoag-base-bioinfo:ubuntu22.04
RUN apt-get update \
    && apt-get install -y --no-install-recommends git make gcc g++ \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /work
ENTRYPOINT ["/bin/bash", "-lc"]
