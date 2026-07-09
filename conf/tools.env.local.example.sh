# Site-specific overrides for conf/tools.env.sh (optional).
# Copy to conf/tools.env.local.sh and edit paths for your deployment.
#
#   cp conf/tools.env.local.example.sh conf/tools.env.local.sh
#
# conf/tools.env.local.sh is gitignored.

# Conda / licensed tools
# export NEOAG_CONDA_BASE="/path/to/miniforge3"
# export POLYSOLVER_HOME="/path/to/polysolver"
# export NOVOALIGN_LICENSE_FILE="/path/to/novoalign.lic"

# External artifact bundle (when tools/ and data/ live outside the lite release tree)
# export NEOAG_TOOLS_ROOT="/path/to/neoag_artifact_bundle"

# Reference data
# export NEOAG_NORMAL_PROTEOME_FASTA="/path/to/Homo_sapiens.GRCh38.pep.all.fa"
# export NEOAG_DBSNP_VCF="/path/to/dbsnp_chr.vcf.gz"
# export NEOAG_SHARED_REF_DIR="/path/to/shared_refs"
#
# VEP cache root (must contain homo_sapiens/105_GRCh38/, not the release dir itself)
# export NEOAG_VEP_CACHE="/path/to/data/vep"
# export NEOAG_VEP_CACHE_VERSION="105"
#
# External tool deployment (LOHHLA / FACETS / ASCAT / Arriba / PRIME):
#   bash scripts/deploy_external_tools.sh
#   bash scripts/verify_external_tools.sh

# Site-local cohort data root, if your deployment uses one
# export NEOAG_COHORT_DATA_ROOT="/path/to/cohort_data"

# RNA FASTQ to TPM (optional; used by neoag-rna-fastq-to-tpm)
# export SALMON_BIN="salmon"
# export SALMON_INDEX="$NEOAG_DATA_ROOT/data/ref/rna/salmon_index"
# export SALMON_TX2GENE="$NEOAG_DATA_ROOT/data/ref/rna/tx2gene.tsv"
# export SALMON_THREADS="8"
# export RSEM_BIN="rsem-calculate-expression"
# export RSEM_REFERENCE="$NEOAG_DATA_ROOT/data/ref/rna/rsem_reference/rsem_ref"
# export RSEM_THREADS="8"
