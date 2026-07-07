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
# VEP cache root (must contain homo_sapiens/105_GRCh38/, not the 105_GRCh38 dir itself)
# export NEOAG_VEP_CACHE="/path/to/data/vep"
# export NEOAG_VEP_CACHE_VERSION="105"

# Patient / cohort data (example: chenxiaoliang)
# export CHENXIAOLIANG_DATA_ROOT="/path/to/chenxiaoliang_data"
