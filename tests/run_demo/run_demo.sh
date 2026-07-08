#DATA_DIR=/public/home/wuwen/00-guoling/00-data/chenxiaoliang_data
WORK_DIR=/home/na/project/neoantigen/test/res
export PATH=$PATH:/home/na/project/neo/bin

source /home/na/project/neo/conf/tools.env.sh

# tools.env.sh puts neoag-tools/bin before neoag-vep/bin on PATH,
# so /usr/bin/env perl resolves to neoag-tools Perl which lacks DBI.pm.
# Prepend neoag-vep/bin so VEP finds a Perl with the DBI module.
export PATH="/home/na/miniforge3/envs/neoag-vep/bin:$PATH"

# VEP offline cache is available on the shared mount;
# tools.env.sh auto-detects it via NEOAG_VEP_CACHE.

neoag-v03 run-full \
  --config run.sliding.private.toml \
  --outdir $WORK_DIR/SAMPLE001_sliding

