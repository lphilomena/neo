# Source before neoag upstream runs: source conf/tools.env.sh
# Optional site overlay: copy conf/tools.env.local.example.sh -> conf/tools.env.local.sh

_NEOAG_TOOLS_ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export NEOAG_PROJECT_ROOT="$(cd "${_NEOAG_TOOLS_ENV_DIR}/.." && pwd)"
export NEOAG_TOOLS_ROOT="${NEOAG_TOOLS_ROOT:-${NEOAG_PROJECT_ROOT}}"
unset _NEOAG_TOOLS_ENV_DIR

export NEOAG_CONDA_ENV="neoag-tools"
export NEOAG_CONDA_BASE="${NEOAG_CONDA_BASE:-$(conda info --base 2>/dev/null || echo ${HOME}/miniconda3)}"

# pVACtools + MHCflurry live in neoag-tools (do not mix neoag-vep onto PATH — Python 3.7)
export LD_LIBRARY_PATH="${NEOAG_CONDA_BASE}/envs/neoag-tools/lib:${LD_LIBRARY_PATH:-}"

# TensorFlow GPU: pip nvidia-*-cu12 libs (see scripts/install_cuda_pip_libs.sh)
NEOAG_NVIDIA_SITE="${NEOAG_CONDA_BASE}/envs/${NEOAG_CONDA_ENV}/lib/python3.11/site-packages/nvidia"
if [[ -d "${NEOAG_NVIDIA_SITE}" ]]; then
  for _neoag_nv_lib in "${NEOAG_NVIDIA_SITE}"/*/lib; do
    if [[ -d "${_neoag_nv_lib}" ]]; then
      export LD_LIBRARY_PATH="${_neoag_nv_lib}:${LD_LIBRARY_PATH}"
    fi
  done
  unset _neoag_nv_lib
fi

# MHCflurry + TensorFlow 2.21: legacy tf.keras via tf_keras (pip install tf_keras)
export TF_USE_LEGACY_KERAS=1
# GPU default on when nvidia libs present; force CPU: export NEOAG_FORCE_CPU=1
if [[ "${NEOAG_FORCE_CPU:-}" == "1" ]]; then
  export CUDA_VISIBLE_DEVICES=""
elif [[ -z "${CUDA_VISIBLE_DEVICES+x}" ]] && [[ ! -d "${NEOAG_NVIDIA_SITE}" ]]; then
  export CUDA_VISIBLE_DEVICES=""
fi

# VEP binary only (separate env neoag-vep)
export NEOAG_VEP_ENV="neoag-vep"
if [[ -x "${NEOAG_TOOLS_ROOT}/bin/vep-neoag" ]]; then
  export NEOAG_VEP_BIN="${NEOAG_TOOLS_ROOT}/bin/vep-neoag"
elif [[ -x "${NEOAG_CONDA_BASE}/envs/${NEOAG_VEP_ENV}/bin/vep" ]]; then
  export NEOAG_VEP_BIN="${NEOAG_CONDA_BASE}/envs/${NEOAG_VEP_ENV}/bin/vep"
  export PATH="${NEOAG_CONDA_BASE}/envs/${NEOAG_VEP_ENV}/bin:${PATH}"
else
  export NEOAG_VEP_BIN="${NEOAG_TOOLS_ROOT}/bin/vep-neoag"
fi
# VEP cache root (must contain homo_sapiens/<version>_GRCh38/, not the release dir itself).
export NEOAG_VEP_CACHE="${NEOAG_TOOLS_ROOT}/data/vep"
if [[ ! -d "${NEOAG_VEP_CACHE}/homo_sapiens" && -d "${NEOAG_TOOLS_ROOT}/../neoag_event_pipeline_artifact_quarantine_20260622_091158/data/vep/homo_sapiens" ]]; then
  export NEOAG_VEP_CACHE="${NEOAG_TOOLS_ROOT}/../neoag_event_pipeline_artifact_quarantine_20260622_091158/data/vep"
fi
export NEOAG_VEP_CACHE_VERSION="105"
export NEOAG_VEP_PLUGINS="${NEOAG_TOOLS_ROOT}/work/vep_plugins"
if [[ ! -f "${NEOAG_VEP_PLUGINS}/Wildtype.pm" && -f "${NEOAG_TOOLS_ROOT}/../neoag_event_pipeline_artifact_quarantine_20260622_091158/work/vep_plugins/Wildtype.pm" ]]; then
  export NEOAG_VEP_PLUGINS="${NEOAG_TOOLS_ROOT}/../neoag_event_pipeline_artifact_quarantine_20260622_091158/work/vep_plugins"
fi
export NEOAG_REFERENCE_FASTA="${NEOAG_TOOLS_ROOT}/data/ref/hg38/Homo_sapiens_assembly38.fasta"

# Ensembl GRCh38 reference proteome (optional; set in conf/tools.env.local.sh)
export NEOAG_NORMAL_PROTEOME_FASTA="${NEOAG_NORMAL_PROTEOME_FASTA:-}"

# DeepImmuno-CNN (optional immunogenicity; 9/10-mer peptide–HLA pairs)
export DEEPIMMUNO_DIR="${NEOAG_TOOLS_ROOT}/tools/DeepImmuno"

NEOAG_TOOL_QUARANTINE="${NEOAG_TOOL_QUARANTINE:-}"

# BigMHC_IM (repo ~5GB incl. models under models/bat*/im/)
export BIGMHC_DIR="${NEOAG_TOOLS_ROOT}/tools/bigmhc"

# Use the neoag-tools Python for BigMHC because other tool envs on PATH may not include torch.
export BIGMHC_PYTHON="${BIGMHC_PYTHON:-${NEOAG_CONDA_BASE}/envs/neoag-tools/bin/python}"
if [[ ! -x "${BIGMHC_PYTHON}" ]]; then
  if [[ -x "${NEOAG_CONDA_BASE}/envs/neoag-tools/bin/python" ]]; then
    export BIGMHC_PYTHON="${NEOAG_CONDA_BASE}/envs/neoag-tools/bin/python"
  elif [[ -x "${HOME}/miniforge3/envs/neoag-tools/bin/python" ]]; then
    export BIGMHC_PYTHON="${HOME}/miniforge3/envs/neoag-tools/bin/python"
  fi
fi
if [[ ! -f "${BIGMHC_DIR}/src/predict.py" && -n "${NEOAG_TOOL_QUARANTINE}" && -f "${NEOAG_TOOL_QUARANTINE}/bigmhc/src/predict.py" ]]; then
  export BIGMHC_DIR="${NEOAG_TOOL_QUARANTINE}/bigmhc"
fi

# PRIME + MixMHCpred (immunogenicity)
export PRIME_HOME="${NEOAG_TOOLS_ROOT}/tools/prime"
export MIXMHCPRED_HOME="${NEOAG_TOOLS_ROOT}/tools/mixMHCpred_install"
if [[ ! -x "${PRIME_HOME}/PRIME" && -n "${NEOAG_TOOL_QUARANTINE}" && -x "${NEOAG_TOOL_QUARANTINE}/prime/PRIME" ]]; then
  export PRIME_HOME="${NEOAG_TOOL_QUARANTINE}/prime"
fi
if [[ ! -x "${MIXMHCPRED_HOME}/MixMHCpred" && -n "${NEOAG_TOOL_QUARANTINE}" && -x "${NEOAG_TOOL_QUARANTINE}/mixMHCpred_install/MixMHCpred" ]]; then
  export MIXMHCPRED_HOME="${NEOAG_TOOL_QUARANTINE}/mixMHCpred_install"
fi
export NEOAG_PRIME_BIN="${PRIME_HOME}/PRIME"
# Parallel PRIME batches (each batch = 1 CPU-bound PRIME.x.bin process)
export NEOAG_PRIME_JOBS="${NEOAG_PRIME_JOBS:-4}"
export MIXMHCPRED_BIN="${MIXMHCPRED_HOME}/MixMHCpred"
export PATH="${PRIME_HOME}:${MIXMHCPRED_HOME}:${PATH}"
export NETMHCPAN_HOME="${NEOAG_NETMHCPAN_HOME:-${NEOAG_TOOLS_ROOT}/tools/netMHCpan}"
if [[ ! -x "${NETMHCPAN_HOME}/netMHCpan" && -n "${NEOAG_TOOL_QUARANTINE}" && -x "${NEOAG_TOOL_QUARANTINE}/netMHCpan/netMHCpan" ]]; then
  export NETMHCPAN_HOME="${NEOAG_TOOL_QUARANTINE}/netMHCpan"
fi
# DTU NetMHCpan reads $NETMHCpan (not NETMHCPAN_HOME) for data/ and bin/ paths.
export NETMHCpan="${NETMHCPAN_HOME}"
export NEOAG_NETMHCPAN_BIN="${NETMHCPAN_HOME}/netMHCpan"
mkdir -p "${NETMHCPAN_HOME}/tmp"
export NEOAG_NETMHCPAN_TMPDIR="${NEOAG_NETMHCPAN_TMPDIR:-${NETMHCPAN_HOME}/tmp}"
export PATH="${NETMHCPAN_HOME}:${PATH}"
# local = project NetMHCpan 4.2; iedb = force IEDB API only (no local attempt)
export NEOAG_NETMHCPAN_BACKEND="${NEOAG_NETMHCPAN_BACKEND:-local}"

# NetMHCstabpan — install: bash scripts/install_netmhcstabpan.sh [--iedb]
export NETMHCSTABPAN_HOME="${NEOAG_TOOLS_ROOT}/tools/netMHCstabpan"
if [[ ! -x "${NETMHCSTABPAN_HOME}/netMHCstabpan" && -n "${NEOAG_TOOL_QUARANTINE}" && -x "${NEOAG_TOOL_QUARANTINE}/netMHCstabpan/netMHCstabpan" ]]; then
  export NETMHCSTABPAN_HOME="${NEOAG_TOOL_QUARANTINE}/netMHCstabpan"
fi
if [[ -x "${NETMHCSTABPAN_HOME}/netMHCstabpan" ]]; then
  export PATH="${NETMHCSTABPAN_HOME}:${PATH}"
fi

# LOHHLA / FACETS / Nextflow
export LOHHLA_HOME="${LOHHLA_HOME:-${NEOAG_TOOLS_ROOT}/tools/lohhla}"
if [[ ! -f "${LOHHLA_HOME}/LOHHLAscript.R" && -f "${NEOAG_TOOLS_ROOT}/../neoag_event_pipeline_artifact_quarantine_20260622_091158/tools/lohhla/LOHHLAscript.R" ]]; then
  export LOHHLA_HOME="${NEOAG_TOOLS_ROOT}/../neoag_event_pipeline_artifact_quarantine_20260622_091158/tools/lohhla"
fi
export POLYSOLVER_HOME="${POLYSOLVER_HOME:-}"
export NOVOALIGN_LICENSE_FILE="${NOVOALIGN_LICENSE_FILE:-}"
export FACETS_HOME="${NEOAG_TOOLS_ROOT}/bin"
if [[ ! -x "${FACETS_HOME}/runFACETS.R" ]]; then
  export FACETS_HOME="${NEOAG_TOOLS_ROOT}/tools/facets"
fi
# ASCAT / PyClone-VI — install: bash scripts/install_ascat_pyclone.sh
export NEOAG_ASCAT_ENV="neoag-ascat"
export ASCAT_HOME="${NEOAG_TOOLS_ROOT}/bin"
if [[ ! -x "${ASCAT_HOME}/ascat.R" ]]; then
  export ASCAT_HOME="${NEOAG_TOOLS_ROOT}/tools/ascat"
fi
export NEOAG_ASCAT_V3_ENV="${NEOAG_ASCAT_V3_ENV:-neoag-ascat-v3}"
export NEOAG_ASCAT_V3_BIN="${NEOAG_TOOLS_ROOT}/bin/ascat-v3"
export NEOAG_PYCLONE_ENV="neoag-pyclone"
export NEOAG_PYCLONE_BIN="${NEOAG_CONDA_BASE}/envs/${NEOAG_PYCLONE_ENV}/bin/pyclone-vi"
# FACETS snp-pileup common SNP VCF (optional; set in conf/tools.env.local.sh)
export NEOAG_DBSNP_VCF="${NEOAG_DBSNP_VCF:-}"
export NEOAG_SV_ENV="neoag-sv"
export NEOAG_MANTA_ENV="neoag-manta"
export PATH="${NEOAG_TOOLS_ROOT}/bin:${LOHHLA_HOME}:${FACETS_HOME}:${ASCAT_HOME}:${PATH}"
if [[ -d "${POLYSOLVER_HOME}/scripts" ]]; then
  export PATH="${POLYSOLVER_HOME}/binaries:${POLYSOLVER_HOME}/scripts:${PATH}"
fi
if [[ -d "${NEOAG_CONDA_BASE}/envs/${NEOAG_ASCAT_ENV}/bin" ]]; then
  export PATH="${NEOAG_CONDA_BASE}/envs/${NEOAG_ASCAT_ENV}/bin:${PATH}"
fi
if [[ -d "${NEOAG_CONDA_BASE}/envs/${NEOAG_SV_ENV}/bin" ]]; then
  export PATH="${NEOAG_CONDA_BASE}/envs/${NEOAG_SV_ENV}/bin:${PATH}"
fi
if [[ -d "${NEOAG_CONDA_BASE}/envs/${NEOAG_MANTA_ENV}/bin" ]]; then
  export PATH="${NEOAG_CONDA_BASE}/envs/${NEOAG_MANTA_ENV}/bin:${PATH}"
fi

# GATK4 (Mutect2 / FilterMutectCalls) — install: bash scripts/install_gatk.sh
export NEOAG_GATK_ENV="neoag-gatk"
if [[ -d "${NEOAG_CONDA_BASE}/envs/${NEOAG_GATK_ENV}/bin" ]]; then
  export PATH="${NEOAG_CONDA_BASE}/envs/${NEOAG_GATK_ENV}/bin:${PATH}"
fi

# RNA fusion callers + EasyFuse — install: bash scripts/install_fusion_tools.sh
# Shared refs (CTAT + EasyFuse tarballs): bash scripts/setup_fusion_refs_from_shared.sh
export NEOAG_FUSION_ENV="neoag-fusion"
export NEOAG_STAR_FUSION_HOME="${NEOAG_TOOLS_ROOT}/tools/STAR-Fusion"
export NEOAG_CTAT_LIB_DIR="${NEOAG_TOOLS_ROOT}/data/ref/ctat"
export NEOAG_FUSIONCATCHER_HOME="${NEOAG_TOOLS_ROOT}/tools/fusioncatcher"
export NEOAG_STAR_FUSION_BIN="${NEOAG_TOOLS_ROOT}/bin/star-fusion-neoag"
export NEOAG_FUSIONCATCHER_BIN="${NEOAG_TOOLS_ROOT}/bin/fusioncatcher-neoag"
export NEOAG_EASYFUSE_HOME="${NEOAG_TOOLS_ROOT}/tools/EasyFuse"
export NEOAG_SHARED_REF_DIR="${NEOAG_SHARED_REF_DIR:-}"
export NEOAG_CTAT_ARCHIVE="${NEOAG_CTAT_ARCHIVE:-${NEOAG_SHARED_REF_DIR:+${NEOAG_SHARED_REF_DIR}/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz}}"
export NEOAG_EASYFUSE_ARCHIVE="${NEOAG_EASYFUSE_ARCHIVE:-${NEOAG_SHARED_REF_DIR:+${NEOAG_SHARED_REF_DIR}/easyfuse_ref_v4.tar.gz}}"
_NEOAG_CTAT_SHARED="${NEOAG_SHARED_REF_DIR}/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play"
_NEOAG_EASYFUSE_SHARED="${NEOAG_SHARED_REF_DIR}/easyfuse_ref_v4"
if [[ -d "${_NEOAG_CTAT_SHARED}" ]]; then
  export CTAT_GENOME_LIB="$(cd "${_NEOAG_CTAT_SHARED}" && pwd)"
elif [[ -L "${NEOAG_CTAT_LIB_DIR}/current" || -d "${NEOAG_CTAT_LIB_DIR}/current" ]]; then
  export CTAT_GENOME_LIB="$(cd "${NEOAG_CTAT_LIB_DIR}/current" && pwd)"
fi
if [[ -d "${_NEOAG_EASYFUSE_SHARED}" ]]; then
  export NEOAG_EASYFUSE_REF="$(cd "${_NEOAG_EASYFUSE_SHARED}" && pwd)"
fi
unset _NEOAG_CTAT_SHARED _NEOAG_EASYFUSE_SHARED
if [[ -d "${NEOAG_CONDA_BASE}/envs/${NEOAG_FUSION_ENV}/bin" ]]; then
  export PATH="${NEOAG_CONDA_BASE}/envs/${NEOAG_FUSION_ENV}/bin:${PATH}"
fi
if [[ -d "${NEOAG_FUSIONCATCHER_HOME}/bin" ]]; then
  export PATH="${NEOAG_FUSIONCATCHER_HOME}/bin:${PATH}"
fi
if [[ -d "${NEOAG_STAR_FUSION_HOME}" ]]; then
  export PATH="${NEOAG_STAR_FUSION_HOME}:${PATH}"
fi

# neoag-tools stays ahead of gatk/sv/manta python shims, while this checkout
# keeps priority for neoag wrapper scripts such as bin/neoag.
export PATH="${NEOAG_TOOLS_ROOT}/bin:${NEOAG_CONDA_BASE}/envs/neoag-tools/bin:${PATH}"

if [[ -f "${NEOAG_PROJECT_ROOT}/conf/tools.env.local.sh" ]]; then
  # shellcheck source=/dev/null
  source "${NEOAG_PROJECT_ROOT}/conf/tools.env.local.sh"
fi


# OptiType (HLA-I typing from DNA/RNA FASTQ or BAM)
export NEOAG_OPTITYPE_ENV="${NEOAG_OPTITYPE_ENV:-neoag-optitype}"
if [[ -z "${OPTITYPE_ENV:-}" ]]; then
  if [[ -d "${NEOAG_CONDA_BASE}/envs/${NEOAG_OPTITYPE_ENV}" ]]; then
    export OPTITYPE_ENV="${NEOAG_CONDA_BASE}/envs/${NEOAG_OPTITYPE_ENV}"
  fi
fi
if [[ -n "${OPTITYPE_ENV:-}" && -x "${OPTITYPE_ENV}/bin/optitype" ]]; then
  export OPTITYPE_BIN="${OPTITYPE_ENV}/bin/optitype"
  export OPTITYPE_REFERENCE="${OPTITYPE_ENV}/share/optitype/data"
  export PATH="${OPTITYPE_ENV}/bin:${PATH}"
fi
