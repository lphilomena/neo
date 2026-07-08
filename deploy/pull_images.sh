#!/bin/bash
# ============================================================
# 项目工具 Docker 镜像拉取脚本
#
# 拉取所有官方维护的 Docker 镜像和 Biocontainers 镜像。
# 所有镜像均指定明确的版本标签 (tag)，不使用 latest。
#
# 用法:
#   bash deploy/pull_images.sh          # 拉取所有镜像
#   bash deploy/pull_images.sh --dry-run  # 仅打印命令，不执行
# ============================================================

set -euo pipefail

DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
fi

pull() {
    local image="$1"
    local desc="$2"
    echo ">>> [$desc]"
    echo "    $image"
    if $DRY_RUN; then
        echo "    (dry-run, skip)"
    else
        docker pull "$image"
    fi
    echo ""
}

echo ""
echo "============================================================"
echo "  一、官方 Docker 镜像"
echo "============================================================"
echo ""

# ---- 变异注释 ----
# VEP 115 - Ensembl Variant Effect Predictor
# 官方镜像: https://hub.docker.com/r/ensemblorg/ensembl-vep
pull "ensemblorg/ensembl-vep:release_115.2" \
     "VEP 115 — 变异功能注释"

# ---- 新抗原预测 ----
# pVACtools (含 pVACseq / pVACfuse / pVACsplice)
# 官方镜像: https://hub.docker.com/r/griffithlab/pvactools
# 版本 6.1.1 对应 conda env 中 "pvactools>=6.1,<7"
pull "griffithlab/pvactools:6.1.1" \
     "pVACtools 6.1.1 — SNV/InDel/融合/剪接位点 新抗原预测"

# ---- HLA 分型 ----
# OptiType - HLA 分型
# 官方镜像: https://hub.docker.com/r/fred2/optitype
pull "fred2/optitype:release-v1.3.1" \
     "OptiType 1.3.1 — HLA 分型"

# ---- 变异检测 ----
# GATK4 - Mutect2 体细胞变异检测
# 官方镜像: https://hub.docker.com/r/broadinstitute/gatk
# 版本 4.6.2.0 对应 conda env 中 gatk4=4.6.2
pull "broadinstitute/gatk:4.6.2.0" \
     "GATK4 4.6.2.0 — Mutect2 体细胞变异检测"

# ---- 融合检测 ----
# Arriba - RNA-seq 融合检测
# 官方镜像: https://hub.docker.com/r/uhrigs/arriba
# 版本 2.5.1 对应 conda env 中 arriba=2.5.1
pull "uhrigs/arriba:2.5.1" \
     "Arriba 2.5.1 — RNA-seq 融合检测"

# EasyFuse - 融合 meta-caller (整合 STAR-Fusion / FusionCatcher / MapSplice2 / InFusion / SOAPfuse)
# 官方镜像: https://hub.docker.com/r/tronbioinformatics/easyfuse
pull "tronbioinformatics/easyfuse:1.3.7" \
     "EasyFuse 1.3.7 — 融合 meta-caller"

# ---- 结构变异检测 ----
# GRIDSS2 - Genomic Rearrangement IDentification Software Suite
# 官方镜像: https://hub.docker.com/r/gridss/gridss
pull "gridss/gridss:2.13.2" \
     "GRIDSS2 2.13.2 — SV 检测 (break-end assembly)"

# DELLY - 结构变异检测 (integrated paired-end + split-read)
# 官方镜像: https://hub.docker.com/r/dellytools/delly
pull "dellytools/delly:v2.3.0" \
     "DELLY 2.3.0 — SV 检测 (paired-end + split-read)"

echo ""
echo "============================================================"
echo "  二、Biocontainers 镜像 (quay.io)"
echo "============================================================"
echo ""

# ---- HLA 分型与 LOH ----
# LOHHLA - HLA 等位基因特异性 LOH
# Biocontainers: https://quay.io/repository/biocontainers/lohhla
pull "quay.io/biocontainers/lohhla:20171108--hdfd78af_3" \
     "LOHHLA — HLA 等位基因特异性 LOH"

# ---- 肿瘤纯度与拷贝数 ----
# ASCAT v2.5.2 - 等位基因特异性拷贝数与纯度 (R 4.0)
# Biocontainers: https://quay.io/repository/biocontainers/ascat
# r40 对应 conda env 中 r-base=4.0
pull "quay.io/biocontainers/ascat:2.5.2--r40hdfd78af_3" \
     "ASCAT 2.5.2 (R 4.0) — 等位基因特异性拷贝数与纯度"

# ASCAT v3.2.0 - 等位基因特异性拷贝数与纯度 (R 4.4)
# r44 对应 conda env 中 r-base=4.4
pull "quay.io/biocontainers/ascat:3.2.0--r44hdfd78af_1" \
     "ASCAT 3.2.0 (R 4.4) — 等位基因特异性拷贝数与纯度 v3"

# ---- 融合检测 ----
# STAR-Fusion - 嵌合转录本检测
# Biocontainers: https://quay.io/repository/biocontainers/star-fusion
# 注: 官方 Docker 镜像 trinityctat/ctatfusion 可能已不可用，改用 biocontainers
pull "quay.io/biocontainers/star-fusion:1.15.1--hdfd78af_1" \
     "STAR-Fusion 1.15.1 — 嵌合转录本检测"

# FusionCatcher - 体细胞融合检测
# Biocontainers: https://quay.io/repository/biocontainers/fusioncatcher
pull "quay.io/biocontainers/fusioncatcher:1.33b--hdfd78af_0" \
     "FusionCatcher 1.33 — 体细胞融合检测"

# ---- 结构变异检测 ----
# Manta - SV/Indel 检测
# Biocontainers: https://quay.io/repository/biocontainers/manta
pull "quay.io/biocontainers/manta:1.6.0--h9ee0642_3" \
     "Manta 1.6.0 — SV/Indel 检测"

# SvABA - 结构变异和 indel 分析 (assembly-based)
# Biocontainers: https://quay.io/repository/biocontainers/svaba
pull "quay.io/biocontainers/svaba:1.2.0--h69ac913_1" \
     "SvABA 1.2.0 — SV/Indel 检测 (assembly-based)"

# ---- 基础依赖（Biocontainers 备用镜像） ----
# samtools - BAM/VCF 操作
# Biocontainers: https://quay.io/repository/biocontainers/samtools
pull "quay.io/biocontainers/samtools:1.23.1--ha83d96e_0" \
     "samtools 1.23.1 — BAM/VCF 操作"

# tabix - VCF/BCF 索引 (htslib 组件)
# Biocontainers: https://quay.io/repository/biocontainers/tabix
pull "quay.io/biocontainers/tabix:1.11--hdfd78af_0" \
     "tabix 1.11 — VCF 索引"

echo ""
echo "============================================================"
echo "  拉取完成"
echo "============================================================"
echo ""
echo "提示: 如需验证已拉取的镜像，请运行:"
echo "  docker images | grep -E 'ensemblorg|griffithlab|fred2|broadinstitute|uhrigs|tronbioinformatics|gridss|dellytools|biocontainers'"
echo ""
