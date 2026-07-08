#!/bin/bash
# ============================================================
# Docker 容器测试: GRIDSS2 (Genomic Rearrangement IDentification)
# 镜像: gridss/gridss:2.13.2
# 参考: tests/test_sv_phase1.py (SV pipeline with GRIDSS2 caller)
#       tests/test_sv_wes_phase1_5.py (SV WES mode)
# 数据: data/fixtures_sv/mini_sv.vcf, mini_ref.fa
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="gridss/gridss:2.13.2"

echo "=== GRIDSS2 容器测试 ==="
echo "镜像: $IMAGE"
echo ""

# ---- 准备测试数据 ----
WORKDIR="${SCRIPT_DIR}/test_output/gridss2"
mkdir -p "${WORKDIR}"

FIXTURE_REF="${PROJECT_ROOT}/data/fixtures_sv/mini_ref.fa"
FIXTURE_VCF="${PROJECT_ROOT}/data/fixtures_sv/mini_sv.vcf"

echo ">>> 检查测试数据:"
echo "    REF: ${FIXTURE_REF}"
echo "    VCF: ${FIXTURE_VCF}"
echo ""

# ---- 验证命令可用性 ----
echo ">>> gridss --help..."
docker run --rm "${IMAGE}" gridss --help 2>&1 | head -25

echo ""
echo ">>> gridsstools --help..."
docker run --rm "${IMAGE}" gridsstools 2>&1 | head -10 || true

# ---- 演示完整调用格式 ----
echo ""
echo "  # GRIDSS2 预处理 (需 BAM 输入):"
echo "  docker run --rm \\"
echo "    -v \$(pwd)/input:/input \\"
echo "    -v \$(pwd)/ref:/ref \\"
echo "    -v \$(pwd)/output:/output \\"
echo "    ${IMAGE} \\"
echo "    gridss -s preprocess \\"
echo "      --reference /ref/hg38.fa \\"
echo "      --output /output/preprocessed.bam \\"
echo "      /input/tumor.bam /input/normal.bam"

echo ""
echo "  # GRIDSS2 体细胞 SV 过滤:"
echo "  docker run --rm \\"
echo "    -v \$(pwd):/data \\"
echo "    ${IMAGE} \\"
echo "    gridss_somatic_filter \\"
echo "      --input /data/gridss.vcf \\"
echo "      --output /data/gridss.somatic.vcf \\"
echo "      --fulloutput /data/gridss.somatic.full.vcf"

echo ""
echo "=== GRIDSS2 容器测试完成 ==="
