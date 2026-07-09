#!/bin/bash
# ============================================================
# Docker 容器测试: GATK4 (Mutect2 体细胞变异检测)
# 镜像: broadinstitute/gatk:4.6.2.0
# 参考: tests/test_snv_phase1_wes.py (snv_wes_full / Mutect2)
#       tests/test_tools.py (run_tool "gatk_mutect2")
# 镜像内置: samtools, bcftools, tabix, Python3
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="broadinstitute/gatk:4.6.2.0"

echo "=== GATK4 容器测试 ==="
echo "镜像: $IMAGE"
echo ""

# ---- 准备测试数据 ----
WORKDIR="${SCRIPT_DIR}/test_output/gatk4"
mkdir -p "${WORKDIR}"

FIXTURE_REF="${PROJECT_ROOT}/data/fixtures_snv/mini_ref.fa"
FIXTURE_VCF="${PROJECT_ROOT}/data/fixtures_snv/mini_somatic.vcf"

echo ">>> 检查测试数据:"
echo "    REF: ${FIXTURE_REF}"
echo "    VCF: ${FIXTURE_VCF}"
echo ""

# ---- 验证 GATK 命令可用性 ----
echo ">>> gatk --help..."
docker run --rm "${IMAGE}" gatk --help 2>&1 | head -15

echo ""
echo ">>> gatk Mutect2 --help..."
docker run --rm "${IMAGE}" gatk Mutect2 --help 2>&1 | head -20

echo ""
echo ">>> 内置工具版本检测:"
echo -n "    samtools: "
docker run --rm "${IMAGE}" samtools --version 2>&1 | head -1
echo -n "    bcftools: "
docker run --rm "${IMAGE}" bcftools --version 2>&1 | head -1
echo -n "    tabix: "
docker run --rm "${IMAGE}" tabix --version 2>&1 | head -1

# ---- 演示完整调用格式 ----
echo ""
echo "  # GATK4 Mutect2 完整调用示例 (需挂载 BAM + REF):"
echo "  docker run --rm \\"
echo "    -v \$(pwd)/input:/input \\"
echo "    -v \$(pwd)/ref:/ref \\"
echo "    -v \$(pwd)/output:/output \\"
echo "    ${IMAGE} \\"
echo "    gatk Mutect2 \\"
echo "      -R /ref/hg38.fa \\"
echo "      -I /input/tumor.bam \\"
echo "      -I /input/normal.bam \\"
echo "      -L /input/wes_capture.bed \\"
echo "      -O /output/somatic.vcf"

echo ""
echo "  # 也可直接使用容器内的 samtools:"
echo "  docker run --rm \\"
echo "    -v \$(pwd):/data \\"
echo "    ${IMAGE} \\"
echo "    samtools index /data/input.bam"

echo ""
echo "=== GATK4 容器测试完成 ==="
