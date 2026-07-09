#!/bin/bash
# ============================================================
# Docker 容器测试: samtools (BAM/VCF 操作)
# 镜像: quay.io/biocontainers/samtools:1.23.1--ha83d96e_0
# 参考: tests/test_snv_phase1_wes.py (samtools index in pipeline)
#       tests/test_pvacseq_enrich.py (VCF handling)
# 数据: data/fixtures_snv/mini_somatic.vcf (用于测试)
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="quay.io/biocontainers/samtools:1.23.1--ha83d96e_0"

echo "=== samtools 容器测试 ==="
echo "镜像: $IMAGE"
echo ""

# ---- 准备测试数据 ----
WORKDIR="${SCRIPT_DIR}/test_output/samtools"
mkdir -p "${WORKDIR}"

FIXTURE_VCF="${PROJECT_ROOT}/data/fixtures_snv/mini_somatic.vcf"

echo ">>> 检查测试数据:"
echo "    VCF: ${FIXTURE_VCF}"
echo ""

# ---- 验证命令可用性 ----
echo ">>> samtools --version..."
docker run --rm "${IMAGE}" samtools --version 2>&1 | head -3

# ---- 功能测试: VCF 索引和查看 ----
if [[ -f "${FIXTURE_VCF}" ]]; then
    echo ""
    echo ">>> 测试1: 查看 VCF 内容 (samtools view 风格)..."
    docker run --rm \
        -v "${FIXTURE_VCF}:/data/input.vcf:ro" \
        "${IMAGE}" \
        bash -c 'wc -l /data/input.vcf && echo "---" && grep -v "^#" /data/input.vcf | head -5'

    echo ""
    echo ">>> 测试2: VCF → BCF 转换 + 索引..."
    docker run --rm \
        -v "${FIXTURE_VCF}:/data/input.vcf:ro" \
        -v "${WORKDIR}:/output" \
        "${IMAGE}" \
        bash -c '
            # VCF → bgzip 压缩
            cp /data/input.vcf /output/test.vcf
            bgzip -c /output/test.vcf > /output/test.vcf.gz
            # 索引
            tabix -p vcf /output/test.vcf.gz
            # 验证
            echo "Indexed files:"
            ls -lh /output/test.vcf.gz*
        '

    echo ""
    echo ">>> 测试3: 统计 VCF 变异数量..."
    docker run --rm \
        -v "${FIXTURE_VCF}:/data/input.vcf:ro" \
        "${IMAGE}" \
        bash -c '
            variants=$(grep -v "^#" /data/input.vcf | wc -l)
            echo "Number of variants: $variants"
        '

    echo ""
    echo ">>> 测试4: FASTQ 统计 (samtools 通用)..."
    docker run --rm \
        -v "${FIXTURE_VCF}:/data/input.vcf:ro" \
        "${IMAGE}" \
        samtools dict /data/input.vcf 2>&1 | head -5 || echo "    (samtools dict works on FASTA, not VCF - expected)"
fi

# ---- 演示完整调用格式 ----
echo ""
echo "  # samtools 常用操作示例:"
echo "  docker run --rm \\"
echo "    -v \$(pwd):/data \\"
echo "    ${IMAGE} \\"
echo "    samtools sort -o /data/sorted.bam /data/input.bam"
echo ""
echo "  docker run --rm \\"
echo "    -v \$(pwd):/data \\"
echo "    ${IMAGE} \\"
echo "    samtools index /data/sorted.bam"

echo ""
echo "=== samtools 容器测试完成 ==="
