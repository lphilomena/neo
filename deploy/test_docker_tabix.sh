#!/bin/bash
# ============================================================
# Docker 容器测试: tabix (VCF/BCF 索引, htslib 组件)
# 镜像: quay.io/biocontainers/tabix:1.11--hdfd78af_0
# 参考: tests/test_snv_phase1_wes.py (VCF handling)
#       tests/test_pvacseq_enrich.py (VCF annotation)
# 数据: data/fixtures_snv/mini_somatic.vcf (VCF 索引测试)
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="quay.io/biocontainers/tabix:1.11--hdfd78af_0"

echo "=== tabix 容器测试 ==="
echo "镜像: $IMAGE"
echo ""

# ---- 准备测试数据 ----
WORKDIR="${SCRIPT_DIR}/test_output/tabix"
mkdir -p "${WORKDIR}"

FIXTURE_VCF="${PROJECT_ROOT}/data/fixtures_snv/mini_somatic.vcf"

echo ">>> 检查测试数据:"
echo "    VCF: ${FIXTURE_VCF}"
echo ""

# ---- 验证命令可用性 ----
echo ">>> tabix --version..."
docker run --rm "${IMAGE}" tabix --version 2>&1 | head -2

echo ""
echo ">>> bgzip --version..."
docker run --rm "${IMAGE}" bgzip --version 2>&1 | head -2

# ---- 功能测试: VCF 压缩 + 索引 + 查询 ----
if [[ -f "${FIXTURE_VCF}" ]]; then
    echo ""
    echo ">>> 测试1: VCF → bgzip 压缩..."
    docker run --rm \
        -v "${FIXTURE_VCF}:/data/input.vcf:ro" \
        -v "${WORKDIR}:/output" \
        "${IMAGE}" \
        bash -c '
            cp /data/input.vcf /output/test.vcf
            bgzip -c /output/test.vcf > /output/test.vcf.gz
            echo "Compressed size:"
            ls -lh /output/test.vcf.gz
        '

    echo ""
    echo ">>> 测试2: tabix 索引 VCF..."
    docker run --rm \
        -v "${WORKDIR}:/data" \
        "${IMAGE}" \
        tabix -p vcf /data/test.vcf.gz

    echo "    Index file created:"
    ls -la "${WORKDIR}/test.vcf.gz.tbi"

    echo ""
    echo ">>> 测试3: tabix 区间查询..."
    docker run --rm \
        -v "${WORKDIR}:/data" \
        "${IMAGE}" \
        bash -c '
            # 从 VCF 中提取第一条变异的染色质和位置
            CHR=$(zgrep -m1 -v "^#" /data/test.vcf.gz | cut -f1)
            POS=$(zgrep -m1 -v "^#" /data/test.vcf.gz | cut -f2)
            echo "Query region: ${CHR}:${POS}-$((POS+1000))"
            tabix /data/test.vcf.gz "${CHR}:${POS}-$((POS+1000))" 2>&1
        '

    echo ""
    echo ">>> 测试4: bgzip 解压验证..."
    docker run --rm \
        -v "${WORKDIR}:/data" \
        "${IMAGE}" \
        bash -c '
            zcat /data/test.vcf.gz | head -10
            echo "..."
            echo "Total lines in compressed VCF: $(zcat /data/test.vcf.gz | wc -l)"
        '
fi

# ---- 演示完整调用格式 ----
echo ""
echo "  # tabix 常用操作示例:"
echo "  # 1. 压缩 VCF"
echo "  docker run --rm \\"
echo "    -v \$(pwd):/data \\"
echo "    ${IMAGE} \\"
echo "    bgzip /data/variants.vcf"
echo ""
echo "  # 2. 索引 VCF"
echo "  docker run --rm \\"
echo "    -v \$(pwd):/data \\"
echo "    ${IMAGE} \\"
echo "    tabix -p vcf /data/variants.vcf.gz"
echo ""
echo "  # 3. 区间查询"
echo "  docker run --rm \\"
echo "    -v \$(pwd):/data \\"
echo "    ${IMAGE} \\"
echo "    tabix /data/variants.vcf.gz chr1:1000000-2000000"

echo ""
echo "=== tabix 容器测试完成 ==="
