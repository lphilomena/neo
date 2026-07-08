#!/bin/bash
# ============================================================
# Docker 容器测试: pVACtools (pVACseq/pVACfuse/pVACsplice)
# 镜像: griffithlab/pvactools:6.1.1
# 参考: tests/test_pvacseq_enrich.py (pVACseq VCF 处理)
#       tests/test_tools.py (run_tool "pvacseq")
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="griffithlab/pvactools:6.1.1"

echo "=== pVACtools 容器测试 ==="
echo "镜像: $IMAGE"
echo ""

# ---- 准备测试数据 ----
WORKDIR="${SCRIPT_DIR}/test_output/pvactools"
mkdir -p "${WORKDIR}"

# 使用项目 fixtures 中的 mini VCF
FIXTURE_VCF="${PROJECT_ROOT}/data/fixtures_snv/mini_somatic.vcf"
FIXTURE_HLA="${PROJECT_ROOT}/data/fixtures_snv/hla.txt"

echo ">>> 检查测试数据:"
if [[ -f "${FIXTURE_VCF}" ]]; then
    echo "    VCF: ${FIXTURE_VCF} (存在)"
else
    echo "    VCF: ${FIXTURE_VCF} (不存在，跳过数据处理测试)"
fi
echo ""

# ---- 验证 pVACtools 命令可用性 ----
echo ">>> pvacseq --help..."
docker run --rm "${IMAGE}" pvacseq --help 2>&1 | head -30

echo ""
echo ">>> pvacfuse --help..."
docker run --rm "${IMAGE}" pvacfuse --help 2>&1 | head -10

echo ""
echo ">>> pvacsplice --help..."
docker run --rm "${IMAGE}" pvacsplice --help 2>&1 | head -10

# ---- 演示完整调用格式 ----
echo ""
echo "  # pVACseq 完整调用示例 (需挂载输入 VCF + HLA 等位基因):"
echo "  docker run --rm \\"
echo "    -v \$(pwd)/input:/input \\"
echo "    -v \$(pwd)/output:/output \\"
echo "    ${IMAGE} \\"
echo "    pvacseq run /input/somatic.vcf \\"
echo "      SAMPLE_NAME \\"
echo "      'HLA-A*02:01,HLA-B*07:02' \\"
echo "      NetMHCpan /output/pvacseq \\"
echo "      -e1 8,9,10,11 \\"
echo "      --iedb-install-directory /opt/iedb"

echo ""
echo "  # pVACfuse 完整调用示例:"
echo "  docker run --rm \\"
echo "    -v \$(pwd)/input:/input \\"
echo "    -v \$(pwd)/output:/output \\"
echo "    ${IMAGE} \\"
echo "    pvacfuse run /input/fusions.tsv \\"
echo "      SAMPLE_NAME \\"
echo "      'HLA-A*02:01' \\"
echo "      NetMHCpan /output/pvacfuse \\"
echo "      --iedb-install-directory /opt/iedb"

echo ""
echo "=== pVACtools 容器测试完成 ==="
