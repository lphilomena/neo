#!/bin/bash
# ============================================================
# Docker 容器测试: OptiType (HLA 分型)
# 镜像: fred2/optitype:release-v1.3.1
# 参考: tests/test_tools.py (TOOL_REGISTRY optitype)
#       tests/test_snv_phase1_wes.py (HLA typing in SNV pipeline)
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="fred2/optitype:release-v1.3.1"

echo "=== OptiType 容器测试 ==="
echo "镜像: $IMAGE"
echo ""

# ---- 验证命令可用性 ----
echo ">>> OptiTypePipeline.py --help..."
docker run --rm "${IMAGE}" OptiTypePipeline.py --help 2>&1 | head -30

echo ""
echo ">>> razers3 --help..."
docker run --rm --entrypoint="razers3" "${IMAGE}" --help 2>&1 | head -10

# ---- 演示完整调用格式 ----
echo ""
echo "  # OptiType 完整调用示例 (DNA 模式, 需 FASTQ 输入):"
echo "  docker run --rm \\"
echo "    -v \$(pwd)/input:/data/input \\"
echo "    -v \$(pwd)/output:/data/output \\"
echo "    ${IMAGE} \\"
echo "    -i /data/input/sample_R1.fastq /data/input/sample_R2.fastq \\"
echo "    --dna -v -o /data/output/"

echo ""
echo "  # OptiType 完整调用示例 (RNA 模式):"
echo "  docker run --rm \\"
echo "    -v \$(pwd)/input:/data/input \\"
echo "    -v \$(pwd)/output:/data/output \\"
echo "    ${IMAGE} \\"
echo "    -i /data/input/sample_R1.fastq /data/input/sample_R2.fastq \\"
echo "    --rna -v -o /data/output/"

echo ""
echo "=== OptiType 容器测试完成 ==="
