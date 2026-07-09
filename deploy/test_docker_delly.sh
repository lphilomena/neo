#!/bin/bash
# ============================================================
# Docker 容器测试: DELLY (结构变异检测)
# 镜像: dellytools/delly:v2.3.0
# 参考: tests/test_sv_phase1.py (SV pipeline)
#       tests/test_sv_wes_phase1_5.py (SV WES mode)
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="dellytools/delly:v2.3.0"

echo "=== DELLY 容器测试 ==="
echo "镜像: $IMAGE"
echo ""

# ---- 验证命令可用性 ----
echo ">>> delly --help..."
docker run --rm "${IMAGE}" delly --help 2>&1 | head -25

echo ""
echo ">>> delly sr --help..."
docker run --rm "${IMAGE}" delly sr --help 2>&1 | head -15 || true

# ---- 演示完整调用格式 ----
echo ""
echo "  # DELLY 胚系 SV 检测 (需 BAM + REF + 排除区域 BED):"
echo "  docker run --rm \\"
echo "    -v \$(pwd)/input:/input \\"
echo "    -v \$(pwd)/ref:/ref \\"
echo "    -v \$(pwd)/output:/output \\"
echo "    ${IMAGE} \\"
echo "    delly sr \\"
echo "      -x /ref/exclude.bed \\"
echo "      -o /output/delly.bcf \\"
echo "      -g /ref/hg38.fa \\"
echo "      /input/sample.bam"

echo ""
echo "  # DELLY 体细胞 SV 检测:"
echo "  docker run --rm \\"
echo "    -v \$(pwd)/input:/input \\"
echo "    -v \$(pwd)/ref:/ref \\"
echo "    -v \$(pwd)/output:/output \\"
echo "    ${IMAGE} \\"
echo "    delly sr \\"
echo "      -x /ref/exclude.bed \\"
echo "      -o /output/somatic.bcf \\"
echo "      -g /ref/hg38.fa \\"
echo "      /input/tumor.bam /input/normal.bam"

echo ""
echo "=== DELLY 容器测试完成 ==="
