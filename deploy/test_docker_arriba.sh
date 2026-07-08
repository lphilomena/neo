#!/bin/bash
# ============================================================
# Docker 容器测试: Arriba (RNA-seq 融合检测)
# 镜像: uhrigs/arriba:2.5.1
# 参考: tests/test_tools.py (TOOL_REGISTRY arriba)
#       tests/test_diagnostic_fusion_rescue.py (fusion detection)
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="uhrigs/arriba:2.5.1"

echo "=== Arriba 容器测试 ==="
echo "镜像: $IMAGE"
echo ""

# ---- 验证命令可用性 ----
echo ">>> arriba --help..."
docker run --rm "${IMAGE}" arriba -h 2>&1 | head -20

echo ""
echo ">>> extract_fastq 脚本..."
docker run --rm "${IMAGE}" run_arriba.sh 2>&1 | head -10 || true

echo ""
echo ">>> 容器内 STAR 版本:"
docker run --rm "${IMAGE}" STAR --version 2>&1 | head -3

# ---- 演示完整调用格式 ----
echo ""
echo "  # Arriba 完整调用示例 (需 STAR 比对后的 BAM):"
echo "  docker run --rm \\"
echo "    -v \$(pwd)/input:/input \\"
echo "    -v \$(pwd)/ref:/ref \\"
echo "    -v \$(pwd)/output:/output \\"
echo "    ${IMAGE} \\"
echo "    arriba \\"
echo "      -x /input/Aligned.out.bam \\"
echo "      -a /ref/assembly.fa \\"
echo "      -g /ref/annotation.gtf \\"
echo "      -b /ref/blacklist.tsv.gz \\"
echo "      -o /output/fusions.tsv \\"
echo "      -f /ref/known_fusions.tsv"

echo ""
echo "=== Arriba 容器测试完成 ==="
