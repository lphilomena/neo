#!/bin/bash
# ============================================================
# Docker 容器测试: Manta (结构变异/Indel 检测)
# 镜像: quay.io/biocontainers/manta:1.6.0--h9ee0642_3
# 参考: tests/test_sv_phase1.py (SV pipeline)
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="quay.io/biocontainers/manta:1.6.0--h9ee0642_3"

echo "=== Manta 容器测试 ==="
echo "镜像: $IMAGE"
echo ""

# ---- 验证命令可用性 ----
echo ">>> configManta.py --help..."
docker run --rm "${IMAGE}" configManta.py --help 2>&1 | head -20 || true

echo ""
echo ">>> 容器内 Python 环境:"
docker run --rm "${IMAGE}" python --version 2>&1 || true

# ---- 演示完整调用格式 ----
echo ""
echo "  # Manta 体细胞 SV 检测完整调用 (需 Tumor/Normal BAM + REF):"
echo "  # 第1步: 配置运行"
echo "  docker run --rm \\"
echo "    -v \$(pwd)/input:/input \\"
echo "    -v \$(pwd)/ref:/ref \\"
echo "    -v \$(pwd)/output:/output \\"
echo "    ${IMAGE} \\"
echo "    configManta.py \\"
echo "      --tumorBam /input/tumor.bam \\"
echo "      --normalBam /input/normal.bam \\"
echo "      --referenceFasta /ref/hg38.fa \\"
echo "      --runDir /output/manta_run"
echo ""
echo "  # 第2步: 执行分析"
echo "  docker run --rm \\"
echo "    -v \$(pwd)/output:/output \\"
echo "    ${IMAGE} \\"
echo "    /output/manta_run/runWorkflow.py -m local -j 8"

echo ""
echo "=== Manta 容器测试完成 ==="
