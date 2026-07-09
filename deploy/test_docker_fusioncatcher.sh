#!/bin/bash
# ============================================================
# Docker 容器测试: FusionCatcher (体细胞融合检测)
# 镜像: quay.io/biocontainers/fusioncatcher:1.33b--hdfd78af_0
# 参考: tests/test_tools.py (TOOL_REGISTRY fusioncatcher)
#       tests/test_diagnostic_fusion_rescue.py (fusion detection)
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="quay.io/biocontainers/fusioncatcher:1.33b--hdfd78af_0"

echo "=== FusionCatcher 容器测试 ==="
echo "镜像: $IMAGE"
echo ""

# ---- 验证命令可用性 ----
echo ">>> fusioncatcher --help..."
docker run --rm "${IMAGE}" fusioncatcher --help 2>&1 | head -25 || true

echo ""
echo ">>> fusioncatcher 版本..."
docker run --rm "${IMAGE}" bash -c 'fusioncatcher --version 2>&1 || fusioncatcher -v 2>&1 || true' | head -5

echo ""
echo ">>> 容器内依赖工具版本:"
docker run --rm "${IMAGE}" bash -c '
  echo -n "  Bowtie: "; bowtie --version 2>&1 | head -1
  echo -n "  Bowtie2: "; bowtie2 --version 2>&1 | head -1
  echo -n "  STAR: "; STAR --version 2>&1 | head -1
' 2>&1 || true

# ---- 演示完整调用格式 ----
echo ""
echo "  # FusionCatcher 完整调用示例 (需参考数据库 + FASTQ):"
echo "  docker run --rm \\"
echo "    -v /path/to/fusioncatcher_data:/data \\"
echo "    -v \$(pwd)/fastq:/fastq \\"
echo "    -v \$(pwd)/output:/output \\"
echo "    ${IMAGE} \\"
echo "    fusioncatcher \\"
echo "      -d /data/human_v98 \\"
echo "      -i /fastq \\"
echo "      -o /output \\"
echo "      -p 8"

echo ""
echo "  # 首次运行前需下载人类参考数据库:"
echo "  docker run --rm \\"
echo "    -v /path/to/fusioncatcher_data:/data \\"
echo "    ${IMAGE} \\"
echo "    download-human-db.sh"

echo ""
echo "=== FusionCatcher 容器测试完成 ==="
