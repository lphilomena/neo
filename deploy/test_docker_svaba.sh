#!/bin/bash
# ============================================================
# Docker 容器测试: SvABA (结构变异/Indel 分析 by assembly)
# 镜像: quay.io/biocontainers/svaba:1.2.0--h69ac913_1
# 参考: tests/test_sv_phase1.py (SV pipeline)
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="quay.io/biocontainers/svaba:1.2.0--h69ac913_1"

echo "=== SvABA 容器测试 ==="
echo "镜像: $IMAGE"
echo ""

# ---- 验证命令可用性 ----
echo ">>> svaba --help..."
docker run --rm "${IMAGE}" svaba --help 2>&1 | head -20 || true

echo ""
echo ">>> svaba 版本..."
docker run --rm "${IMAGE}" svaba 2>&1 | head -5 || true

# ---- 演示完整调用格式 ----
echo ""
echo "  # SvABA 体细胞 SV/Indel 检测 (需 Tumor/Normal BAM + REF):"
echo "  docker run --rm \\"
echo "    -v \$(pwd)/input:/input \\"
echo "    -v \$(pwd)/ref:/ref \\"
echo "    -v \$(pwd)/output:/output \\"
echo "    ${IMAGE} \\"
echo "    svaba run \\"
echo "      -t /input/tumor.bam \\"
echo "      -n /input/normal.bam \\"
echo "      -a SAMPLE_ID \\"
echo "      -G /ref/hg38.fa \\"
echo "      -p 8 \\"
echo "      -k /ref/dbsnp_common.vcf.gz"

echo ""
echo "  # SvABA 仅 germline SV 检测:"
echo "  docker run --rm \\"
echo "    -v \$(pwd)/input:/input \\"
echo "    -v \$(pwd)/ref:/ref \\"
echo "    -v \$(pwd)/output:/output \\"
echo "    ${IMAGE} \\"
echo "    svaba run \\"
echo "      -t /input/normal.bam \\"
echo "      -a SAMPLE_ID \\"
echo "      -G /ref/hg38.fa \\"
echo "      -p 8 \\"
echo "      --germline"

echo ""
echo "=== SvABA 容器测试完成 ==="
