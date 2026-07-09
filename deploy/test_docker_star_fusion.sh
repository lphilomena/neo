#!/bin/bash
# ============================================================
# Docker 容器测试: STAR-Fusion (嵌合转录本检测)
# 镜像: quay.io/biocontainers/star-fusion:1.15.1--hdfd78af_1
# 参考: tests/test_tools.py (TOOL_REGISTRY star_fusion)
#       tests/test_diagnostic_fusion_rescue.py (fusion detection)
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="quay.io/biocontainers/star-fusion:1.15.1--hdfd78af_1"

echo "=== STAR-Fusion 容器测试 ==="
echo "镜像: $IMAGE"
echo ""

# ---- 验证命令可用性 ----
echo ">>> STAR-Fusion --help..."
docker run --rm "${IMAGE}" STAR-Fusion --help 2>&1 | head -25

echo ""
echo ">>> STAR 版本:"
docker run --rm "${IMAGE}" STAR --version 2>&1 | head -3

# ---- 演示完整调用格式 ----
echo ""
echo "  # STAR-Fusion 完整调用示例 (需 CTAT genome lib + FASTQ):"
echo "  docker run --rm \\"
echo "    -v /path/to/ctat_lib:/ctat \\"
echo "    -v \$(pwd)/fastq:/data \\"
echo "    -v \$(pwd)/output:/output \\"
echo "    ${IMAGE} \\"
echo "    STAR-Fusion \\"
echo "      --genome_lib_dir /ctat/GRCh38_gencode_v37_CTAT_lib \\"
echo "      --left_fq /data/sample_R1.fastq.gz \\"
echo "      --right_fq /data/sample_R2.fastq.gz \\"
echo "      --output_dir /output/star_fusion_out \\"
echo "      --FusionInspector validate"

echo ""
echo "  # STAR-Fusion 仅预测 (已有 chimeric junction):"
echo "  docker run --rm \\"
echo "    -v /path/to/ctat_lib:/ctat \\"
echo "    -v \$(pwd)/output:/output \\"
echo "    ${IMAGE} \\"
echo "    STAR-Fusion \\"
echo "      --genome_lib_dir /ctat/GRCh38_gencode_v37_CTAT_lib \\"
echo "      --chimeric_junction /output/Chimeric.out.junction \\"
echo "      --output_dir /output/star_fusion_out"

echo ""
echo "=== STAR-Fusion 容器测试完成 ==="
