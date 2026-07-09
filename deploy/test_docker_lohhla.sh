#!/bin/bash
# ============================================================
# Docker 容器测试: LOHHLA (HLA 等位基因特异性 LOH)
# 镜像: quay.io/biocontainers/lohhla:20171108--hdfd78af_3
# 参考: tests/test_facets_lohhla.py (parse_lohhla_prediction)
# 数据: tools/lohhla/example-file/correct-example-out/
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="quay.io/biocontainers/lohhla:20171108--hdfd78af_3"
LOHHLA_EXAMPLE="${PROJECT_ROOT}/tools/lohhla/example-file/correct-example-out/example.10.DNA.HLAlossPrediction_CI.xls"

echo "=== LOHHLA 容器测试 ==="
echo "镜像: $IMAGE"
echo ""

# ---- 验证命令可用性 ----
# 注意: biocontainers lohhla 镜像的 activate 脚本有语法问题,
# 需 --entrypoint 跳过 conda activate
echo ">>> 容器内 R 环境 (跳过 conda activate)..."
docker run --rm --entrypoint="" "${IMAGE}" R --version 2>&1 | head -5

echo ""
echo ">>> LOHHLA 关键脚本..."
docker run --rm --entrypoint="" "${IMAGE}" bash -c 'ls /usr/local/bin/LOHHLAscript.R 2>/dev/null || find / -name "LOHHLAscript.R" -type f 2>/dev/null | head -5 || echo "LOHHLAscript.R not found in container"'

echo ""
echo ">>> 检查示例输出文件:"
if [[ -f "${LOHHLA_EXAMPLE}" ]]; then
    echo "    ${LOHHLA_EXAMPLE}"
    echo "    前5行:"
    head -5 "${LOHHLA_EXAMPLE}"
else
    echo "    示例文件不存在: ${LOHHLA_EXAMPLE}"
fi

# ---- 演示完整调用格式 ----
echo ""
echo "  # LOHHLA 完整调用示例 (需 Tumor/Normal BAM + HLA 等位基因):"
echo "  docker run --rm \\"
echo "    -v \$(pwd)/input:/input \\"
echo "    -v \$(pwd)/output:/output \\"
echo "    -v \$(pwd)/ref:/ref \\"
echo "    ${IMAGE} \\"
echo "    LOHHLAscript.R \\"
echo "      --patientId PATIENT1 \\"
echo "      --outputDir /output \\"
echo "      --normalBAMfile /input/normal.bam \\"
echo "      --tumorBAMfile /input/tumor.bam \\"
echo "      --HLAfastaLoc /ref/hla_gen.fasta \\"
echo "      --HLAalleles 'HLA-A*01:01,HLA-A*02:01,HLA-B*07:02' \\"
echo "      --gatkDir /usr/local/bin"

echo ""
echo "=== LOHHLA 容器测试完成 ==="
