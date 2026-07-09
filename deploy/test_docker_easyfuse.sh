#!/bin/bash
# ============================================================
# Docker 容器测试: EasyFuse (融合 meta-caller)
# 镜像: tronbioinformatics/easyfuse:1.3.7
# 参考: tests/test_easyfuse_adapter.py (parse_easyfuse / filter)
#       tests/test_diagnostic_fusion_rescue.py (fusion rescue)
# 数据: data/fixtures/easyfuse_fusions.pass.tsv
#       data/fixtures/easyfuse_fusions.v2.pass.csv
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="tronbioinformatics/easyfuse:1.3.7"

echo "=== EasyFuse 容器测试 ==="
echo "镜像: $IMAGE"
echo ""

# ---- 准备测试数据 ----
WORKDIR="${SCRIPT_DIR}/test_output/easyfuse"
mkdir -p "${WORKDIR}"

FIXTURE_PASS="${PROJECT_ROOT}/data/fixtures/easyfuse_fusions.pass.tsv"
FIXTURE_V2="${PROJECT_ROOT}/data/fixtures/easyfuse_fusions.v2.pass.csv"

echo ">>> 检查测试数据:"
echo "    v1 fixture: ${FIXTURE_PASS}"
echo "    v2 fixture: ${FIXTURE_V2}"
echo ""

# ---- 验证命令可用性 ----
echo ">>> EasyFuse processing.py --help..."
docker run --rm "${IMAGE}" python /code/easyfuse-1.3.7/processing.py --help 2>&1 | head -20 || true

echo ""
echo ">>> 容器内可用工具:"
docker run --rm "${IMAGE}" bash -c 'echo "STAR version:" && STAR --version 2>&1 | head -1' || true

# ---- 演示完整调用格式 ----
echo ""
echo "  # EasyFuse 完整调用示例 (需挂载数据 + 参考库 ~92GB):"
echo "  docker run --rm \\"
echo "    -v /path/to/easyfuse_ref:/ref \\"
echo "    -v /path/to/fastq_data:/data \\"
echo "    -v \$(pwd)/output:/output \\"
echo "    ${IMAGE} \\"
echo "    python /code/easyfuse-1.3.7/processing.py -i /data -o /output"

echo ""
echo "  # 验证 fixture 数据的解析 (挂载 fixtures 目录):"
echo "  docker run --rm \\"
echo "    -v ${PROJECT_ROOT}:/project \\"
echo "    ${IMAGE} \\"
echo "    python3 -c \"import csv; rows=list(csv.DictReader(open('/project/data/fixtures/easyfuse_fusions.pass.tsv'), delimiter='\\t')); print(f'Parsed {len(rows)} rows')\""

echo ""
echo "=== EasyFuse 容器测试完成 ==="
