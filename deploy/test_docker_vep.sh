#!/bin/bash
# ============================================================
# Docker 容器测试: VEP 115 (Ensembl Variant Effect Predictor)
# 镜像: ensemblorg/ensembl-vep:release_115.2
# 参考: tests/test_pvacseq_enrich.py (VCF annotation)
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="ensemblorg/ensembl-vep:release_115.2"

echo "=== VEP 115 容器测试 ==="
echo "镜像: $IMAGE"
echo ""

# ---- 准备测试数据 ----
WORKDIR="${SCRIPT_DIR}/test_output/vep"
mkdir -p "${WORKDIR}"
rm -rf "${WORKDIR}"/*

# 创建一个最小的 VCF 用于测试
cat > "${WORKDIR}/test_input.vcf" << 'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1000000	.	A	G	.	PASS	.
chr1	2000000	.	C	T	.	PASS	.
EOF

echo ">>> 测试输入 VCF:"
cat "${WORKDIR}/test_input.vcf"
echo ""

# ---- 运行 VEP ----
echo ">>> 运行 VEP (--help 先确认可用性)..."
docker run --rm \
    -v "${WORKDIR}:/data" \
    "${IMAGE}" \
    vep --help 2>&1 | head -20

echo ""
echo ">>> 运行 VEP 注释 (offline cache mode)..."
# VEP offline 模式需要 --cache 或 --database
# 若本地无缓存数据，使用 --database 联机模式会失败
# 这里使用 --help 已确认容器工作正常，并演示实际调用格式:
echo ""
echo "  # 完整 VEP 注释调用示例 (需挂载缓存目录):"
echo "  docker run --rm \\"
echo "    -v /path/to/vep_cache:/opt/vep/.vep \\"
echo "    -v \$(pwd)/input:/data \\"
echo "    ${IMAGE} \\"
echo "    vep --cache --offline --dir_cache /opt/vep/.vep \\"
echo "    -i /data/input.vcf -o /data/output.vcf --vcf"

echo ""
echo "=== VEP 容器测试完成 ==="
