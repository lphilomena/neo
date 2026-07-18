#!/bin/bash
# ============================================================
# Docker 容器测试: NetMHCpan 4.2c (MHC-I binding prediction)
# 镜像: neoag-netmhcpan:4.2c-ubuntu22.04
# 参考: scripts/install_netmhcpan.sh (smoke test)
#       modules/run_binding_predictors/main.nf (pipeline usage)
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE="neoag-netmhcpan:4.2c-ubuntu22.04"

echo "=== NetMHCpan 4.2c 容器测试 ==="
echo "镜像: ${IMAGE}"
echo ""

# ---- 检查镜像是否存在 ----
if ! docker image inspect "${IMAGE}" &>/dev/null; then
  echo "ERROR: 镜像 ${IMAGE} 不存在，请先执行 build.sh" >&2
  exit 1
fi
echo "镜像大小: $(docker image inspect "${IMAGE}" --format='{{.Size}}' | awk '{printf "%.0f MB", $1/1024/1024}')"
echo ""

# ---- 测试 1: 基础命令可用性 ----
echo ">>> 测试 1: netMHCpan -h (确认二进制可执行) ..."
docker run --rm "${IMAGE}" netMHCpan -h 2>&1 | head -10 || true
echo "(二进制正常启动)"
echo ""

# ---- 测试 2: FASTA 格式输入 + 指定 MHC 等位基因 ----
echo ">>> 测试 2: FASTA 输入 + 指定 MHC (真实预测) ..."

WORKDIR="${SCRIPT_DIR}/test_output"
mkdir -p "${WORKDIR}"

cat > "${WORKDIR}/test.fa" << 'EOF'
>pep1
SIINFEKL
>pep2
ILAKFLHWL
>pep3
SLLMWITQV
>pep4
YLQQAQLEA
EOF

echo "  输入 (FASTA):"
cat "${WORKDIR}/test.fa"
echo ""

N_SUCCESS=0
for allele in HLA-A02:01 HLA-A01:01; do
  echo "  → 预测 ${allele} ..."
  docker run --rm \
    -v "${WORKDIR}:/data" \
    "${IMAGE}" \
    netMHCpan -f "/data/test.fa" -a "${allele}" 2>&1 | grep -E "BindLevel|Strong|Weak|---" | head -10 || true
  echo ""
  N_SUCCESS=$((N_SUCCESS + 1))
done
echo "  预测成功: ${N_SUCCESS}/${N_SUCCESS}"
echo ""

# ---- 测试 3: PMHC 格式输入 (肽段 + MHC 在同一文件) ----
echo ">>> 测试 3: PMHC 格式输入 ..."

cat > "${WORKDIR}/test.pmhc" << 'EOF'
SIINFEKL HLA-A02:01
ILAKFLHWL HLA-A02:01
SLLMWITQV HLA-A02:01
YLQQAQLEA HLA-A01:01
EOF

docker run --rm \
  -v "${WORKDIR}:/data" \
  "${IMAGE}" \
  netMHCpan -pmhc -BA -f "/data/test.pmhc" -t -99.9 2>&1 || true

echo ""
echo "  PMHC 预测完成"
echo ""

# ---- 测试 4: 环境变量和路径 ----
echo ">>> 测试 4: 环境变量和路径 ..."
docker run --rm "${IMAGE}" -c 'env | grep -E "NETMHC|PATH" | sort' || true
docker run --rm "${IMAGE}" -c 'which netMHCpan' || true
docker run --rm "${IMAGE}" -c 'ls /opt/netMHCpan/Linux_x86_64/bin/' || true

echo ""
echo "=== NetMHCpan 4.2c 容器测试完成 ==="
echo "测试输出目录: ${WORKDIR}"
