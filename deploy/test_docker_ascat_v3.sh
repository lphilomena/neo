#!/bin/bash
# ============================================================
# Docker 容器测试: ASCAT 3.2.0 (等位基因特异性拷贝数与纯度, R 4.4)
# 镜像: quay.io/biocontainers/ascat:3.2.0--r44hdfd78af_1
# 参考: tests/test_facets_lohhla.py (CNV/purity parsing)
#       环境: conda_envs/env.neoag-ascat-v3.yml (ascat=3.2.0)
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="quay.io/biocontainers/ascat:3.2.0--r44hdfd78af_1"

echo "=== ASCAT 3.2.0 容器测试 ==="
echo "镜像: $IMAGE"
echo ""

# ---- 验证命令可用性 ----
echo ">>> R 环境..."
docker run --rm "${IMAGE}" R --version 2>&1 | head -5

echo ""
echo ">>> 加载 ASCAT R 包..."
docker run --rm "${IMAGE}" R -e 'if (require("ASCAT")) { cat("ASCAT loaded, version=", as.character(packageVersion("ASCAT")), "\n") } else { cat("ASCAT NOT found\n") }' 2>&1

echo ""
echo ">>> 与 ASCAT v2.5.2 对比 — v3 新增功能检查:"
docker run --rm "${IMAGE}" R -e '
  library(ASCAT);
  cat("Available functions (sample):\n");
  funcs <- ls("package:ASCAT");
  v3_funcs <- grep("ascat.run|ascat.aspc|ascat.predict|ascat.mult", funcs, value=TRUE);
  for(f in v3_funcs) cat("  ", f, "\n");
' 2>&1

# ---- 演示完整调用格式 ----
echo ""
echo "  # ASCAT v3 完整调用示例 (需 Tumor/Normal LogR + BAF 数据):"
echo "  docker run --rm \\"
echo "    -v \$(pwd)/input:/input \\"
echo "    -v \$(pwd)/output:/output \\"
echo "    ${IMAGE} \\"
echo "    R -e '"
echo "      library(ASCAT);"
echo "      ascat.bc <- ascat.loadData(\\"/input/tumor.LogR\\", \\"/input/tumor.BAF\\", \\"/input/normal.LogR\\", \\"/input/normal.BAF\\");"
echo "      ascat.gg <- ascat.loadAlleleCounts(\\"/input/tumor.alleleCounts\\", \\"/input/normal.alleleCounts\\");"
echo "      # v3 特有: 多样本模式与预测功能"
echo "      ascat.aspcf <- ascat.aspc(ascat.bc, ascat.gg);"
echo "      ascat.output <- ascat.runAscat(ascat.aspcf, gamma=1);"
echo "      write.table(ascat.output[[\"segments\"]], \\"/output/ascat_v3_segments.tsv\\", sep=\\\"\\\\t\\\", quote=F, row.names=F)"
echo "    '"

echo ""
echo "=== ASCAT 3.2.0 容器测试完成 ==="
