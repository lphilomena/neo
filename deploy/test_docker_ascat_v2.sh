#!/bin/bash
# ============================================================
# Docker 容器测试: ASCAT 2.5.2 (等位基因特异性拷贝数与纯度, R 4.0)
# 镜像: quay.io/biocontainers/ascat:2.5.2--r40hdfd78af_3
# 参考: tests/test_facets_lohhla.py (CNV/purity parsing)
#       环境: conda_envs/env.neoag-ascat.yml (ascat=2.5.2)
# ============================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
IMAGE="quay.io/biocontainers/ascat:2.5.2--r40hdfd78af_3"

echo "=== ASCAT 2.5.2 容器测试 ==="
echo "镜像: $IMAGE"
echo ""

# ---- 验证命令可用性 ----
echo ">>> R 环境..."
docker run --rm "${IMAGE}" R --version 2>&1 | head -5

echo ""
echo ">>> 加载 ASCAT R 包..."
docker run --rm "${IMAGE}" R -e 'if (require("ASCAT")) { cat("ASCAT loaded, version=", as.character(packageVersion("ASCAT")), "\n") } else { cat("ASCAT NOT found\n") }' 2>&1

echo ""
echo ">>> 检查 ASCAT 及依赖 R 包..."
docker run --rm "${IMAGE}" R -e '
  cat("R packages available:\n");
  for(p in c("ASCAT","data.table","foreach","doParallel","RColorBrewer")) {
    v <- tryCatch(as.character(packageVersion(p)), error=function(e) "NOT_AVAILABLE");
    cat(sprintf("  %-20s -> %s\n", p, v));
  }
' 2>&1

# ---- Full invocation example ----
cat <<'EXAMPLE'

  # ASCAT full invocation example (requires Tumor/Normal LogR + BAF data):
  docker run --rm \
    -v $(pwd)/input:/input \
    -v $(pwd)/output:/output \
    ${IMAGE} \
    R -e '
      library(ASCAT);
      ascat.aspcf <- ascat.aspc(
        ascat.bc = ascat.loadData("/input/tumor.LogR", "/input/tumor.BAF", "/input/normal.LogR", "/input/normal.BAF"),
        ascat.gg = ascat.loadAlleleCounts("/input/tumor.alleleCounts", "/input/normal.alleleCounts")
      );
      ascat.output <- ascat.runAscat(ascat.aspcf);
      write.table(ascat.output[["segments"]], "/output/ascat_segments.tsv", sep="\t", quote=F, row.names=F)
    '
EXAMPLE

echo ""
echo "=== ASCAT 2.5.2 容器测试完成 ==="
