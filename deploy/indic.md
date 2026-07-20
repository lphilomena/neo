  脚本设计要点

  每个脚本遵循统一结构：
  1. 镜像版本声明 — 与 pull_images.sh 保持一致
  2. 测试数据检查 — 引用 data/fixtures*/ 中的现有测试数据
  3. 命令可用性验证 — docker run --rm <image> <tool> --help
  4. 功能测试 (如数据可用) — samtools/tabix 直接处理 mini_somatic.vcf
  5. 完整调用示例 — 注释化的 docker run 生产级调用格式
  6. 测试数据安全 — 输出写入 deploy/test_output/<tool>/，不污染 fixtures

# docker 测试运行方式
``` sh
  bash deploy/test_docker_samtools.sh     # 单个
  for f in deploy/test_docker_*.sh; do bash "$f"; done  # 全部
```

全流程2种使用方式：

  # Conda 模式（默认，行为不变）
  neoag run-full --config run.toml --outdir output/

  # Docker 模式（有镜像的用容器，无镜像的自动用 conda）
  NEOAG_RUNNER_MODE=docker neoag run-full --config run.toml --outdir output/
