# Release Boundary And Online Release Guide

本文件合并了原先的 `docs/ONLINE_RELEASE.md` 与 `docs/RELEASE_BOUNDARY_V04.md`（内容高度重叠，
两者都在描述"轻量在线发布包里包含/排除什么"），原两份文档已删除。

## 打包边界

**包含在轻量在线发布包中：**

- `src/` 下的Python源码（`src/neoag_v03/`）
- `bin/` 下的CLI/包装脚本
- `profiles/`、`conf/`（示例配置）
- `modules/`、`workflows/`（Nextflow模块/工作流）
- `tests/`
- 小型fixture：`data/fixtures`、`data/fixtures_snv`、`data/fixtures_sv`、`data/improve`、`assets/`、`resources/`
- Markdown文档（`docs/`）、license、notice、citation、release notes、manifest
- `scripts/` 下的打包/复现脚本
- conda环境manifest（如`conda/env.neoag-ascat-v3.yml`）
- `.agents/skills/`（AI可路由的skill定义）

**排除在轻量在线发布包外：**

- `.git`、虚拟环境（`.venv/`、`.venv.local/`）、缓存、字节码、`.pytest_cache/`
- `/tools/`（外部工具二进制；**注意`.gitignore`里这条必须锚定根目录写成`/tools/`，不能写成`tools/`，
  否则会连带排除`src/neoag_v03/tools/`这个真实Python包**）
- `results/`、`work/`、`dist/`、`conda_packs/`
- `.nextflow*`元数据和日志
- 大型参考数据：`data/ref`、`data/vep`、`data/external`、`data/examples`
- 本地/私有部署文件：`conf/tools.env.local.sh`、`conf/site.config`、`conf/private/*`、
  `conf/*.private.toml`、`conf/*.local.toml`
- 患者专属脚本、样本标识符，以及**站点本地绝对路径**（含真实用户名/内部集群挂载路径的文档或脚本不应进入发布包，
  参见下方"站点本地路径"说明）
- 可编辑的office产物，如`.docx`/`.pptx`

## 站点本地路径（重要）

此前 `docs/PROJECT_DATA_PATHS.md` 里直接写了内部集群的绝对路径和用户名（如
`/mnt/.../peixunban/gl/liup/neodata`、`/home/na/ref/...`），这违反了上面"排除站点本地绝对路径"的边界，
本轮已删除该文件。真实的数据/参考路径请写在**不提交到仓库**的本地文件里
（`conf/tools.env.local.sh`、`conf/site.config`，参见 [`SITE_CONFIG_BOUNDARY.md`](SITE_CONFIG_BOUNDARY.md)），
仓库里只保留形如`$NEOAG_TOOLS_ROOT`、`$NEOAG_REFERENCE_FASTA`这样的环境变量占位。

## 可选/外部依赖

部分工具集成对轻量发布来说是可选的：代码应该在这些工具缺失时记录为"缺失"而不是让发布测试失败。
当前视workflow模式而定的可选/外部工具包括NetMHCpan、PRIME、STAR-Fusion、FusionCatcher等。

## 打包

```bash
python scripts/package_online_release.py --outdir work/releases
```

产出：

- `neoag_event_pipeline_v043_online_<date>.tar.gz`
- `neoag_event_pipeline_v043_online_<date>.tar.gz.sha256`（外部校验和）
- `neoag_event_pipeline_v043_online_<date>.manifest.json`

## 冒烟测试

```bash
tmpdir=$(mktemp -d)
tar -xzf work/releases/neoag_event_pipeline_v043_online_<date>.tar.gz -C "$tmpdir"
cd "$tmpdir"/neoag_event_pipeline_v043_online_<date>
python -m pip install -e '.[test]'
pytest -q
neoag-v03 run-demo --entry-mode snv_indel --outdir work/demo_v043 --sample-id DEMO001
```

Nextflow冒烟（在线包含轻量`bin/nextflow`launcher，不含运行时依赖缓存，首次使用可能需要联网下载到`NXF_HOME`）：

```bash
bin/neoag-nextflow run workflows/main.nf \
  -w /tmp/neoag_nf_work \
  --pvac_files data/fixtures/pvacseq_aggregated.tsv \
  --outdir /tmp/neoag_nf_demo \
  --sample_id NF_DEMO
```

## 部署提示

完整环境/工具/参考数据配置见 [`INSTALL_AND_DATA.md`](INSTALL_AND_DATA.md)。真实数据运行前设置
`NEOAG_TOOLS_ROOT`指向共享工具安装目录。fixture demo和默认测试套件不需要licensed/重量级工具。

```bash
export NEOAG_TOOLS_ROOT=/path/to/neoag_artifacts
source conf/tools.env.sh
neoag-v03 check-tools          # 全量核对，部署/CI阶段用
neoag-v03 run-demo --entry-mode snv_indel --outdir /tmp/demo --sample-id DEMO   # 按入口验收，日常开发用
```

## 本次发布的验证项

在线发布至少应通过：

- `pytest -q`
- `run-demo --entry-mode <每个入口> --outdir ... --sample-id ...`（`<每个入口>`依次替换为
  `snv_indel`/`fusion`/`splice_junction`/`sv_wgs`/`sv_wes`/`peptide_only`，例如
  `neoag-v03 run-demo --entry-mode snv_indel --outdir work/demo --sample-id DEMO`）
- `bin/neoag-nextflow run workflows/main.nf`（轻量fixture输入）
- `scripts/check_release_boundary.sh`

## 临床使用边界

本发布仍然是计算原型：可以对新抗原候选、安全性标志、免疫逃逸假设进行排序和注释，但不产出临床治疗建议。
