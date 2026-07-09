# utility相关
 assets/ 目录的作用

  这个目录下的 9 个文件都是空的 TSV 占位文件，它们仅仅包含了表头（列名），没有实际数据行。这些文件在 NeoAg（肿瘤新抗原）分析流水线中充当可选输入文件的默认占位符。

  各文件含义

  ┌──────────────────────────────┬────────────────────┬───────────────────────────────────────────────────────────────────┐
  │             文件             │    对应数据类型    │                               列名                                │
  ├──────────────────────────────┼────────────────────┼───────────────────────────────────────────────────────────────────┤
  │ empty_vep.tsv                │ VEP 注释结果       │ gene, consequence                                                 │
  ├──────────────────────────────┼────────────────────┼───────────────────────────────────────────────────────────────────┤
  │ empty_cnv.tsv                │ 拷贝数变异         │ chrom, start, end, total_cn                                       │
  ├──────────────────────────────┼────────────────────┼───────────────────────────────────────────────────────────────────┤
  │ empty_expression.tsv         │ 基因表达量         │ gene, TPM                                                         │
  ├──────────────────────────────┼────────────────────┼───────────────────────────────────────────────────────────────────┤
  │ empty_purity.tsv             │ 肿瘤纯度           │ sample_id, purity                                                 │
  ├──────────────────────────────┼────────────────────┼───────────────────────────────────────────────────────────────────┤
  │ empty_normal_hla_ligands.tsv │ 正常组织 HLA 配体  │ peptide, source_tissue, hla_allele                                │
  ├──────────────────────────────┼────────────────────┼───────────────────────────────────────────────────────────────────┤
  │ empty_mhcflurry.tsv          │ MHCflurry 结合预测 │ peptide, allele, affinity, ...                                    │
  ├──────────────────────────────┼────────────────────┼───────────────────────────────────────────────────────────────────┤
  │ empty_normal_expression.tsv  │ 正常组织表达       │ gene, normal_tissue_max_tpm, normal_hspc_tpm, critical_tissue_hit │
  ├──────────────────────────────┼────────────────────┼───────────────────────────────────────────────────────────────────┤
  │ empty_netmhcpan.tsv          │ NetMHCpan 结合预测 │ Peptide, HLA, Score_EL, %Rank_EL, ...                             │
  ├──────────────────────────────┼────────────────────┼───────────────────────────────────────────────────────────────────┤
  │ empty_hla_loh.tsv            │ HLA 杂合性缺失     │ hla_allele, loh_status                                            │
  └──────────────────────────────┴────────────────────┴───────────────────────────────────────────────────────────────────┘

  使用机制

  在 src/neoag_v03/cli.py:518-522 等多处可以看到这样的代码：

  vep_appm=_sv_optional_path(args.vep_appm) or ROOT / "assets" / "empty_vep.tsv"
  hla_loh=_sv_optional_path(args.hla_loh) or ROOT / "assets" / "empty_hla_loh.tsv"
  purity=_sv_optional_path(args.purity) or ROOT / "assets" / "empty_purity.tsv"
  cnv=_sv_optional_path(args.cnv) or ROOT / "assets" / "empty_cnv.tsv"

  当用户未提供某个可选的输入文件时（命令行参数为空），流水线不会崩溃报错，而是自动使用这些仅含表头的空文件作为 fallback。这样下游代码可以正常读取文件（格式正确），但由于没有数据行，实际不会产生任何匹配结果，实现了"跳过该分析"的效果。


# 流程相关 

与 Nextflow 工作流存在紧密的配合关系，进一步关联到modules/。

三层架构总结：
  profiles/           ← 参数配置层 (TOML): 权重、阈值、策略
      │                 定义"怎么评分"
      │
  workflows/          ← 流程编排层 (Nextflow): 步骤顺序、并行、依赖
      │                 定义"先做什么、后做什么"
      │
  modules/            ← 执行模块层 (Nextflow + Python): 调用具体工具
                        定义"每个步骤怎么做"

  profiles/ 和 workflows/ 是正交设计，profiles/ 目录中的 TOML 配置文件是整个流水线的"参数大脑"，而 workflows/ 中的 Nextflow 文件是"执行骨架"，二者通过 profile_name 这个关键参数串联起来：
  - Profile 控制评分逻辑的参数（权重、阈值、免疫原性工具选择等），不改变流程结构
  - Workflow 控制执行步骤的顺序和依赖关系（先 Mutect2 再 pVAC 再评分），不关心具体参数值
  - 二者通过 --profile_name 解耦，同一个 workflow 可以搭配不同 profile 适应不同癌种

  ---
  配合机制：
  用户指定 --profile_name leukemia
           │
           ▼
  ┌─────────────────────────────────────────────┐
  │  Nextflow workflow (执行编排层)              │
  │  workflows/main.nf                          │
  │    params.profile_name = 'leukemia'         │
  │    → 传给 NEOAG_V03_RC                     │
  │         │                                    │
  │         ▼                                    │
  │  neoag_v03_rc.nf (子工作流)                 │
  │    → 将 profile_name 分发给所有模块          │
  │      PARSE_PVAC(sample_id, profile_name)     │
  │      APPM_2(sample_id, profile_name)         │
  │      CCF_2(events, profile_name)             │
  │      PEPTIDE_SAFETY(..., profile_name)       │
  │      SCORE_V041(..., profile_name)           │
  │      IMMUNE_ESCAPE(..., profile_name)        │
  │         │                                    │
  │         ▼                                    │
  │  ┌──────────────────────────────────────┐    │
  │  │  Python 模块层 (执行层)              │    │
  │  │  --profile 'leukemia'                │    │
  │  │  → load_profile("leukemia")          │    │
  │  │  → profiles/leukemia.toml            │    │
  │  │  → 返回 dict 控制评分/权重/阈值      │    │
  │  └──────────────────────────────────────┘    │
  └─────────────────────────────────────────────┘

  Nextflow 本身不解析 TOML，它只负责：
  1. 接收用户传入的 --profile_name
  2. 将其作为字符串参数透传给每个 Python CLI 模块
  3. Python 模块内部调用 load_profile(profile_name) 从 profiles/ 目录加载实际配置

  ---
  工作流与 Profile 的具体对应关系：
  ┌────────────────────┬─────────────────┬───────────────────────────────────────────┐
  │  Nextflow 工作流   │  默认 profile   │                   场景                    │
  ├────────────────────┼─────────────────┼───────────────────────────────────────────┤
  │ main.nf            │ default         │ 通用实体瘤 SNV+InDel+SV 分析              │
  ├────────────────────┼─────────────────┼───────────────────────────────────────────┤
  │ main_full.nf       │ leukemia        │ 白血病端到端全流程                        │
  ├────────────────────┼─────────────────┼───────────────────────────────────────────┤
  │ snv_phase1_wes.nf  │ default         │ WES SNV Phase1: BAM→Mutect2→pVAC→评分     │
  ├────────────────────┼─────────────────┼───────────────────────────────────────────┤
  │ sv_phase1_wgs.nf   │ sv_wgs_phase1   │ WGS SV 检测+评分                          │
  ├────────────────────┼─────────────────┼───────────────────────────────────────────┤
  │ sv_phase1_5_wes.nf │ sv_wes_phase1_5 │ WES SV 捕获感知评分                       │
  ├────────────────────┼─────────────────┼───────────────────────────────────────────┤
  │ neoag_v03_rc.nf    │ (由上层传入)    │ 核心评分子流程（被其他 workflow include） │
  └────────────────────┴─────────────────┴───────────────────────────────────────────┘

## profiles 详解
 profiles/ 目录存放的是预定义的参数配置文件（profiles），使用 TOML 格式。每个 profile 对应一种肿瘤类型或分析场景，控制 NeoAg 新抗原预测流水线的各个评分环节的权重、阈值和策略。配置文件由 src/neoag_v03/config.py 中的 load_profile() 加载，在 CLI 命令层、pipeline
  层以及各个评分模块中被广泛使用。

  加载机制

  # src/neoag_v03/config.py:8-18
  def load_profile(profile: str | Path) -> dict:
      # 1. 如果是完整路径，直接读取
      # 2. 否则在 profiles/ 目录下查找 "<profile>.toml"

  用户通过命令行参数 --profile <name> 指定使用哪个 profile（如 --profile sarcoma），系统自动从 profiles/ 目录加载对应的 .toml 文件。

  ---
  各 Profile 文件详解

  1. default.toml — 默认通用配置

  最完整的配置模板，适用于通用实体瘤 WGS/WES 分析场景。

  关键设定：
  - 突变来源：SNV、InDel、SV、CNV、Fusion、Splice、TAA（全覆盖），但不强制执行（enforce_enabled_sources = false，保留向后兼容）
  - L3 统一免疫学评分权重：这是最完整的权重体系，HLA 呈递 (hla_presentation=0.25) 和 HLA 结合 (hla_binding=0.20) 占比最高，免疫原性占 0.15
  - 免疫原性：三源加权组合 — BigMHC_IM (0.45) + PRIME (0.35) + IEDB (0.20)
  - Gate 过滤：要求呈递等级 A 或 B，%EL Rank ≤ 2.0，StabPan ≥ 1.4
  - 安全阈值：关键组织 TPM > 10 硬排除，HSPC TPM > 8 失败

  2. leukemia.toml — 白血病专用

  针对血液肿瘤（白血病） 的优化配置。

  与 default 的关键差异：
  - 移除 CNV 和 SV：白血病中这些变异类型不适用，只启用 SNV、InDel、Fusion、Splice、TAA
  - Fusion 优先级大幅提升：1.35（vs default 的 1.15），因为融合基因在白血病中更为常见和重要
  - TAA 优先级降低：0.25（vs default 0.50）
  - HSPC 安全阈值大幅收紧：HSPC（造血干祖细胞）失败阈值从 8.0 TPM 降到 2.0，警告阈值从 2.0 降到 0.5 — 这是因为白血病治疗中保护正常造血干细胞至关重要
  - 高表达阈值降低：10.0 TPM（vs default 20.0）
  - 免疫原性：仅用 PRIME + BigMHC_IM，均权 mean 组合，不使用 IEDB
  - 无 L3 统一权重、无 APM 相关配置（精简）

  3. sarcoma.toml — 肉瘤专用

  针对肉瘤的优化配置。

  与 default 的关键差异：
  - Fusion 优先级大幅提升：1.55（所有 profile 中最高），因为肉瘤常由特异性融合基因驱动
  - TAA 优先级降低：0.35
  - 高表达阈值：15.0 TPM
  - 免疫原性：与 leukemia 一致，仅用 PRIME + BigMHC_IM 均权组合
  - 无 L3 统一权重（精简）

  4. sv_wgs_phase1.toml — SV WGS Phase 1

  针对结构变异（SV）全基因组测序（WGS） 的分析场景。

  关键设定：
  - 仅启用 SV 子类型：SV_Fusion、SV_Frameshift、SV_Junction、SV_Insertion
  - SV Phase 1 特有的过滤器（[sv_phase1]）：
    - 断点合并距离 200bp
    - 允许 Tier2 事件
    - 最小肿瘤支持 reads：SR ≥ 2, PE ≥ 4
    - 正常样本严格过滤：SR = 0, PE ≤ 1
    - 要求编码区影响
  - 事件权重偏向置信度和表达：event_confidence=0.28, event_expression=0.22
  - Gate 放宽：允许 C_BINDING_ONLY 等级通过，StabPan 阈值设为 0.0（不限制）
  - 免疫原性关闭：enabled = false，这是纯结构变异检测阶段，不需要免疫原性评分

  5. sv_wes_phase1_5.toml — SV WES Phase 1.5

  针对结构变异（SV）外显子组测序（WES） 的分析场景。

  与 sv_wgs_phase1 的关键差异：
  - WES 特有的置信度上限（[wes_confidence_caps]）：根据捕获覆盖范围对事件分级设上限
    - Tier1 最高为 B，Tier2 为 B_CAUTION，Tier3 为 C
    - DNA-only 事件上限 C，有 RNA 支持可提升到 B
  - WES 过滤器（[wes_filter]）：
    - 要求断点接近捕获靶区（250bp 内）
    - 调用区域扩展 1000bp
    - Tier1 必须有 RNA 支持（require_rna_for_tier1 = true）
    - RNA 连接 reads ≥ 3 才作为 Tier1
  - [sv_wes_phase1_5] 特有段：min_rna_junction_reads_tier1 = 3
  - Source 优先级整体下调：SV_Fusion 从 1.20 降为 0.85，SV_Frameshift 从 1.15 降为 0.90（反映外显子组捕获对 SV 检测的不确定性）

  6. immunogenicity_extended.toml — 扩展免疫原性

  这是一个部分覆盖 profile，仅覆盖免疫原性和评分相关的 section，其余从 default 继承。

  独有设定：
  - 三种免疫原性工具：PRIME + BigMHC_IM + DeepImmuno
  - 权重偏向 BigMHC_IM：BigMHC_IM (0.55) > PRIME (0.35) > DeepImmuno (0.10)
  - DeepImmuno（基于 CNN 的免疫原性预测器）作为额外的补充信号

  ---
  设计模式总结

  ┌─────────────────────────┬────────────┬─────────────────────────────────────────────┐
  │         Profile         │ 癌种/场景  │                  核心特点                   │
  ├─────────────────────────┼────────────┼─────────────────────────────────────────────┤
  │ default                 │ 通用实体瘤 │ 最完整配置，7种突变来源全覆盖，三源免疫原性 │
  ├─────────────────────────┼────────────┼─────────────────────────────────────────────┤
  │ leukemia                │ 白血病     │ 去CNV/SV，升Fusion，收紧HSPC阈值保护造血    │
  ├─────────────────────────┼────────────┼─────────────────────────────────────────────┤
  │ sarcoma                 │ 肉瘤       │ Fusion权重最高(1.55)，反映融合驱动特征      │
  ├─────────────────────────┼────────────┼─────────────────────────────────────────────┤
  │ sv_wgs_phase1           │ SV WGS     │ 仅SV子类型，严格过滤，关闭免疫原性          │
  ├─────────────────────────┼────────────┼─────────────────────────────────────────────┤
  │ sv_wes_phase1_5         │ SV WES     │ WES捕获感知，RNA支持分级，优先级下调        │
  ├─────────────────────────┼────────────┼─────────────────────────────────────────────┤
  │ immunogenicity_extended │ 扩展免疫   │ 三工具免疫原性，DeepImmuno补充信号          │
  └─────────────────────────┴────────────┴─────────────────────────────────────────────┘

  这是一种典型的策略模式：通过切换 profile 文件，同一套流水线代码可以适配不同肿瘤类型和分析场景，而无需修改任何 Python 代码。Profile 之间采用浅层合并（TOML section 级别覆盖），未在子 profile 中定义的参数自然继承 default.toml 的值。


# 流程执行程序
## 主脚本
runner.py  （src/neoag_v03/tools/runner.py）的结构如下：

  总体架构

  ┌─────────────────────────────────────────────────┐
  │  RUNNERS (字典, L1086-1099)                       │
  │  "tool_name" → lambda(ctx, path): run_xxx(...)   │
  │  供外部通过 run_tool(name, ctx, output) 调度      │
  └─────────────────────────────────────────────────┘
            │
            │ 每个 runner 内部调用 ctx.exe("tool_name")
            │ 得到可执行文件路径，再用 subprocess.run() 执行
            ▼
  ┌─────────────────────────────────────────────────┐
  │  核心函数 _run_cmd(cmd, workdir)                  │
  │  L173: subprocess.run(cmd, cwd=workdir, ...)     │
  │  所有工具的底层执行都经过这里                      │
  └─────────────────────────────────────────────────┘

每个 run_xxx() 函数都遵循同样的三步模式：

  def run_xxx(ctx: RunContext, out_path: Path) -> Path:
      # 1. stub 模式 → 直接复制 fixture 文件
      if ctx.stub:
          _stub_copy(spec.fixture_outputs[...], out_path)
          return out_path

      # 2. 验证输入数据
      if not ctx.xxx:
          raise ValueError("xxx requires xxx")

      # 3. 组装命令 + subprocess 执行
      cmd = [ctx.exe("xxx"), ...args...]
      _run_cmd(cmd, workdir)
      return out_path

  当前所有工具的调用路径是：
  run_tool(name, ctx, output)
    → ctx.exe(name)        # 返回 conda 环境的可执行文件路径
    → _run_cmd(cmd, work)  # subprocess.run() 本地执行


