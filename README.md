# 鸡卵巢颗粒细胞单细胞测序服务器流程

本目录用于把论文复现脚本整理成“正式项目可复用”的服务器流程。

这套流程的方法学来源仍然是论文 *Single-cell RNA sequencing uncovers cellular heterogeneity of granulosa cells and provides a signature for follicular development in chicken*，但工程目标已经改成：

- 论文数据可以作为验证入口
- 你自己的新项目可以从原始数据开始直接复用
- 固定部分和可变部分分开
- 环境拆开，避免主流程、velocity、SCENIC 互相污染

## 0. 当前规范入口

从这一版开始，workflow 的规范执行顺序已经重构为：

- `workflow/00_init_project.sh`
- `workflow/01_register_delivery.sh`
- `workflow/02_validate_metadata.sh`
- `workflow/03_audit_inputs.sh`
- `workflow/04_standardize_inputs.sh`
- `workflow/05_input_summary.sh`
- `workflow/10_run_cellranger_from_fastq.sh`
- `workflow/20_run_main_pipeline.sh`
- `workflow/30_run_velocyto.sh`
- `workflow/31_run_scvelo.sh`
- `workflow/40_build_ortholog_cache.sh`
- `workflow/41_download_scenic_resources.sh`
- `workflow/42_run_scenic.sh`

旧编号脚本仍然保留，但仅作为兼容 wrapper，不再是推荐入口。

## 1. 两种运行模式

### 1.1 论文验证模式

适用于作者已经提供了聚合后的 `10X matrix`，例如：

- `HF_GCs`
- `PHF_GCs`

这时主流程直接从 `data/<sample>` 读矩阵，RNA velocity 仍然单独回到原始 `cellranger_out/<SRR>`。

这就是为什么论文验证模式里会出现：

- `SAMPLE_NAMES=HF_GCs,PHF_GCs`
- `RAW_SAMPLES=SRR...`

二者不相同。

### 1.2 正式项目模式

适用于你自己的真实项目。

推荐起点：

- `FASTQ -> Cell Ranger -> matrix -> 主流程`

也支持：

- 已有 `cellranger_out`
- 已有 `10X matrix`

正式项目模式下，通常应当让：

- `SAMPLE_NAMES`
- `RAW_SAMPLES`

保持一致，都是你自己的样本名，例如：

- `Ctrl_1`
- `Ctrl_2`
- `Treat_1`
- `Treat_2`

## 2. 正式项目的标准入口

新增入口脚本：

- `workflow/00_init_project.sh`

它的作用是：

1. 新建一个项目目录
2. 自动生成标准目录结构
3. 把外部 `FASTQ`、参考文件、已有 `cellranger_out`、已有 `matrix` 用软链接接进来
4. 自动生成项目专属配置 `config/project_config.sh`
5. 自动把全局 `config/server_config.sh` 切换到这个项目

也就是说，以后新项目不应该从手改 `config/server_config.sh` 开始，而应该先跑 `00_init_project.sh`。

## 3. 项目目录结构

初始化完成后，一个正式项目目录会包含：

- `fastq/`
- `reference/`
- `cellranger_out/`
- `data/`
- `results/`
- `logs/`
- `envs/`
- `resources/`
- `config/project_config.sh`

其中：

- `fastq/`：原始 FASTQ 输入
- `reference/`：参考基因组和 GTF
- `cellranger_out/`：Cell Ranger 输出
- `data/`：主流程输入矩阵目录
- `results/`：主流程、velocity、SCENIC 输出
- `config/project_config.sh`：该项目真正的配置文件

## 4. 变量对应关系

这是最容易混淆的部分。

### 4.1 主流程变量

- `SAMPLE_NAMES`
  - 主流程读取的样本目录名
  - 对应 `data/<sample>`

- `ANALYSIS_GROUP_1_NAME`
- `ANALYSIS_GROUP_1_SAMPLES`
- `ANALYSIS_GROUP_2_NAME`
- `ANALYSIS_GROUP_2_SAMPLES`
  - 主流程分组定义
  - 成员必须写 `SAMPLE_NAMES` 中的样本名

- `DEG_IDENT_1`
- `DEG_IDENT_2`
  - 差异分析对比组
  - 应填写 `ANALYSIS_GROUP_*_NAME`

### 4.2 RNA velocity 变量

- `RAW_SAMPLES`
  - velocity 使用的原始样本名
  - 对应 `cellranger_out/<sample>`

- `RAW_GROUP_1_NAME`
- `RAW_GROUP_1_SAMPLES`
- `RAW_GROUP_2_NAME`
- `RAW_GROUP_2_SAMPLES`
  - velocity 参考对象的分组

### 4.3 正式项目中的推荐关系

正式项目里，一般这样写：

- `SAMPLE_NAMES=Ctrl_1,Ctrl_2,Treat_1,Treat_2`
- `RAW_SAMPLES=Ctrl_1,Ctrl_2,Treat_1,Treat_2`
- `ANALYSIS_GROUP_1_NAME=Ctrl`
- `ANALYSIS_GROUP_1_SAMPLES=Ctrl_1,Ctrl_2`
- `ANALYSIS_GROUP_2_NAME=Treat`
- `ANALYSIS_GROUP_2_SAMPLES=Treat_1,Treat_2`
- `DEG_IDENT_1=Ctrl`
- `DEG_IDENT_2=Treat`
- `RAW_GROUP_1_NAME=Ctrl`
- `RAW_GROUP_1_SAMPLES=Ctrl_1,Ctrl_2`
- `RAW_GROUP_2_NAME=Treat`
- `RAW_GROUP_2_SAMPLES=Treat_1,Treat_2`

## 5. 流程入口脚本

### 5.1 项目初始化

- `workflow/00_init_project.sh`
  - 正式项目入口
  - 生成项目框架并自动切换全局配置

### 5.2 收件登记

- `workflow/01_register_delivery.sh`
  - 冻结原始交付物
  - 写入 delivery manifest 和 received files manifest

### 5.3 metadata 校验

- `workflow/02_validate_metadata.sh`
  - 校验 `metadata/samples.tsv` 和 `metadata/comparisons.tsv`
  - 生成 canonical sample sheet

### 5.4 输入审计

- `workflow/03_audit_inputs.sh`
  - 识别 `fastq / cellranger_out / matrix`
  - 生成 input inventory 和 branch readiness

### 5.5 输入标准化

- `workflow/04_standardize_inputs.sh`
  - 把 workflow 可用主输入统一挂到 `data/<sample>`
  - 不再区分旧的 `00_prepare_*` 入口

### 5.6 intake summary

- `workflow/05_input_summary.sh`
  - 生成 intake summary
  - 更新 `status/workflow_status.json`

### 5.7 环境安装

- `workflow/01_install_envs.sh`
  - 安装：
    - `r_main`
    - `r_legacy`
    - `r_scenic`
    - `velocity`
    - `pyscenic`
  - 并下载 SCENIC 资源

### 5.8 从 FASTQ 开始跑 Cell Ranger

- `workflow/10_run_cellranger_from_fastq.sh`
  - 从 `FASTQ` 开始
  - 生成 `cellranger_out`

### 5.9 主流程

- `workflow/20_run_main_pipeline.sh`
  - 分阶段执行：
    - `01_build_raw_objects.R`
    - `01a_pre_qc_eda.R`
    - `01b_qc_doublet.R`
    - `01c_post_qc_eda.R`
    - `02_build_reductions.R`
    - `02a_integration_eda.R`
    - `02b_finalize_clustering.R`
    - `03_annotation.R`
    - `03a_annotation_eda.R`
    - `04_deg_enrichment.R`
    - `05_trajectory.R`
  - 会在 `pre_qc / post_qc / integration / annotation` 四个 gate 自动停住
  - 审阅对应 `reports/eda/<stage>/report.md` 后，在 `config/eda_gates.tsv` 中把该 gate 设为 `approved`，再重跑同一个 `20_run_main_pipeline.sh`

### 5.10 RNA velocity

- `workflow/30_run_velocyto.sh`
  - 需要：
    - `outs/possorted_genome_bam.bam`
    - `outs/filtered_feature_bc_matrix`
  - 生成 `.loom`

- `workflow/31_run_scvelo.sh`
  - 先准备参考对象
  - 再跑 `scVelo`

### 5.11 SCENIC

- `workflow/40_build_ortholog_cache.sh`
  - 构建同源映射缓存

- `workflow/41_download_scenic_resources.sh`
  - 下载或检查 SCENIC 资源

- `workflow/42_run_scenic.sh`
  - 导出同源基因表达矩阵
  - 跑 pySCENIC
  - 跑 RSS / CSI 下游

## 6. 当你拿到原始数据后，实际怎么跑

### 6.1 情况 A：你拿到的是 FASTQ

先初始化新项目：

```bash
cd /home/user_test/scRNA_pipeline_chicken

bash workflow/00_init_project.sh \
  --project-root /home/user_test/projects/chicken_case01 \
  --fastq-source-dir /path/to/fastq \
  --genome-fasta-gz /path/to/genome.fa.gz \
  --reference-gtf /path/to/genes.gtf \
  --sample-names Ctrl_1,Ctrl_2,Treat_1,Treat_2 \
  --group1-name Ctrl \
  --group1-samples Ctrl_1,Ctrl_2 \
  --group2-name Treat \
  --group2-samples Treat_1,Treat_2
```

然后按顺序跑：

```bash
bash workflow/01_register_delivery.sh --source-path /path/to/fastq
bash workflow/02_validate_metadata.sh
bash workflow/03_audit_inputs.sh
bash workflow/10_run_cellranger_from_fastq.sh
bash workflow/04_standardize_inputs.sh
bash workflow/05_input_summary.sh
bash workflow/01_install_envs.sh
bash workflow/20_run_main_pipeline.sh
```

注意：

- `20_run_main_pipeline.sh` 第一次不会直接跑到底
- 它会在 `reports/eda/pre_qc`、`reports/eda/post_qc`、`reports/eda/integration`、`reports/eda/annotation` 依次停住
- 每次审阅完成后，把 `config/eda_gates.tsv` 对应行改成 `approved`，再重跑同一个命令

如果要 RNA velocity：

```bash
bash workflow/30_run_velocyto.sh
bash workflow/31_run_scvelo.sh
```

如果要 SCENIC：

```bash
bash workflow/42_run_scenic.sh
```

### 6.2 情况 B：你已经有 cellranger_out

初始化时再加：

```bash
--cellranger-out-source-dir /path/to/cellranger_out
```

然后跑：

```bash
bash workflow/01_register_delivery.sh --source-path /path/to/cellranger_out
bash workflow/02_validate_metadata.sh
bash workflow/03_audit_inputs.sh
bash workflow/04_standardize_inputs.sh
bash workflow/05_input_summary.sh
bash workflow/01_install_envs.sh
bash workflow/20_run_main_pipeline.sh
```

### 6.3 情况 C：你已经有 10X matrix

初始化时再加：

```bash
--matrix-source-dir /path/to/matrix_root
```

然后直接跑：

```bash
bash workflow/01_register_delivery.sh --source-path /path/to/matrix_root
bash workflow/02_validate_metadata.sh
bash workflow/03_audit_inputs.sh
bash workflow/04_standardize_inputs.sh
bash workflow/05_input_summary.sh
bash workflow/01_install_envs.sh
bash workflow/20_run_main_pipeline.sh
```

注意：

- 这种情况下主流程可以直接跑
- 但 velocity 仍然需要 `cellranger_out` 和 `BAM`

## 7. 主流程每一步做什么

### 7.1 QC 和去双细胞

- 输入：`data/<sample>`
- 输出：
  - `results/checkpoints/01_after_qc_doublet.rds`
  - `results/tables/qc_doublet_summary.csv`

### 7.2 聚类

- 标准化、高变基因、PCA、Harmony、UMAP、聚类
- 输出：
  - `results/checkpoints/02_after_clustering.rds`
  - `results/tables/cluster_summary.csv`

### 7.3 注释

- 按论文里的颗粒细胞 marker 做注释
- 输出：
  - `results/checkpoints/03_after_annotation.rds`
  - marker 表
  - cluster 到 cell type 的映射表

### 7.4 差异分析和富集

- 使用 `ANALYSIS_GROUP_*` 定义的组
- 输出：
  - DEG 表
  - GO 结果
  - KEGG 结果

### 7.5 拟时序

- `Slingshot + tradeSeq`
- 输出：
  - pseudotime 表
  - 动态基因结果
  - 拟时序热图

## 8. 哪些部分是固定的，哪些需要你改

### 8.1 固定部分

- 主流程顺序
- checkpoint 结构
- velocity 和 SCENIC 独立分支
- 论文来源的方法框架

### 8.2 每个新项目都要改的部分

- 项目根目录
- FASTQ 路径
- 参考文件路径
- 样本名
- 分组定义
- DEG 对比组

### 8.3 可能需要跟着项目变化调整的部分

- `TARGET_CLUSTERS`
- `TRAJECTORY_START`
- QC 阈值
- `celltype_marker_list()`

最后这一点非常重要：

当前注释 marker 和 `pGC/eGC/rgGC/lGC` 命名是按论文的鸡颗粒细胞写的。

- 如果你的新数据仍然是鸡卵巢颗粒细胞，可以先沿用
- 如果你换了物种或组织，必须改 marker、富集数据库和 SCENIC 同源映射逻辑

## 9. 环境拆分原则

### 9.1 r_main

只负责主流程。

当前固定策略：

- 运行时使用独立 conda 环境
- Seurat、Harmony、slingshot、tradeSeq、clusterProfiler 等优先由 conda 安装
- `DoubletFinder` 作为少数 conda 中没有的包，单独从 GitHub 固定提交安装
- `Monocle3` 作为正式轨迹验证/扩展后端进入共享环境，但默认主流程仍然先走 `Slingshot + tradeSeq`

负责内容：

- QC
- DoubletFinder
- Harmony
- clustering
- annotation
- DEG
- enrichment
- Slingshot
- tradeSeq
- Monocle3

### 9.2 velocity

只负责：

- `velocyto`
- `scVelo`
- `scanorama`
- `scrublet`

### 9.3 pyscenic

只负责：

- pySCENIC 网络推断
- AUC

### 9.4 r_scenic

只负责：

当前固定策略：

- 运行时使用独立 conda 环境
- AUCell、RcisTarget、GENIE3、ComplexHeatmap 等优先由 conda 安装
- `SCENIC` 本体从官方 GitHub release `v1.3.0` 安装
- `scFunctions` 作为 CSI 增强组件单独从 GitHub 固定提交安装
- `r-tidyverse`、`r-philentropy`、`r-svmisc`、`r-getopt`、`r-upsetr`、`r-dt` 预先放进 conda 环境，避免 `scFunctions` 二次装包时缺依赖

负责内容：

- SCENIC 下游整合
- RSS
- CSI

### 9.5 扩展环境

- `r_legacy`
  - `Monocle2`
  - `PHATE`
  - `topGO`
- `r_interaction`
  - `CellChat`
  - `NicheNet`
  - `decoupleR`
  - `gprofiler2`
- `r_spatial`
  - `BayesSpace`
  - `RCTD/spacexr`
  - `CARD`
  - `SPARK`
- `py_spatial`
  - `Squidpy`
  - `stLearn`
- `py_spatial_legacy`
  - `SpaGCN`
  - `STAGATE_pyG`
- `py_cell2location`
  - `cell2location`
- `r_validation`
  - `scDesign3`

### 9.6 复用原则

- 已有的 `r_main / velocity / r_scenic / pyscenic` 直接复用，不另外复制一套
- `r_legacy` 逻辑层已纳入共享体系，迁移期默认仍复用现有 shared prefix `r_heart_legacy`
- 共享环境根通过 `SHARED_ENV_DIR` 指向 `syf_f5/01shared_resources/envs/scRNA_pipeline_chicken`
- `workflow/01_install_envs.sh` 对已存在 prefix 走原地 `conda env update`
- 只有 prefix 不存在时才会新建
- 新增的扩展环境只在真正需要时再建，不默认一次性全装

## 10. 环境安装日志

运行 `workflow/01_install_envs.sh` 后，安装日志会自动写到项目的 `logs/` 目录：

- `install_r_main_env.log`
- `install_r_legacy_env.log`
- `install_r_scenic_env.log`
- `install_velocity_env.log`
- `install_pyscenic_env.log`
- `install_r_main_pkgs.log`
- `install_r_legacy_pkgs.log`
- `install_r_scenic_pkgs.log`
- `install_scenic_resources.log`

推荐排查顺序：

1. 先看 conda 环境是否创建成功
2. 再看 R 包安装日志里是否缺依赖
3. 最后看 SCENIC 资源下载日志

按需更新示例：

```bash
INSTALL_TARGETS=r_main,velocity INSTALL_SCENIC_RESOURCES=no bash workflow/01_install_envs.sh
```

这条命令只会原地更新 `r_main` 和 `velocity`，不会碰 `r_legacy / r_scenic / pyscenic`。

如果要把当前 5 套正式共享环境都同步到最新定义：

```bash
INSTALL_TARGETS=all bash workflow/01_install_envs.sh
```

如果要把扩展环境也一起建好：

```bash
INSTALL_TARGETS=r_main,r_legacy,velocity,r_scenic,pyscenic,r_interaction,r_spatial,py_spatial,py_spatial_legacy,py_cell2location,r_validation bash workflow/01_install_envs.sh
```

## 11. 当前实际状态

截至目前：

- `velocity` 环境已经可用
- `pyscenic` 环境已经可用，运行时会出现 `pkg_resources` 弃用警告，但不影响命令执行
- `r_legacy` 当前共享前缀路径名仍是 `r_heart_legacy`，但 smoke test 已通过
- `r_main` 当前未完全对齐重整方案：远端 smoke test 还缺 `gprofiler2`
- `r_scenic` 当前未完全对齐重整方案：远端 smoke test 还缺 `SCENIC`
- `r_interaction` 当前未完全对齐重整方案：远端 prefix 存在但 `Rscript` 不可用，需要重建或原地修复
- `r_spatial / py_spatial / py_spatial_legacy / py_cell2location / r_validation` 仍是已定义未装全状态

补充说明：

- `Monocle3 / scanorama / scrublet / gprofiler2 / SCENIC / CellChat / NicheNet` 已经在共享定义或安装脚本中声明，但是否真正可用仍以远端 smoke test 为准

建议查看日志：

```bash
tail -f /home/user_test/projects/chicken_sichuan_formal_demo/logs/install_r_main_pkgs.log
tail -f /home/user_test/projects/chicken_sichuan_formal_demo/logs/install_r_scenic_pkgs.log
tail -f /home/user_test/projects/chicken_sichuan_formal_demo/logs/install_pyscenic_env.log
```
