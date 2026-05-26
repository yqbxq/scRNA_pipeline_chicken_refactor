# 鸡卵巢颗粒细胞单细胞测序服务器流程

本目录用于把论文复现脚本整理成“正式项目可复用”的服务器流程。

这套流程的方法学来源仍然是论文 *Single-cell RNA sequencing uncovers cellular heterogeneity of granulosa cells and provides a signature for follicular development in chicken*，但工程目标已经改成：

- 论文数据可以作为验证入口
- 你自己的新项目可以从原始数据开始直接复用
- 固定部分和可变部分分开
- 环境拆开，避免主流程、velocity、SCENIC 互相污染

## all-repo 分支说明

`all-repo` 是当前单 repo 内部布局整合分支。旧的 `v00`-`v07`
里程碑提交保留在本分支祖先历史中，当前本地 `v08` 调控模块 WIP 已作为单独提交导入。

布局和迁移索引见：

- `docs/monorepo_layout.md`
- `docs/legacy_branch_migration.md`

## 0. 当前规范入口

`workflow/` 是唯一的流程目录，包含 stage shell、R 单模块脚本、Python helper 和工具脚本。推荐入口是：

- `workflow/03stages/init_project.sh`
- `workflow/03stages/register_delivery.sh`
- `workflow/03stages/validate_metadata.sh`
- `workflow/03stages/audit_inputs.sh`
- `workflow/03stages/standardize_inputs.sh`
- `workflow/03stages/input_summary.sh`
- `workflow/03stages/install_envs.sh`
- `workflow/03stages/alignment.sh`
- `workflow/03stages/00_ortholog.sh`
- `workflow/01run.sh <stage>`
- `workflow/03stages/10_velocity.sh`
- `workflow/03stages/08_regulation.sh`

`cell_count_inventory` 是项目级单一事实源：04c 生成
`cluster_eligibility_tsv`，05b/05a 用它决定 pseudobulk 是否可跑，07a/07c
用它过滤 CellChat/NicheNet 的输入 cluster。不要在下游模块里重新计算一套
细胞数阈值。

旧的顶层 wrapper 和旧 `workflow/r` 结构已在 `all-repo` 中移除，避免和新的单目录 workflow 重复。

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

- `workflow/03stages/init_project.sh`

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

- `workflow/03stages/init_project.sh`
  - 正式项目入口
  - 生成项目框架并自动切换全局配置

### 5.2 收件登记

- `workflow/03stages/register_delivery.sh`
  - 冻结原始交付物
  - 写入 delivery manifest 和 received files manifest

### 5.3 metadata 校验

- `workflow/03stages/validate_metadata.sh`
  - 校验 `metadata/samples.tsv` 和 `metadata/comparisons.tsv`
  - `comparisons.tsv` 支持可选 `subset_column` / `subset_value`，两列必须同时填写或同时留空
  - 生成 canonical sample sheet

### 5.3a 分析意图驱动 metadata

`metadata/analysis_questions.tsv` 是分析问题的 Tier 1 单一真相。05/06/07
stage 不直接维护模块输入表，而是在运行前调用 `ensure_metadata_fresh`：

- 当 `analysis_questions.tsv` 比生成表更新，或任一 Tier 2 表缺失时，自动运行 `workflow/03stages/95_run_metadata_generator.sh`。
- 随后运行 `workflow/03stages/96_validate_metadata.sh`，保证生成表和真实对象 metadata 的当前一致性。
- 生成表包括 `comparisons.tsv`、`communication_pairs.tsv`、`trajectory_pairs.tsv`、`scenic_targets.tsv`、`enrichment_targets.tsv`、`gene_program_targets.tsv`、`deconv_pairs.tsv` 和 `spatial_pairs.tsv`。

Tier 2 表为自动生成文件，不应手工编辑。后续新增 ST 或联合分析 stage 时，按 `docs/st_stage_metadata_hook_template.md` 在 `source common.sh` 后接入同一行 hook。

### 5.3b ST intake scaffold

ST 作为第二个 modality 接入共享 intake：

- `metadata/samples.tsv` 新增 `modality / section_id / chip_id / bundle_layout / image_path / run_spatial / run_deconv / run_joint`。
- `metadata/sections.tsv` 登记物理切片，`metadata/spatial_reference_inventory.tsv` 登记联合分析使用的 scRNA reference。
- `config/spatial_qc_thresholds.tsv.template`、`config/spatial_object_layers.tsv.template`、`config/spatial_intake_contract.tsv.template` 和 `config/marker_panels/spatial_region_panel.tsv.template` 是 ST 配置模板。
- `workflow/04python/spatial_intake.py` 负责 ST bundle 文件发现、ready 状态和标准化链接/复制；SAW 在 M2 只做存在性检查，不解析 GEF。
- `workflow/03stages/50_run_spatial_pipeline.sh`、`60_run_spatial_extensions.sh`、`70_run_joint_scrna_spatial.sh` 是非 R 入口脚手架；对应 R stage 尚未在本轮实现。

### 5.4 输入审计

- `workflow/03stages/audit_inputs.sh`
  - 识别 `fastq / cellranger_out / matrix`
  - 生成 input inventory 和 branch readiness

### 5.5 输入标准化

- `workflow/03stages/standardize_inputs.sh`
  - 把 workflow 可用主输入统一挂到 `data/<sample>`
  - 不再区分旧的 `00_prepare_*` 入口

### 5.6 intake summary

- `workflow/03stages/input_summary.sh`
  - 生成 intake summary
  - 更新 `status/workflow_status.json`

### 5.7 环境安装

- `workflow/03stages/install_envs.sh`
  - 安装：
    - `r_main`
    - `r_legacy`
    - `r_scenic`
    - `r_decoupler`
    - `velocity`
    - `pyscenic`
  - ST/联合分析相关环境可按需启用：
    - `--with-spatial` 安装 `r_spatial,py_spatial`
    - `--with-spatial-legacy` 安装 `py_spatial_legacy`
    - `--with-cell2location` 安装 `py_cell2location`
    - `--with-validation` 安装 `r_validation`
    - `--with-saw` 安装 `py_saw`
  - 并下载 SCENIC 资源

### 5.8 从 FASTQ 开始跑 Cell Ranger

- `workflow/03stages/alignment.sh`
  - 从 `FASTQ` 开始
  - 生成 `cellranger_out`

### 5.9 主流程

- `workflow/01run.sh <stage>`
  - 分阶段调度 `workflow/03stages/*.sh`
  - 例如：
    - `workflow/01run.sh 00_ortholog`
    - `workflow/01run.sh 01_build_raw`
    - `workflow/01run.sh 02_qc`
    - `workflow/01run.sh 03_panorama`
    - `workflow/01run.sh 04_subcluster`
    - `workflow/01run.sh 04d_cluster_robustness`
    - `workflow/01run.sh 04_robustness`
    - `workflow/01run.sh 05_deg`
    - `workflow/01run.sh 06_enrichment`
    - `workflow/01run.sh 07_communication`
    - `workflow/01run.sh 08_regulation`
  - 会在 `pre_qc / post_qc / integration / annotation` 四个 gate 自动停住
  - 审阅对应 `reports/eda/<stage>/report.md` 后，在 `config/eda_gates.tsv` 中把该 gate 设为 `approved`，再重跑对应 stage

### 5.10 04 子聚类模块（v04 分支）

04 子聚类模块走 workflow stage 入口：

- `workflow/03stages/04_subcluster.sh`
  - `04a_subcluster_build.R`：按 `object_layers.tsv` 生成 subcluster 层；默认继承 panorama 选定的 normalization/integration，只有 layer 显式配置多候选时才进入候选模式。
  - `04a_review.R`：汇总 candidate diagnostics，写 `subcluster_review_summary.tsv` 和每层 `selected_integration.txt`。
  - `04b_subcluster_annotate.R`：对已经 finalize 的 `clustered_<layer_id>` 复用 annotation helper，写 `annotated_<layer_id>` 和跨层注释汇总，并把子层注释回填到 panorama 的 `cell_subtype` metadata。
  - `04c_subcluster_eda.R`：生成 subcluster review 报告，汇总 cluster count、注释置信度、condition split、panorama-vs-subcluster 对照。
- `workflow/03stages/04d_cluster_robustness.sh`
  - `04d_cluster_robustness.R`：注册 question-driven scDesign3 gate、target、threshold 和 preflight 状态；cluster metrics 用 `pending_engine` 占位，避免被误读为真实 ARI/NMI。
- `workflow/03stages/04_robustness.sh`
  - `04e_scdesign3_engine.R`：在 `scdesign3_targets` gate 批准后，对 M1 `cluster_robustness` target 运行真实 scDesign3 fit/simulate/recluster/score 引擎；默认 `SCDESIGN3_N_SIM=5`。
  - `04f_scdesign3_finalize.R`：把 04e 的 ARI/NMI/Jaccard 结果回填到 target/question gate，并覆盖 root `results/tables/scdesign3_all_questions_status.tsv`。

gate 规则采用 Option B：

- 运行 `04_subcluster.sh` 前必须已经通过 `annotation` gate，且 `03_panorama_completed=true`。
- 如果所有 subcluster 都是 inherited 模式，`04a_review` 不会强制设置 `subcluster` gate，流程直接进入 04b/04c。
- 如果存在 candidate 模式 layer，`04a_review` 会把 `subcluster` gate 设为 `pending`；审阅 `subcluster_review_summary.tsv` 和各层 `selected_integration.txt` 后，把 `subcluster` gate 改为 `approved`，再重跑 `04_subcluster.sh`。
- `04c` 完成后不会再次重置 `subcluster` gate；是否运行 04d 由你在审阅 04c 报告后决定。
- 独立运行 `04d_cluster_robustness.sh` 时会同时校验 `subcluster_gate_passed=true` 和 `04_subcluster_completed=true`。
- `04d` 完成后会把 `scdesign3_targets` gate 置为 `pending`；审阅 `results/tables/04d_cluster_robustness/target_gate_status.tsv` 后批准该 gate，才能运行 `04_robustness.sh`。
- `04_robustness.sh` 完成后会把 `scdesign3_validated` gate 置为 `pending`；07 communication 会等待该 gate 批准后再进入核心解释。

`comparisons.tsv` 在 04c 中支持可选子集列：

- `subset_column` 和 `subset_value` 必须同时填写或同时留空。
- `subset_value` 支持 CSV，语义是 OR。
- 典型用法是先按 `panorama_cell_type` 或 `<layer_id>_cell_type` 过滤细胞，再做 condition split 或后续比较预览。

04 模块相关 status 字段：

| status 字段 | 写入时机 |
| --- | --- |
| `status.04a_subcluster_build_completed` | `04_subcluster.sh` 完成 04a/04a_review/finalize 后随 04c 完成态一起写入 |
| `status.04b_subcluster_annotated` | `04b_subcluster_annotate.R` 完成后随 04c 完成态一起写入 |
| `status.04c_subcluster_eda_completed` | `04c_subcluster_eda.R` 完成后写入 |
| `status.04_subcluster_completed` | 04a→04b→04c 全链完成后写入 |
| `status.04d_cluster_robustness_completed` | 独立 04d scDesign3 gate 注册 stage 完成后写入 |
| `status.04e_scdesign3_engine_completed` | 04e scDesign3 engine 完成后随 04_robustness 写入 |
| `status.04f_scdesign3_finalize_completed` | 04f gate finalize 完成后随 04_robustness 写入 |
| `status.04_robustness_completed` | 04e→04f 全链完成后写入 |

### 5.11 05 DEG 模块（workflow standalone）

05 DEG 走当前 workflow-native stage，底层脚本在 `workflow/05single_script/05*`：

- `workflow/03stages/05_deg.sh`
  - `05a_marker_discovery.R`：按 `analysis_mode` 运行 cell-level marker/DEG；`composition` 行跳过并交给 05c。
  - `05b_pseudobulk_de.R`：仅对 `condition_within_type` 行按 `aggregation_group_var` 运行 muscat pseudobulk DE；复制门失败时只写探索性状态和 05a fallback 路径。
  - `05c_composition.R`：仅对 `composition` 行按 `composition_group_var` 运行 sample-level proportion / propeller。
  - `05d_deg_eda.R`：汇总 marker / pseudobulk / composition manifest，输出 `reports/eda/deg/report.md` 和 `results/tables/deg/gene_program_registry.tsv`。

`comparisons.tsv` 支持 05 专用可选列：

- `subset_column` / `subset_value`：先过滤细胞，再做组间比较；`subset_value` 支持 CSV。
- `analysis_mode`：区分 `subtype_marker`、`subtype_pairwise`、`condition_within_type` 和 `composition`。
- `aggregation_group_var`：05b pseudobulk 的样本内聚合单位，例如 `cell_type`、`cell_subtype` 或 `all_cells`。
- `composition_group_var`：05c proportion 的分类单位，例如 `cell_subtype` 或 `cell_type`。
- `produces_gene_program` / `gene_program_role`：控制 05d gene-program registry 和 06/07/08 下游可用性。
- `force_exploratory`：默认 `no`；设为 `yes` 时 N=1 等复制门失败结果会标记为 `exploratory_forced`，仅供探索和下游富集预览。
- `min_cells_per_group`：默认 `3`，控制 05a 的 FindMarkers 最小细胞数。
- `logfc_threshold`：默认 `0`，控制 05a 的 FindMarkers logFC 阈值。

05 完成后会把 `deg` gate 置为 `pending`。审阅 `reports/eda/deg/report.md` 后，在 `config/eda_gates.tsv` 中批准 `deg`，后续 06 enrichment 才能继续。

### 5.12 06 Enrichment 模块（workflow standalone）

06 enrichment 依赖 05 DEG 完成且 `deg` gate 已批准：

- `workflow/03stages/06_enrichment.sh`
  - `06a_go_enrichment.R`：只消费 `metadata/enrichment_targets.tsv` 中 `database=GO` 的目标，再通过 `results/tables/deg/gene_program_registry.tsv` 找对应 gene program。
  - `06b_kegg_enrichment.R`：只消费 `metadata/enrichment_targets.tsv` 中 `database=KEGG` 的目标，再通过同一个 registry 找对应 gene program。
  - `06c_enrichment_eda.R`：汇总 GO/KEGG manifest，输出 `reports/eda/enrichment/report.md`、跨 cluster 热图和 shared pathways。

06 不再从 05 manifest 自动膨胀所有 DEG。没有出现在 `enrichment_targets.tsv` 的 DEG/marker 不会自动富集。每个 layer 会尽量生成 `results/tables/enrichment/background/background_genes_<layer>.tsv`，clusterProfiler GO/KEGG 使用该 layer-specific expressed universe；如果对象不可用，manifest 会保留背景为空/不可用的状态信息。

物种策略默认是 `ENRICHMENT_SPECIES_STRATEGY=chicken_primary`。如果设置为包含 `human`、`dual`、`both` 或 `mapped`，06 会额外使用 00 ortholog cache 将鸡基因映射到人类符号后跑辅助通道。运行后端优先使用可加载的 `clusterProfiler`；如果当前环境中 `clusterProfiler` 因 `DOSE` 等依赖不可加载，会回退到 `gprofiler2` 并在 manifest 的 `status` / `reason` 中记录。KEGG 在线查询默认 `ENRICHMENT_KEGG_TIMEOUT_SEC=60`，超时会记录为 `timeout` 状态。

06 manifest 和图标题会携带 `formal_status` / `result_level` / `biological_replicates`，使 `formal`、`exploratory_only`、`exploratory_forced` 不会在下游报告中丢失。

主要输出：

- `results/tables/enrichment/go/go_enrichment_manifest.tsv`
- `results/tables/enrichment/kegg/kegg_enrichment_manifest.tsv`
- `results/tables/enrichment/enrichment_summary.tsv`
- `results/figures/enrichment/cross_cluster_go_bp_heatmap.png`
- `results/figures/enrichment/cross_cluster_kegg_heatmap.png`
- `reports/eda/enrichment/report.md`

### 5.13 07 Communication 模块（workflow standalone）

07 communication 使用 `metadata/communication_pairs.tsv` 驱动，不复用 `comparisons.tsv`。每行是一个 `layer_scope + sender + receiver + condition_split` 通讯任务，适合 TC→GC、GC→TC、GC 亚群发育流等有方向的问题。

- `workflow/03stages/07_communication.sh`
  - `07a_cellchat.R`：按 layer / pair / condition 在 `r_interaction` 中运行 CellChat，并用 00 ortholog cache 将鸡表达矩阵映射到人类符号。
  - `07b_liana_consensus.py`：读取 H5AD 镜像并运行 LIANA+ 多算法共识；依赖缺失时写 schema-valid 空输出和 manifest。
  - `07c_nichenet.R`：按 sender→receiver 方向运行 NicheNet；通过 `communication_pairs.tsv` 的显式 `receiver_deg_comparison_id` / `baseline_marker_comparison_id` 到 `gene_program_registry.tsv` 查 gene program。
  - `07d_communication_consensus.R`：通信证据链 consensus 占位，完整 evidence tier 判定见 P-R03-E。
  - `07e_communication_eda.R`：按 `pair_id + layer + condition` 汇总 CellChat / NicheNet 共识，并按 `direction_filter` 决定是否保留全网络。

`communication_pairs.tsv` 支持 `sender` / `receiver` 精确 cell type、CSV、`*` 和 `prefix_*`。`condition_split_var` / `condition_split_values` 用于把 syf、f5 等阶段拆开分别跑。`requires_cell_subtype=yes` 时必须存在 `cell_subtype` 且 sender/receiver 必须能在该列解析，07 不允许 fallback 到 `cell_type` / `annotation_label` / `cluster_id`。NicheNet 缺失 gene program 会记录失败或跳过，不再临时跑 receiver FindMarkers，也不再用 `pair_id == comparison_id` 猜 DEG。NicheNet 三件套资源可用 `workflow/06tools/download_nichenet_resources.sh` 准备。

主要输出：

- `results/tables/communication/cellchat/cellchat_index.tsv`
- `results/tables/communication/nichenet/nichenet_index.tsv`
- `results/tables/communication/cross_validation.tsv`
- `results/tables/communication/consensus_lr.tsv`
- `reports/eda/communication/report.md`

07 完成后会把 `communication` gate 置为 `pending`，该 gate 目前只作为审阅记录，不阻塞下游模块。

### 5.14 RNA velocity

- `workflow/03stages/09_10_preflight.sh`
  - 在重型运行前检查 `trajectory_pairs.tsv` schema、layer_scope、split values、root/terminal 标签、velocity loom/BAM/H5/GTF 和 Seurat/loom cell id 命名空间
  - 输出 `results/tables/trajectory_velocity_preflight/preflight_checks.tsv` 和 `reports/eda/trajectory_velocity_preflight/report.md`
- `workflow/03stages/10_velocity.sh`
  - `10a_run_velocyto.sh`：从 dnbc4tools/Cell Ranger BAM 和 `filtered_feature_bc_matrix.h5` 生成 per-sample `.loom`；如果已存在匹配 `.loom`，会直接复用并写入 index
  - `10b_prepare_velocity_reference.R`：按 `trajectory_pairs.tsv` 中 `method=velocity` 的 H05/H06 行导出 UMAP 和 metadata reference，cell id 规范为 scVelo 使用的 `sample:barcode`
  - 完成 10a/10b 后会把 `velocity_inputs` gate 置为 `pending`
  - `10c_scvelo_dynamical.py`：按 pooled/split unit 运行 scVelo dynamical，并记录 stochastic 旁路指标
  - `10d_velocyto_steady_state.R`：运行 velocyto.R steady-state 交叉验证
  - `10e_scvelo_drivers.py`：导出 scVelo driver gene top list 和 split overlap
  - `10f_cellrank_fate.py`：默认运行 CellRank fate/macrostates，可用 `methods_extra=-cellrank` 禁用
  - `10h_velocity_root_terminal.R`：输出 09c 可反向读取的 `velocity_root_terminal_<unit>.tsv`
- `workflow/03stages/09_10_trajectory_velocity.sh`
  - 按计划顺序串联 `10_velocity.sh -> 09_trajectory.sh -> 10_velocity_finalize.sh`
- `workflow/03stages/10_velocity_finalize.sh`
  - `10g_velocity_consistency.R`：汇总 scVelo/stochastic/velocyto/Slingshot/CellRank 一致性和 split 比较
  - `10i_velocity_eda.R`：生成 `velocity_module_status.tsv`、`velocity_triage.tsv` 和 velocity finalize 报告

### 5.15 Regulation / SCENIC / decoupleR

- `workflow/03stages/00_ortholog.sh`
  - 构建同源映射缓存

- `workflow/03stages/08_regulation.sh`
  - `08a-08d`：SCENIC export、pySCENIC GRN、regulon/AUCell、RSS/CSI 下游
  - `08e`：decoupleR TF activity 和 pathway activity
  - `08f`：综合调控 EDA 报告，允许 SCENIC-only 或 decoupleR-only partial report

SCENIC 和 decoupleR 是互补证据，不是完全等价的方法。SCENIC 把 chicken expression matrix 映射到 human symbol 后使用 human cisTarget/motif 资源；decoupleR 则把 human DoRothEA/PROGENy prior network 映射到 chicken symbol，表达矩阵保留 chicken gene。08 的 TF/pathway 活性不能替代 05 DEG、06 enrichment 或 07 communication。

推荐先跑 panorama：

```bash
REGULATION_LAYERS=panorama bash workflow/03stages/08_regulation.sh
```

扩展到子层时显式指定：

```bash
REGULATION_LAYERS=panorama,GC_subcluster bash workflow/03stages/08_regulation.sh
```

详细 schema 见：

- `docs/module_08_regulation.md`
- `docs/regulation_schema.md`

## 6. 当你拿到原始数据后，实际怎么跑

### 6.1 情况 A：你拿到的是 FASTQ

先初始化新项目：

```bash
cd /home/user_test/scRNA_pipeline_chicken

bash workflow/03stages/init_project.sh \
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
bash workflow/03stages/register_delivery.sh --source-path /path/to/fastq
bash workflow/03stages/validate_metadata.sh
bash workflow/03stages/audit_inputs.sh
bash workflow/03stages/alignment.sh
bash workflow/03stages/standardize_inputs.sh
bash workflow/03stages/input_summary.sh
bash workflow/03stages/install_envs.sh
bash workflow/01run.sh <stage>
```

注意：

- `20_run_main_pipeline.sh` 第一次不会直接跑到底
- 它会在 `reports/eda/pre_qc`、`reports/eda/post_qc`、`reports/eda/integration`、`reports/eda/annotation` 依次停住
- 每次审阅完成后，把 `config/eda_gates.tsv` 对应行改成 `approved`，再重跑同一个命令

如果要 RNA velocity：

```bash
bash workflow/03stages/09_10_preflight.sh
bash workflow/03stages/09_10_trajectory_velocity.sh
```

等价的显式顺序是：

```bash
bash workflow/03stages/10_velocity.sh
bash workflow/03stages/09_trajectory.sh
bash workflow/03stages/10_velocity_finalize.sh
```

如果要 SCENIC：

```bash
bash workflow/03stages/08_regulation.sh
```

### 6.2 情况 B：你已经有 cellranger_out

初始化时再加：

```bash
--cellranger-out-source-dir /path/to/cellranger_out
```

然后跑：

```bash
bash workflow/03stages/register_delivery.sh --source-path /path/to/cellranger_out
bash workflow/03stages/validate_metadata.sh
bash workflow/03stages/audit_inputs.sh
bash workflow/03stages/standardize_inputs.sh
bash workflow/03stages/input_summary.sh
bash workflow/03stages/install_envs.sh
bash workflow/01run.sh <stage>
```

### 6.3 情况 C：你已经有 10X matrix

初始化时再加：

```bash
--matrix-source-dir /path/to/matrix_root
```

然后直接跑：

```bash
bash workflow/03stages/register_delivery.sh --source-path /path/to/matrix_root
bash workflow/03stages/validate_metadata.sh
bash workflow/03stages/audit_inputs.sh
bash workflow/03stages/standardize_inputs.sh
bash workflow/03stages/input_summary.sh
bash workflow/03stages/install_envs.sh
bash workflow/01run.sh <stage>
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
- `workflow/03stages/install_envs.sh` 对已存在 prefix 走原地 `conda env update`
- 只有 prefix 不存在时才会新建
- 新增的扩展环境只在真正需要时再建，不默认一次性全装

## 10. 环境安装日志

运行 `workflow/03stages/install_envs.sh` 后，安装日志会自动写到项目的 `logs/` 目录：

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
INSTALL_TARGETS=r_main,velocity INSTALL_SCENIC_RESOURCES=no bash workflow/03stages/install_envs.sh
```

这条命令只会原地更新 `r_main` 和 `velocity`，不会碰 `r_legacy / r_scenic / pyscenic`。

如果要把当前 5 套正式共享环境都同步到最新定义：

```bash
INSTALL_TARGETS=all bash workflow/03stages/install_envs.sh
```

如果要把扩展环境也一起建好：

```bash
INSTALL_TARGETS=r_main,r_legacy,velocity,r_scenic,pyscenic,r_interaction,r_spatial,py_spatial,py_spatial_legacy,py_cell2location,r_validation bash workflow/03stages/install_envs.sh
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
