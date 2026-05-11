# 鸡颗粒细胞单细胞流程正式 SOP

## 1. 目的

这套 SOP 的目标不是把论文做到像素级复刻，而是把这篇论文里**真正能复用的方法学**整理成一套以后可以直接在服务器上跑你自己项目的数据流程。

因此我们明确分成两层：

- **标准流程**
  - 面向以后自己的真实项目
  - 只保留可复用、可迁移的通用方法
- **论文模式**
  - 面向这篇鸡颗粒细胞论文的对照、出图和特化验证
  - 单独放置，不污染标准流程

## 1.1 当前规范入口

从这一版开始，正式推荐的 shell 入口顺序是：

- `workflow/03stages/init_project.sh`
- `workflow/03stages/register_delivery.sh`
- `workflow/03stages/validate_metadata.sh`
- `workflow/03stages/audit_inputs.sh`
- `workflow/03stages/standardize_inputs.sh`
- `workflow/03stages/input_summary.sh`
- `workflow/03stages/alignment.sh`
- `workflow/01run.sh <stage>`
- `workflow/03stages/10_velocity.sh`
- `workflow/03stages/10_velocity.sh`
- `workflow/03stages/00_ortholog.sh`
- `workflow/03stages/08_regulation.sh`
- `workflow/03stages/08_regulation.sh`
- `workflow/06tools/download_nichenet_resources.sh`

旧编号脚本仍保留兼容，但不再是推荐的主入口。

---

## 1.2 工程输出语言和布尔字段规范

- CI、shell/R/Python 运行日志、`stop()`/exception 错误信息统一使用英文，便于自动化匹配和远端排查。
- 面向人工审阅的分诊表、报告 `recommended_action` 可以使用中文，优先保证项目使用者能直接理解下一步动作。
- metadata 和运行索引中的布尔开关继续使用 `yes/no` 字符串；脚本入口必须通过 `normalize_flag()` 或同等逻辑标准化，不直接依赖大小写或 TRUE/FALSE 类型。

---

## 2. 当前标准流程进度

截至目前，标准流程已经具备下列能力：

### 2.1 已完成并进入标准流程的部分

- 服务器环境拆分完成
  - `r_main`
  - `r_scenic`
  - `r_decoupler`
  - `velocity`
  - `pyscenic`
- 主流程可运行
  - `QC`
  - `DoubletFinder`
  - `Harmony`
  - `聚类`
  - `细胞类型注释`
  - `DEG`
  - `GO/KEGG`
  - `Slingshot`
  - `tradeSeq`
- RNA velocity 已进入标准流程
  - `velocyto`
  - `scVelo`
- SCENIC 已进入标准流程
  - `pySCENIC`
  - `AUCell`
  - `RSS`
  - `CSI`
- decoupleR 已进入 08 调控证据层
  - DoRothEA TF activity
  - PROGENy pathway activity
  - SCENIC/decoupleR partial report
- 轨迹推断已经固定为**双输出**
  - 粗粒度：`cell_type`
  - 细粒度：`seurat_clusters`
- 通用出图已经进入标准流程
  - `Figure 2`
  - `Figure 3`
  - `Figure 5F`
  - `Figure 6A-6D`
- 同源映射缓存构建已经补入标准流程
  - 默认 `strict`
  - 可选 `relaxed`
- 细胞通讯 07 模块已经进入标准流程
  - `communication_pairs.tsv` 四元组驱动
  - CellChat / NicheNet
  - syf/f5 condition split

### 2.2 已完成但仍属于论文模式的部分

- 论文式 `Figure 4` 复现热图
- 论文式 `Figure 6` 特化排版
- 论文补充表对照逻辑

### 2.3 当前仍属于试验/验证中的部分

- `relaxed ortholog -> SCENIC -> Figure 6` 这一条仍在试验线
- 论文 `Figure 6` 与原图的高拟合对齐仍未完全完成

---

## 3. 这次复现用到了哪些数据

### 3.1 论文作者提供的数据

- 作者提供的聚合矩阵
  - `HF_GCs`
  - `PHF_GCs`
- 原始样本级 `cellranger_out`
  - `SRR32621891`
  - `SRR32621892`
  - `SRR32621893`
  - `SRR32621894`
  - `SRR32621895`

### 3.2 作者补充表

位于：

- [鸡卵巢(四川)](/D:/code/r/鸡卵巢(四川))

其中可以直接对照标准流程的有：

- [mmc2.xlsx](/D:/code/r/鸡卵巢(四川/1-s2.0-S0032579125006741-mmc2.xlsx)
  - 测序统计
- [mmc3.xlsx](/D:/code/r/鸡卵巢(四川/1-s2.0-S0032579125006741-mmc3.xlsx)
  - cluster population
- [mmc4.xlsx](/D:/code/r/鸡卵巢(四川/1-s2.0-S0032579125006741-mmc4.xlsx)
  - 各 cell type top 10 marker
- [mmc5.xlsx](/D:/code/r/鸡卵巢(四川/1-s2.0-S0032579125006741-mmc5.xlsx)
  - 各 cell type GO
- [mmc6.xlsx](/D:/code/r/鸡卵巢(四川/1-s2.0-S0032579125006741-mmc6.xlsx)
  - 各 cell type KEGG
- [mmc7.xlsx](/D:/code/r/鸡卵巢(四川/1-s2.0-S0032579125006741-mmc7.xlsx)
  - HF vs PHF 的 DEG
- [mmc8.xlsx](/D:/code/r/鸡卵巢(四川/1-s2.0-S0032579125006741-mmc8.xlsx)
  - DEG 的 GO
- [mmc9.xlsx](/D:/code/r/鸡卵巢(四川/1-s2.0-S0032579125006741-mmc9.xlsx)
  - DEG 的 KEGG

注意：

- 这些补充表**没有作者公开的 SCENIC 原始表**
- 也没有现成的 `RSS/CSI/regulon` xlsx

因此：

- `mmc2-9` 适合验证标准流程主分析
- `Figure 6 / SCENIC` 不能只靠 `xlsx` 做完全对照

---

## 4. 标准流程的正式结构

标准流程目录：

- [server_final_pipeline](/D:/code/r/scRNA_liucheng/server_final_pipeline)

核心入口：

- `workflow/03stages/init_project.sh`
- `workflow/03stages/register_delivery.sh`
- `workflow/03stages/validate_metadata.sh`
- `workflow/03stages/audit_inputs.sh`
- `workflow/03stages/standardize_inputs.sh`
- `workflow/03stages/input_summary.sh`
- `workflow/03stages/alignment.sh`
- `workflow/03stages/install_envs.sh`
- `workflow/01run.sh <stage>`
- `workflow/03stages/10_velocity.sh`
- `workflow/03stages/10_velocity.sh`
- `workflow/03stages/00_ortholog.sh`
- `workflow/03stages/08_regulation.sh`
- `workflow/03stages/08_regulation.sh`
- `workflow/06tools/download_nichenet_resources.sh`

---

## 5. 标准流程每一步做什么

### 5.1 项目初始化

入口：

- `workflow/03stages/init_project.sh`

作用：

- 新建项目根目录
- 自动生成标准目录框架
- 自动链接外部 `FASTQ / reference / cellranger_out / matrix`
- 自动生成项目专属配置文件

### 5.2 环境安装

入口：

- `workflow/03stages/install_envs.sh`

作用：

- 安装 `r_main`
- 安装 `r_scenic`
- 安装 `velocity`
- 安装 `pyscenic`
- 下载 SCENIC 数据库

### 5.3 主流程

入口：

- `workflow/01run.sh <stage>`

内部步骤：

1. `01_build_raw_objects.R`
2. `01a_pre_qc_eda.R`
3. `01b_qc_doublet.R`
4. `01c_post_qc_eda.R`
5. `02_build_reductions.R`
6. `02a_integration_eda.R`
7. `02b_finalize_clustering.R`
8. `03_annotation.R`
9. `03a_annotation_eda.R`
10. `04_deg_enrichment.R`
11. `05_trajectory.R`

门控规则：

- `20_run_main_pipeline.sh` 会在 `pre_qc / post_qc / integration / annotation` 四个 EDA 节点停住
- 操作者需要审阅 `reports/eda/<stage>/report.md`
- 审阅通过后，把 `config/eda_gates.tsv` 中对应 `gate_id` 的 `status` 改为 `approved`
- 然后重跑同一个 `20_run_main_pipeline.sh`

输出：

- `results/checkpoints`
- `results/tables`
- `results/figures`

### 5.4 RNA velocity

入口：

- `workflow/03stages/10_velocity.sh`
- `workflow/03stages/10_velocity.sh`

要求：

- 必须有 `cellranger_out`
- 必须有 `BAM`
- 必须有 `filtered_feature_bc_matrix`

### 5.5 同源映射

入口：

- `workflow/03stages/00_ortholog.sh`

作用：

- 从**标准项目目录**下的 `${RESULTS_DIR}/checkpoints/03_after_annotation.rds` 提取鸡基因
- 用 `biomaRt` 建鸡到人的同源映射缓存
- 输出供 SCENIC 复用的本地映射表
- 后续 SCENIC 直接读取本地缓存，不重复在线查

输出目录：

- `${RESULTS_DIR}/ortholog_cache`

关键输出：

- `chicken_human_orthologs.csv`
- `chicken_human_orthologs.rds`
- `chicken_human_orthologs_for_pipeline.csv`
- `ortholog_coverage.csv`
- `ortholog_coverage.txt`
- `ortholog_key_gene_hits.csv`

### 5.6 Regulation / SCENIC / decoupleR

入口：

- `workflow/03stages/08_regulation.sh`

运行顺序：

1. 下载/检查 SCENIC 数据库
2. 构建或读取同源映射缓存
3. 导出人类表达矩阵
4. 运行 `pySCENIC grnboost2`
5. 运行 regulon / AUCell
6. 运行 `RSS / CSI`
7. 运行 decoupleR TF activity / pathway activity
8. 输出 `Figure 6A-6D` 和 `reports/eda/regulation/report.md`

关键原则：

- `08a-08d` 是 SCENIC。
- `08e` 是 decoupleR。
- `08f` 是综合报告层，只读取 08a-08e 与 05/06/07 已有 summary。
- SCENIC 把 chicken expression matrix 映射到 human symbol。
- decoupleR 把 human prior network 映射到 chicken symbol，表达矩阵保留 chicken gene。
- 08 是调控证据，不能替代 05 DEG、06 enrichment 或 07 communication。

推荐先跑：

```bash
REGULATION_LAYERS=panorama bash workflow/03stages/08_regulation.sh
```

扩展子层：

```bash
REGULATION_LAYERS=panorama,GC_subcluster bash workflow/03stages/08_regulation.sh
```

---

## 6. 拟时序为什么必须双输出

现在标准流程已经固定成双输出：

- 粗粒度：`cell_type`
- 细粒度：`seurat_clusters`

原因不是重复，而是它们解决的问题不同。

### 6.1 `cell_type` 粗版

适合：

- 你自己项目的常规解释
- 稳定展示主发育方向
- `tradeSeq`
- 跨项目复用

优点：

- 生物学解释直接
- 对分群编号不敏感
- 更稳

### 6.2 `seurat_clusters` 细版

适合：

- 看更细的分叉
- 对照论文里的细粒度轨迹
- 发现大类内部的小状态差异

优点：

- 分辨率更高
- 更接近论文 `Figure 5F` 的展示逻辑

### 6.3 当前标准流程产物

标准流程现在会同时产出：

- `slingshot_pseudotime_cell_type.csv`
- `slingshot_pseudotime_seurat_clusters.csv`
- `Figure_5F_Lineages_cell_type.png`
- `Figure_5F_Lineages_seurat_clusters.png`
- `Figure_5F_Lineages.png`

其中：

- `Figure_5F_Lineages.png` 默认给细粒度 top-N 版本
- `tradeSeq` 默认走粗粒度

---

## 7. 同源映射为什么要保留 strict / relaxed 两种模式

### 7.1 默认模式：strict

标准流程默认应使用：

- `SCENIC_ORTHOLOG_MODE=strict`

含义：

- 只保留 one-to-one ortholog

适合：

- 标准流程默认版
- 自己项目的正式结果
- 更稳、更可复用的分析

### 7.2 可选模式：relaxed

可选：

- `SCENIC_ORTHOLOG_MODE=relaxed`

含义：

- 允许 `one-to-many / many-to-many`
- 但仍然只为每个鸡基因选一个最佳人同源

适合：

- `strict` 覆盖率明显偏低
- 关键 TF 明显缺失
- 你想探索更完整的调控网络

### 7.3 当前建议

- 标准流程默认：`strict`
- 论文模式或探索模式：必要时试 `relaxed`

不要反过来：

- 不要把 `relaxed` 直接当成标准流程默认值

---

## 8. 通用出图哪些已经进标准流程

### 8.1 已经纳入标准流程的通用图

- `Figure 2`
  - cluster split UMAP
  - cluster proportion
- `Figure 3`
  - cell type split UMAP
  - feature plot
  - dot plot
  - cell type proportion
- `Figure 5F`
  - coarse / fine 双轨迹
- `Figure 6A-6D`
  - regulon heatmap
  - RSS rank
  - CSI clustering
  - CSI module activity

这些图的通用函数集中在：

- `workflow/05single_script/helpers/plotting_utils.R`

### 8.2 仍然留在论文模式的图

- 论文式 `Figure 4`
- 论文式 `Figure 6` 固定基因/固定 panel 排版

原因：

- 这些图带有明显的论文特化逻辑
- 不适合直接进标准流程默认版

---

## 9. 当前对照结果该如何理解

### 9.1 可以认为标准流程已经基本过关的证据

- `DEG` 对照最好
  - 作者 `mmc7` 中的 DEG 当前都能落在我们的显著 DEG 集里
- 主流程可以从服务器独立跑通
- velocity 可以跑通
- SCENIC 可以跑通

### 9.2 不能简单拿来判死刑的差异

- `cluster 编号不同`
  - 这是正常现象
  - cluster 数字只是标签，不是固定生物学编号
- `marker/GO/KEGG` 不完全一致
  - 这块容易受注释规则、背景集和排序逻辑影响

### 9.3 当前最不稳定的部分

- `SCENIC / Figure 6`

原因：

- 同源映射规则
- TF 列表版本
- motif 数据库版本
- 论文资源版本不可见

因此：

- `Figure 6` 不应当作为标准流程是否可用的唯一判据

---

## 10. 你以后拿到自己的数据后，应该怎么跑

### 10.1 如果拿到的是 FASTQ

顺序：

1. `00_init_project.sh`
2. `01_register_delivery.sh`
3. `02_validate_metadata.sh`
4. `03_audit_inputs.sh`
5. `10_run_cellranger_from_fastq.sh`
6. `04_standardize_inputs.sh`
7. `05_input_summary.sh`
8. `01_install_envs.sh`
9. `20_run_main_pipeline.sh`
10. `30_run_velocyto.sh`
11. `31_run_scvelo.sh`
12. `42_run_scenic.sh`

### 10.2 如果拿到的是 cellranger_out

顺序：

1. `00_init_project.sh`
2. `01_register_delivery.sh`
3. `02_validate_metadata.sh`
4. `03_audit_inputs.sh`
5. `04_standardize_inputs.sh`
6. `05_input_summary.sh`
7. `01_install_envs.sh`
8. `20_run_main_pipeline.sh`
9. `30_run_velocyto.sh`
10. `31_run_scvelo.sh`
11. `42_run_scenic.sh`

### 10.3 如果拿到的是 10X matrix

顺序：

1. `00_init_project.sh`
2. `01_register_delivery.sh`
3. `02_validate_metadata.sh`
4. `03_audit_inputs.sh`
5. `04_standardize_inputs.sh`
6. `05_input_summary.sh`
7. `01_install_envs.sh`
8. `20_run_main_pipeline.sh`
9. 如果以后要 velocity，再补 `cellranger_out + BAM`
10. `42_run_scenic.sh`

---

## 11. 最终建议

当前最合理的工作方式是：

1. **冻结标准流程 v1**
   - 主流程
   - velocity
   - SCENIC
   - 双轨迹
   - 同源映射缓存
   - 通用出图

2. **保留论文模式**
   - 单独对照原文
   - 单独做高拟合图
   - 不反向污染标准流程

3. **等你自己的数据下来时**
   - 优先跑标准流程
   - 再根据问题决定是否启用论文特化逻辑

这才最符合你一开始的目标：

- 不是为了复现而复现
- 而是通过复现，把以后项目的阻力降到最低
