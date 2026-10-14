# Codex 执行任务：SYF/F5 历史代码完整取证与两层 GitHub 归档

> 2026-10-14。**这是可执行任务书，当前尚未在本地 WSL/服务器执行。** 凡 GitHub 启动包中标为 historical report 的信息，必须在现场核验。

## 0. 任务目标

沿用 F6 的治理实践：先独立归档全部历史来源代码与证据，完成双向「目的 → 代码 → 输入 → 对象 → 运行 → 表/图」溯源，保留全部方法版本、重复版本、图版本；经用户选版和实际测试后，再建立正式可复用模块。

建议正式成果为两个**独立**仓库/工作区：
1. `syf-f5-legacy-analysis-audit`：只读代码证据、文件清点、日志、对象/表/图登记、重复族、人工选版画廊、Git 状态/证据收据（不得提交生物数据大文件）。
2. `syf-f5-analysis`：被批准且测试通过的正式可复用代码、参数、模块、对象契约、测试和文档。

本分支 `governance/syf_f5/` 是启动包和现场核验前的索引入口，**不预设两份独立仓库已经创建**。不得直接覆盖 F6 两仓库，也不得把现有通用流水线当作唯一事实来源。

## 1. Skills 与环境核实（必须真实执行）

在用户的 WSL/Codex 工作区：
1. 查找当前项目 `AGENTS.md`，定位 `~/.codex/skills/`、`.agents/skills/` 及其它明确登记的 Skill 安装目录。逐个打开并阅读以下 `SKILL.md`：
   - `research-agent-orchestrator`
   - `research-repo-audit`
   - `research-provenance`
   - `reuse-impact-control`
   - `se-method-router`
   - `safe-git-sync`（只在提交前应用）
2. 有 `governance_tools/research_guard.py`、`governance_tools/knowledge_route.py` 时，先运行 `--help` 核实实际参数；按 Skill 要求依次 `inventory / impact / capture / verify`，按需路由知识库方法，不得虚构 CLI 调用。
3. 需要 `tools/apply_knowledge.py` 时，先核对知识工程目录并运行对应的 `check_knowledge.py` 与 `test_repairs.py`。知识库未批准候选按 candidate 处理，不把候选方法卡视作现场真实性证据。
4. 输出 `skill_execution_receipt.md`：每个 Skill 的实际路径、读取结果、调用命令、stdout/退出码、有效范围；缺失/不兼容时标记 unavailable 并明确降级方案。

## 2. 现场保护与身份核查（最先执行）

- 以用户当前正式分析的 SYF/F5 任务为主：`~/.codex/tasks/*-syf-f5/` 与服务器 `05projects/syf_f5` 下的 `results/`、`workflow/`、项目配置、日志、Git 工作树、远程部署副本。查找所有旧任务和历史路径，不能只读这两个目录就认定覆盖完成。
- 记录 WSL、服务器、GitHub 的独立快照：`git status`、branch/HEAD、tracked/untracked 状态、remotes、工作树修改、仓库边界和源文件 SHA。扫描以只读为主，原始文件不改、不移动、不删。
- 在**私有/本地**清单中保留真实绝对路径；公开 GitHub 仅保留脱敏路径、来源 ID、hash 和非敏感的小体积源码。禁止提交 FASTQ、RDS、h5ad、凭据、SSH key、token、个人信息和非公开样本元数据。
- 检查 `DATA_VERSION_LEDGER.md`：严禁把不同 cellset 的 21,556/18,550 与 24,395/16,277 混同。通过 barcode 哈希、assay/layer、QC lineage、annotation authority 及运行日志确认关系；最终 GROUP vs JOINT 是否已有用户授权，必须单独取证。
- 对 F5_2、SYF_2、F6_1 的交叉参与代码做旁支登记，不让疑难样本或 pilot 成为 CORE4 biological replicate。

## 3. 第一轮范围：完整历史原样入库，先不要重写/去重

从所有可访问的 `.R`、`.Rmd`、`.py`、`.sh`、Notebook、config、helper、Git 历史和临时任务脚本提取历史版本。对每份小体积源文件保留**逐字节**快照与 SHA，路径冲突时采用 source ID 而非覆盖；每份重复代码也要留原位置及来源记录。

### 分组（发现额外功能时扩充）

| 代码域 | 涵盖内容 |
|---|---|
| A 输入准备 | FASTQ/计数矩阵/dnbc4tools/Cell Ranger、gene ID、矩阵转换、样本元数据、资源引用 |
| B 污染与质量控制 | RAW、pre-QC EDA、SoupX/DecontX、scaled-MAD 与 fixed QC、scDblFinder、doublet 比率、post-QC EDA、失败样本检查 |
| C 单样本分析 | SYF_1/SYF_3/F5_1/F5_3 各自归一化、HVG、PCA、UMAP、聚类、人工注释、marker、DEG、边界判定、图版本 |
| D 阶段内整合 | SYF 两重复与 F5 两重复各自 NONE/RPCA/Harmony/FastMNN、参数网格、稳定性、标签迁移和审核 |
| E 四样本共同分析 | CORE4 GROUP/JOINT、多批次整合/注释一致性、同阶段重复性、跨阶段比较、pseudobulk、composition、DA 等 |
| F GC/TC 亚群 | 父对象来源、GC 与 TC subset、每次聚类、marker/DEG、边界、注释、GC/TC cross-stage 映射、所有图 |
| G 富集与程序 | GO/KEGG/Reactome/ORA/GSEA、鸡基因/ortholog、背景集和来源 DEG、cNMF/模块评分（若存在） |
| H 轨迹与动力学 | Slingshot/Monocle/PAGA/tradeSeq/scVelo 等（以实际运行发现为准） |
| I 细胞通讯 | CellChat、LIANA、LR 来源、GC↔TC / 其它细胞方向、NicheNet、下游支持和不同图版本 |
| J 调控与专项 | SCENIC、AUCell、decoupler、TF、空间推断及其它实际存在的任务 |
| K 可视化/报告 | 专门汇图脚本、DotPlot/UMAP/Heatmap、制图主题/颜色/排序、图表合成和 HTML 报告 |

**每个领域都要同时登记已完成、失败、跳过、仅写代码、仅有方案、执行状态未知的资产。** 并行方法、参数敏感性、校正版本、已选图版本都要留下原始文件。F6 审计用的 N/I/E/C 功能登记口径可以复用为索引视角，但不得因此遗漏 SYF/F5 的任何领域。

## 4. 溯源、版本族与人工选版

- 产生 `REGISTRY_SCHEMA.md` 列出的全部登记表，含 `source_manifest.csv`、`script_inventory.csv`、`object_registry.csv`、`run_registry.csv`、`table_registry.csv`、`figure_registry.csv`、`trace_graph.csv`、`analysis_families.csv`、`duplicate_implementations.csv`、`coverage_matrix.csv`、`approval_queue.csv`、`excluded_side_branches.csv`。
- **必须查真实运行证据**：生成命令/脚本、调度日志、stderr、历史任务报告、sessionInfo、文件哈希、产物内容。仅代码里出现 `saveRDS()` 或 `ggsave()` 不等于实际执行；仅存在图片不等于已锁定生成代码。
- 同一分析目的整理成版本族：完全重复、逻辑重复但参数不同、输入对象不同、修订、平行科学方法、用户偏好图。分别标明差异产生于对象、算法、阈值、注释、颜色/主题还是绘制数据。
- 搜寻用户过去**已经明确选定的图片版本**的原始出处；若只有历史对话的选择、没有可核验源文件，则标为 `user_preference_reported_needs_source`。凡未确认的统一加入 `review_gallery.html`，绝不以 mtime 或文件名 final 自动选版。
- 每个模块都要可做双向反查：给定任意 `figure_id` 能得到 `script_id/run_id/object_id/table_id`；给定 `script_id` 能定位每个执行版本和输出。失败的追溯也要写入 `open_gaps.md`。

## 5. GitHub 同步与后续模块重构

- 第一阶段归档完成后，对历史代码快照做原件哈希一致性、敏感信息排除、覆盖率核验；合格的小体积审计资产可以通过 `safe-git-sync` 推送到**独立历史审计仓库**的任务分支，保留完整推送 SHA/URL 收据。无需把未批准候选写成正式代码。
- 第二阶段由用户在 `approval_queue.csv` 审核每个分析目的、对象/图版本。用户确定的版本才移植到**独立正式重构仓库**，优先共用函数与依赖注入，保留必要的并行科学方法为 adapter；不要提前删除重复实现。
- 变更前用 `reuse-impact-control` 展开上下游影响，添加可替换输入契约/样本防混逻辑/无数据的静态检查/小规模测试。真实全量执行需另行获得计算授权，记录成本；回归同输入和人工批准后才 Git 提交。
- 不得改动 `yqbxq/f6-analysis`、`yqbxq/f6-legacy-analysis-audit`、已存在的正式对象或任一历史出图文件。与 F6/疑难样本相关的跨样本候选只列为旁支或明确另行批准的研究。

## 6. 本次 Codex 第一次实际交付（不能只交计划）

1. 列出**实际访问**的本地/服务器目录、Git 仓库状态和 Skill 安装核查收据，缺失权限必须具体写明。
2. 输出 A–K 各自源码、对象、表格、图、日志**可证实数量**；先建完整清单、不可只做 GC/TC 或 DotPlot。
3. 将全部可获取的小体积历史源代码归档为只读素材，并登记精确 SHA；重复版本仍保留。
4. 建立至少一条 RAW → SoupX → QC → doublet → 阶段整合 → 注释 → 具体图/表的**实证**链。无法闭环时标为 unresolved，附所缺日志或产物。
5. 生成版本差异清单、已有用户指定图版本核对记录和 `review_gallery.html`。审计仓库安全推送后给出 commit SHA；推送失败列出真实错误。
6. 回复准确的完成/未完成矩阵；不能把仅在 GitHub 看到的通用工作流说成 SYF/F5 真实生产代码。

以上步骤由能访问 WSL 和服务器的 Codex 按已安装 Skill 执行。当前 ChatGPT GitHub 启动分支本身不含完整服务器历史资产。
