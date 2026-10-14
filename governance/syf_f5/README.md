# SYF–F5 历史代码治理 · 启动记录

**创建日期**：2026-10-14　**状态**：GitHub 侧启动包（尚未连接服务器/WSL 现场）

此目录承接 F6 的方法，分为「完整不可变历史素材」→「代码—输入对象—运行—图表证据关联」→「用户选版」→「可复用模块与回归测试」。当前提交**只含操作契约和已知版本线索**，不代表已收集 SYF/F5 服务器历史代码、不代表任何分析重跑或测试通过。

## 所用治理技能

在具备历史源码访问权限的 Codex 工作区，先实际阅读：
- `research-agent-orchestrator`
- `research-repo-audit`
- `research-provenance`
- `reuse-impact-control`
- `se-method-router`
- 经人工批准进入发布阶段时才用 `safe-git-sync`

逐一检查 `SKILL.md`、工作区 `AGENTS.md` 以及 `governance_tools/research_guard.py` 是否存在；未安装或不可用者须在运行记录中登记。知识库原文为 candidate/anchored_unreviewed 的规则继续按候选约束处理，不能冒充验证结论。

## 界限与阶段

- **正式分析对象**：SYF_1、SYF_3、F5_1、F5_3；阶段内重复与跨阶段比较必须记录真实样本单位。
- **历史发现范围**：这四份样本的全部单样本、阶段内及四样本方案，含失败、修订、可视化多版本；发现 F5_2、SYF_2、F6 参与的脚本要登记为独立 **excluded/side branch**，不得自动并入正式四样本分析。
- **当前操作**：先只读清点和证据归档，完整保留原始代码、helpers、运行命令、参数、所有图版本；允许标注完全相同 SHA 与近似功能，暂不裁剪、不重写、不删除。
- **发布**：用户明确选定版本并经同输入验证后才能进入新的正式仓库和 GitHub 任务分支。

## 已可核实的 GitHub 线索

本仓库 `main` 可见通用 `workflow/r/01_build_raw_objects.R`、`01ab_ambient_branch.R`、`01b_qc_doublet.R`、`02a_integration_eda.R`、`02b_finalize_clustering.R`、`03_annotation.R`、`04_marker_discovery.R`、`04a_pseudobulk_ds.R`、`04b_composition.R`、`05_trajectory.R` 以及调控/velocity 工作流。这些仅为 **source candidates**；尚无证据说明它们是 SYF/F5 实际执行版本。

F6 的历史证据仓库为 `yqbxq/f6-legacy-analysis-audit`，模块化正式代码仓库为 `yqbxq/f6-analysis`。SYF/F5 应建立独立历史审计和独立可复用分析工作区，不要把 F6 旧分支当作正式四样本运行链。

入口：[`CODEX_TASK.md`](CODEX_TASK.md) · [`DATA_VERSION_LEDGER.md`](DATA_VERSION_LEDGER.md) · [`REGISTRY_SCHEMA.md`](REGISTRY_SCHEMA.md)。
