# 审计产物与核验契约

第一阶段应生成**实际有数据行**的登记表及只读源码快照，完整覆盖单样本、阶段内、CORE4、GC/TC 和所有发现的专项。推荐目录与字段如下；所有 stable ID 跨阶段保持不变，允许一个脚本对应多个运行/输出。

| 文件 | 必须字段 |
|---|---|
| `source_manifest.csv` | source_id, original_path, archived_path, sha256, size, git_ref, sample_scope, module, scan_status |
| `script_inventory.csv` | script_id, source_id, language, functions, source_imports, purpose, static_inputs, static_outputs, execution_status, evidence_id |
| `run_registry.csv` | run_id, script_id, sample_ids, start_time, command, parameter_manifest, environment, log_path, exit_code, run_status, evidence_id |
| `object_registry.csv` | object_id, source_path, sha256, sample_scope, cellset_id, cell_count, assay, layers, cluster_key, label_authority_id, verification |
| `object_lineage.csv` | parent_object_id, operation_id, child_object_id, run_id, transformation, evidence_id, confidence |
| `table_registry.csv` | table_id, path, sha256, table_type, object_id, producer_script_id, run_id, evidence_id, confidence |
| `figure_registry.csv` | figure_id, path, sha256, plot_type, data_object_ids, table_ids, producer_script_id, run_id, theme_version, selection_status, evidence_id |
| `analysis_families.csv` | family_id, module, question, candidate_script_ids, candidate_object_ids, candidate_output_ids, difference_summary, approval_status |
| `duplicate_implementations.csv` | duplicate_group_id, script_ids, same_sha, same_purpose, changed_inputs_or_params, risk_note, recommended_reuse |
| `coverage_matrix.csv` | module, inventory_count, traced_input_count, traced_run_count, traced_output_count, unresolved_count, user_review_status |
| `excluded_side_branches.csv` | source_id, sample_scope, reason, legacy_result_path, allowed_reuse, evidence_id |
| `approval_queue.csv` | family_id, candidate_version_ids, figure_preview_refs, differences, historical_user_preference, decision, decision_evidence |
| `open_gaps.md` | unresolved paths, permissions, version conflicts, missing logs, uncertain links, next evidence query |

额外必备：`trace_graph.csv` （`producer → consumed input → run → object → table/figure`）、`review_gallery.html`（每个图片族**全部版本**预览）、`module_status.md`、`validation.md`、`git_snapshot.md`、可复现的验证脚本、`git_sha_receipts.md`。

- **evidence confidence**：`confirmed_by_execution` / `reproduced_under_control` / `supported_by_indirect_evidence` / `unresolved`。
- **execution status**：`completed` / `failed` / `skipped` / `not_run` / `unverified`。
- **selection status**：`user_selected`（需明确历史选择凭据）/ `candidate` / `alternative` / `superseded` / `pending_review`。历史“终版/final”文件名本身不构成选择凭据。

## 审计与发布门槛

1. 基于 SHA 保留每份原始小体积源码，重复版本也逐件登记。对大型 RDS、h5ad、FASTQ 等保留路径、校验值和只读验证结果；默认不入 GitHub。
2. 通过 `source()`、`readRDS()`、`saveRDS()`、`ggsave()`、`write.table()`、shell/job 历史和日志取证。静态解析只是候选，不足以确定真实运行。
3. 证据图双向可查询；随机抽查一个图 → 生成代码/参数/表/对象，随机抽查一个代码 → 历史各运行/输出，失败或找不到的在 gaps 中记录。
4. 候选代码选版和图版本由用户确认；正式模块做可替换输入契约、科学方法差异 adapter、单元/回归测试和跨样本防混门禁。
5. 仅在 Git 工作树状态可核查、无凭据/大对象泄露且测试**实际运行**有收据后，经 `safe-git-sync` 向任务分支推送。禁止修改 F6 两仓库、禁止重写历史服务器原件、禁止对未知格式对象做破坏性操作。
