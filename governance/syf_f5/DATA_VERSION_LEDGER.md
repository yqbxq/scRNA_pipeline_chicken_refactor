# SYF–F5 输入版本与细胞集基线（待服务器核对）

此表只记录**已有历史报告的描述**。运行前必须用对象元数据、barcode 清单、脚本及日志核验，不得据此宣告最终联合版本已批准。

| 版本/用途 | SYF_1 | SYF_3 | F5_1 | F5_3 | 说明 |
|---|---:|---:|---:|---:|---|
| 2026-09-02 统一单样本注释表来源 | 10,396 | 11,160 | 8,582 | 9,968 | 混合不同 QC 层级；SYF_1 为旧 post-QC singlets，其余三份可追溯到完整 called-cell SoupX unit |
| `stage_integration_2plus2_v1` 阶段内整合输入 | 14,553 | 9,842 | 7,838 | 8,439 | 四份样本使用 RAW → scaled-MAD QC → RAW scDblFinder → 按 final singlet barcode 提取历史 SoupX corrected counts |

阶段内整合输入合计：SYF **24,395**；F5 **16,277**；CORE4 合计 **40,672**。两行版本是不同 cell-set universe，不能直接按“整合前后细胞变化”比较。

已知历史方案：每阶段分别运行 NONE、RPCA、Harmony、FastMNN；统一 HVG=2,000、LogNormalize scale.factor=10,000、dims=30、k=20、seed=42，并比较多个 resolution。Harmony res0.8 的注释/DEG历史结果需分别定位生成代码、日志、authority 与对象。四样本最终 JOINT vs GROUP、最终对象的批准/冻结状态仍需进一步核查。

## 注释 authority 线索

- 历史正式四样本联合整合注释工作簿：`SYF_F5_四样本整合统一三级细胞注释_终版_更新_20260908.xlsx`。
- 历史 MD5 记录：`35d55affd794e288af48ffa411e78689`。使用前重新计算原件校验值；检查 `annotation_authority_fix` 路线及 cluster/cell label crosswalk。
- 单样本另有 `SYF_F5_统一三级细胞注释_含判断依据图片.xlsx`，不可默认为与整合版同一套细胞/注释权威。
- 既往 `GC_TC_boundary` 大类 taxonomy 与修复后的 `corrected_fixed_labels` 应作为不同标签版本保存，不可静默混淆。

## 严格区分的旁支

`F5_2`、`SYF_2` 属疑难样本；`F6_1` 是 pilot/challenge。其 rescue、mapping、add-in、stress/sensitivity、stage-similarity 相关代码保留在侧支来源目录和排除清单内，正式 CORE4 不将其作为正常重复。

## 本地现场检索线索（不公开用户名和认证配置）

- 工作站 Codex 任务目录名匹配 `*-syf-f5`；远端项目相对目录 `05projects/syf_f5`，尤其 `results/`、`workflow/`。
- 关联审计：`stage_integration_2plus2_v1`、`decontamination_comparison`、`CORE4`、`annotation_authority_fix`、`F6_CORE4_bad_samples_reaudit_20260925`。
- 复用的通用工作流与当前项目现场要分别记 SHA；GitHub 文件相同名字不等于实际运行脚本。
