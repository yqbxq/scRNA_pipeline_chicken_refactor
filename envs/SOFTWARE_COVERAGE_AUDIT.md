# 项目流程软件覆盖审计

这份审计只对照 [项目流程-完整版.md](/home/yqb/code/项目流程-完整版.md) 里明确提到的软件，不掺别的仓库需求。

审计对象：

- 共享环境定义：
  - [environment_r_main.yml](/home/yqb/code/chicken_project/pipelines/scRNA_chicken/envs/environment_r_main.yml)
  - [environment_r_legacy.yml](/home/yqb/code/chicken_project/pipelines/scRNA_chicken/envs/environment_r_legacy.yml)
  - [environment_velocity.yml](/home/yqb/code/chicken_project/pipelines/scRNA_chicken/envs/environment_velocity.yml)
  - [environment_r_scenic.yml](/home/yqb/code/chicken_project/pipelines/scRNA_chicken/envs/environment_r_scenic.yml)
  - [environment_pyscenic.yml](/home/yqb/code/chicken_project/pipelines/scRNA_chicken/envs/environment_pyscenic.yml)

## 结论

不是都包含在内。

当前可以分成 3 类：

1. 已经纳入共享环境或鸡心脏兼容环境。
2. 已经补成扩展环境定义，但还需要远端安装和验证。
3. 本来就不该塞进 conda 共享环境，应该单独安装。

## 1. 已纳入

| 软件 | 当前状态 | 建议落点 | 说明 |
| --- | --- | --- | --- |
| Seurat | 已纳入 | `r_main` | `SCTransform`、标准化、聚类、ST 基础分析都归这里 |
| Harmony | 已纳入 | `r_main` | 主流程整合 |
| DoubletFinder | 已纳入 | `r_main` | 由 GitHub 固定提交安装 |
| Scanpy | 已纳入 | `velocity` | 当前 Python 主单细胞环境 |
| scVelo | 已纳入 | `velocity` | RNA velocity |
| velocyto | 已纳入 | `velocity` | RNA velocity 上游接口 |
| Scanorama | 已补入共享定义 | `velocity` | 之前只在鸡心脏 overlay 里，现在提升到共享 Python 环境 |
| Scrublet | 已补入共享定义 | `velocity` | 作为 DoubletFinder 的 Python 备选 |
| Slingshot | 已纳入 | `r_main` | 当前正式默认轨迹 |
| tradeSeq | 已纳入 | `r_main` | 当前正式默认轨迹统计 |
| Monocle3 | 已补入共享定义 | `r_main` | 作为正式轨迹验证/扩展后端 |
| Monocle2 | 已纳入 | `r_legacy` | 仅历史兼容和鸡心脏论文复现，迁移期共享 prefix 路径名仍是 `r_heart_legacy` |
| clusterProfiler | 已纳入 | `r_main` | scRNA/ST 富集主工具 |
| gprofiler2 | 已补入共享定义 | `r_main` | 作为 clusterProfiler 的补充路线 |
| SCENIC | 已纳入 | `r_scenic` | R 侧下游整合 |
| pySCENIC | 已纳入 | `pyscenic` | Python 侧 regulon 推断 |

## 2. 已补成扩展环境定义，待安装验证

| 软件 | 当前状态 | 建议落点 | 说明 |
| --- | --- | --- | --- |
| CellChat | 已定义 | `r_interaction` | 细胞通讯分析，不该塞进 `r_main` |
| NicheNet | 已定义 | `r_interaction` | 和 CellChat 一样，独立更稳 |
| decoupleR | 已定义 | `r_interaction` / `r_spatial` | scRNA/ST 调控活性共用 |
| BayesSpace | 已定义 | `r_spatial` | 空间 domain 增强 |
| RCTD | 已定义 | `r_spatial` | 通过 `spacexr` 提供 |
| CARD | 已定义 | `r_spatial` | 去卷积备选 |
| SPARK | 已定义 | `r_spatial` | 空间变异基因 |
| Squidpy | 已定义 | `py_spatial` | 空间邻域与图分析 |
| cell2location | 已定义 | `py_cell2location` | 概率去卷积单独分层 |
| SpatialDE | 已定义 | `r_spatial` | 通过 Bioconductor `spatialDE` |
| SpaGCN | 已定义 | `py_spatial_legacy` | 老方法，和现代 spatial stack 分开 |
| STAGATE | 已定义 | `py_spatial_legacy` | 通过 `STAGATE_pyG` 分层 |
| stLearn | 已定义 | `py_spatial` | 现代空间 Python 栈 |
| scDesign3 | 已定义 | `r_validation` | 验证/模拟单独分层，不污染主环境 |

## 3. 不建议放进共享 conda 环境

| 软件 | 当前状态 | 建议落点 | 说明 |
| --- | --- | --- | --- |
| DNBelab C Series pipeline | 不进共享环境 | 独立上游工具目录 | 上游定量工具，和分析环境分离 |
| STARsolo | 不进共享环境 | 独立上游工具目录 | 上游比对/计数工具，版本和参考强绑定 |
| SAW / STOmics | 不进共享环境 | 独立上游工具目录 | ST 上游工具链，不适合并入分析环境 |

## 4. 映射关系说明

`项目流程-完整版.md` 里的一些写法其实不是独立软件：

- `SCTransform`、`Seurat integration`、`Seurat NormalizeData/ScaleData/FindVariableFeatures`
  - 都归 `Seurat`
- `Monocle`
  - 正式新项目建议理解为 `Monocle3`
  - 鸡心脏复现里的论文兼容线仍然是 `Monocle2`
- `SCENIC / pySCENIC`
  - 分成 `r_scenic + pyscenic`

## 5. 当前判断

如果按“鸡心脏复现 + 未来正式项目”这个目标来衡量：

- 主干环境已经有了
- 但还没有做到“项目流程-完整版.md 里提到的软件全部可落地”

现在的缺口分成两类：

1. 结构性缺口。
2. 远端安装和 smoke test 缺口。

其中：

- `r_legacy` 之前属于结构性缺口，已经补进共享定义和安装入口。
- `r_main / r_scenic / r_interaction` 目前仍然是远端 smoke test 缺口。

接下来要做的是：

- 在共享根下原地更新 `r_main / r_legacy / velocity / r_scenic / pyscenic`
- 新建并验证 `r_interaction / r_spatial / py_spatial / py_spatial_legacy / py_cell2location / r_validation`
