# 共享环境重整方案

这份文档面向两个目标：

1. 修复当前服务器上因目录迁移导致的环境债务。
2. 为后续真实项目提供一套可长期复用的共享环境基线。

它的设计依据主要来自：

- [项目流程-完整版.md](/home/yqb/code/项目流程-完整版.md)
- [README.md](/home/yqb/code/chicken_project/pipelines/scRNA_chicken/README.md)
- [SOP_标准流程与论文模式.md](/home/yqb/code/chicken_project/pipelines/scRNA_chicken/SOP_标准流程与论文模式.md)
- [README.md](/home/yqb/code/chicken_heart_repro/envs/README.md)

## 1. 先说明白：为什么这次看起来只有 R 环境炸了

这次实际暴露出来的是 `Rscript` 在运行时硬编码了旧前缀：

- 旧前缀：`/home/user_test/scRNA_pipeline_chicken/envs/conda/r_main`
- 新前缀：`/home/user_test/syf_f5/01shared_resources/envs/scRNA_pipeline_chicken/conda/r_main`

R 生态里这类问题很常见：

- `Rscript`
- `R`
- 部分包的内部脚本 / 配置
- 部分编译产物

会把安装时的绝对路径写进二进制或启动逻辑里。目录迁移后，文件虽然还在，但运行时会继续去找旧路径。

Python 环境这次没有首先爆出来，不代表它们“天然没问题”，只是当前验证到的几类命令还可运行：

- `scanorama`
- `scvelo`
- `pyscenic`

所以：

- 从“先救当前流程”角度，不一定每个 conda 前缀都要立刻重建。
- 从“给未来项目做长期基线”角度，建议把共享环境整体重整，而不只修 R。

## 2. 不建议的做法

不要把所有分析软件塞进一个巨型环境。原因很明确：

- `Seurat / Bioconductor / monocle2` 的 R 依赖链和 `SCENIC`、`pySCENIC` 不是一个时代的东西。
- `scVelo / velocyto / scanorama / scanpy / squidpy / cell2location` 的 Python 依赖跨度大，混到一个环境里会让升级和回归很难控。
- `Monocle2`、`PHATE` 这种历史方法应该留在“兼容层”，不应该污染正式主环境。

最佳实践不是“一个最大全家桶”，而是“一个最稳的分层环境体系”。

## 3. 依据未来项目流程抽出来的软件需求

从 [项目流程-完整版.md](/home/yqb/code/项目流程-完整版.md) 可以把软件需求分成 4 层。

### 3.1 主流程必需

- `Seurat`
- `Harmony`
- `DoubletFinder`
- `Slingshot`
- `tradeSeq`
- `Monocle3`
- `clusterProfiler`
- `gprofiler2`
- `biomaRt`
- `org.Gg.eg.db`
- 基础绘图与数据整理包

这是未来多数 scRNA 正式项目都会长期用到的核心层。

### 3.2 Python 主流程增强

- `scanpy`
- `scvelo`
- `velocyto`
- `scanorama`
- `scrublet`
- `loompy`
- `h5py`
- `numpy / scipy / pandas / scikit-learn`

这层既覆盖鸡心脏当前 `scanorama`，也覆盖后续 velocity 线路。

### 3.3 SCENIC 线路

- `pySCENIC`
- `ctxcore`
- `arboreto`
- R 侧的 `SCENIC / AUCell / RcisTarget / GENIE3 / ComplexHeatmap / scFunctions`

这条线不适合和主流程强行合并。

### 3.4 历史兼容线

- `Monocle2`
- `PHATE`
- `topGO`

这条线目前主要服务鸡心脏和少数历史分析，不该进主环境。

### 3.5 通讯与调控扩展

- `CellChat`
- `NicheNet`
- `decoupleR`
- `gprofiler2`

这层不是每个项目都要跑，但 [项目流程-完整版.md](/home/yqb/code/项目流程-完整版.md) 已经把它列成正式分析模块，不能再视为“项目外功能”。

### 3.6 空间扩展

- `BayesSpace`
- `SpaGCN`
- `STAGATE`
- `RCTD`
- `CARD`
- `cell2location`
- `Squidpy`
- `SpatialDE`
- `SPARK`
- `stLearn`

这层覆盖 ST 聚类增强、去卷积、邻域分析、空间差异与空间轨迹。

## 4. 推荐的共享环境结构

正式建议保留 5 套共享环境，外加 3 套可选扩展环境。

### 4.1 P0：必须长期维护

#### `r_main`

用途：

- scRNA 主流程
- Seurat-based ST 主流程
- Harmony
- annotation
- DEG / enrichment
- Slingshot / tradeSeq
- Monocle3
- 常规绘图

应该包含：

- `Seurat`
- `SeuratObject`
- `Harmony`
- `SingleCellExperiment`
- `slingshot`
- `tradeSeq`
- `Monocle3`
- `clusterProfiler`
- `gprofiler2`
- `org.Gg.eg.db`
- `biomaRt`
- `hdf5r`
- `reticulate`
- `yaml`
- `patchwork`
- `cowplot`
- `ggrepel`
- `pheatmap`
- `circlize`
- `data.table`

说明：

- 这是最核心的正式环境。
- 鸡心脏里的 ST transfer、spatial pseudotime 主流程，也应优先落在这里。

#### `py_main`

用途：

- `scanorama`
- `scanpy`
- `scvelo`
- `scrublet`
- `velocyto`
- 通用 Python 单细胞矩阵处理

建议来源：

- 以当前 `environment_velocity.yml` 为主干
- 把鸡心脏 `environment_scovel_heart.yml` 里的 `scanorama` 吸收进来

说明：

- 当前 `velocity` 这个名字已经太窄了。
- 长期建议概念上把它看成 `py_main`，即使迁移期路径名暂时还叫 `velocity`。

#### `r_scenic`

用途：

- SCENIC 下游
- AUCell
- RSS
- CSI

说明：

- 保持独立，避免把 `SCENIC` 依赖污染 `r_main`。

#### `py_scenic`

用途：

- `pySCENIC`
- `ctxcore`
- `arboreto`

说明：

- 和 `r_scenic` 配套，但不要和 `py_main` 合并。

#### `r_legacy`

用途：

- `Monocle2`
- `PHATE`
- `topGO`
- 历史兼容分析

说明：

- 这套环境只为兼容老方法存在。
- 未来新正式项目的轨迹分析应优先用 `Slingshot / Monocle3` 或更新方法，不再默认落到这里。

### 4.2 P1：建议后续补齐的扩展层

#### `r_spatial`

用途：

- `BayesSpace`
- `RCTD`
- `CARD`
- `SPARK`
- `decoupleR`

说明：

- 这层和 `r_main` 可以共享很多依赖，但不建议一开始硬塞进去。
- 等你正式把空间转录组作为长期主线后再独立维护，更稳。

#### `py_spatial`

用途：

- `Squidpy`
- `stLearn`
- 需要时再评估 `SpaGCN / STAGATE`

配套拆分：

- `cell2location` 单独放进 `py_cell2location`
- `SpatialDE` 当前按共享实现放进 `r_spatial`

说明：

- 这层目前在你的正式流程中还没有固定下来，先作为扩展层。

#### `py_cell2location`

用途：

- `cell2location`

说明：

- 单独拆层比混进 `py_spatial` 更稳。
- 这样 GPU / `pyro` 依赖不会污染通用空间 Python 环境。

#### `r_interaction`

用途：

- `CellChat`
- `NicheNet`
- `decoupleR`
- `gprofiler2`

说明：

- 这层负责细胞通讯、配体受体和调控活性扩展。
- 不建议把 `CellChat / NicheNet` 直接塞进 `r_main`，否则主环境会继续膨胀。

## 5. 什么应该放在 conda 外面

这些不建议混进共享分析环境：

- `Cell Ranger`
- `DNBelab C Series pipeline`
- `STARsolo`
- `SAW / STOmics`
- 可能需要 root / system deps 的大型上游工具

它们更适合：

- 单独安装
- 单独版本管理
- 在项目配置里显式引用路径

## 6. 对当前服务器的具体建议

### 6.1 不是“只重建 R”，而是“分两阶段重整”

#### 阶段 A：先把正式共享环境重新建干净

优先级：

1. `r_main`
2. `py_main`（迁移期路径名可暂保留 `velocity`）
3. `r_scenic`
4. `py_scenic`
5. `r_legacy`

理由：

- 这 5 套正式环境加 `r_interaction / r_spatial / py_spatial / py_cell2location`，已经覆盖鸡心脏复现和未来正式项目的大头。
- 一次重整后，后续新项目只消费共享环境，不再在实验仓里“借环境”。
- 对已经存在且能工作的 `r_main / velocity / r_scenic / pyscenic / r_legacy`，优先原地更新，不重复复制安装。
- `r_legacy` 迁移期可继续复用现有 shared prefix 路径名 `r_heart_legacy`，不用再复制出第二套 legacy 环境。

#### 阶段 B：把鸡心脏从“借旧环境”切到“消费正式共享环境”

目标：

- `chicken_heart_repro` 不再依赖旧前缀兼容链接
- 只依赖 `syf_f5/01shared_resources/envs/...`
- 兼容链接全部下线

### 6.2 当前不建议做的事

- 不要现在把所有空间方法也一次性打进主环境
- 不要把 `Monocle2` 并回 `r_main`
- 不要把 `SCENIC` 并回主流程环境
- 不要把 `CellChat / NicheNet` 并回 `r_main`

## 7. 推荐的重建验收标准

每个共享环境至少通过下面的 smoke tests。

### `r_main`

```bash
Rscript -e "library(Seurat); library(harmony); library(slingshot); library(tradeSeq); library(monocle3); library(clusterProfiler); library(gprofiler2)"
```

### `py_main`

```bash
python - <<'PY'
import scanpy, scvelo, scanorama
print(scanpy.__version__)
print(scvelo.__version__)
print(scanorama.__version__)
PY
velocyto --help >/dev/null
```

### `r_scenic`

```bash
Rscript -e "library(SCENIC); library(AUCell); library(RcisTarget); library(GENIE3); library(ComplexHeatmap)"
```

### `py_scenic`

```bash
pyscenic --help >/dev/null
python - <<'PY'
import ctxcore, arboreto
print('ok')
PY
```

### `r_legacy`

```bash
Rscript -e "library(monocle); library(phateR); library(topGO)"
```

### `r_interaction`

```bash
Rscript -e "library(CellChat); library(decoupleR); library(gprofiler2)"
```

## 8. 当前结论

如果只是修鸡心脏眼前问题，确实只重建 `r_main / r_heart_legacy` 就够。

但如果目标是支撑你后续 [项目流程-完整版.md](/home/yqb/code/项目流程-完整版.md) 里的真实项目，那么更合理的做法是：

- 不是“继续修修补补当前环境”
- 而是“把共享环境体系整体重整为 5 套正式环境 + 3 套可选扩展环境”

这才是长期最干净、最稳的版本。
