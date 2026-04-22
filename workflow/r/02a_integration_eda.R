source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))
source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "object_layer_helpers.R"))
source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "plotting_helpers.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

layer_specs <- read_object_layer_specs(cfg)
panorama_id <- panorama_layer_id(layer_specs)
panorama_spec <- layer_specs[[panorama_id]]
panorama_paths <- object_layer_paths(cfg, panorama_id)
ensure_object_layer_dirs(panorama_paths)

input_rds <- panorama_paths$candidate_rds
if (!file.exists(input_rds)) {
  input_rds <- file.path(cfg$checkpoint_dir, "02_reduction_candidates.rds")
}
if (!file.exists(input_rds)) {
  stop("缺少输入检查点: panorama reduction candidates", call. = FALSE)
}

compat_stage_dir <- ensure_eda_stage_dir(cfg, "integration")
stage_dir <- panorama_paths$report_dir
dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)
obj <- maybe_join_layers(obj)

has_harmony <- "harmony" %in% Reductions(obj)
selected_reduction <- selected_layer_reduction(obj, panorama_spec)
resolution_dims <- usable_reduction_dims(obj, panorama_spec$pca_dims, selected_reduction)

cell_cycle_genes <- intersect(c("TOP2A", "MKI67", "PCNA", "UBE2C", "CENPF", "TYMS"), rownames(obj))
stress_genes <- intersect(c("FOS", "JUN", "JUNB", "HSPA1A", "HSP90AA1", "ATF3"), rownames(obj))

loadings_mat <- Loadings(obj, reduction = "pca")
pc_ids <- seq_len(min(5, ncol(loadings_mat)))
driver_rows <- bind_rows(lapply(pc_ids, function(pc_idx) {
  pc_name <- colnames(loadings_mat)[pc_idx]
  ord <- order(abs(loadings_mat[, pc_idx]), decreasing = TRUE)
  top_genes <- rownames(loadings_mat)[ord][seq_len(min(20, nrow(loadings_mat)))]
  data.frame(
    pc = pc_name,
    top_genes = paste(top_genes, collapse = ","),
    cell_cycle_hits = paste(intersect(top_genes, cell_cycle_genes), collapse = ","),
    stress_hits = paste(intersect(top_genes, stress_genes), collapse = ","),
    n_cell_cycle_hits = length(intersect(top_genes, cell_cycle_genes)),
    n_stress_hits = length(intersect(top_genes, stress_genes)),
    stringsAsFactors = FALSE
  )
}))

r2_by_factor <- function(values, group) {
  if (length(unique(stats::na.omit(group))) < 2) {
    return(NA_real_)
  }
  fit <- stats::lm(values ~ as.factor(group))
  summary(fit)$r.squared
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  max(x)
}

pca_embed <- Embeddings(obj, reduction = "pca")
pc_association <- bind_rows(lapply(pc_ids, function(pc_idx) {
  pc_name <- colnames(pca_embed)[pc_idx]
  data.frame(
    reduction = "pca",
    component = pc_name,
    sample_r2 = r2_by_factor(pca_embed[, pc_idx], obj$orig.ident),
    group_r2 = if ("analysis_group" %in% colnames(obj@meta.data)) r2_by_factor(pca_embed[, pc_idx], obj$analysis_group) else NA_real_,
    stringsAsFactors = FALSE
  )
}))

if (has_harmony) {
  harmony_embed <- Embeddings(obj, reduction = "harmony")
  harmony_ids <- seq_len(min(5, ncol(harmony_embed)))
  pc_association <- bind_rows(
    pc_association,
    bind_rows(lapply(harmony_ids, function(comp_idx) {
      comp_name <- colnames(harmony_embed)[comp_idx]
      data.frame(
        reduction = "harmony",
        component = comp_name,
        sample_r2 = r2_by_factor(harmony_embed[, comp_idx], obj$orig.ident),
        group_r2 = if ("analysis_group" %in% colnames(obj@meta.data)) r2_by_factor(harmony_embed[, comp_idx], obj$analysis_group) else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  )
}

evaluate_resolution <- function(seu, reduction_name, dims, res) {
  tmp <- FindNeighbors(
    seu,
    reduction = reduction_name,
    dims = dims,
    verbose = FALSE
  )
  tmp <- FindClusters(
    tmp,
    resolution = res,
    random.seed = cfg$random_seed,
    verbose = FALSE
  )
  data.frame(
    reduction = reduction_name,
    resolution = res,
    n_clusters = length(levels(tmp$seurat_clusters)),
    stringsAsFactors = FALSE
  )
}

resolution_rows <- bind_rows(lapply(panorama_spec$res_range, function(res) {
  evaluate_resolution(obj, selected_reduction, resolution_dims, res)
}))

triage_rows <- list()
append_triage <- function(severity, signal_id, evidence, recommended_action) {
  triage_rows[[length(triage_rows) + 1]] <<- data.frame(
    severity = severity,
    signal_id = signal_id,
    evidence = evidence,
    recommended_action = recommended_action,
    manual_review_required = "yes",
    stringsAsFactors = FALSE
  )
}

for (i in seq_len(nrow(driver_rows))) {
  row <- driver_rows[i, , drop = FALSE]
  if (row$n_cell_cycle_hits >= 2) {
    append_triage(
      "medium",
      paste0("cell_cycle_axis_", row$pc),
      sprintf("%s 包含 %s 个细胞周期高载荷基因", row$pc, row$n_cell_cycle_hits),
      "如果 cell cycle 遮蔽了你的核心科学问题，再考虑回归，而不是默认直接回归。"
    )
  }
  if (row$n_stress_hits >= 2) {
    append_triage(
      "medium",
      paste0("stress_axis_", row$pc),
      sprintf("%s 包含 %s 个 stress/IEG 高载荷基因", row$pc, row$n_stress_hits),
      "优先怀疑 dissociation stress，而不是立刻把它解释成新状态。"
    )
  }
}

pca_sample_dom <- pc_association %>%
  filter(reduction == "pca") %>%
  summarise(max_sample_r2 = safe_max(sample_r2), max_group_r2 = safe_max(group_r2))

if (nrow(pca_sample_dom) > 0 && is.finite(pca_sample_dom$max_sample_r2) && pca_sample_dom$max_sample_r2 >= 0.30 &&
    (is.na(pca_sample_dom$max_group_r2) || pca_sample_dom$max_group_r2 < pca_sample_dom$max_sample_r2)) {
  append_triage(
    "medium",
    "batch_dominant_pca_axes",
    sprintf("PCA 中 sample_r2 最大值=%.3f, group_r2 最大值=%.3f", pca_sample_dom$max_sample_r2, pca_sample_dom$max_group_r2),
    "优先检查是否需要 batch correction / integration，再决定是否沿用未整合空间。"
  )
}

if (has_harmony) {
  pca_group <- pc_association %>% filter(reduction == "pca") %>% summarise(max_group_r2 = safe_max(group_r2))
  harmony_group <- pc_association %>% filter(reduction == "harmony") %>% summarise(max_group_r2 = safe_max(group_r2))
  if (nrow(pca_group) > 0 && nrow(harmony_group) > 0 &&
      is.finite(pca_group$max_group_r2) && is.finite(harmony_group$max_group_r2) &&
      pca_group$max_group_r2 >= 0.15 && harmony_group$max_group_r2 < pca_group$max_group_r2 * 0.5) {
    append_triage(
      "medium",
      "possible_overcorrection",
      sprintf("group_r2 从 PCA %.3f 降到 Harmony %.3f", pca_group$max_group_r2, harmony_group$max_group_r2),
      "如果整合后条件结构被明显抹平，优先改用更保守策略或退回未整合分析。"
    )
  }
}

triage_df <- if (length(triage_rows) > 0) bind_rows(triage_rows) else data.frame(
  severity = character(0),
  signal_id = character(0),
  evidence = character(0),
  recommended_action = character(0),
  manual_review_required = character(0),
  stringsAsFactors = FALSE
)

if ("umap_rna" %in% Reductions(obj)) {
  umap_rna <- as.data.frame(Embeddings(obj, reduction = "umap_rna"))
  colnames(umap_rna)[1:2] <- c("UMAP_1", "UMAP_2")
  umap_rna$orig.ident <- obj$orig.ident
  umap_sample_plot <- ggplot(umap_rna, aes(x = UMAP_1, y = UMAP_2, color = orig.ident)) +
    geom_point(size = 0.25, alpha = 0.7) +
    theme_classic(base_size = 11) +
    labs(title = "Panorama RNA UMAP by sample", x = "UMAP_1", y = "UMAP_2", color = NULL)
  save_plot_dual(umap_sample_plot, panorama_paths$umap_rna_by_sample_png, width = 7, height = 6)
}

if (has_harmony) {
  umap_harmony <- as.data.frame(Embeddings(obj, reduction = "umap_harmony"))
  colnames(umap_harmony)[1:2] <- c("UMAP_1", "UMAP_2")
  umap_harmony$orig.ident <- obj$orig.ident
  harmony_plot <- ggplot(umap_harmony, aes(x = UMAP_1, y = UMAP_2, color = orig.ident)) +
    geom_point(size = 0.25, alpha = 0.7) +
    theme_classic(base_size = 11) +
    labs(title = "Panorama Harmony UMAP by sample", x = "UMAP_1", y = "UMAP_2", color = NULL)
  save_plot_dual(harmony_plot, panorama_paths$umap_harmony_by_sample_png, width = 7, height = 6)
}

elbow_plot <- ElbowPlot(obj, ndims = max(panorama_spec$pca_dims)) +
  ggtitle("Panorama PCA elbow plot")
save_plot_dual(elbow_plot, panorama_paths$pca_elbow_png, width = 7, height = 5)

resolution_plot <- ggplot(resolution_rows, aes(x = resolution, y = n_clusters)) +
  geom_line(color = "#2F6CB3") +
  geom_point(color = "#2F6CB3", size = 2) +
  geom_hline(yintercept = panorama_spec$target_clusters, linetype = "dashed", color = "#D62728") +
  theme_classic(base_size = 11) +
  labs(
    title = sprintf("Panorama resolution candidates on %s", selected_reduction),
    x = "Resolution",
    y = "Cluster count"
  )
save_plot_dual(resolution_plot, panorama_paths$resolution_plot_png, width = 7, height = 5)

write_tsv(driver_rows, panorama_paths$pca_driver_summary_tsv)
write_tsv(pc_association, panorama_paths$pc_association_tsv)
write_tsv(resolution_rows, panorama_paths$resolution_candidates_tsv)
write_tsv(triage_df, panorama_paths$triage_tsv)

report_lines <- c(
  "# Panorama Integration / Resolution EDA Report",
  "",
  sprintf("- layer_id: `%s`", panorama_spec$layer_id),
  sprintf("- object layer config: `%s`", cfg$object_layer_config_file),
  sprintf("- 配置的 integration_mode: `%s`", panorama_spec$integration_mode),
  sprintf("- 实际候选 reduction: `%s`", selected_reduction),
  sprintf("- target_clusters: `%s`", panorama_spec$target_clusters),
  sprintf("- pca_dims: `%s`", format_index_spec(panorama_spec$pca_dims)),
  sprintf("- Harmony 可用: `%s`", ifelse(has_harmony, "yes", "no")),
  sprintf("- triage 信号条数: `%s`", nrow(triage_df)),
  "",
  "## Key Files",
  sprintf("- `pca_driver_summary.tsv`: `%s`", panorama_paths$pca_driver_summary_tsv),
  sprintf("- `pc_association.tsv`: `%s`", panorama_paths$pc_association_tsv),
  sprintf("- `resolution_candidates.tsv`: `%s`", panorama_paths$resolution_candidates_tsv),
  sprintf("- `integration_triage.tsv`: `%s`", panorama_paths$triage_tsv),
  "",
  "## Review Focus",
  "- 先判断主导变异轴是 biology 还是 nuisance，再决定是否继续沿用当前 integration_mode。",
  "- 当前 gate 审的是 panorama；后续子对象层会在 panorama finalize 之后按 object_layers.tsv 继续拆分重建。",
  "- resolution 候选曲线只提供证据，不自动改写其它 layer 的 target_clusters。"
)

if (nrow(triage_df) > 0) {
  report_lines <- c(report_lines, "", "## Triage Summary")
  for (i in seq_len(nrow(triage_df))) {
    row <- triage_df[i, , drop = FALSE]
    report_lines <- c(report_lines, sprintf("- [%s] `%s`: %s", row$severity, row$signal_id, row$recommended_action))
  }
}

write_markdown(report_lines, panorama_paths$report_md)
write_markdown(report_lines, file.path(compat_stage_dir, "report.md"))

message("panorama integration EDA 已输出到: ", panorama_paths$report_md)
message("兼容 gate 报告已更新: ", file.path(compat_stage_dir, "report.md"))
