source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(tradeSeq)
  library(pheatmap)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

input_rds <- file.path(cfg$checkpoint_dir, "03_after_annotation.rds")
if (!file.exists(input_rds)) {
  stop("缺少输入检查点: 03_after_annotation.rds", call. = FALSE)
}

obj <- readRDS(input_rds)
obj <- maybe_join_layers(obj)
sce_base <- as.SingleCellExperiment(obj)
umap_all <- Embeddings(obj, "umap")
lineage_cols <- colorRampPalette(c(
  "#2C7BB6", "#00A6CA", "#00CCBC", "#90EB9D",
  "#FFFF8C", "#F9D057", "#F29E2E", "#E76818", "#D7191C"
))(100)

get_label_vector <- function(seu, label_var) {
  if (!label_var %in% colnames(seu@meta.data)) {
    stop(sprintf("轨迹标签列不存在: %s", label_var), call. = FALSE)
  }
  as.character(seu@meta.data[[label_var]])
}

find_fine_start_cluster <- function(seu, fine_label_var, coarse_label_var, coarse_start, explicit_cluster = "") {
  fine_labels <- get_label_vector(seu, fine_label_var)
  fine_labels <- fine_labels[!is.na(fine_labels) & nzchar(fine_labels)]

  if (nzchar(explicit_cluster) && explicit_cluster %in% fine_labels) {
    return(explicit_cluster)
  }

  if (coarse_start %in% fine_labels) {
    return(coarse_start)
  }

  if (coarse_label_var %in% colnames(seu@meta.data)) {
    coarse_labels <- get_label_vector(seu, coarse_label_var)
    idx <- which(coarse_labels == coarse_start & !is.na(seu@meta.data[[fine_label_var]]))
    if (length(idx) > 0) {
      candidate_tab <- sort(table(as.character(seu@meta.data[[fine_label_var]][idx])), decreasing = TRUE)
      if (length(candidate_tab) > 0) {
        return(names(candidate_tab)[1])
      }
    }
  }

  overall_tab <- sort(table(fine_labels), decreasing = TRUE)
  names(overall_tab)[1]
}

run_slingshot_with_labels <- function(sce_obj, label_vec, start_label) {
  slingshot(
    sce_obj,
    clusterLabels = label_vec,
    reducedDim = "UMAP",
    start.clus = start_label,
    stretch = 0
  )
}

write_pseudotime_csv <- function(sce_obj, out_csv) {
  pseudotime <- as.data.frame(slingPseudotime(sce_obj))
  pseudotime$cell_id <- rownames(pseudotime)
  write.csv(pseudotime, out_csv, row.names = FALSE)
}

pick_top_lineages <- function(sce_obj, top_n = 4) {
  pt_mat <- slingPseudotime(sce_obj)
  if (is.null(dim(pt_mat)) || ncol(pt_mat) == 0) {
    return(integer(0))
  }
  lineage_support <- vapply(seq_len(ncol(pt_mat)), function(i) {
    sum(!is.na(pt_mat[, i]))
  }, numeric(1))
  order(lineage_support, decreasing = TRUE)[seq_len(min(top_n, length(lineage_support)))]
}

plot_lineage_panels <- function(sce_obj, seu, out_png, lineage_indices = NULL, title_prefix = "Lineage") {
  pt_mat <- slingPseudotime(sce_obj)
  curve_list <- slingCurves(sce_obj)

  if (is.null(dim(pt_mat)) || ncol(pt_mat) == 0) {
    warning(sprintf("未检测到可绘制的 lineage: %s", out_png))
    return(invisible(NULL))
  }

  if (is.null(lineage_indices)) {
    lineage_indices <- seq_len(ncol(pt_mat))
  }
  lineage_indices <- lineage_indices[lineage_indices >= 1 & lineage_indices <= ncol(pt_mat)]
  if (length(lineage_indices) == 0) {
    warning(sprintf("lineage_indices 为空，跳过作图: %s", out_png))
    return(invisible(NULL))
  }

  png(
    out_png,
    width = max(420 * length(lineage_indices), 900),
    height = 430,
    res = 150
  )
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(mfrow = c(1, length(lineage_indices)), mar = c(4, 4, 3, 1))

  for (panel_idx in seq_along(lineage_indices)) {
    lineage_idx <- lineage_indices[panel_idx]
    plot(
      umap_all,
      col = "#D8E0F2",
      pch = 16,
      cex = 0.45,
      xlab = "UMAP_1",
      ylab = "UMAP_2",
      main = paste(title_prefix, panel_idx)
    )

    pt_vec <- pt_mat[, lineage_idx]
    in_lineage <- !is.na(pt_vec)
    if (any(in_lineage)) {
      pt_rescaled <- pt_vec[in_lineage]
      pt_rescaled <- (pt_rescaled - min(pt_rescaled)) / (max(pt_rescaled) - min(pt_rescaled) + 1e-8)
      pt_cols <- lineage_cols[pmax(1, pmin(100, floor(pt_rescaled * 99) + 1))]
      points(umap_all[in_lineage, , drop = FALSE], col = pt_cols, pch = 16, cex = 0.55)
    }
    lines(curve_list[[lineage_idx]], lwd = 2, col = "black")
  }
}

# 1. 粗粒度轨迹：cell_type
coarse_label_var <- cfg$trajectory_coarse_label
coarse_start <- cfg$trajectory_start
coarse_labels <- get_label_vector(obj, coarse_label_var)
sce_coarse <- run_slingshot_with_labels(sce_base, coarse_labels, coarse_start)
write_pseudotime_csv(
  sce_coarse,
  file.path(cfg$table_dir, sprintf("slingshot_pseudotime_%s.csv", coarse_label_var))
)
plot_lineage_panels(
  sce_coarse,
  obj,
  file.path(cfg$figure_dir, sprintf("Figure_5F_Lineages_%s.png", coarse_label_var))
)

# 保留旧文件名，兼容已有流程；默认指向粗粒度 pseudotime
write_pseudotime_csv(sce_coarse, file.path(cfg$table_dir, "slingshot_pseudotime.csv"))

# 2. 细粒度轨迹：seurat_clusters
fine_label_var <- cfg$trajectory_fine_label
fine_labels <- get_label_vector(obj, fine_label_var)
fine_start <- find_fine_start_cluster(
  seu = obj,
  fine_label_var = fine_label_var,
  coarse_label_var = coarse_label_var,
  coarse_start = coarse_start,
  explicit_cluster = cfg$trajectory_fine_start_cluster
)
sce_fine <- run_slingshot_with_labels(sce_base, fine_labels, fine_start)
write_pseudotime_csv(
  sce_fine,
  file.path(cfg$table_dir, sprintf("slingshot_pseudotime_%s.csv", fine_label_var))
)
fine_top_lineages <- pick_top_lineages(sce_fine, top_n = cfg$trajectory_fine_top_n)
plot_lineage_panels(
  sce_fine,
  obj,
  file.path(cfg$figure_dir, sprintf("Figure_5F_Lineages_%s.png", fine_label_var)),
  lineage_indices = seq_len(ncol(slingPseudotime(sce_fine)))
)

# 旧文件名 Figure_5F_Lineages.png 默认输出细粒度 top-N 版本，更贴近论文展示习惯
plot_lineage_panels(
  sce_fine,
  obj,
  file.path(cfg$figure_dir, "Figure_5F_Lineages.png"),
  lineage_indices = fine_top_lineages
)

# 3. tradeSeq 默认走粗粒度，更稳
tradeseq_label_var <- cfg$trajectory_tradeseq_label
sce_for_gam <- if (identical(tradeseq_label_var, fine_label_var)) sce_fine else sce_coarse

counts_matrix <- as.matrix(counts(sce_for_gam))
hvg_genes <- intersect(VariableFeatures(obj), rownames(counts_matrix))
counts_matrix <- counts_matrix[hvg_genes, , drop = FALSE]

sce_gam <- fitGAM(
  counts = counts_matrix,
  pseudotime = slingPseudotime(sce_for_gam, na = FALSE),
  cellWeights = slingCurveWeights(sce_for_gam),
  nknots = cfg$tradeseq_knots,
  verbose = FALSE
)

asso_res <- associationTest(sce_gam)
asso_res <- asso_res[order(asso_res$pvalue), , drop = FALSE]
write.csv(
  cbind(gene = rownames(asso_res), as.data.frame(asso_res)),
  file.path(cfg$table_dir, "tradeseq_association.csv"),
  row.names = FALSE
)

top100_genes <- rownames(asso_res)[seq_len(min(100, nrow(asso_res)))]
pt_vec <- slingPseudotime(sce_for_gam, na = FALSE)[, 1]
cell_order <- order(pt_vec)
data_matrix <- as.matrix(get_assay_matrix(obj, type = "data"))
heatmap_matrix <- data_matrix[top100_genes, cell_order, drop = FALSE]
heatmap_matrix <- t(scale(t(heatmap_matrix)))
heatmap_matrix[heatmap_matrix > 3] <- 3
heatmap_matrix[heatmap_matrix < -3] <- -3

png(file.path(cfg$figure_dir, "trajectory_top100_heatmap.png"), width = 1200, height = 1200, res = 150)
pheatmap(
  heatmap_matrix,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "Top 100 genes along pseudotime"
)
dev.off()

saveRDS(sce_coarse, file.path(cfg$checkpoint_dir, sprintf("slingshot_%s.rds", coarse_label_var)))
saveRDS(sce_fine, file.path(cfg$checkpoint_dir, sprintf("slingshot_%s.rds", fine_label_var)))
saveRDS(sce_gam, file.path(cfg$checkpoint_dir, "tradeseq_fitgam.rds"))

message(
  sprintf(
    "轨迹分析完成：粗粒度=%s (%d 条), 细粒度=%s (%d 条), tradeSeq=%s",
    coarse_label_var,
    ncol(slingPseudotime(sce_coarse)),
    fine_label_var,
    ncol(slingPseudotime(sce_fine)),
    tradeseq_label_var
  )
)
