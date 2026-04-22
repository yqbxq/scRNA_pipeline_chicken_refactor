source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))
source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "plotting_helpers.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

input_rds <- file.path(cfg$checkpoint_dir, "01_after_qc_doublet.rds")
if (!file.exists(input_rds)) {
  stop("缺少输入检查点: 01_after_qc_doublet.rds", call. = FALSE)
}

obj <- readRDS(input_rds)
obj <- maybe_join_layers(obj)

DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = cfg$hvg_nfeatures, verbose = FALSE)
obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)
obj <- RunHarmony(obj, group.by.vars = "orig.ident", theta = 3, lambda = 1, plot_convergence = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = cfg$pca_dims, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = cfg$pca_dims, verbose = FALSE)

evaluate_resolution <- function(seu, res) {
  tmp <- FindClusters(
    seu,
    resolution = res,
    random.seed = cfg$random_seed,
    verbose = FALSE
  )
  list(
    object = tmp,
    n_clusters = length(levels(tmp$seurat_clusters))
  )
}

search_records <- list()
evaluated_res <- numeric(0)
found_res <- NA_real_
selected_obj <- NULL
selected_cluster_count <- NA_integer_

for (res in cfg$res_range) {
  eval_res <- evaluate_resolution(obj, res)
  search_records[[length(search_records) + 1]] <- data.frame(
    stage = "coarse",
    resolution = res,
    n_clusters = eval_res$n_clusters
  )
  evaluated_res <- c(evaluated_res, res)

  if (eval_res$n_clusters == cfg$target_clusters) {
    found_res <- res
    selected_obj <- eval_res$object
    selected_cluster_count <- eval_res$n_clusters
    break
  }
}

if (is.na(found_res)) {
  coarse_df <- do.call(rbind, search_records)
  interval_ids <- integer(0)

  if (nrow(coarse_df) >= 2) {
    for (i in seq_len(nrow(coarse_df) - 1)) {
      left_n <- coarse_df$n_clusters[i]
      right_n <- coarse_df$n_clusters[i + 1]
      lower_n <- min(left_n, right_n)
      upper_n <- max(left_n, right_n)
      if (cfg$target_clusters >= lower_n && cfg$target_clusters <= upper_n) {
        interval_ids <- c(interval_ids, i)
      }
    }
  }

  if (length(interval_ids) == 0 && nrow(coarse_df) >= 2) {
    closest_idx <- which.min(abs(coarse_df$n_clusters - cfg$target_clusters))
    if (closest_idx == nrow(coarse_df)) {
      interval_ids <- closest_idx - 1L
    } else {
      interval_ids <- closest_idx
    }
  }

  for (interval_idx in unique(interval_ids)) {
    left_res <- coarse_df$resolution[interval_idx]
    right_res <- coarse_df$resolution[interval_idx + 1]
    fine_grid <- seq(left_res, right_res, by = cfg$res_fine_step)
    fine_grid <- sort(unique(round(fine_grid, 6)))
    fine_grid <- fine_grid[!fine_grid %in% evaluated_res]

    for (res in fine_grid) {
      eval_res <- evaluate_resolution(obj, res)
      search_records[[length(search_records) + 1]] <- data.frame(
        stage = "fine",
        resolution = res,
        n_clusters = eval_res$n_clusters
      )
      evaluated_res <- c(evaluated_res, res)

      if (eval_res$n_clusters == cfg$target_clusters) {
        found_res <- res
        selected_obj <- eval_res$object
        selected_cluster_count <- eval_res$n_clusters
        break
      }
    }

    if (!is.na(found_res)) {
      break
    }
  }
}

if (is.na(found_res)) {
  search_df <- do.call(rbind, search_records)
  best_idx <- which.min(abs(search_df$n_clusters - cfg$target_clusters))
  found_res <- search_df$resolution[best_idx]
  fallback_eval <- evaluate_resolution(obj, found_res)
  selected_obj <- fallback_eval$object
  selected_cluster_count <- fallback_eval$n_clusters
  search_records[[length(search_records) + 1]] <- data.frame(
    stage = "fallback",
    resolution = found_res,
    n_clusters = selected_cluster_count
  )
}

obj <- selected_obj
current_ids <- levels(obj$seurat_clusters)
new_ids <- setNames(as.character(seq_along(current_ids)), current_ids)
obj <- RenameIdents(obj, new_ids)
obj$seurat_clusters <- factor(Idents(obj), levels = as.character(seq_along(current_ids)))

group_var <- resolve_plot_group_var(obj)
cluster_summary <- obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(plot_group = .data[[group_var]]) %>%
  count(plot_group, seurat_clusters, name = "cell_number") %>%
  group_by(plot_group) %>%
  mutate(proportion = round(cell_number / sum(cell_number) * 100, 5)) %>%
  ungroup()

search_df <- do.call(rbind, search_records) %>%
  distinct(stage, resolution, .keep_all = TRUE) %>%
  arrange(resolution, stage)

write.csv(cluster_summary, file.path(cfg$table_dir, "cluster_summary.csv"), row.names = FALSE)
write.csv(search_df, file.path(cfg$table_dir, "resolution_search.csv"), row.names = FALSE)
writeLines(
  c(
    sprintf("selected_resolution=%s", found_res),
    sprintf("selected_cluster_count=%s", selected_cluster_count),
    sprintf("target_clusters=%s", cfg$target_clusters)
  ),
  file.path(cfg$table_dir, "selected_resolution.txt")
)
saveRDS(obj, file.path(cfg$checkpoint_dir, "02_after_clustering.rds"))

tryCatch({
  figure_2a <- plot_cluster_split_umap_paper(obj, group_var = group_var, cluster_var = "seurat_clusters")
  figure_2b <- plot_cluster_proportion_paper(obj, group_var = group_var, cluster_var = "seurat_clusters")
  figure_2 <- (figure_2a | figure_2b) +
    plot_layout(widths = c(2.4, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 26, face = "bold"))

  save_plot_dual(
    plot_obj = figure_2a,
    png_path = file.path(cfg$figure_dir, "Figure_2A.png"),
    width = 14,
    height = 6
  )
  save_plot_dual(
    plot_obj = figure_2b,
    png_path = file.path(cfg$figure_dir, "Figure_2B.png"),
    width = 6.5,
    height = 6
  )
  save_plot_dual(
    plot_obj = figure_2,
    png_path = file.path(cfg$figure_dir, "Figure_2.png"),
    width = 16,
    height = 6.5
  )

  pca_elbow <- ElbowPlot(obj, ndims = max(cfg$pca_dims)) +
    ggtitle("PCA elbow plot")
  save_plot_dual(
    plot_obj = pca_elbow,
    png_path = file.path(cfg$figure_dir, "pca_elbow.png"),
    width = 7,
    height = 5
  )
}, error = function(e) {
  warning(sprintf("02_clustering 绘图失败，但检查点与表格已保存: %s", conditionMessage(e)))
})

message("已保存 02_after_clustering.rds，并输出 Figure_2")
