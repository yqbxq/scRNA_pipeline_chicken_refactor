source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))
source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "annotation_stats_helpers.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

layer_status_df <- built_layer_status(cfg)
comparison_df <- read_comparison_sheet(cfg)

if (nrow(comparison_df) == 0) {
  stop("comparisons.tsv 中没有启用的比较设计。", call. = FALSE)
}

marker_root_dir <- file.path(cfg$table_dir, "marker_discovery")
dir.create(marker_root_dir, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(marker_root_dir, "marker_discovery_manifest.tsv")
stage_dir <- ensure_eda_stage_dir(cfg, "marker_discovery")

choose_annotated_object_path <- function(layer_row) {
  ann_paths <- layer_annotation_paths(cfg, layer_row$layer_id[[1]])
  if (file.exists(ann_paths$annotated_rds)) {
    return(ann_paths$annotated_rds)
  }
  layer_row$clustered_rds[[1]]
}

run_exploratory_findmarkers <- function(seu, cluster_id, comparison_row, group_var) {
  hit <- as.character(seu$cluster_id) == cluster_id
  cluster_obj <- subset(seu, cells = colnames(seu)[hit])
  cluster_obj <- maybe_join_layers(cluster_obj)
  cluster_meta <- cluster_obj@meta.data
  groups <- unique(as.character(cluster_meta[[group_var]]))
  groups <- groups[nzchar(groups)]
  if (!all(c(comparison_row$ident_1[[1]], comparison_row$ident_2[[1]]) %in% groups)) {
    return(NULL)
  }

  group_counts <- table(as.character(cluster_meta[[group_var]]))
  if (any(group_counts[c(comparison_row$ident_1[[1]], comparison_row$ident_2[[1]])] < 3)) {
    return(NULL)
  }

  Idents(cluster_obj) <- group_var
  res <- tryCatch(
    FindMarkers(
      cluster_obj,
      ident.1 = comparison_row$ident_1[[1]],
      ident.2 = comparison_row$ident_2[[1]],
      logfc.threshold = 0,
      test.use = "wilcox",
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(res) || nrow(res) == 0) {
    return(NULL)
  }
  tibble::rownames_to_column(res, "gene")
}

manifest_rows <- list()

for (idx in seq_len(nrow(layer_status_df))) {
  layer_row <- layer_status_df[idx, , drop = FALSE]
  layer_id <- layer_row$layer_id[[1]]
  obj_path <- choose_annotated_object_path(layer_row)
  obj <- readRDS(obj_path)
  obj <- maybe_join_layers(obj)
  obj <- standardize_design_metadata(obj)

  cluster_var <- if ("cluster_id" %in% colnames(obj@meta.data)) {
    "cluster_id"
  } else {
    normalize_scalar_value(layer_row$cluster_column[[1]], "seurat_clusters")
  }
  if (!cluster_var %in% colnames(obj@meta.data)) {
    cluster_var <- "seurat_clusters"
  }
  obj$cluster_id <- as.character(obj@meta.data[[cluster_var]])

  if (!"annotation_label" %in% colnames(obj@meta.data)) {
    obj$annotation_label <- if ("cell_type" %in% colnames(obj@meta.data)) as.character(obj$cell_type) else obj$cluster_id
  }

  for (cmp_idx in seq_len(nrow(comparison_df))) {
    comparison_row <- comparison_df[cmp_idx, , drop = FALSE]
    if (!comparison_applies_to_layer(comparison_row, layer_id)) {
      next
    }

    group_var <- normalize_scalar_value(comparison_row$group_var[[1]], "group_id")
    if (!group_var %in% colnames(obj@meta.data)) {
      group_var <- "group_id"
    }

    obj_sub <- subset(obj, cells = colnames(obj)[obj@meta.data[[group_var]] %in% c(comparison_row$ident_1[[1]], comparison_row$ident_2[[1]])])
    obj_sub <- maybe_join_layers(obj_sub)

    out_paths <- marker_discovery_paths(cfg, layer_id, comparison_row$comparison_id[[1]])
    ensure_parent_dirs(out_paths)

    cluster_ids <- sort(unique(as.character(obj_sub$cluster_id)))
    result_rows <- list()
    summary_rows <- list()

    for (cluster_id in cluster_ids) {
      res <- run_exploratory_findmarkers(obj_sub, cluster_id, comparison_row, group_var)
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        comparison_id = comparison_row$comparison_id[[1]],
        layer_id = layer_id,
        cluster_id = cluster_id,
        annotation_label = unique(as.character(obj_sub$annotation_label[obj_sub$cluster_id == cluster_id]))[1],
        result_available = ifelse(is.null(res), "no", "yes"),
        inference_scope = "exploratory_cell_level",
        stringsAsFactors = FALSE
      )
      if (!is.null(res)) {
        res$comparison_id <- comparison_row$comparison_id[[1]]
        res$layer_id <- layer_id
        res$cluster_id <- cluster_id
        res$annotation_label <- unique(as.character(obj_sub$annotation_label[obj_sub$cluster_id == cluster_id]))[1]
        res$inference_scope <- "exploratory_cell_level"
        result_rows[[length(result_rows) + 1]] <- res
      }
    }

    if (length(result_rows) > 0) {
      write_tsv(bind_rows(result_rows), out_paths$exploratory_tsv)
    } else {
      write_tsv(data.frame(stringsAsFactors = FALSE), out_paths$exploratory_tsv)
    }
    write_tsv(bind_rows(summary_rows), out_paths$exploratory_summary_tsv)

    manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
      comparison_id = comparison_row$comparison_id[[1]],
      layer_id = layer_id,
      group_var = group_var,
      exploratory_results_tsv = out_paths$exploratory_tsv,
      exploratory_summary_tsv = out_paths$exploratory_summary_tsv,
      stringsAsFactors = FALSE
    )
  }
}

write_tsv(bind_rows(manifest_rows), manifest_path)
write_markdown(
  c(
    "# Marker Discovery Report",
    "",
    "- inference scope: `exploratory_cell_level`",
    "- `FindMarkers/FindAllMarkers` 仅用于 marker discovery、cluster characterization 和 fallback。",
    sprintf("- manifest: `%s`", manifest_path)
  ),
  file.path(stage_dir, "report.md")
)
message("exploratory marker discovery 完成，manifest: ", manifest_path)
