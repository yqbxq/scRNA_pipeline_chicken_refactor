#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(.script_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_common.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "project_paths_spatial.R"), encoding = "UTF-8")

load_required_packages(c("Seurat", "Matrix", "ggplot2", "jsonlite"))

cfg <- get_spatial_script_config()
module_name <- "spatial_04a_subcluster_build"
prepare_dirs_spatial(cfg)
set.seed(cfg$random_seed)

parse_cli <- function(args) {
  out <- list(layer_file = cfg$spatial_object_layer_file)
  positional <- character()
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--layer-file") && i < length(args)) {
      out$layer_file <- args[[i + 1L]]
      i <- i + 2L
    } else {
      positional <- c(positional, arg)
      i <- i + 1L
    }
  }
  if (length(positional) > 0 && file.exists(positional[[1]])) {
    out$layer_file <- positional[[1]]
  }
  out
}

index_schema <- function() {
  data.frame(
    layer_id = character(),
    parent_layer = character(),
    layer_role = character(),
    enabled = character(),
    selection_column = character(),
    selection_values = character(),
    n_spots = integer(),
    normalization_method = character(),
    hvg_nfeatures = integer(),
    pca_dims = character(),
    target_clusters = integer(),
    picked_resolution = character(),
    n_clusters_picked = integer(),
    status = character(),
    reason = character(),
    out_rds = character(),
    runtime_sec = numeric(),
    stringsAsFactors = FALSE
  )
}

index_row <- function(layer, n_spots = NA_integer_, picked_resolution = "", n_clusters_picked = NA_integer_, status, reason = "", out_rds = "", runtime_sec = NA_real_) {
  data.frame(
    layer_id = spatial_cell(layer, "layer_id", ""),
    parent_layer = spatial_cell(layer, "parent_layer", ""),
    layer_role = spatial_cell(layer, "layer_role", "region_subset"),
    enabled = spatial_cell(layer, "enabled", ""),
    selection_column = spatial_cell(layer, "selection_column", ""),
    selection_values = spatial_cell(layer, "selection_values", ""),
    n_spots = as.integer(n_spots),
    normalization_method = spatial_layer_normalization_method(as.list(layer[1, , drop = FALSE])),
    hvg_nfeatures = as.integer(spatial_numeric_or(layer$hvg_nfeatures, NA_real_)),
    pca_dims = spatial_cell(layer, "pca_dims", ""),
    target_clusters = as.integer(spatial_numeric_or(layer$target_clusters, NA_real_)),
    picked_resolution = picked_resolution,
    n_clusters_picked = as.integer(n_clusters_picked),
    status = status,
    reason = reason,
    out_rds = out_rds,
    runtime_sec = as.numeric(runtime_sec),
    stringsAsFactors = FALSE
  )
}

write_index_manifest <- function(index_df, inputs) {
  index_path <- file.path(cfg$spatial_subcluster_build_table_dir, "region_subset_index.tsv")
  spatial_write_tsv(index_df, index_path)
  write_manifest_local(
    manifest_path = cfg$module_04a_subcluster_build_manifest_path,
    new_outputs = list(
      region_subset_index = build_output_entry(index_path, "tsv", module_name, "region subset build status and output index", base_dir = cfg$project_root, schema = infer_schema_from_df(index_df)),
      region_subset_objects_dir = build_output_entry(cfg$spatial_region_dir, "directory", module_name, "per-layer subclustered region subset RDS files", base_dir = cfg$project_root),
      tables_dir = build_output_entry(cfg$spatial_subcluster_build_table_dir, "directory", module_name, "per-layer subcluster metrics", base_dir = cfg$project_root),
      figures_dir = build_output_entry(cfg$spatial_subcluster_build_figure_dir, "directory", module_name, "per-layer subcluster figures", base_dir = cfg$project_root)
    ),
    module_name = module_name,
    base_dir = cfg$project_root,
    inputs = inputs,
    version = cfg$module_04_version,
    depends_on = list(spatial_03a_region_annotation_eda = cfg$module_03a_region_annotation_eda_manifest_path)
  )
}

plot_umap_cluster <- function(obj, path) {
  ensure_dir(dirname(path))
  emb <- Seurat::Embeddings(obj, reduction = "umap_subcluster")
  plot_df <- data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2], cluster = as.character(obj@meta.data[rownames(emb), "cluster_default", drop = TRUE]))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
    ggplot2::geom_point(size = 0.9, alpha = 0.9, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(color = "cluster")
  ggplot2::ggsave(path, p, width = 5.8, height = 4.8, dpi = 180, bg = "white")
  path
}

plot_spatial_cluster <- function(obj, path) {
  ensure_dir(dirname(path))
  coords <- spatial_filter_coords(obj)
  coords <- coords[intersect(rownames(coords), rownames(obj@meta.data)), , drop = FALSE]
  coords$cluster <- as.character(obj@meta.data[rownames(coords), "cluster_default", drop = TRUE])
  p <- ggplot2::ggplot(coords, ggplot2::aes(x = col, y = row, color = cluster)) +
    ggplot2::geom_point(size = 0.9, alpha = 0.9, na.rm = TRUE) +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(color = "cluster")
  ggplot2::ggsave(path, p, width = 5.8, height = 5.2, dpi = 180, bg = "white")
  path
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
if (!file.exists(cfg$spatial_panorama_annotated_rds)) {
  stop(sprintf("missing annotated panorama: %s", cfg$spatial_panorama_annotated_rds), call. = FALSE)
}
panorama <- readRDS(cfg$spatial_panorama_annotated_rds)
layers <- spatial_region_subset_layers(args$layer_file)
enabled_flags <- if (nrow(layers) > 0 && "enabled" %in% colnames(layers)) vapply(layers$enabled, spatial_bool, logical(1), default = FALSE) else logical(0)
enabled <- if (nrow(layers) > 0) layers[enabled_flags, , drop = FALSE] else layers

if (nrow(enabled) == 0) {
  empty_layer <- data.frame(layer_id = "", parent_layer = "", layer_role = "region_subset", enabled = "no", selection_column = "", selection_values = "", hvg_nfeatures = NA, pca_dims = "", target_clusters = NA, normalization_methods = "", stringsAsFactors = FALSE)
  index_df <- rbind(index_schema(), index_row(empty_layer, status = "skipped_no_enabled_layers", reason = "no enabled region_subset layers"))
  write_index_manifest(index_df, inputs = list(panorama_annotated = cfg$spatial_panorama_annotated_rds, spatial_object_layer_file = args$layer_file))
  message("spatial subcluster build skipped: no enabled region_subset layers")
  quit(status = 0)
}

index_rows <- list()
for (i in seq_len(nrow(enabled))) {
  layer <- enabled[i, , drop = FALSE]
  started <- proc.time()[["elapsed"]]
  layer_id <- spatial_cell(layer, "layer_id", sprintf("layer_%d", i))
  settings <- as.list(layer[1, , drop = FALSE])
  selection_column <- spatial_cell(layer, "selection_column", "")
  selection_values <- spatial_cell(layer, "selection_values", "")
  selected_values <- spatial_tokenize(selection_values)
  selected_cells <- if (nzchar(selection_column) && selection_column %in% colnames(panorama@meta.data) && length(selected_values) > 0) {
    rownames(panorama@meta.data)[as.character(panorama@meta.data[[selection_column]]) %in% selected_values]
  } else {
    character(0)
  }
  n_spots <- length(selected_cells)
  if (n_spots < cfg$subcluster_min_spots) {
    index_rows[[length(index_rows) + 1L]] <- index_row(
      layer,
      n_spots = n_spots,
      status = "skipped_too_few_spots",
      reason = sprintf("selected %d spots, below SPATIAL_SUBCLUSTER_MIN_SPOTS=%d", n_spots, cfg$subcluster_min_spots),
      runtime_sec = round(proc.time()[["elapsed"]] - started, 3)
    )
    next
  }
  result <- tryCatch({
    ignored_backends <- spatial_cell(layer, "clustering_backends", "")
    if (nzchar(ignored_backends) && !identical(ignored_backends, "b1_seurat_snn")) {
      warning(sprintf("ignoring clustering_backends=%s; forcing b1_seurat_snn for region subset %s", ignored_backends, layer_id), call. = FALSE)
    }
    subset_obj <- spatial_subset_panorama_by_region(panorama, selection_column, selection_values)
    subset_obj <- rebuild_layer_normalization_st(subset_obj, settings)
    clustered <- cluster_subset_snn_st(subset_obj, settings)
    subset_obj <- merge_small_subclusters(clustered$obj, "cluster_default", cfg$subcluster_small_cluster_frac)
    metrics_df <- clustered$metrics_df
    metrics_df$layer_id <- layer_id
    metrics_df <- metrics_df[, c("layer_id", setdiff(colnames(metrics_df), "layer_id")), drop = FALSE]
    layer_table_dir <- file.path(cfg$spatial_subcluster_build_table_dir, layer_id)
    layer_figure_dir <- file.path(cfg$spatial_subcluster_build_figure_dir, layer_id)
    ensure_dir(layer_table_dir)
    ensure_dir(layer_figure_dir)
    metrics_path <- file.path(layer_table_dir, "cluster_metrics.tsv")
    spatial_write_tsv(metrics_df, metrics_path)
    plot_umap_cluster(subset_obj, file.path(layer_figure_dir, "umap_cluster.png"))
    plot_spatial_cluster(subset_obj, file.path(layer_figure_dir, "spatial_cluster.png"))
    out_rds <- file.path(cfg$spatial_region_dir, sprintf("%s.rds", spatial_safe_id(layer_id)))
    saveRDS(subset_obj, out_rds)
    n_clusters <- length(unique(as.character(subset_obj@meta.data$cluster_default)))
    index_row(
      layer,
      n_spots = ncol(subset_obj),
      picked_resolution = clustered$picked_resolution,
      n_clusters_picked = n_clusters,
      status = "ok",
      reason = "",
      out_rds = relative_path_local(out_rds, cfg$project_root),
      runtime_sec = round(proc.time()[["elapsed"]] - started, 3)
    )
  }, error = function(e) {
    index_row(
      layer,
      n_spots = n_spots,
      status = "failed",
      reason = conditionMessage(e),
      runtime_sec = round(proc.time()[["elapsed"]] - started, 3)
    )
  })
  index_rows[[length(index_rows) + 1L]] <- result
}

index_df <- if (length(index_rows) > 0) do.call(rbind, index_rows) else index_schema()
write_index_manifest(index_df, inputs = list(panorama_annotated = cfg$spatial_panorama_annotated_rds, spatial_object_layer_file = args$layer_file))
message(sprintf("spatial subcluster build complete: %d ok layer(s)", sum(index_df$status == "ok", na.rm = TRUE)))
