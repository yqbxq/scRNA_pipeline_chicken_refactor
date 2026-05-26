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

load_required_packages(c("Seurat", "ggplot2", "jsonlite"))

cfg <- get_spatial_script_config()
module_name <- "spatial_02a2_compute_umap"
prepare_dirs_spatial(cfg)

umap_seed <- as.integer(Sys.getenv("UMAP_SEED", as.character(cfg$random_seed)))
if (!is.finite(umap_seed)) {
  umap_seed <- cfg$random_seed
}
set.seed(umap_seed)

plot_reduction <- function(obj, reduction, color_col, path, title) {
  emb <- Seurat::Embeddings(obj, reduction = reduction)
  color <- if (color_col %in% colnames(obj@meta.data)) obj@meta.data[rownames(emb), color_col, drop = TRUE] else "all"
  plot_df <- data.frame(x = emb[, 1], y = emb[, 2], color = color, stringsAsFactors = FALSE)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, color = color)) +
    ggplot2::geom_point(size = 0.7, alpha = 0.85) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(x = paste0(reduction, "_1"), y = paste0(reduction, "_2"), color = color_col, title = title)
  ggplot2::ggsave(path, p, width = 5.8, height = 4.8, dpi = 180, bg = "white")
  path
}

panorama <- readRDS(cfg$spatial_panorama_integrated_rds)
reduction_specs <- data.frame(
  integration_mode = c("none", "harmony", "cca"),
  source_reduction = c("pca_none", "harmony", "cca_integrated"),
  umap_name = c("umap_none", "umap_harmony", "umap_cca"),
  stringsAsFactors = FALSE
)

rows <- list()
figures <- list()
for (i in seq_len(nrow(reduction_specs))) {
  spec <- reduction_specs[i, , drop = FALSE]
  source_reduction <- spec$source_reduction[[1]]
  umap_name <- spec$umap_name[[1]]
  if (!source_reduction %in% names(panorama@reductions)) {
    rows[[length(rows) + 1L]] <- data.frame(
      integration_mode = spec$integration_mode[[1]],
      source_reduction = source_reduction,
      umap_name = umap_name,
      status = "skipped_missing_reduction",
      dims = "",
      figure = "",
      message = sprintf("missing reduction: %s", source_reduction),
      stringsAsFactors = FALSE
    )
    next
  }

  emb <- Seurat::Embeddings(panorama, reduction = source_reduction)
  dims <- seq_len(min(ncol(emb), max(as.integer(spatial_parse_dims(Sys.getenv("PCA_DIMS_PANORAMA", "1:30"), fallback = 1:30)))))
  started <- proc.time()[["elapsed"]]
  panorama <- spatial_run_umap_or_fallback(panorama, source_reduction, umap_name, dims)
  figure_path <- file.path(cfg$spatial_umap_figure_dir, sprintf("%s_by_section.png", umap_name))
  figure_path <- plot_reduction(panorama, umap_name, "section_id", figure_path, sprintf("UMAP %s by section", spec$integration_mode[[1]]))
  figures[[umap_name]] <- figure_path
  rows[[length(rows) + 1L]] <- data.frame(
    integration_mode = spec$integration_mode[[1]],
    source_reduction = source_reduction,
    umap_name = umap_name,
    status = "ok",
    dims = paste(dims, collapse = ","),
    figure = normalizePath(figure_path, winslash = "/", mustWork = FALSE),
    message = "",
    runtime_sec = round(proc.time()[["elapsed"]] - started, 3),
    umap_n_neighbors = as.integer(Sys.getenv("UMAP_N_NEIGHBORS", "30")),
    umap_min_dist = as.numeric(Sys.getenv("UMAP_MIN_DIST", "0.3")),
    umap_spread = as.numeric(Sys.getenv("UMAP_SPREAD", "1.0")),
    umap_metric = Sys.getenv("UMAP_METRIC", "cosine"),
    umap_seed = umap_seed,
    stringsAsFactors = FALSE
  )
}

saveRDS(panorama, cfg$spatial_panorama_umap_rds)
umap_index <- do.call(rbind, rows)
umap_index_tsv <- file.path(cfg$spatial_umap_table_dir, "spatial_umap_index.tsv")
spatial_write_tsv(umap_index, umap_index_tsv)

report_path <- file.path(cfg$spatial_umap_table_dir, "report.md")
report_lines <- c(
  "# Spatial UMAP",
  "",
  sprintf("- Source object: `%s`", relative_path_local(cfg$spatial_panorama_integrated_rds, cfg$project_root)),
  sprintf("- Output object: `%s`", relative_path_local(cfg$spatial_panorama_umap_rds, cfg$project_root)),
  "",
  "## UMAP Summary",
  render_markdown_table_local(umap_index)
)
write_markdown_local(report_lines, report_path)

write_manifest_local(
  manifest_path = cfg$module_02a2_umap_manifest_path,
  new_outputs = list(
    panorama_umap = build_output_entry(cfg$spatial_panorama_umap_rds, "rds", module_name, "spatial panorama object with UMAP reductions", base_dir = cfg$project_root),
    umap_index = build_output_entry(umap_index_tsv, "tsv", module_name, "one row per spatial UMAP embedding", base_dir = cfg$project_root, schema = infer_schema_from_df(umap_index)),
    report = build_output_entry(report_path, "md", module_name, "spatial UMAP summary report", base_dir = cfg$project_root),
    figures_dir = build_output_entry(cfg$spatial_umap_figure_dir, "directory", module_name, "spatial UMAP figures", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(module_02a_integration_manifest = cfg$module_02a_integration_manifest_path),
  version = cfg$module_02a2_version,
  depends_on = list(spatial_02a_integration_eda = cfg$module_02a_integration_manifest_path)
)

message("spatial UMAP complete: ", umap_index_tsv)
