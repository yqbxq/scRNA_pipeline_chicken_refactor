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
module_name <- "spatial_03_region_annotation"
prepare_dirs_spatial(cfg)
set.seed(cfg$random_seed)

parse_cli <- function(args) {
  out <- list(
    marker_panel_dir = file.path(cfg$config_dir, "marker_panels"),
    override_file = cfg$selected_region_annotation_override_file,
    cluster_col = "cluster_default"
  )
  positional <- character()
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--marker-panel-dir") && i < length(args)) {
      out$marker_panel_dir <- args[[i + 1L]]
      i <- i + 2L
    } else if (identical(arg, "--override-file") && i < length(args)) {
      out$override_file <- args[[i + 1L]]
      i <- i + 2L
    } else if (identical(arg, "--cluster-col") && i < length(args)) {
      out$cluster_col <- args[[i + 1L]]
      i <- i + 2L
    } else {
      positional <- c(positional, arg)
      i <- i + 1L
    }
  }
  if (length(positional) > 0 && dir.exists(positional[[1]])) {
    out$marker_panel_dir <- positional[[1]]
  }
  out
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
if (!file.exists(cfg$spatial_panorama_clustered_rds)) {
  stop(sprintf("missing clustered panorama: %s", cfg$spatial_panorama_clustered_rds), call. = FALSE)
}
panorama <- readRDS(cfg$spatial_panorama_clustered_rds)
if (!args$cluster_col %in% colnames(panorama@meta.data)) {
  stop(sprintf("cluster column is missing in clustered panorama: %s", args$cluster_col), call. = FALSE)
}
panel <- read_spatial_region_panel(args$marker_panel_dir, layer_id = "panorama_st")

annotation <- region_three_evidence_chain(
  panorama = panorama,
  cluster_col = args$cluster_col,
  panel = panel,
  override_file = args$override_file
)
panorama <- assign_region_labels(annotation$panorama, annotation$evidence_df)
panorama@misc$region_annotation$panel_dir <- args$marker_panel_dir
panorama@misc$region_annotation$override_file <- args$override_file
saveRDS(panorama, cfg$spatial_panorama_annotated_rds)

region_assignment <- annotation$evidence_df
cluster_markers <- annotation$marker_table
marker_cols <- c("cluster", "gene", "avg_log2FC", "p_val_adj", "pct.1", "pct.2")
for (col in marker_cols) {
  if (!col %in% colnames(cluster_markers)) {
    cluster_markers[[col]] <- NA
  }
}
cluster_markers <- cluster_markers[, marker_cols, drop = FALSE]
module_score_matrix <- annotation$module_score_matrix
module_score_matrix$cluster <- rownames(module_score_matrix)
module_score_matrix <- module_score_matrix[, c("cluster", setdiff(colnames(module_score_matrix), "cluster")), drop = FALSE]

region_assignment_path <- file.path(cfg$spatial_region_annotation_table_dir, "region_assignment.tsv")
cluster_markers_path <- file.path(cfg$spatial_region_annotation_table_dir, "cluster_markers.tsv")
module_score_matrix_path <- file.path(cfg$spatial_region_annotation_table_dir, "module_score_matrix.tsv")
spatial_write_tsv(region_assignment, region_assignment_path)
spatial_write_tsv(cluster_markers, cluster_markers_path)
spatial_write_tsv(module_score_matrix, module_score_matrix_path)

write_manifest_local(
  manifest_path = cfg$module_03_region_annotation_manifest_path,
  new_outputs = list(
    panorama_annotated = build_output_entry(cfg$spatial_panorama_annotated_rds, "rds", module_name, "spatial panorama with cluster-level region labels", base_dir = cfg$project_root),
    region_assignment = build_output_entry(region_assignment_path, "tsv", module_name, "cluster to region assignment and evidence", base_dir = cfg$project_root, schema = infer_schema_from_df(region_assignment)),
    cluster_markers = build_output_entry(cluster_markers_path, "tsv", module_name, "FindAllMarkers evidence used by region annotation", base_dir = cfg$project_root, schema = infer_schema_from_df(cluster_markers)),
    module_score_matrix = build_output_entry(module_score_matrix_path, "tsv", module_name, "cluster by region module score matrix", base_dir = cfg$project_root, schema = infer_schema_from_df(module_score_matrix))
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(panorama_clustered = cfg$spatial_panorama_clustered_rds, marker_panel_dir = args$marker_panel_dir, region_override_file = args$override_file),
  version = cfg$module_03_version,
  depends_on = list(spatial_02b_finalize_clustering = cfg$module_02b_clustering_manifest_path)
)

message(sprintf("spatial region annotation complete: %d clusters", nrow(region_assignment)))
