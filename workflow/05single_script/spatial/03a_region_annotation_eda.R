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
module_name <- "spatial_03a_region_annotation_eda"
prepare_dirs_spatial(cfg)

plot_spatial_region <- function(obj, path) {
  coords <- spatial_filter_coords(obj)
  coords$region <- as.character(obj@meta.data[rownames(coords), "region", drop = TRUE])
  p <- ggplot2::ggplot(coords, ggplot2::aes(x = col, y = row, color = region)) +
    ggplot2::geom_point(size = 0.9, alpha = 0.9, na.rm = TRUE) +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(color = "region")
  ggplot2::ggsave(path, p, width = 5.8, height = 5.2, dpi = 180, bg = "white")
  path
}

plot_region_bar <- function(summary_df, path) {
  p <- ggplot2::ggplot(summary_df, ggplot2::aes(x = region, y = n_spots, fill = region)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none") +
    ggplot2::labs(x = NULL, y = "spots")
  ggplot2::ggsave(path, p, width = 5.8, height = 4.2, dpi = 180, bg = "white")
  path
}

plot_cluster_region <- function(region_assignment, path) {
  p <- ggplot2::ggplot(region_assignment, ggplot2::aes(x = cluster, y = region, fill = confidence)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(x = "cluster", y = "region", fill = "confidence")
  ggplot2::ggsave(path, p, width = 6.2, height = 4.2, dpi = 180, bg = "white")
  path
}

plot_module_heatmap <- function(module_score_matrix, path) {
  value_cols <- setdiff(colnames(module_score_matrix), "cluster")
  if (length(value_cols) == 0) {
    plot_df <- data.frame(cluster = "none", region = "none", score = 0)
  } else {
    plot_df <- do.call(rbind, lapply(value_cols, function(col) {
      data.frame(cluster = module_score_matrix$cluster, region = col, score = as.numeric(module_score_matrix[[col]]), stringsAsFactors = FALSE)
    }))
  }
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = region, y = cluster, fill = score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
    ggplot2::labs(x = "region module", y = "cluster", fill = "score")
  ggplot2::ggsave(path, p, width = 6.2, height = 4.6, dpi = 180, bg = "white")
  path
}

if (!file.exists(cfg$spatial_panorama_annotated_rds)) {
  stop(sprintf("missing annotated panorama: %s", cfg$spatial_panorama_annotated_rds), call. = FALSE)
}
panorama <- readRDS(cfg$spatial_panorama_annotated_rds)
if (!"region" %in% colnames(panorama@meta.data)) {
  stop("annotated panorama is missing region metadata", call. = FALSE)
}

region_assignment_path <- file.path(cfg$spatial_region_annotation_table_dir, "region_assignment.tsv")
cluster_markers_path <- file.path(cfg$spatial_region_annotation_table_dir, "cluster_markers.tsv")
module_score_matrix_path <- file.path(cfg$spatial_region_annotation_table_dir, "module_score_matrix.tsv")
region_assignment <- spatial_read_tsv(region_assignment_path)
cluster_markers <- spatial_read_tsv(cluster_markers_path)
module_score_matrix <- spatial_read_tsv(module_score_matrix_path)

region_counts <- as.data.frame(table(region = panorama$region), stringsAsFactors = FALSE)
colnames(region_counts) <- c("region", "n_spots")
region_counts$fraction <- region_counts$n_spots / sum(region_counts$n_spots)
section_col <- if ("section_id" %in% colnames(panorama@meta.data)) "section_id" else if ("spatial_section_id" %in% colnames(panorama@meta.data)) "spatial_section_id" else ""
if (nzchar(section_col)) {
  section_tab <- as.data.frame.matrix(table(panorama$region, panorama@meta.data[[section_col]]))
  section_tab$region <- rownames(section_tab)
  region_summary <- merge(region_counts, section_tab, by = "region", all.x = TRUE)
} else {
  region_summary <- region_counts
}
triage_df <- compute_region_triage(panorama, region_assignment)

region_summary_path <- file.path(cfg$spatial_region_annotation_dir, "region_summary.tsv")
region_triage_path <- file.path(cfg$spatial_region_annotation_dir, "region_triage.tsv")
spatial_write_tsv(region_summary, region_summary_path)
spatial_write_tsv(triage_df, region_triage_path)

figures <- c(
  spatial_region_assignment = plot_spatial_region(panorama, file.path(cfg$spatial_region_annotation_figure_dir, "spatial_region_assignment.png")),
  region_marker_heatmap = plot_module_heatmap(module_score_matrix, file.path(cfg$spatial_region_annotation_figure_dir, "region_marker_heatmap.png")),
  cluster_region_sankey = plot_cluster_region(region_assignment, file.path(cfg$spatial_region_annotation_figure_dir, "cluster_region_sankey.png")),
  region_spot_barchart = plot_region_bar(region_summary, file.path(cfg$spatial_region_annotation_figure_dir, "region_spot_barchart.png"))
)

report_lines <- c(
  "# Spatial Region Annotation EDA",
  "",
  sprintf("- Annotated panorama: `%s`", relative_path_local(cfg$spatial_panorama_annotated_rds, cfg$project_root)),
  sprintf("- Region summary: `%s`", relative_path_local(region_summary_path, cfg$project_root)),
  sprintf("- Triage table: `%s`", relative_path_local(region_triage_path, cfg$project_root)),
  "",
  "## Triage Signals",
  render_markdown_table_local(triage_df),
  "",
  "## Region Summary",
  render_markdown_table_local(region_summary),
  "",
  "## Cluster Region Mapping",
  render_markdown_table_local(region_assignment),
  "",
  "## Marker Evidence Preview",
  render_markdown_table_local(utils::head(cluster_markers, 50))
)
write_markdown_local(report_lines, file.path(cfg$spatial_region_annotation_dir, "report.md"))

write_manifest_local(
  manifest_path = cfg$module_03a_region_annotation_eda_manifest_path,
  new_outputs = list(
    region_summary = build_output_entry(region_summary_path, "tsv", module_name, "spot counts and section balance by region", base_dir = cfg$project_root, schema = infer_schema_from_df(region_summary)),
    region_triage = build_output_entry(region_triage_path, "tsv", module_name, "region annotation triage signals", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df)),
    report = build_output_entry(file.path(cfg$spatial_region_annotation_dir, "report.md"), "md", module_name, "spatial region annotation EDA report", base_dir = cfg$project_root),
    figures_dir = build_output_entry(cfg$spatial_region_annotation_figure_dir, "directory", module_name, "region annotation EDA figures", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(panorama_annotated = cfg$spatial_panorama_annotated_rds, region_assignment = region_assignment_path, cluster_markers = cluster_markers_path),
  version = cfg$module_03_version,
  depends_on = list(spatial_03_region_annotation = cfg$module_03_region_annotation_manifest_path)
)

message(sprintf("spatial region annotation EDA complete: %d regions", nrow(region_summary)))
