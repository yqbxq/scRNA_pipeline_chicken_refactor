#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
      getwd()
    }
  }
)

source_utf8 <- function(path) {
  source(path, encoding = "UTF-8")
}

source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "config.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_02.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_03.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "triage_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "clustering_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "annotation_utils.R"))

load_required_packages(c("Seurat", "dplyr", "ggplot2", "tibble", "tidyr", "jsonlite"))

cfg <- get_single_script_config_03()
module_name <- "03d_annotate"
prepare_dirs_03(cfg)
set.seed(cfg$random_seed)

ensure_placeholder_png <- function(path, title) {
  if (file.exists(path)) {
    return(invisible(TRUE))
  }
  plot_obj <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::labs(title = title)
  save_plot_local(plot_obj, path, width = 6, height = 4)
}

panorama_spec <- panorama_layer_spec(cfg)
manifest_03c <- read_manifest_local(cfg$module_03c_manifest_path)
clustered_rds <- resolve_output_local(manifest_03c, "clustered_object")
seu <- readRDS(clustered_rds)
seu <- maybe_join_layers(seu)

cluster_col <- paste0(panorama_spec$layer_id, "_cluster")
if (!cluster_col %in% colnames(seu@meta.data)) {
  stop(sprintf("clustered object 缺少 cluster 列: %s", cluster_col), call. = FALSE)
}

tissue_values <- if ("tissue" %in% colnames(seu@meta.data)) unique(normalize_flag(seu$tissue, "*")) else character(0)
tissue_values <- tissue_values[nzchar(tissue_values) & tissue_values != "*"]
panel_df <- read_marker_panel_rows(
  cfg,
  layer_id = panorama_spec$layer_id,
  tissue = if (length(tissue_values) == 1) tissue_values[[1]] else NULL
)

paths_module <- list(
  cluster_var = cluster_col,
  table_dir = cfg$annotation_table_dir_layer,
  figure_dir = cfg$annotation_report_dir_layer,
  cluster_markers_tsv = file.path(cfg$annotation_table_dir_layer, "cluster_markers.tsv"),
  annotation_table_tsv = file.path(cfg$annotation_table_dir_layer, "annotation_table.tsv"),
  annotation_evidence_tsv = file.path(cfg$annotation_table_dir_layer, "annotation_evidence.tsv"),
  module_score_summary_tsv = file.path(cfg$annotation_table_dir_layer, "module_score_summary.tsv"),
  top_markers_dotplot_png = file.path(cfg$annotation_report_dir_layer, "top_markers_dotplot.png"),
  panel_validation_dotplot_png = file.path(cfg$annotation_report_dir_layer, "panel_validation_dotplot.png"),
  module_score_heatmap_png = file.path(cfg$annotation_report_dir_layer, "module_score_heatmap.png")
)

result <- annotate_one_layer(seu, panorama_spec$layer_id, panel_df, cfg, paths_module)
seu <- result$seu
ensure_placeholder_png(paths_module$top_markers_dotplot_png, "Top markers unavailable")
ensure_placeholder_png(paths_module$panel_validation_dotplot_png, "Panel validation unavailable")
ensure_placeholder_png(paths_module$module_score_heatmap_png, "Module score unavailable")

ensure_dir(dirname(cfg$panorama_annotated_rds))
saveRDS(seu, cfg$panorama_annotated_rds)
ensure_dir(dirname(cfg$compat_annotated_rds))
if (!identical(normalizePath(cfg$panorama_annotated_rds, winslash = "/", mustWork = FALSE), normalizePath(cfg$compat_annotated_rds, winslash = "/", mustWork = FALSE))) {
  invisible(file.copy(cfg$panorama_annotated_rds, cfg$compat_annotated_rds, overwrite = TRUE))
}

status_df <- read_tsv_optional(cfg$layer_status_file)
status_row <- if (nrow(status_df) > 0 && any(status_df$layer_id == panorama_spec$layer_id)) {
  status_df[status_df$layer_id == panorama_spec$layer_id, , drop = FALSE][1, , drop = FALSE]
} else {
  data.frame(layer_id = panorama_spec$layer_id, stringsAsFactors = FALSE)
}
status_row$status <- "annotated"
status_row$annotated_rds <- normalizePath(cfg$panorama_annotated_rds, winslash = "/", mustWork = FALSE)
status_row$cluster_column <- cluster_col
invisible(upsert_layer_status(cfg$layer_status_file, status_row))

annotation_table <- result$annotation_table
evidence_table <- result$evidence_table
triage_rows <- list()
determined_fraction <- if (nrow(annotation_table) > 0) mean(annotation_table$confidence == "确定") else 0
if (determined_fraction < 0.5) {
  triage_rows[[length(triage_rows) + 1]] <- make_triage_row(
    sample_id = panorama_spec$layer_id,
    severity = "medium",
    signal_id = "annotation_low_confidence",
    evidence = sprintf("determined_fraction=%.3f; clusters=%s", determined_fraction, nrow(annotation_table))
  )
}
if (nrow(evidence_table) > 0) {
  conflict_clusters <- evidence_table %>%
    dplyr::filter(overlap_n > 0) %>%
    dplyr::count(cluster_id, name = "hit_panels") %>%
    dplyr::filter(hit_panels >= 2)
  if (nrow(conflict_clusters) > 0) {
    triage_rows[[length(triage_rows) + 1]] <- make_triage_row(
      sample_id = panorama_spec$layer_id,
      severity = "medium",
      signal_id = "annotation_panel_conflict",
      evidence = sprintf("conflict_clusters=%s", paste(conflict_clusters$cluster_id, collapse = ","))
    )
  }
}
triage_df <- if (length(triage_rows) > 0) dplyr::bind_rows(triage_rows) else empty_triage_df(include_sample = TRUE)
triage_tsv <- file.path(cfg$annotation_table_dir_layer, "annotation_triage.tsv")
write_tsv_local(triage_df, triage_tsv)

write_manifest_local(
  manifest_path = cfg$module_03d_manifest_path,
  new_outputs = list(
    annotated_object = build_output_entry(cfg$panorama_annotated_rds, "rds", module_name, "panorama annotated Seurat object", base_dir = cfg$project_root),
    compatibility_annotated_object = build_output_entry(cfg$compat_annotated_rds, "rds", module_name, "compatibility checkpoint for downstream legacy modules", base_dir = cfg$project_root),
    annotation_table_tsv = build_output_entry(paths_module$annotation_table_tsv, "tsv", module_name, "one row per cluster annotation decision", base_dir = cfg$project_root, schema = infer_schema_from_df(annotation_table)),
    annotation_evidence_tsv = build_output_entry(paths_module$annotation_evidence_tsv, "tsv", module_name, "one row per cluster/panel overlap", base_dir = cfg$project_root, schema = infer_schema_from_df(evidence_table)),
    module_score_summary_tsv = build_output_entry(paths_module$module_score_summary_tsv, "tsv", module_name, "module score validation summary by cluster", base_dir = cfg$project_root, schema = infer_schema_from_df(result$module_score_summary)),
    top_markers_dotplot_png = build_output_entry(paths_module$top_markers_dotplot_png, "png", module_name, "top marker dotplot", base_dir = cfg$project_root),
    panel_validation_dotplot_png = build_output_entry(paths_module$panel_validation_dotplot_png, "png", module_name, "panel validation dotplot", base_dir = cfg$project_root),
    module_score_heatmap_png = build_output_entry(paths_module$module_score_heatmap_png, "png", module_name, "module score heatmap", base_dir = cfg$project_root),
    annotation_triage_tsv = build_output_entry(triage_tsv, "tsv", module_name, "annotation triage signals", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df)),
    layer_status_tsv = build_output_entry(cfg$layer_status_file, "tsv", module_name, "one row per built/annotated object layer", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    clustered_object = clustered_rds,
    marker_panel_dir = cfg$marker_panel_dir,
    layer_status_tsv = cfg$layer_status_file
  ),
  version = cfg$module_version,
  depends_on = list(module_03c = cfg$module_03c_manifest_path)
)

message("03d 完成。annotated_object: ", cfg$panorama_annotated_rds)
