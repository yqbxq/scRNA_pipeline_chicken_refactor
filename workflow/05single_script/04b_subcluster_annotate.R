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
source_utf8(file.path(.script_dir, "helpers", "project_paths_04.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "triage_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "clustering_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "annotation_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "cell_subtype_backfill_utils.R"))

load_required_packages(c("Seurat", "dplyr", "ggplot2", "tibble", "tidyr", "jsonlite"))

cfg <- get_single_script_config_04()
module_name <- "04b_subcluster_annotate"
prepare_dirs_04(cfg)
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

empty_annotation_summary_04 <- function() {
  data.frame(
    layer_id = character(0),
    cluster_count = integer(0),
    determined_n = integer(0),
    tentative_n = integer(0),
    undetermined_n = integer(0),
    determined_fraction = numeric(0),
    annotated_rds = character(0),
    annotation_table_tsv = character(0),
    annotation_evidence_tsv = character(0),
    module_score_summary_tsv = character(0),
    stringsAsFactors = FALSE
  )
}

manifest_clustered_keys_04 <- function(manifest) {
  output_names <- names(manifest$outputs %||% list())
  if (is.null(output_names)) {
    return(character(0))
  }
  output_names[startsWith(output_names, "clustered_")]
}

layer_id_from_clustered_key_04 <- function(key) {
  sub("^clustered_", "", key)
}

manifest_04a <- read_manifest_local(cfg$module_04a_manifest_path)
clustered_keys <- manifest_clustered_keys_04(manifest_04a)

summary_rows <- list()
triage_rows <- list()
output_entries <- list()

for (clustered_key in clustered_keys) {
  layer_id <- layer_id_from_clustered_key_04(clustered_key)
  clustered_rds <- resolve_output_local(manifest_04a, clustered_key)
  if (!file.exists(clustered_rds)) {
    warning(sprintf("跳过缺失 clustered RDS: %s -> %s", clustered_key, clustered_rds), call. = FALSE)
    next
  }

  message("04b annotate layer: ", layer_id)
  seu <- readRDS(clustered_rds)
  seu <- maybe_join_layers(seu)
  cluster_col <- paste0(layer_id, "_cluster")
  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop(sprintf("clustered object 缺少 cluster 列: %s", cluster_col), call. = FALSE)
  }

  tissue_values <- if ("tissue" %in% colnames(seu@meta.data)) unique(normalize_flag(seu$tissue, "*")) else character(0)
  tissue_values <- tissue_values[nzchar(tissue_values) & tissue_values != "*"]
  panel_df <- read_marker_panel_rows(
    cfg,
    layer_id = layer_id,
    tissue = if (length(tissue_values) == 1) tissue_values[[1]] else NULL
  )

  table_dir <- layer_annotation_table_dir_04(cfg, layer_id)
  report_dir <- layer_annotation_report_dir_04(cfg, layer_id)
  checkpoint_dir <- layer_checkpoint_dir_04(cfg, layer_id)
  ensure_dir(table_dir)
  ensure_dir(report_dir)
  ensure_dir(checkpoint_dir)

  paths_module <- list(
    cluster_var = cluster_col,
    table_dir = table_dir,
    figure_dir = report_dir,
    cluster_markers_tsv = file.path(table_dir, sprintf("cluster_markers_%s.tsv", layer_id)),
    annotation_table_tsv = file.path(table_dir, sprintf("annotation_table_%s.tsv", layer_id)),
    annotation_evidence_tsv = file.path(table_dir, sprintf("annotation_evidence_%s.tsv", layer_id)),
    module_score_summary_tsv = file.path(table_dir, sprintf("module_score_summary_%s.tsv", layer_id)),
    top_markers_dotplot_png = file.path(report_dir, sprintf("top_markers_dotplot_%s.png", layer_id)),
    panel_validation_dotplot_png = file.path(report_dir, sprintf("panel_validation_dotplot_%s.png", layer_id)),
    module_score_heatmap_png = file.path(report_dir, sprintf("module_score_heatmap_%s.png", layer_id))
  )

  result <- annotate_one_layer(seu, layer_id, panel_df, cfg, paths_module)
  seu <- result$seu
  ensure_placeholder_png(paths_module$top_markers_dotplot_png, "Top markers unavailable")
  ensure_placeholder_png(paths_module$panel_validation_dotplot_png, "Panel validation unavailable")
  ensure_placeholder_png(paths_module$module_score_heatmap_png, "Module score unavailable")

  seu[[paste0(layer_id, "_cell_type")]] <- seu$annotation_label
  seu[[paste0(layer_id, "_cell_type_confidence")]] <- seu$annotation_confidence
  seu[[paste0(layer_id, "_annotation_relation")]] <- seu$annotation_evidence_relation
  seu$cell_subtype <- seu$annotation_label

  annotated_rds <- file.path(checkpoint_dir, sprintf("%s_after_annotation.rds", layer_id))
  saveRDS(seu, annotated_rds)

  annotation_table <- result$annotation_table
  evidence_table <- result$evidence_table
  determined_fraction <- if (nrow(annotation_table) > 0) mean(annotation_table$confidence == "确定") else 0
  determined_n <- if (nrow(annotation_table) > 0) sum(annotation_table$confidence == "确定", na.rm = TRUE) else 0L
  tentative_n <- if (nrow(annotation_table) > 0) sum(annotation_table$confidence == "暂定", na.rm = TRUE) else 0L
  undetermined_n <- if (nrow(annotation_table) > 0) sum(annotation_table$confidence == "未定", na.rm = TRUE) else 0L

  if (determined_fraction < 0.5) {
    triage_rows[[length(triage_rows) + 1]] <- make_triage_row(
      sample_id = layer_id,
      severity = "medium",
      signal_id = "subcluster_low_confidence",
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
        sample_id = layer_id,
        severity = "medium",
        signal_id = "annotation_panel_conflict",
        evidence = sprintf("conflict_clusters=%s", paste(conflict_clusters$cluster_id, collapse = ","))
      )
    }
  }

  status_df <- read_tsv_optional(cfg$layer_status_file)
  status_row <- if (nrow(status_df) > 0 && any(status_df$layer_id == layer_id)) {
    status_df[status_df$layer_id == layer_id, , drop = FALSE][1, , drop = FALSE]
  } else {
    data.frame(layer_id = layer_id, stringsAsFactors = FALSE)
  }
  status_row$status <- "annotated"
  status_row$annotated_rds <- normalizePath(annotated_rds, winslash = "/", mustWork = FALSE)
  status_row$annotation_status <- "annotated"
  status_row$determined_fraction <- sprintf("%.6f", determined_fraction)
  status_row$cluster_column <- cluster_col
  invisible(upsert_layer_status(cfg$layer_status_file, status_row))

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    layer_id = layer_id,
    cluster_count = nrow(annotation_table),
    determined_n = determined_n,
    tentative_n = tentative_n,
    undetermined_n = undetermined_n,
    determined_fraction = determined_fraction,
    annotated_rds = normalizePath(annotated_rds, winslash = "/", mustWork = FALSE),
    annotation_table_tsv = normalizePath(paths_module$annotation_table_tsv, winslash = "/", mustWork = FALSE),
    annotation_evidence_tsv = normalizePath(paths_module$annotation_evidence_tsv, winslash = "/", mustWork = FALSE),
    module_score_summary_tsv = normalizePath(paths_module$module_score_summary_tsv, winslash = "/", mustWork = FALSE),
    stringsAsFactors = FALSE
  )

  output_entries[[paste0("annotated_", layer_id)]] <- build_output_entry(annotated_rds, "rds", module_name, sprintf("annotated subcluster layer %s", layer_id), base_dir = cfg$project_root)
  output_entries[[paste0("annotation_table_tsv_", layer_id)]] <- build_output_entry(paths_module$annotation_table_tsv, "tsv", module_name, sprintf("annotation table for %s", layer_id), base_dir = cfg$project_root, schema = infer_schema_from_df(annotation_table))
  output_entries[[paste0("annotation_evidence_tsv_", layer_id)]] <- build_output_entry(paths_module$annotation_evidence_tsv, "tsv", module_name, sprintf("annotation evidence for %s", layer_id), base_dir = cfg$project_root, schema = infer_schema_from_df(evidence_table))
  output_entries[[paste0("module_score_summary_tsv_", layer_id)]] <- build_output_entry(paths_module$module_score_summary_tsv, "tsv", module_name, sprintf("module score summary for %s", layer_id), base_dir = cfg$project_root, schema = infer_schema_from_df(result$module_score_summary))
  output_entries[[paste0("top_markers_dotplot_png_", layer_id)]] <- build_output_entry(paths_module$top_markers_dotplot_png, "png", module_name, sprintf("top marker dotplot for %s", layer_id), base_dir = cfg$project_root)
  output_entries[[paste0("panel_validation_dotplot_png_", layer_id)]] <- build_output_entry(paths_module$panel_validation_dotplot_png, "png", module_name, sprintf("panel validation dotplot for %s", layer_id), base_dir = cfg$project_root)
  output_entries[[paste0("module_score_heatmap_png_", layer_id)]] <- build_output_entry(paths_module$module_score_heatmap_png, "png", module_name, sprintf("module score heatmap for %s", layer_id), base_dir = cfg$project_root)
}

summary_df <- if (length(summary_rows) > 0) dplyr::bind_rows(summary_rows) else empty_annotation_summary_04()
summary_tsv <- file.path(cfg$subcluster_table_dir, "subcluster_annotation_summary.tsv")
write_tsv_local(summary_df, summary_tsv)

triage_df <- if (length(triage_rows) > 0) dplyr::bind_rows(triage_rows) else empty_triage_df(include_sample = TRUE)
triage_tsv <- file.path(cfg$subcluster_table_dir, "subcluster_annotation_triage.tsv")
write_tsv_local(triage_df, triage_tsv)

backfill_result <- backfill_panorama_cell_subtype_04(cfg, summary_df)

output_entries$summary_tsv <- build_output_entry(summary_tsv, "tsv", module_name, "one row per annotated subcluster layer", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df))
output_entries$subcluster_annotation_triage_tsv <- build_output_entry(triage_tsv, "tsv", module_name, "one row per subcluster annotation triage signal", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df))
output_entries$layer_status_tsv <- build_output_entry(cfg$layer_status_file, "tsv", module_name, "one row per built/annotated object layer", base_dir = cfg$project_root)
output_entries$panorama_cell_subtype_rds <- build_output_entry(backfill_result$panorama_rds, "rds", module_name, "panorama object with backfilled cell_subtype metadata", base_dir = cfg$project_root)
output_entries$panorama_cell_subtype_backfill_tsv <- build_output_entry(backfill_result$summary_tsv, "tsv", module_name, "cell_subtype backfill summary by source layer", base_dir = cfg$project_root, schema = infer_schema_from_df(backfill_result$summary))
if (nzchar(backfill_result$compat_rds)) {
  output_entries$compatibility_panorama_cell_subtype_rds <- build_output_entry(backfill_result$compat_rds, "rds", module_name, "compatibility panorama object with backfilled cell_subtype metadata", base_dir = cfg$project_root)
}

if (file.exists(cfg$module_04b_manifest_path)) {
  unlink(cfg$module_04b_manifest_path)
}
write_manifest_local(
  manifest_path = cfg$module_04b_manifest_path,
  new_outputs = output_entries,
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    module_04a_manifest = cfg$module_04a_manifest_path,
    module_03d_manifest = cfg$module_03d_manifest_path,
    marker_panel_dir = cfg$marker_panel_dir,
    layer_status_tsv = cfg$layer_status_file
  ),
  version = cfg$module_version,
  depends_on = list(module_03d = cfg$module_03d_manifest_path, module_04a = cfg$module_04a_manifest_path)
)

message("04b 完成。annotated layers: ", length(summary_rows), "; panorama cell_subtype backfill: ", backfill_result$summary_tsv)
