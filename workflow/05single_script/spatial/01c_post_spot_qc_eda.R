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
module_name <- "spatial_01c_post_qc_eda"
prepare_dirs_spatial(cfg)

thresholds <- spatial_thresholds_table(cfg)
qc_summary <- spatial_read_tsv(cfg$qc_filter_summary_tsv)
pre_triage <- spatial_read_tsv(cfg$pre_qc_triage_tsv)
if (nrow(qc_summary) == 0) {
  stop(sprintf("QC filter summary is missing or empty: %s", cfg$qc_filter_summary_tsv), call. = FALSE)
}

for (col in c(
  "sample_id", "section_id", "object_stem", "raw_rds", "post_qc_rds", "n_before",
  "n_after", "retention", "n_drop_high_mito", "n_drop_low_feature",
  "spatial_aware_summary", "spatial_aware_max_cluster_frac", "mask_png"
)) {
  if (!col %in% colnames(qc_summary)) {
    qc_summary[[col]] <- ""
  }
}
if (nrow(pre_triage) > 0 && !"section_id" %in% colnames(pre_triage)) {
  pre_triage$section_id <- ""
}

save_metric_plot <- function(plot_obj, path, width = 7, height = 5) {
  ensure_dir(dirname(path))
  ggplot2::ggsave(path, plot_obj, width = width, height = height, dpi = 180, bg = "white")
}

plot_metric_overlay <- function(raw_metrics, filtered_metrics, metric, section_id) {
  raw_metrics$qc_state <- "raw"
  filtered_metrics$qc_state <- "post_qc"
  plot_df <- rbind(
    raw_metrics[, c(metric, "qc_state"), drop = FALSE],
    filtered_metrics[, c(metric, "qc_state"), drop = FALSE]
  )
  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[metric]], fill = qc_state)) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.55, position = "identity", na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(x = metric, y = "spots", fill = NULL, title = section_id)
}

resolve_project_path <- function(path, base_dir) {
  path <- as.character(path %||% "")
  if (!nzchar(path) || grepl("^/", path)) {
    return(path)
  }
  file.path(base_dir, path)
}

triage_rows <- list()
pre_post_rows <- list()
figure_paths <- character()
report_sections <- list()
retention_values <- suppressWarnings(as.numeric(qc_summary$retention))
unbalanced_global <- (length(retention_values[is.finite(retention_values)]) >= 2L) &&
  (max(retention_values, na.rm = TRUE) - min(retention_values, na.rm = TRUE) > 0.25)

for (idx in seq_len(nrow(qc_summary))) {
  row <- qc_summary[idx, , drop = FALSE]
  section_id <- spatial_safe_id(spatial_cell(row, "section_id", sprintf("section_%s", idx)))
  sample_id <- spatial_safe_id(spatial_cell(row, "sample_id", section_id))
  stem <- spatial_safe_id(spatial_cell(row, "object_stem", sample_id))
  raw_path <- resolve_project_path(spatial_cell(row, "raw_rds"), cfg$project_root)
  post_qc_path <- resolve_project_path(spatial_cell(row, "post_qc_rds"), cfg$project_root)
  if (!file.exists(raw_path) || !file.exists(post_qc_path)) {
    stop(sprintf("missing raw or post-QC RDS for section_id=%s", section_id), call. = FALSE)
  }

  raw_obj <- readRDS(raw_path)
  post_obj <- readRDS(post_qc_path)
  raw_metrics <- compute_spot_qc_metrics(raw_obj)
  post_metrics <- compute_spot_qc_metrics(post_obj)

  section_figures <- character()
  for (metric in c("nFeature_Spatial", "nCount_Spatial", "percent.mito")) {
    fig_path <- file.path(cfg$spatial_post_qc_overlay_dir, sprintf("%s__raw_vs_post_qc_%s.png", stem, metric))
    save_metric_plot(plot_metric_overlay(raw_metrics, post_metrics, metric, section_id), fig_path, width = 6.4, height = 4.4)
    section_figures <- c(section_figures, fig_path)
  }
  mask_rel <- spatial_cell(row, "mask_png")
  mask_path <- if (nzchar(mask_rel)) resolve_project_path(mask_rel, cfg$project_root) else ""
  if (file.exists(mask_path)) {
    section_figures <- c(mask_path, section_figures)
  }
  figure_paths <- c(figure_paths, section_figures)

  threshold <- spatial_threshold_for_section(thresholds, section_id, sample_id, cfg)
  pre_row <- if (nrow(pre_triage) > 0) pre_triage[pre_triage$section_id == section_id, , drop = FALSE] else data.frame(stringsAsFactors = FALSE)
  if (nrow(pre_row) == 0 && nrow(pre_triage) > 0) {
    pre_row <- pre_triage[pre_triage$section_id == sample_id, , drop = FALSE]
  }
  pre_high_mito <- nrow(pre_row) > 0 && spatial_bool(spatial_cell(pre_row[1, , drop = FALSE], "high_mito", "false"), default = FALSE)
  pre_low_feature <- nrow(pre_row) > 0 && spatial_bool(spatial_cell(pre_row[1, , drop = FALSE], "low_feature", "false"), default = FALSE)
  n_drop_high_mito <- spatial_numeric_or(row$n_drop_high_mito, 0)
  n_drop_low_feature <- spatial_numeric_or(row$n_drop_low_feature, 0)
  retention <- spatial_numeric_or(row$retention, NA_real_)
  boundary_drop <- identical(spatial_cell(row, "spatial_aware_summary"), "boundary_drop_detected") ||
    spatial_numeric_or(row$spatial_aware_max_cluster_frac, 0) > 0.1
  drift <- (pre_high_mito && n_drop_high_mito == 0) || (pre_low_feature && n_drop_low_feature == 0)

  triage_rows[[length(triage_rows) + 1L]] <- data.frame(
    section_id = section_id,
    excessive_drop = isTRUE(is.finite(retention) && retention < threshold$excessive_drop_threshold),
    boundary_drop = isTRUE(boundary_drop),
    drift_from_pre_qc = isTRUE(drift),
    unbalanced_sections = isTRUE(unbalanced_global),
    stringsAsFactors = FALSE
  )
  pre_post_rows[[length(pre_post_rows) + 1L]] <- data.frame(
    section_id = section_id,
    pre_high_mito = pre_high_mito,
    pre_low_feature = pre_low_feature,
    n_drop_high_mito = n_drop_high_mito,
    n_drop_low_feature = n_drop_low_feature,
    retention = retention,
    stringsAsFactors = FALSE
  )

  rel_figures <- vapply(section_figures, function(path) relative_path_local(path, cfg$project_root), character(1))
  report_sections[[section_id]] <- c(
    sprintf("### %s", section_id),
    "",
    sprintf("- Retention: %.3f (%s of %s spots).", retention, spatial_cell(row, "n_after"), spatial_cell(row, "n_before")),
    sprintf("- Spatial-aware filter: %s / %s.", spatial_cell(row, "spatial_aware_flag", "false"), spatial_cell(row, "spatial_aware_summary", "disabled")),
    "",
    paste(sprintf("![](%s)", rel_figures), collapse = "\n\n")
  )
}

triage_df <- do.call(rbind, triage_rows)
pre_post_df <- do.call(rbind, pre_post_rows)
triage_df <- triage_df[, c("section_id", "excessive_drop", "boundary_drop", "drift_from_pre_qc", "unbalanced_sections"), drop = FALSE]
spatial_write_tsv(triage_df, cfg$post_filter_triage_tsv)

report_lines <- c(
  "# Spatial Post-QC Report",
  "",
  sprintf("- QC filter summary: `%s`", relative_path_local(cfg$qc_filter_summary_tsv, cfg$project_root)),
  sprintf("- Post-filter triage: `%s`", relative_path_local(cfg$post_filter_triage_tsv, cfg$project_root)),
  "",
  "## Current Thresholds",
  render_markdown_table_local(thresholds),
  "",
  "## Filter Summary",
  render_markdown_table_local(qc_summary),
  "",
  "## Post-Filter Triage",
  render_markdown_table_local(triage_df),
  "",
  "## Pre-QC Alignment",
  render_markdown_table_local(pre_post_df),
  "",
  unlist(report_sections, use.names = FALSE),
  "",
  "## Next Step",
  sprintf("Approve `spatial_post_qc` in `%s`, adjust `%s`, or rerun intake if the retained tissue layout is not acceptable.", relative_path_local(cfg$eda_gate_file, cfg$project_root), relative_path_local(cfg$spatial_qc_threshold_file, cfg$project_root))
)
write_markdown_local(report_lines, cfg$post_qc_report_md)

write_manifest_local(
  manifest_path = cfg$module_02_eda_manifest_path,
  new_outputs = list(
    report = build_output_entry(cfg$post_qc_report_md, "md", module_name, "human-readable spatial post-filter QC report", base_dir = cfg$project_root),
    post_filter_triage = build_output_entry(cfg$post_filter_triage_tsv, "tsv", module_name, "one row per section post-filter triage", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df)),
    figures_dir = build_output_entry(cfg$spatial_post_qc_overlay_dir, "directory", module_name, "raw-vs-filtered QC figures", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    qc_filter_summary = cfg$qc_filter_summary_tsv,
    pre_qc_triage = cfg$pre_qc_triage_tsv,
    spatial_qc_threshold_file = cfg$spatial_qc_threshold_file
  ),
  version = cfg$module_02_version,
  depends_on = list(
    spatial_01b_qc_filter = cfg$module_02_filter_manifest_path,
    spatial_01a_pre_qc_eda = cfg$module_01a_manifest_path
  )
)

message(sprintf("spatial post-QC EDA complete: %s sections", nrow(triage_df)))
