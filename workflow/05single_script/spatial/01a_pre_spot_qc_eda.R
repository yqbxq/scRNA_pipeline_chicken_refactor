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
module_name <- "spatial_01a_pre_qc_eda"
prepare_dirs_spatial(cfg)
set.seed(cfg$random_seed)

thresholds <- spatial_read_tsv(cfg$spatial_qc_threshold_file)
if (nrow(thresholds) == 0) {
  thresholds <- data.frame(
    section_id = "__DEFAULT__",
    qc_min_nfeature = cfg$qc_min_nfeature_default,
    qc_min_ncount = cfg$qc_min_ncount_default,
    qc_max_mito_pct = cfg$qc_max_mito_pct_default,
    stringsAsFactors = FALSE
  )
}
for (col in c("section_id", "qc_min_nfeature", "qc_min_ncount", "qc_max_mito_pct")) {
  if (!col %in% colnames(thresholds)) {
    thresholds[[col]] <- ""
  }
}

threshold_for_section <- function(section_id) {
  hit <- thresholds[thresholds$section_id == section_id, , drop = FALSE]
  if (nrow(hit) == 0) {
    hit <- thresholds[thresholds$section_id == "__DEFAULT__", , drop = FALSE]
  }
  if (nrow(hit) == 0) {
    hit <- thresholds[1, , drop = FALSE]
  }
  min_feature <- suppressWarnings(as.numeric(hit$qc_min_nfeature[[1]]))
  min_count <- suppressWarnings(as.numeric(hit$qc_min_ncount[[1]]))
  max_mito <- suppressWarnings(as.numeric(hit$qc_max_mito_pct[[1]]))
  if (is.na(min_feature)) min_feature <- cfg$qc_min_nfeature_default
  if (is.na(min_count)) min_count <- cfg$qc_min_ncount_default
  if (is.na(max_mito)) max_mito <- cfg$qc_max_mito_pct_default
  list(qc_min_nfeature = min_feature, qc_min_ncount = min_count, qc_max_mito_pct = max_mito)
}

safe_fraction <- function(x) {
  if (length(x) == 0) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

save_metric_plot <- function(plot_obj, path, width = 7, height = 5) {
  ensure_dir(dirname(path))
  ggplot2::ggsave(path, plot_obj, width = width, height = height, dpi = 180, bg = "white")
}

raw_paths <- list.files(cfg$spatial_checkpoint_dir, pattern = "_raw\\.rds$", full.names = TRUE)
if (length(raw_paths) == 0) {
  stop(sprintf("no raw spatial objects found in %s", cfg$spatial_checkpoint_dir), call. = FALSE)
}

triage_rows <- list()
metrics_paths <- character()
figure_paths <- character()
report_sections <- list()

for (raw_path in raw_paths) {
  obj <- readRDS(raw_path)
  metrics <- compute_spot_qc_metrics(obj)
  section_id <- unique(metrics$section_id)[[1]]
  section_id <- spatial_safe_id(section_id)
  threshold <- threshold_for_section(section_id)

  metrics_path <- file.path(cfg$spatial_pre_qc_table_dir, sprintf("spot_metrics_%s.tsv", section_id))
  spatial_write_tsv(metrics, metrics_path)
  metrics_paths <- c(metrics_paths, metrics_path)

  metric_names <- c("nFeature_Spatial", "nCount_Spatial", "log10_nCount", "percent.mito")
  section_figures <- character()
  for (metric in metric_names) {
    violin_path <- file.path(cfg$spatial_pre_qc_figure_dir, sprintf("%s__violin_%s.png", section_id, metric))
    spatial_path <- file.path(cfg$spatial_pre_qc_figure_dir, sprintf("%s__spatial_%s.png", section_id, metric))
    save_metric_plot(spatial_eda_violin(metrics, metric), violin_path, width = 6, height = 4)
    save_metric_plot(spatial_eda_spatial_plot(metrics, metric), spatial_path, width = 5, height = 5)
    section_figures <- c(section_figures, violin_path, spatial_path)
  }
  scatter_path <- file.path(cfg$spatial_pre_qc_figure_dir, sprintf("%s__scatter_ncount_nfeature.png", section_id))
  save_metric_plot(spatial_eda_scatter(metrics), scatter_path, width = 6, height = 5)
  section_figures <- c(section_figures, scatter_path)
  figure_paths <- c(figure_paths, section_figures)

  high_mito <- safe_fraction(metrics$percent.mito > threshold$qc_max_mito_pct) > 0.25
  low_feature <- safe_fraction(metrics$nFeature_Spatial < threshold$qc_min_nfeature) > 0.35
  sparse_tissue <- safe_fraction(metrics$in_tissue %in% TRUE) < 0.5
  sectioning_artifact <- safe_fraction(!(metrics$in_tissue %in% TRUE)) > 0.4
  mito_qc_unreliable <- !any(metrics$mito_qc_valid %in% TRUE, na.rm = TRUE)
  triage_rows[[length(triage_rows) + 1L]] <- data.frame(
    section_id = section_id,
    qc_min_nfeature = threshold$qc_min_nfeature,
    qc_min_ncount = threshold$qc_min_ncount,
    qc_max_mito_pct = threshold$qc_max_mito_pct,
    spots = nrow(metrics),
    median_nFeature = stats::median(metrics$nFeature_Spatial, na.rm = TRUE),
    median_nCount = stats::median(metrics$nCount_Spatial, na.rm = TRUE),
    median_percent_mito = stats::median(metrics$percent.mito, na.rm = TRUE),
    high_mito = isTRUE(high_mito),
    low_feature = isTRUE(low_feature),
    sparse_tissue = isTRUE(sparse_tissue),
    sectioning_artifact = isTRUE(sectioning_artifact),
    mito_qc_unreliable = isTRUE(mito_qc_unreliable),
    stringsAsFactors = FALSE
  )

  rel_figures <- vapply(section_figures, function(path) relative_path_local(path, cfg$project_root), character(1))
  report_sections[[section_id]] <- c(
    sprintf("### %s", section_id),
    "",
    sprintf("- Current thresholds: min features %s; min counts %s; max mito %.2f%%.", threshold$qc_min_nfeature, threshold$qc_min_ncount, threshold$qc_max_mito_pct),
    sprintf("- Spots: %s; median features %.1f; median counts %.1f; median mito %.2f%%.", nrow(metrics), stats::median(metrics$nFeature_Spatial, na.rm = TRUE), stats::median(metrics$nCount_Spatial, na.rm = TRUE), stats::median(metrics$percent.mito, na.rm = TRUE)),
    "",
    paste(sprintf("![](%s)", rel_figures), collapse = "\n\n")
  )
}

triage_df <- do.call(rbind, triage_rows)
spatial_write_tsv(triage_df, cfg$pre_qc_triage_tsv)

report_lines <- c(
  "# Spatial Pre-QC Report",
  "",
  sprintf("- Threshold file: `%s`", relative_path_local(cfg$spatial_qc_threshold_file, cfg$project_root)),
  sprintf("- Triage table: `%s`", relative_path_local(cfg$pre_qc_triage_tsv, cfg$project_root)),
  "",
  "## Current Thresholds",
  render_markdown_table_local(thresholds),
  "",
  "## Triage Signals",
  render_markdown_table_local(triage_df),
  "",
  unlist(report_sections, use.names = FALSE),
  "",
  "## Next Step",
  sprintf("Edit `%s` as needed, then approve `spatial_pre_qc` in `%s` before spatial filtering.", relative_path_local(cfg$spatial_qc_threshold_file, cfg$project_root), relative_path_local(cfg$eda_gate_file, cfg$project_root))
)
write_markdown_local(report_lines, cfg$pre_qc_report_md)

write_manifest_local(
  manifest_path = cfg$module_01a_manifest_path,
  new_outputs = list(
    report = build_output_entry(cfg$pre_qc_report_md, "md", module_name, "human-readable spatial pre-QC report", base_dir = cfg$project_root),
    triage = build_output_entry(cfg$pre_qc_triage_tsv, "tsv", module_name, "one row per section spatial pre-QC triage", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df)),
    spot_metrics_dir = build_output_entry(cfg$spatial_pre_qc_table_dir, "directory", module_name, "per-section spot metrics tables", base_dir = cfg$project_root),
    figures_dir = build_output_entry(cfg$spatial_pre_qc_figure_dir, "directory", module_name, "spatial pre-QC figures", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    raw_object_dir = cfg$spatial_checkpoint_dir,
    spatial_qc_threshold_file = cfg$spatial_qc_threshold_file
  ),
  version = cfg$module_version,
  depends_on = list(module_01 = cfg$module_01_manifest_path)
)

message(sprintf("spatial pre-QC EDA complete: %s sections", nrow(triage_df)))
