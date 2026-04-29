#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

source_utf8 <- function(path) source(path, encoding = "UTF-8")

source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "config.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_02.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_03.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_04.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_05.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "deg_utils.R"))

load_required_packages(c("dplyr", "tibble", "jsonlite"))

cfg <- get_single_script_config_05()
module_name <- "05d_deg_eda"
prepare_dirs_05(cfg)

marker_manifest_tsv <- file.path(cfg$marker_discovery_table_dir, "marker_discovery_manifest.tsv")
pseudobulk_manifest_tsv <- file.path(cfg$pseudobulk_table_dir, "pseudobulk_manifest.tsv")
composition_manifest_tsv <- file.path(cfg$composition_table_dir, "composition_manifest.tsv")

marker_df <- read_tsv_optional(marker_manifest_tsv)
pb_df <- read_tsv_optional(pseudobulk_manifest_tsv)
comp_df <- read_tsv_optional(composition_manifest_tsv)

key_cols <- c("layer_id", "comparison_id")
empty_keyed <- empty_df_05(key_cols)
keys <- dplyr::bind_rows(
  if (all(key_cols %in% colnames(marker_df))) marker_df[, key_cols, drop = FALSE] else empty_keyed,
  if (all(key_cols %in% colnames(pb_df))) pb_df[, key_cols, drop = FALSE] else empty_keyed,
  if (all(key_cols %in% colnames(comp_df))) comp_df[, key_cols, drop = FALSE] else empty_keyed
)
keys <- unique(keys[nzchar(keys$layer_id) & nzchar(keys$comparison_id), , drop = FALSE])

status_rows <- list()
summary_rows <- list()
for (idx in seq_len(nrow(keys))) {
  layer_id <- keys$layer_id[[idx]]
  comparison_id <- keys$comparison_id[[idx]]
  marker_hit <- marker_df[marker_df$layer_id == layer_id & marker_df$comparison_id == comparison_id, , drop = FALSE]
  pb_hit <- pb_df[pb_df$layer_id == layer_id & pb_df$comparison_id == comparison_id, , drop = FALSE]
  comp_hit <- comp_df[comp_df$layer_id == layer_id & comp_df$comparison_id == comparison_id, , drop = FALSE]

  marker_status <- if (nrow(marker_hit) > 0) normalize_scalar_value(marker_hit$status[[1]], "available") else "missing"
  pb_status <- if (nrow(pb_hit) > 0) normalize_scalar_value(pb_hit$inference_status[[1]], "skipped") else "missing"
  comp_status <- if (nrow(comp_hit) > 0) normalize_scalar_value(comp_hit$inference_status[[1]], "skipped") else "missing"
  pb_path <- if (nrow(pb_hit) > 0 && "ds_results_tsv" %in% colnames(pb_hit)) normalize_scalar_value(pb_hit$ds_results_tsv[[1]]) else ""
  comp_path <- if (nrow(comp_hit) > 0 && "formal_results_tsv" %in% colnames(comp_hit)) normalize_scalar_value(comp_hit$formal_results_tsv[[1]]) else ""
  marker_path <- if (nrow(marker_hit) > 0 && "exploratory_results_tsv" %in% colnames(marker_hit)) normalize_scalar_value(marker_hit$exploratory_results_tsv[[1]]) else ""

  marker_counts <- count_significant_rows_05(marker_path, alpha = cfg$deg_alpha)
  pb_counts <- count_significant_rows_05(pb_path, alpha = cfg$deg_alpha)
  comp_counts <- count_significant_rows_05(comp_path, alpha = cfg$deg_alpha)

  status_rows[[length(status_rows) + 1]] <- data.frame(
    layer_id = layer_id,
    comparison_id = comparison_id,
    marker_status = marker_status,
    pseudobulk_status = pb_status,
    composition_status = comp_status,
    marker_results_tsv = marker_path,
    pseudobulk_results_tsv = pb_path,
    composition_results_tsv = comp_path,
    stringsAsFactors = FALSE
  )
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    layer_id = layer_id,
    comparison_id = comparison_id,
    exploratory_marker_total_n = marker_counts$total_n,
    exploratory_marker_significant_n = marker_counts$significant_n,
    formal_de_total_n = pb_counts$total_n,
    formal_de_significant_n = pb_counts$significant_n,
    significant_composition_cluster_n = comp_counts$significant_n,
    alpha = cfg$deg_alpha,
    stringsAsFactors = FALSE
  )
}

status_df <- if (length(status_rows) > 0) dplyr::bind_rows(status_rows) else empty_df_05(c(
  "layer_id", "comparison_id", "marker_status", "pseudobulk_status",
  "composition_status", "marker_results_tsv", "pseudobulk_results_tsv", "composition_results_tsv"
))
summary_df <- if (length(summary_rows) > 0) dplyr::bind_rows(summary_rows) else empty_df_05(c(
  "layer_id", "comparison_id", "exploratory_marker_total_n", "exploratory_marker_significant_n",
  "formal_de_total_n", "formal_de_significant_n", "significant_composition_cluster_n", "alpha"
))

paths <- deg_report_paths_05(cfg)
write_tsv_local(status_df, paths$status_matrix_tsv)
write_tsv_local(summary_df, paths$summary_tsv)

forced_or_exploratory <- status_df[
  status_df$pseudobulk_status %in% c("exploratory_only", "exploratory_forced") |
    status_df$composition_status %in% c("exploratory_only", "exploratory_forced"),
  ,
  drop = FALSE
]

report_lines <- c(
  "# 05d DEG EDA",
  "",
  sprintf("- marker_manifest: `%s`", marker_manifest_tsv),
  sprintf("- pseudobulk_manifest: `%s`", pseudobulk_manifest_tsv),
  sprintf("- composition_manifest: `%s`", composition_manifest_tsv),
  sprintf("- status_matrix: `%s`", paths$status_matrix_tsv),
  sprintf("- summary: `%s`", paths$summary_tsv),
  "",
  "## Inference Status Matrix",
  render_markdown_table_local(status_df),
  "",
  "## DEG Summary",
  render_markdown_table_local(summary_df),
  "",
  "## Exploratory Warnings",
  if (nrow(forced_or_exploratory) == 0) "No exploratory-only or forced exploratory results recorded." else render_markdown_table_local(forced_or_exploratory[, c("layer_id", "comparison_id", "pseudobulk_status", "composition_status", "marker_results_tsv"), drop = FALSE])
)
write_markdown_local(report_lines, paths$report_md)

if (file.exists(cfg$module_05d_manifest_path)) {
  unlink(cfg$module_05d_manifest_path)
}
write_manifest_local(
  manifest_path = cfg$module_05d_manifest_path,
  new_outputs = list(
    report = build_output_entry(paths$report_md, "md", module_name, "DEG EDA report", base_dir = cfg$project_root),
    deg_status_matrix_tsv = build_output_entry(paths$status_matrix_tsv, "tsv", module_name, "one row per layer/comparison DEG status", base_dir = cfg$project_root, schema = infer_schema_from_df(status_df)),
    deg_summary_tsv = build_output_entry(paths$summary_tsv, "tsv", module_name, "one row per layer/comparison DEG summary", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df))
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    marker_manifest_tsv = marker_manifest_tsv,
    pseudobulk_manifest_tsv = pseudobulk_manifest_tsv,
    composition_manifest_tsv = composition_manifest_tsv
  ),
  version = cfg$module_version,
  depends_on = list(
    module_05a = cfg$module_05a_manifest_path,
    module_05b = cfg$module_05b_manifest_path,
    module_05c = cfg$module_05c_manifest_path
  )
)

message("05d completed. report: ", paths$report_md)
