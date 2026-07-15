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
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))

load_required_packages(c("dplyr", "ggplot2", "tibble", "jsonlite"))

cfg <- get_single_script_config_03()
module_name <- "03g_annotation_eda"
prepare_dirs_03(cfg)

panorama_spec <- panorama_layer_spec(cfg)
manifest_03f <- read_manifest_local(cfg$module_03f_manifest_path)
annotation_table_path <- resolve_output_local(manifest_03f, "annotation_table_tsv")
annotation_table <- read_tsv_optional(annotation_table_path)
cluster_marker_qc <- read_tsv_optional(resolve_output_local(read_manifest_local(cfg$module_03c3_manifest_path), "cluster_marker_qc_tsv"))
cluster_marker_risk <- read_tsv_optional(resolve_output_local(read_manifest_local(cfg$module_03d_marker_risk_manifest_path), "cluster_marker_risk_tsv"))

if (nrow(annotation_table) == 0) {
  summary_df <- data.frame(
    layer_id = panorama_spec$layer_id,
    cluster_count = 0L,
    annotated_cluster_n = 0L,
    high_confidence_n = 0L,
    medium_confidence_n = 0L,
    low_confidence_n = 0L,
    stringsAsFactors = FALSE
  )
} else {
  conf <- tolower(annotation_table$confidence)
  summary_df <- data.frame(
    layer_id = panorama_spec$layer_id,
    cluster_count = nrow(annotation_table),
    annotated_cluster_n = sum(nzchar(annotation_table$manual_annotation), na.rm = TRUE),
    high_confidence_n = sum(conf %in% c("high", "确定"), na.rm = TRUE),
    medium_confidence_n = sum(conf %in% c("medium", "暂定"), na.rm = TRUE),
    low_confidence_n = sum(conf %in% c("low", "未定"), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
summary_df$annotated_fraction <- ifelse(summary_df$cluster_count > 0, summary_df$annotated_cluster_n / summary_df$cluster_count, NA_real_)

review_table <- annotation_table
if (nrow(review_table) > 0 && nrow(cluster_marker_qc) > 0) {
  review_table <- dplyr::left_join(
    review_table,
    cluster_marker_qc[, intersect(c("cluster", "strict_marker_n", "deg_quality_flag", "risk_level"), colnames(cluster_marker_qc)), drop = FALSE],
    by = "cluster"
  )
}
if (nrow(review_table) > 0 && nrow(cluster_marker_risk) > 0) {
  review_table <- dplyr::left_join(
    review_table,
    cluster_marker_risk[, intersect(c("cluster", "marker_quality_flag", "manual_review_reason"), colnames(cluster_marker_risk)), drop = FALSE],
    by = "cluster"
  )
}

summary_tsv <- file.path(cfg$annotation_table_dir_layer, "annotation_eda_summary.tsv")
review_tsv <- file.path(cfg$annotation_table_dir_layer, "annotation_review_table.tsv")
report_path <- file.path(cfg$annotation_report_dir_layer, "report.md")
write_tsv_local(summary_df, summary_tsv)
write_tsv_local(review_table, review_tsv)

confidence_counts <- if (nrow(annotation_table) > 0) {
  as.data.frame(table(annotation_table$confidence), stringsAsFactors = FALSE)
} else {
  data.frame()
}
colnames(confidence_counts) <- if (ncol(confidence_counts) == 2) c("confidence", "cluster_n") else colnames(confidence_counts)

report_lines <- c(
  "# 03g Panorama Annotation EDA Report",
  "",
  sprintf("- layer_id: `%s`", panorama_spec$layer_id),
  sprintf("- annotation_table: `%s`", annotation_table_path),
  sprintf("- summary_tsv: `%s`", summary_tsv),
  sprintf("- cluster_count: `%s`", summary_df$cluster_count[[1]]),
  sprintf("- annotated_fraction: `%s`", fmt_num(summary_df$annotated_fraction[[1]], 3)),
  "",
  "## Confidence Counts",
  render_markdown_table_local(confidence_counts),
  "",
  "## Cluster Review",
  render_markdown_table_local(review_table[, intersect(c(
    "cluster", "n_cells", "manual_annotation", "manual_annotation_zh", "confidence",
    "strict_marker_n", "deg_quality_flag", "marker_quality_flag", "reason", "manual_review_reason"
  ), colnames(review_table)), drop = FALSE]),
  "",
  "## Review Action",
  sprintf("Review this report and `%s`; approve the `annotation` gate before running 04_subcluster.", cfg$eda_gate_file)
)
write_markdown_local(report_lines, report_path)

write_manifest_local(
  manifest_path = cfg$module_03g_manifest_path,
  new_outputs = list(
    report = build_output_entry(report_path, "md", module_name, "human-readable panorama annotation EDA report", base_dir = cfg$project_root),
    summary_tsv = build_output_entry(summary_tsv, "tsv", module_name, "one row panorama annotation summary", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df)),
    annotation_review_table_tsv = build_output_entry(review_tsv, "tsv", module_name, "manual annotation review table with marker risk context", base_dir = cfg$project_root, schema = infer_schema_from_df(review_table))
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(module_03f_manifest = cfg$module_03f_manifest_path, annotation_table_tsv = annotation_table_path),
  version = cfg$module_version,
  depends_on = list(module_03f = cfg$module_03f_manifest_path)
)

write_manifest_local(
  manifest_path = cfg$module_03e_manifest_path,
  new_outputs = list(
    report = build_output_entry(report_path, "md", "03g_annotation_eda", "human-readable panorama annotation EDA report", base_dir = cfg$project_root),
    summary_tsv = build_output_entry(summary_tsv, "tsv", "03g_annotation_eda", "one row panorama annotation summary", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df)),
    annotation_review_table_tsv = build_output_entry(review_tsv, "tsv", "03g_annotation_eda", "manual annotation review table with marker risk context", base_dir = cfg$project_root, schema = infer_schema_from_df(review_table))
  ),
  module_name = "03e_annotation_eda",
  base_dir = cfg$project_root,
  inputs = list(module_03g_manifest = cfg$module_03g_manifest_path, annotation_table_tsv = annotation_table_path),
  version = cfg$module_version,
  depends_on = list(module_03g = cfg$module_03g_manifest_path)
)

message("03g complete. Report: ", report_path)
