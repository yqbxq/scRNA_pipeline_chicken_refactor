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
module_name <- "03e_annotation_eda"
prepare_dirs_03(cfg)
set.seed(cfg$random_seed)

panorama_spec <- panorama_layer_spec(cfg)
manifest_03d <- read_manifest_local(cfg$module_03d_manifest_path)
annotation_table_path <- resolve_output_local(manifest_03d, "annotation_table_tsv")
annotation_table <- read_tsv_optional(annotation_table_path)

if (nrow(annotation_table) == 0) {
  summary_df <- data.frame(
    layer_id = panorama_spec$layer_id,
    cluster_count = 0L,
    determined_n = 0L,
    tentative_n = 0L,
    undetermined_n = 0L,
    support_n = 0L,
    conflict_n = 0L,
    unrelated_n = 0L,
    stringsAsFactors = FALSE
  )
} else {
  summary_df <- data.frame(
    layer_id = panorama_spec$layer_id,
    cluster_count = nrow(annotation_table),
    determined_n = sum(annotation_table$confidence == "确定", na.rm = TRUE),
    tentative_n = sum(annotation_table$confidence == "暂定", na.rm = TRUE),
    undetermined_n = sum(annotation_table$confidence == "未定", na.rm = TRUE),
    support_n = sum(annotation_table$relation == "支持", na.rm = TRUE),
    conflict_n = sum(annotation_table$relation == "冲突", na.rm = TRUE),
    unrelated_n = sum(annotation_table$relation == "无关", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
summary_df$determined_fraction <- ifelse(summary_df$cluster_count > 0, summary_df$determined_n / summary_df$cluster_count, NA_real_)

summary_tsv <- file.path(cfg$annotation_table_dir_layer, "annotation_eda_summary.tsv")
report_path <- file.path(cfg$annotation_report_dir_layer, "report.md")
write_tsv_local(summary_df, summary_tsv)

cluster_review <- if (nrow(annotation_table) == 0) {
  "No annotation rows available."
} else {
  render_markdown_table_local(annotation_table[, intersect(c(
    "cluster_id", "candidate_celltype", "confidence", "relation", "n_cells",
    "data_driven_candidate", "literature_candidate", "module_score_check"
  ), colnames(annotation_table)), drop = FALSE])
}

confidence_counts <- if (nrow(annotation_table) > 0) as.data.frame(table(annotation_table$confidence), stringsAsFactors = FALSE) else data.frame()
relation_counts <- if (nrow(annotation_table) > 0) as.data.frame(table(annotation_table$relation), stringsAsFactors = FALSE) else data.frame()

report_lines <- c(
  "# 03e Panorama Annotation EDA Report",
  "",
  sprintf("- layer_id: `%s`", panorama_spec$layer_id),
  sprintf("- annotation_table: `%s`", annotation_table_path),
  sprintf("- summary_tsv: `%s`", summary_tsv),
  sprintf("- cluster_count: `%s`", summary_df$cluster_count[[1]]),
  sprintf("- determined_fraction: `%s`", fmt_num(summary_df$determined_fraction[[1]], 3)),
  "",
  "## Confidence Counts",
  render_markdown_table_local(confidence_counts),
  "",
  "## Relation Counts",
  render_markdown_table_local(relation_counts),
  "",
  "## Cluster Review",
  cluster_review,
  "",
  "## Review Action",
  sprintf("Review this report and `%s`; approve the `annotation` gate before running 04_subcluster.", cfg$eda_gate_file)
)
write_markdown_local(report_lines, report_path)

write_manifest_local(
  manifest_path = cfg$module_03e_manifest_path,
  new_outputs = list(
    report = build_output_entry(report_path, "md", module_name, "human-readable panorama annotation EDA report", base_dir = cfg$project_root),
    summary_tsv = build_output_entry(summary_tsv, "tsv", module_name, "one row panorama annotation summary", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df))
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(module_03d_manifest = cfg$module_03d_manifest_path, annotation_table_tsv = annotation_table_path),
  version = cfg$module_version,
  depends_on = list(module_03d = cfg$module_03d_manifest_path)
)

message("03e 完成。report: ", report_path)
