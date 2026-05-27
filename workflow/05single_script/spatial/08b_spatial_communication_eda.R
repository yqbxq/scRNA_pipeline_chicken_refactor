#!/usr/bin/env Rscript

.script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(.script_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/helpers/manifest_utils.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/helpers/metadata_io.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/helpers/commot_signal_utils.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_common.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "project_paths_spatial.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_deconv_utils.R"), encoding = "UTF-8")

cfg <- get_spatial_script_config()
module_name <- "spatial_08b_spatial_communication_eda"
prepare_dirs_spatial(cfg)

commot_lr_tsv <- file.path(cfg$spatial_commot_table_dir, "commot_lr.tsv")
commot_df <- spatial_read_tsv(commot_lr_tsv)
summary_df <- summarize_commot_signal(commot_df, signal_threshold = cfg$commot_signal_threshold)

summary_tsv <- file.path(cfg$spatial_commot_table_dir, "commot_spatial_summary.tsv")
report_md <- file.path(cfg$spatial_commot_report_dir, "report.md")
spatial_write_tsv(summary_df, summary_tsv)

status_counts <- if (nrow(summary_df) == 0) {
  data.frame(status = "empty", n = 0L, stringsAsFactors = FALSE)
} else {
  as.data.frame(table(status = summary_df$status), stringsAsFactors = FALSE)
}
support_counts <- if (nrow(summary_df) == 0) {
  data.frame(spatial_support = "no", n = 0L, stringsAsFactors = FALSE)
} else {
  as.data.frame(table(spatial_support = summary_df$spatial_support), stringsAsFactors = FALSE)
}

ensure_dir(dirname(report_md))
writeLines(c(
  "# Spatial COMMOT Communication",
  "",
  sprintf("- LR axes evaluated: %d", nrow(summary_df)),
  sprintf("- Signal threshold: %s", cfg$commot_signal_threshold),
  "",
  "## Status Counts",
  "",
  capture.output(print(status_counts, row.names = FALSE)),
  "",
  "## Spatial Support Counts",
  "",
  capture.output(print(support_counts, row.names = FALSE))
), report_md, useBytes = TRUE)

if (file.exists(cfg$module_08b_spatial_communication_eda_manifest_path)) {
  unlink(cfg$module_08b_spatial_communication_eda_manifest_path)
}
st07_write_manifest_local(
  manifest_path = cfg$module_08b_spatial_communication_eda_manifest_path,
  new_outputs = list(
    commot_spatial_summary_tsv = build_output_entry(summary_tsv, "tsv", module_name, "one row per COMMOT LR axis spatial-support decision", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df)),
    report_md = build_output_entry(report_md, "md", module_name, "spatial COMMOT EDA report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(commot_lr_tsv = commot_lr_tsv, commot_signal_threshold = cfg$commot_signal_threshold),
  version = cfg$module_07_version,
  depends_on = list(spatial_08f_commot_spatial = cfg$module_08f_commot_spatial_manifest_path)
)

message(sprintf("08b COMMOT EDA completed. summary=%s", summary_tsv))
