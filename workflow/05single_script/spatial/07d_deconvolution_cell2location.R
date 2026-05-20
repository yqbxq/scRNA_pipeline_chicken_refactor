#!/usr/bin/env Rscript

.script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(.script_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_common.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "project_paths_spatial.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_deconv_utils.R"), encoding = "UTF-8")

cfg <- get_spatial_script_config()
module_name <- "spatial_07d_deconvolution_cell2location"
prepare_dirs_spatial(cfg)

pairs <- read_deconv_pairs_st(cfg)
rows <- list()
if (nrow(pairs) == 0) {
  pair <- data.frame(deconv_id = "", enabled = "yes", tool = "cell2location", stringsAsFactors = FALSE)
  outputs <- write_empty_deconv_outputs_st(cfg, cfg$spatial_c2l_table_dir, pair, "all", "cell2location", "skipped_no_deconv_pairs", "deconv_pairs.tsv has no enabled rows")
  rows[[1L]] <- st07_deconv_row(pair, "all", "cell2location", "skipped_no_deconv_pairs", "deconv_pairs.tsv has no enabled rows", outputs = outputs)
} else {
  sidecar <- file.path(cfg$pipeline_root, "workflow", "04python", "cell2location_pipeline.py")
  use_gpu <- tolower(cfg$spatial_c2l_use_gpu) %in% c("yes", "true", "1", "on")
  args <- c(sidecar, "--check-env", if (use_gpu) "--use-gpu" else character(0))
  check <- tryCatch(system2(cfg$py_cell2location_bin, args = args, stdout = TRUE, stderr = TRUE), error = function(e) structure(conditionMessage(e), status = 127))
  code <- attr(check, "status")
  if (is.null(code)) code <- 0L
  check_text <- paste(as.character(check), collapse = " ")
  status <- if (identical(as.integer(code), 0L)) "failed_sidecar" else if (grepl("skipped_no_gpu", check_text)) "skipped_no_gpu" else if (grepl("skipped_no_python_env", check_text) || identical(as.integer(code), 127L)) "skipped_no_python_env" else "failed_sidecar"
  reason <- if (identical(as.integer(code), 0L)) "cell2location environment exists; formal training is gated to server runtime execution" else check_text
  for (idx in seq_len(nrow(pairs))) {
    pair <- pairs[idx, , drop = FALSE]
    section <- st07_scalar(pair$section_filter, "all")
    row_status <- status
    row_reason <- reason
    if (!deconv_pair_enabled_st(pair) || !deconv_pair_allows_tool_st(pair, "cell2location")) {
      row_status <- "skipped_disabled"
      row_reason <- "pair disabled or tool filter does not include cell2location"
    }
    outputs <- write_empty_deconv_outputs_st(cfg, cfg$spatial_c2l_table_dir, pair, section, "cell2location", row_status, row_reason)
    rows[[length(rows) + 1L]] <- st07_deconv_row(pair, section, "cell2location", row_status, row_reason, outputs = outputs)
  }
}

manifest_df <- do.call(rbind, rows)
manifest_tsv <- file.path(cfg$spatial_c2l_table_dir, "cell2location_manifest.tsv")
report_path <- file.path(cfg$spatial_c2l_figure_dir, "cell2location_report.md")
write_deconv_module_manifest_st(
  cfg,
  manifest_df,
  manifest_tsv,
  report_path,
  module_name,
  cfg$module_07d_deconvolution_cell2location_manifest_path,
  "cell2location_manifest",
  depends_on = list(spatial_07a_deconvolution_rctd = cfg$module_07a_deconvolution_rctd_manifest_path)
)

message(sprintf("spatial cell2location deconvolution complete: %d manifest rows", nrow(manifest_df)))
