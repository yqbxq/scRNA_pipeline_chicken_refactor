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
module_name <- "spatial_07c_deconvolution_card"
manifest_df <- run_deconv_method_scaffold_st(cfg, "card", "CARD", cfg$spatial_card_table_dir, "failed_card_deconvolution")
manifest_tsv <- file.path(cfg$spatial_card_table_dir, "card_manifest.tsv")
report_path <- file.path(cfg$spatial_card_figure_dir, "card_report.md")

write_deconv_module_manifest_st(
  cfg,
  manifest_df,
  manifest_tsv,
  report_path,
  module_name,
  cfg$module_07c_deconvolution_card_manifest_path,
  "card_manifest",
  depends_on = list(spatial_07a_deconvolution_rctd = cfg$module_07a_deconvolution_rctd_manifest_path)
)

message(sprintf("spatial CARD deconvolution complete: %d manifest rows", nrow(manifest_df)))
