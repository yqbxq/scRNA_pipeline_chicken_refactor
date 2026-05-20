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
module_name <- "spatial_07f_deconvolution_validation"
prepare_dirs_spatial(cfg)

reference <- load_spatial_reference_inventory_st(cfg)
if (!identical(reference$status, "ok")) {
  status <- "skipped_no_reference"
  reason <- reference$reason
} else if (!requireNamespace("scDesign3", quietly = TRUE)) {
  status <- "skipped_no_packages"
  reason <- "scDesign3 is not installed in this R environment"
} else {
  status <- "failed_synthesis"
  reason <- "scDesign3 package is available, but formal synthetic-spot validation is gated to the validation runtime"
}

validation_id <- "deconv_validation_default"
out_dir <- file.path(cfg$spatial_deconv_validation_table_dir, validation_id)
ensure_dir(out_dir)
truth_tsv <- file.path(out_dir, "synthetic_truth.tsv")
summary_tsv <- file.path(out_dir, "method_summary.tsv")
st07_write_tsv(st07_empty_df(c("spot_id", "cell_type", "true_proportion")), truth_tsv)
summary_df <- data.frame(method = character(), mean_rmse = numeric(), mean_pearson = numeric(), runtime_sec = numeric(), stringsAsFactors = FALSE)
st07_write_tsv(summary_df, summary_tsv)
manifest_tsv <- file.path(cfg$spatial_deconv_validation_table_dir, "validation_manifest.tsv")
report_path <- file.path(cfg$spatial_deconv_validation_report_dir, "report.md")
manifest_df <- data.frame(
  validation_id = validation_id,
  synthetic_n_spots = cfg$spatial_validation_n_spots,
  methods_validated = "",
  mean_rmse_best_method = NA_real_,
  status = status,
  reason = reason,
  synthetic_truth_tsv = truth_tsv,
  method_summary_tsv = summary_tsv,
  stringsAsFactors = FALSE
)
st07_write_tsv(manifest_df, manifest_tsv)

report_lines <- c(
  "# Spatial Deconvolution Validation",
  "",
  sprintf("- status: `%s`", status),
  sprintf("- synthetic_n_spots: `%s`", cfg$spatial_validation_n_spots),
  sprintf("- dirichlet_alpha: `%s`", cfg$spatial_validation_dirichlet_alpha),
  "",
  "## Manifest",
  render_markdown_table_local(manifest_df)
)
write_markdown_local(report_lines, report_path)

st07_write_manifest_local(
  manifest_path = cfg$module_07f_deconvolution_validation_manifest_path,
  new_outputs = list(
    validation_manifest = build_output_entry(manifest_tsv, "tsv", module_name, "deconvolution validation manifest", base_dir = cfg$project_root, schema = infer_schema_from_df(manifest_df)),
    synthetic_truth = build_output_entry(truth_tsv, "tsv", module_name, "synthetic spot truth proportions", base_dir = cfg$project_root),
    method_summary = build_output_entry(summary_tsv, "tsv", module_name, "per-method validation metrics", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df)),
    report = build_output_entry(report_path, "md", module_name, "spatial deconvolution validation report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(spatial_reference_inventory = cfg$spatial_reference_inventory_file),
  version = cfg$module_07_version,
  depends_on = list(spatial_reference_inventory = cfg$spatial_reference_inventory_file)
)

message(sprintf("spatial deconvolution validation complete: %s", status))
