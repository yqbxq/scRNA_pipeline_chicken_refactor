#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(utils)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
.script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else ""
.script_dir <- if (nzchar(.script_path)) dirname(normalizePath(.script_path, winslash = "/", mustWork = FALSE)) else getwd()

source_utf8 <- function(path) source(path, encoding = "UTF-8")
source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "h5ad_contract_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "h5ad_export_utils.R"))

upstream_rds <- Sys.getenv("H5AD_UPSTREAM_RDS", unset = "")
target_dir <- Sys.getenv("H5AD_TARGET_DIR", unset = "")
target_manifest <- Sys.getenv("H5AD_TARGET_MANIFEST", unset = file.path(target_dir, "_manifest.json"))
source_module <- Sys.getenv("H5AD_SOURCE_MODULE", unset = "upstream")
output_key <- Sys.getenv("H5AD_OUTPUT_KEY", unset = "object")

if (!nzchar(upstream_rds)) stop("H5AD_UPSTREAM_RDS is required", call. = FALSE)
if (!nzchar(target_dir)) stop("H5AD_TARGET_DIR is required", call. = FALSE)

ensure_dir(target_dir)
obj <- h5ad_read_upstream_object(upstream_rds)
contract <- load_h5ad_contract(Sys.getenv("H5AD_CONTRACT_FILE", unset = "metadata/h5ad_export_contract.tsv"), modality = "scrna")
fail_on <- h5ad_env_or_default("H5AD_EXPORT_CONTRACT_FAIL_ON", h5ad_env_or_default("H5AD_CONTRACT_FAIL_ON", "warn"))

h5ad_path <- file.path(target_dir, sprintf("%s_%s_scrna.h5ad", source_module, output_key))
result <- h5ad_export_one(obj, h5ad_path, "scrna", contract, fail_on = if (identical(fail_on, "error")) "warn" else fail_on)
if (!isTRUE(result$contract_check$passed)) {
  h5ad_write_skip_json(result$skipped_json, result)
}

summary <- h5ad_summary_row("scrna", upstream_rds, result)
summary_tsv <- file.path(target_dir, "h5ad_export_summary.tsv")
write_tsv_local(summary, summary_tsv)

outputs <- list(
  export_summary_tsv = build_output_entry(summary_tsv, "tsv", "90a_export_h5ad_scrna", "H5AD export status summary", base_dir = Sys.getenv("PROJECT_ROOT", unset = getwd()), schema = infer_schema_from_df(summary))
)
if (isTRUE(result$contract_check$passed)) {
  outputs$h5ad_file <- build_output_entry(result$h5ad_path, "h5ad", "90a_export_h5ad_scrna", "scRNA H5AD mirror", base_dir = Sys.getenv("PROJECT_ROOT", unset = getwd()))
} else {
  outputs$h5ad_skipped_json <- build_output_entry(result$skipped_json, "json", "90a_export_h5ad_scrna", "H5AD export skip reason", base_dir = Sys.getenv("PROJECT_ROOT", unset = getwd()))
}

h5ad_write_manifest_compat(
  manifest_path = target_manifest,
  new_outputs = outputs,
  module_name = "90a_export_h5ad_scrna",
  base_dir = Sys.getenv("PROJECT_ROOT", unset = getwd()),
  inputs = list(upstream_rds = upstream_rds, upstream_manifest = Sys.getenv("H5AD_UPSTREAM_MANIFEST", unset = "")),
  version = Sys.getenv("MODULE_90_EXPORT_H5AD_VERSION", unset = "1.0"),
  depends_on = list(upstream_manifest = Sys.getenv("H5AD_UPSTREAM_MANIFEST", unset = ""))
)

message("90a H5AD export status: ", result$status)
