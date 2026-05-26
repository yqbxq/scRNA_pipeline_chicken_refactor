#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(utils)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
.script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else ""
.script_dir <- if (nzchar(.script_path)) dirname(normalizePath(.script_path, winslash = "/", mustWork = FALSE)) else getwd()
.root_script_dir <- normalizePath(file.path(.script_dir, ".."), winslash = "/", mustWork = TRUE)

source_utf8 <- function(path) source(path, encoding = "UTF-8")
source_utf8(file.path(.root_script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.root_script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.root_script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.root_script_dir, "helpers", "h5ad_contract_utils.R"))
source_utf8(file.path(.root_script_dir, "helpers", "h5ad_export_utils.R"))

upstream_rds <- Sys.getenv("H5AD_UPSTREAM_RDS", unset = "")
target_dir <- Sys.getenv("H5AD_TARGET_DIR", unset = "")
target_manifest <- Sys.getenv("H5AD_TARGET_MANIFEST", unset = file.path(target_dir, "_manifest.json"))
source_module <- Sys.getenv("H5AD_SOURCE_MODULE", unset = "upstream")
output_key <- Sys.getenv("H5AD_OUTPUT_KEY", unset = "object")

if (!nzchar(upstream_rds)) stop("H5AD_UPSTREAM_RDS is required", call. = FALSE)
if (!nzchar(target_dir)) stop("H5AD_TARGET_DIR is required", call. = FALSE)

ensure_dir(target_dir)
obj <- h5ad_read_upstream_object(upstream_rds)
contract <- load_h5ad_contract(Sys.getenv("H5AD_CONTRACT_FILE", unset = "metadata/h5ad_export_contract.tsv"), modality = "spatial")
fail_on <- h5ad_env_or_default("H5AD_EXPORT_CONTRACT_FAIL_ON", h5ad_env_or_default("H5AD_CONTRACT_FAIL_ON", "warn"))
per_section <- h5ad_env_or_default("H5AD_EXPORT_SPATIAL_PER_SECTION", "yes")

section_ids <- ""
if (h5ad_is_mock_object(obj) && !is.null(obj$obs$section_id)) {
  section_ids <- unique(as.character(obj$obs$section_id))
}
if (!identical(per_section, "yes") || length(section_ids) == 0L || !nzchar(section_ids[[1]])) {
  section_ids <- ""
}

rows <- list()
outputs <- list()
for (section_id in section_ids) {
  section_obj <- obj
  suffix <- if (nzchar(section_id)) paste0("_", section_id) else ""
  h5ad_path <- file.path(target_dir, sprintf("%s_%s_spatial%s.h5ad", source_module, output_key, suffix))
  result <- h5ad_export_one(section_obj, h5ad_path, "spatial", contract, fail_on = if (identical(fail_on, "error")) "warn" else fail_on)
  if (!isTRUE(result$contract_check$passed)) {
    h5ad_write_skip_json(result$skipped_json, result)
  }
  rows[[length(rows) + 1L]] <- h5ad_summary_row("spatial", upstream_rds, result, section_id = section_id)
  key <- if (nzchar(section_id)) paste0("h5ad_file_", section_id) else "h5ad_file"
  skip_key <- if (nzchar(section_id)) paste0("h5ad_skipped_json_", section_id) else "h5ad_skipped_json"
  if (isTRUE(result$contract_check$passed)) {
    outputs[[key]] <- build_output_entry(result$h5ad_path, "h5ad", "90b_export_h5ad_spatial", "spatial H5AD mirror", base_dir = Sys.getenv("PROJECT_ROOT", unset = getwd()))
  } else {
    outputs[[skip_key]] <- build_output_entry(result$skipped_json, "json", "90b_export_h5ad_spatial", "H5AD export skip reason", base_dir = Sys.getenv("PROJECT_ROOT", unset = getwd()))
  }
}

summary <- do.call(rbind, rows)
summary_tsv <- file.path(target_dir, "h5ad_export_summary.tsv")
write_tsv_local(summary, summary_tsv)
outputs$export_summary_tsv <- build_output_entry(summary_tsv, "tsv", "90b_export_h5ad_spatial", "H5AD export status summary", base_dir = Sys.getenv("PROJECT_ROOT", unset = getwd()), schema = infer_schema_from_df(summary))

h5ad_write_manifest_compat(
  manifest_path = target_manifest,
  new_outputs = outputs,
  module_name = "90b_export_h5ad_spatial",
  base_dir = Sys.getenv("PROJECT_ROOT", unset = getwd()),
  inputs = list(upstream_rds = upstream_rds, upstream_manifest = Sys.getenv("H5AD_UPSTREAM_MANIFEST", unset = "")),
  version = Sys.getenv("MODULE_90_EXPORT_H5AD_VERSION", unset = "1.0"),
  depends_on = list(upstream_manifest = Sys.getenv("H5AD_UPSTREAM_MANIFEST", unset = ""))
)

message("90b H5AD export rows: ", nrow(summary))
