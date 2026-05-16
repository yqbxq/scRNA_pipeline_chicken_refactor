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

load_required_packages(c("Seurat", "Matrix", "jsonlite"))

cfg <- get_spatial_script_config()
module_name <- "spatial_01_build_objects"
prepare_dirs_spatial(cfg)
set.seed(cfg$random_seed)

combine_rows <- function(primary, extra) {
  out <- primary
  if (!is.null(extra) && nrow(extra) > 0) {
    for (col in colnames(extra)) {
      out[[col]] <- extra[[col]][[1]]
    }
  }
  out
}

find_contract_row <- function(contracts, platform, layout) {
  if (nrow(contracts) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  for (col in c("platform", "bundle_layout")) {
    if (!col %in% colnames(contracts)) {
      contracts[[col]] <- ""
    }
  }
  platform <- tolower(platform)
  layout <- tolower(layout)
  hit <- contracts[tolower(contracts$platform) == platform & tolower(contracts$bundle_layout) == layout, , drop = FALSE]
  if (nrow(hit) == 0) {
    hit <- contracts[tolower(contracts$platform) == platform, , drop = FALSE]
  }
  if (nrow(hit) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  hit[1, , drop = FALSE]
}

samples <- spatial_read_tsv(if (file.exists(cfg$canonical_sample_sheet)) cfg$canonical_sample_sheet else cfg$sample_sheet)
sections <- spatial_read_tsv(cfg$section_sheet)
inventory <- spatial_read_tsv(cfg$spatial_input_inventory_file)
contracts <- spatial_read_tsv(cfg$spatial_intake_contract_file)

if (nrow(samples) == 0) {
  stop(sprintf("sample sheet is missing or empty: %s", cfg$canonical_sample_sheet), call. = FALSE)
}
if (nrow(inventory) == 0) {
  stop(sprintf("spatial input inventory is missing or empty: %s", cfg$spatial_input_inventory_file), call. = FALSE)
}
for (col in c("sample_id", "section_id", "run_spatial", "modality", "bundle_layout", "platform")) {
  if (!col %in% colnames(samples)) {
    samples[[col]] <- ""
  }
}
for (col in c("sample_id", "section_id", "ready", "ready_status", "standardized_outs")) {
  if (!col %in% colnames(inventory)) {
    inventory[[col]] <- ""
  }
}
for (col in c("section_id", "enabled")) {
  if (!col %in% colnames(sections)) {
    sections[[col]] <- ""
  }
}

spatial_samples <- samples[
  tolower(samples$modality) == "spatial" &
    tolower(samples$run_spatial) != "no",
  ,
  drop = FALSE
]
if (nrow(spatial_samples) == 0) {
  stop("no spatial samples are enabled in samples.canonical.tsv", call. = FALSE)
}

diagnostics <- list()
objects <- list()
for (idx in seq_len(nrow(spatial_samples))) {
  sample_row <- spatial_samples[idx, , drop = FALSE]
  sample_id <- spatial_cell(sample_row, "sample_id")
  section_id <- spatial_cell(sample_row, "section_id", sample_id)
  started <- proc.time()[["elapsed"]]
  output_rds <- file.path(cfg$spatial_checkpoint_dir, sprintf("%s_raw.rds", spatial_safe_id(sample_id)))

  diag_row <- data.frame(
    sample_id = sample_id,
    section_id = section_id,
    status = "failed",
    reason = "",
    elapsed_seconds = NA_real_,
    output_rds = output_rds,
    n_spots = NA_integer_,
    stringsAsFactors = FALSE
  )

  section_row <- sections[sections$section_id == section_id, , drop = FALSE]
  if (nrow(section_row) == 0) {
    diag_row$reason <- "section_not_registered"
    diagnostics[[length(diagnostics) + 1L]] <- diag_row
    next
  }
  section_row <- section_row[1, , drop = FALSE]
  if (!spatial_bool(spatial_cell(section_row, "enabled", "yes"))) {
    diag_row$status <- "skipped"
    diag_row$reason <- "section_disabled"
    diagnostics[[length(diagnostics) + 1L]] <- diag_row
    next
  }

  inv_row <- inventory[inventory$sample_id == sample_id | inventory$section_id == section_id, , drop = FALSE]
  if (nrow(inv_row) == 0) {
    diag_row$reason <- "missing_spatial_inventory_row"
    diagnostics[[length(diagnostics) + 1L]] <- diag_row
    next
  }
  inv_row <- inv_row[1, , drop = FALSE]
  if (tolower(spatial_cell(inv_row, "ready")) != "true") {
    diag_row$reason <- sprintf("spatial_inventory_not_ready:%s", spatial_cell(inv_row, "ready_status"))
    diagnostics[[length(diagnostics) + 1L]] <- diag_row
    next
  }

  combined_sample <- combine_rows(sample_row, inv_row)
  contract_row <- find_contract_row(
    contracts,
    spatial_cell(combined_sample, "platform", spatial_cell(section_row, "platform", "visium")),
    spatial_cell(combined_sample, "bundle_layout", "outs_visium")
  )

  result <- tryCatch({
    obj <- load_spatial_object(combined_sample, section_row, contract_row)
    mito_override <- spatial_cell(section_row, "mito_set_override")
    mito_detection <- detect_spatial_mito_features(obj, gtf_path = cfg$clean_gtf, manual_override = mito_override)
    obj <- inject_spatial_mito_qc(obj, mito_detection)
    saveRDS(obj, output_rds)
    obj
  }, error = function(e) e)

  diag_row$elapsed_seconds <- round(proc.time()[["elapsed"]] - started, 3)
  if (inherits(result, "error")) {
    reason <- conditionMessage(result)
    diag_row$status <- if (grepl("deferred to M5", reason, fixed = TRUE)) "pending_m5" else "failed"
    diag_row$reason <- reason
  } else {
    diag_row$status <- "ok"
    diag_row$reason <- ""
    diag_row$n_spots <- ncol(result)
    objects[[sample_id]] <- result
  }
  diagnostics[[length(diagnostics) + 1L]] <- diag_row
}

diagnostics_df <- do.call(rbind, diagnostics)
spatial_write_tsv(diagnostics_df, cfg$load_diagnostics_tsv)

if (length(objects) == 0) {
  stop(sprintf("all spatial sections failed; see %s", cfg$load_diagnostics_tsv), call. = FALSE)
}

load_summary <- summarize_load_metrics(objects)
write.csv(load_summary, cfg$load_summary_csv, row.names = FALSE)

write_manifest_local(
  manifest_path = cfg$module_01_manifest_path,
  new_outputs = list(
    load_summary = build_output_entry(cfg$load_summary_csv, "csv", module_name, "one row per loaded ST section", base_dir = cfg$project_root, schema = infer_schema_from_df(load_summary)),
    load_diagnostics = build_output_entry(cfg$load_diagnostics_tsv, "tsv", module_name, "one row per ST section load attempt", base_dir = cfg$project_root, schema = infer_schema_from_df(diagnostics_df)),
    raw_object_dir = build_output_entry(cfg$spatial_checkpoint_dir, "directory", module_name, "per-sample raw spatial Seurat objects", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    spatial_input_inventory = cfg$spatial_input_inventory_file,
    section_sheet = cfg$section_sheet,
    sample_sheet = cfg$canonical_sample_sheet,
    spatial_intake_contract = cfg$spatial_intake_contract_file
  ),
  version = cfg$module_version
)

message(sprintf("spatial raw object build complete: %s objects", length(objects)))
