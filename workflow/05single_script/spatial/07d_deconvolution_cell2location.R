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

c2l_latest_h5ad <- function(cfg, module, section = "") {
  explicit <- if (identical(module, "03d_panorama")) Sys.getenv("SPATIAL_C2L_REF_H5AD", unset = "") else Sys.getenv("SPATIAL_C2L_ST_H5AD", unset = "")
  if (nzchar(explicit)) {
    return(normalizePath(explicit, winslash = "/", mustWork = FALSE))
  }
  h5ad_dir <- file.path(cfg$results_dir, "90a_export_h5ad", module)
  if (!dir.exists(h5ad_dir)) {
    return("")
  }
  files <- list.files(h5ad_dir, pattern = "\\.h5ad$", full.names = TRUE)
  if (length(files) == 0) {
    return("")
  }
  section <- st07_scalar(section, "")
  if (nzchar(section) && !section %in% c("all", "*", "ALL", "__ALL__")) {
    section_hits <- files[grepl(section, basename(files), fixed = TRUE)]
    if (length(section_hits) > 0) {
      files <- section_hits
    }
  }
  files <- files[order(file.info(files)$mtime)]
  normalizePath(files[[length(files)]], winslash = "/", mustWork = FALSE)
}

c2l_output_paths <- function(cfg, pair, section) {
  out_dir <- file.path(cfg$spatial_c2l_table_dir, spatial_safe_id(st07_scalar(pair$deconv_id, "default_deconv")), spatial_safe_id(section))
  list(
    out_dir = out_dir,
    proportion_tsv = file.path(out_dir, "spot_celltype_proportions.tsv"),
    proportion_wide_tsv = file.path(out_dir, "spot_celltype_proportions_wide.tsv"),
    spot_metadata_tsv = file.path(out_dir, "spot_metadata.tsv"),
    summary_tsv = file.path(out_dir, "method_summary.tsv"),
    method_object_rds = ""
  )
}

c2l_sidecar_status <- function(code, text) {
  code <- as.integer(code %||% 0L)
  if (identical(code, 0L)) {
    return(list(status = "ok", reason = ""))
  }
  if (grepl("skipped_no_gpu", text, fixed = TRUE) || identical(code, 21L)) {
    return(list(status = "skipped_no_gpu", reason = text))
  }
  if (grepl("skipped_no_python_env", text, fixed = TRUE) || identical(code, 20L) || identical(code, 127L)) {
    return(list(status = "skipped_no_python_env", reason = text))
  }
  list(status = "failed_sidecar", reason = text)
}

c2l_read_summary <- function(path) {
  df <- st07_read_tsv(path)
  if (nrow(df) == 0) {
    return(list(n_spots = 0L, n_celltypes = 0L, runtime_sec = NA_real_))
  }
  list(
    n_spots = suppressWarnings(as.integer(df$n_spots[[1]] %||% 0L)),
    n_celltypes = suppressWarnings(as.integer(df$n_celltypes[[1]] %||% 0L)),
    runtime_sec = suppressWarnings(as.numeric(df$runtime_sec[[1]] %||% NA_real_))
  )
}

pairs <- read_deconv_pairs_st(cfg)
rows <- list()
if (nrow(pairs) == 0) {
  pair <- data.frame(deconv_id = "", enabled = "yes", tool = "cell2location", stringsAsFactors = FALSE)
  outputs <- write_empty_deconv_outputs_st(cfg, cfg$spatial_c2l_table_dir, pair, "all", "cell2location", "skipped_no_deconv_pairs", "deconv_pairs.tsv has no enabled rows")
  rows[[1L]] <- st07_deconv_row(pair, "all", "cell2location", "skipped_no_deconv_pairs", "deconv_pairs.tsv has no enabled rows", outputs = outputs)
} else {
  sidecar <- file.path(cfg$pipeline_root, "workflow", "04python", "cell2location_pipeline.py")
  use_gpu <- tolower(cfg$spatial_c2l_use_gpu) %in% c("yes", "true", "1", "on")
  reference <- load_spatial_reference_inventory_st(cfg)
  ref_h5ad <- c2l_latest_h5ad(cfg, Sys.getenv("SPATIAL_C2L_REF_MODULE", unset = "03d_panorama"))
  for (idx in seq_len(nrow(pairs))) {
    pair <- pairs[idx, , drop = FALSE]
    section <- st07_scalar(pair$section_filter, "all")
    if (!deconv_pair_enabled_st(pair) || !deconv_pair_allows_tool_st(pair, "cell2location")) {
      outputs <- write_empty_deconv_outputs_st(cfg, cfg$spatial_c2l_table_dir, pair, section, "cell2location", "skipped_disabled", "pair disabled or tool filter does not include cell2location")
      rows[[length(rows) + 1L]] <- st07_deconv_row(pair, section, "cell2location", "skipped_disabled", "pair disabled or tool filter does not include cell2location", outputs = outputs)
      next
    }
    if (!identical(reference$status, "ok")) {
      outputs <- write_empty_deconv_outputs_st(cfg, cfg$spatial_c2l_table_dir, pair, section, "cell2location", reference$status, reference$reason)
      rows[[length(rows) + 1L]] <- st07_deconv_row(pair, section, "cell2location", reference$status, reference$reason, outputs = outputs)
      next
    }
    outputs <- c2l_output_paths(cfg, pair, section)
    ensure_dir(outputs$out_dir)
    st_h5ad <- c2l_latest_h5ad(cfg, Sys.getenv("SPATIAL_C2L_ST_MODULE", unset = "spatial_03_region"), section)
    args <- c(
      sidecar,
      "--annotation-col", reference$annotation_col,
      "--out-dir", outputs$out_dir,
      "--ref-epochs", as.character(cfg$spatial_c2l_ref_epochs),
      "--st-epochs", as.character(cfg$spatial_c2l_st_epochs),
      "--detection-alpha", as.character(cfg$spatial_c2l_detection_alpha),
      if (use_gpu) "--use-gpu" else character(0)
    )
    if (nzchar(ref_h5ad)) {
      args <- c(args, "--ref-h5ad", ref_h5ad)
    }
    if (nzchar(st_h5ad)) {
      args <- c(args, "--st-h5ad", st_h5ad)
    }
    if (nzchar(section) && !section %in% c("all", "*", "ALL", "__ALL__")) {
      args <- c(args, "--section-id", section)
    }
    timing <- system.time(run <- tryCatch(system2(cfg$py_cell2location_bin, args = args, stdout = TRUE, stderr = TRUE), error = function(e) structure(conditionMessage(e), status = 127)))
    code <- attr(run, "status")
    if (is.null(code)) code <- 0L
    text <- paste(as.character(run), collapse = " ")
    outcome <- c2l_sidecar_status(code, text)
    if (!identical(outcome$status, "ok") || !file.exists(outputs$proportion_tsv) || file.info(outputs$proportion_tsv)$size == 0) {
      empty <- write_empty_deconv_outputs_st(cfg, cfg$spatial_c2l_table_dir, pair, section, "cell2location", outcome$status, outcome$reason)
      rows[[length(rows) + 1L]] <- st07_deconv_row(pair, section, "cell2location", outcome$status, outcome$reason, outputs = empty)
      next
    }
    summary <- c2l_read_summary(outputs$summary_tsv)
    outputs$out_dir <- NULL
    rows[[length(rows) + 1L]] <- st07_deconv_row(pair, section, "cell2location", "ok", "", n_spots = summary$n_spots, n_celltypes = summary$n_celltypes, runtime_sec = unname(timing[["elapsed"]]), outputs = outputs)
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
