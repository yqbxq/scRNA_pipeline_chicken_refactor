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

load_required_packages(c("Seurat", "Matrix", "ggplot2", "jsonlite"))

cfg <- get_spatial_script_config()
module_name <- "spatial_01b_qc_filter"
prepare_dirs_spatial(cfg)
set.seed(cfg$random_seed)

thresholds <- spatial_thresholds_table(cfg)
raw_paths <- list.files(cfg$spatial_checkpoint_dir, pattern = "_raw\\.rds$", full.names = TRUE)
if (length(raw_paths) == 0) {
  stop(sprintf("no raw spatial objects found in %s", cfg$spatial_checkpoint_dir), call. = FALSE)
}

summary_rows <- list()
mask_paths <- character()
post_qc_paths <- character()

for (raw_path in raw_paths) {
  obj <- readRDS(raw_path)
  ids <- spatial_object_ids(obj, fallback = sub("_raw\\.rds$", "", basename(raw_path)))
  stem <- sub("_raw\\.rds$", "", basename(raw_path))
  threshold <- spatial_threshold_for_section(thresholds, ids$section_id, ids$sample_id, cfg)
  filter_result <- apply_qc_filter(obj, threshold)
  keep <- filter_result$keep
  if (!any(keep)) {
    stop(sprintf("QC filter removed all spots for section_id=%s sample_id=%s", ids$section_id, ids$sample_id), call. = FALSE)
  }

  drop_cluster <- list(n_clusters = NA_integer_, max_cluster_frac = NA_real_, summary = "disabled", status = "disabled")
  if (isTRUE(threshold$spatial_aware_filter)) {
    drop_cluster <- detect_spatial_drop_cluster(obj, keep, eps_factor = cfg$spatial_dbscan_eps_factor)
  }

  keep_cells <- names(keep)[keep]
  filt <- subset(obj, cells = keep_cells)
  filt@misc$qc_filter <- list(
    thresholds = threshold,
    n_before = ncol(obj),
    n_after = ncol(filt),
    retention = ncol(filt) / ncol(obj),
    timestamp = as.character(Sys.time()),
    spatial_aware_summary = drop_cluster
  )
  if (isTRUE(threshold$spatial_aware_filter)) {
    filt@misc$spatial_drop_mask <- drop_cluster
  }

  post_qc_path <- file.path(cfg$spatial_checkpoint_dir, sprintf("%s_post_qc.rds", stem))
  saveRDS(filt, post_qc_path)
  post_qc_paths <- c(post_qc_paths, post_qc_path)

  mask_path <- file.path(cfg$spatial_post_qc_mask_dir, sprintf("qc_filter_mask_%s.png", stem))
  write_qc_filter_mask_plot(obj, keep, mask_path)
  mask_paths <- c(mask_paths, mask_path)

  reason_summary <- filter_result$summary[1, , drop = FALSE]
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    sample_id = ids$sample_id,
    section_id = ids$section_id,
    object_stem = stem,
    raw_rds = relative_path_local(raw_path, cfg$project_root),
    post_qc_rds = relative_path_local(post_qc_path, cfg$project_root),
    n_before = ncol(obj),
    n_after = ncol(filt),
    retention = ncol(filt) / ncol(obj),
    n_drop_high_mito = reason_summary$n_drop_high_mito,
    n_drop_low_feature = reason_summary$n_drop_low_feature,
    n_drop_high_feature = reason_summary$n_drop_high_feature,
    n_drop_low_count = reason_summary$n_drop_low_count,
    n_drop_high_count = reason_summary$n_drop_high_count,
    spatial_aware_flag = isTRUE(threshold$spatial_aware_filter),
    spatial_aware_status = as.character(drop_cluster$status %||% ""),
    spatial_aware_summary = as.character(drop_cluster$summary %||% ""),
    spatial_aware_n_clusters = as.integer(drop_cluster$n_clusters %||% NA_integer_),
    spatial_aware_max_cluster_frac = as.numeric(drop_cluster$max_cluster_frac %||% NA_real_),
    mask_png = relative_path_local(mask_path, cfg$project_root),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
spatial_write_tsv(summary_df, cfg$qc_filter_summary_tsv)

write_manifest_local(
  manifest_path = cfg$module_02_filter_manifest_path,
  new_outputs = list(
    qc_filter_summary = build_output_entry(cfg$qc_filter_summary_tsv, "tsv", module_name, "one row per spatial object after spot QC filtering", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df)),
    post_qc_object_dir = build_output_entry(cfg$spatial_checkpoint_dir, "directory", module_name, "per-section post-QC Seurat objects", base_dir = cfg$project_root),
    qc_filter_mask_dir = build_output_entry(cfg$spatial_post_qc_mask_dir, "directory", module_name, "per-section spatial keep/drop mask figures", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    raw_object_dir = cfg$spatial_checkpoint_dir,
    spatial_qc_threshold_file = cfg$spatial_qc_threshold_file
  ),
  version = cfg$module_02_version,
  depends_on = list(
    spatial_01_build_objects = cfg$module_01_manifest_path,
    spatial_01a_pre_qc_eda = cfg$module_01a_manifest_path
  )
)

message(sprintf("spatial spot QC filtering complete: %s sections", nrow(summary_df)))
