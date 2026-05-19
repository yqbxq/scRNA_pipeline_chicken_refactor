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
source(file.path(.script_dir, "helpers", "spatial_de_helpers.R"), encoding = "UTF-8")

load_required_packages(c("Seurat", "Matrix", "ggplot2", "jsonlite"))

cfg <- get_spatial_script_config()
module_name <- "spatial_05_region_marker"
prepare_dirs_spatial(cfg)

comparisons <- read_spatial_comparisons_05(cfg)
manifest_rows <- list()
panorama_input <- ""

if (nrow(comparisons) == 0) {
  manifest_df <- spatial_manifest_row_05(
    vars = list(comparison_id = "", analysis_mode = "", analysis_unit = "", gene_program_role = "", group_var = "", subset_column = "", subset_value = "", force_exploratory = FALSE, min_cells_per_group = 0L, min_biological_replicates = 0L, logfc_threshold = cfg$spatial_marker_logfc_threshold),
    layer_id = "", group_by = "region", status = "skipped_no_spatial_comparisons", reason = "comparisons.tsv has no enabled analysis_modality=spatial rows", cfg = cfg
  )
} else {
  loaded <- load_panorama_for_05(cfg)
  panorama <- loaded$object
  panorama_input <- loaded$path
  include_undetermined <- tolower(cfg$spatial_de_include_undetermined) %in% c("yes", "true", "1", "on")
  if (!"region" %in% colnames(panorama@meta.data)) {
    stop("spatial panorama is missing required region metadata", call. = FALSE)
  }

  for (idx in seq_len(nrow(comparisons))) {
    vars <- spatial_comparison_vars_05(comparisons[idx, , drop = FALSE], cfg)
    subset_res <- subset_panorama_for_comparison_05(panorama, vars)
    for (layer_id in spatial_layer_ids_for_comparison(vars)) {
      for (group_by in spatial_group_bys_05(panorama)) {
        out_dir <- file.path(cfg$spatial_marker_table_dir, group_by, spatial_safe_id(vars$comparison_id), spatial_safe_id(layer_id))
        fig_dir <- file.path(cfg$spatial_marker_figure_dir, group_by, spatial_safe_id(vars$comparison_id), spatial_safe_id(layer_id))
        if (!identical(subset_res$status, "ok")) {
          summary <- data.frame(group_value = character(), n_spots = integer(), status = subset_res$status, reason = subset_res$reason, stringsAsFactors = FALSE)
          paths <- write_spatial_marker_outputs(out_dir, fig_dir, data.frame(), summary, vars)
          manifest_rows[[length(manifest_rows) + 1L]] <- spatial_manifest_row_05(vars, layer_id, group_by, inference_status = "skipped", status = subset_res$status, reason = subset_res$reason, exploratory_results_tsv = paths$markers, exploratory_summary_tsv = paths$summary, cfg = cfg)
          next
        }
        run <- compute_region_markers_st(subset_res$object, group_by, vars, include_undetermined)
        summary <- run$summary
        if (nrow(summary) > 0) {
          summary$comparison_id <- vars$comparison_id
          summary$layer_id <- layer_id
          summary$group_by <- group_by
        }
        markers <- run$markers
        if (nrow(markers) > 0) {
          markers$comparison_id <- vars$comparison_id
          markers$layer_id <- layer_id
          markers$group_by <- group_by
          markers$inference_status <- "exploratory_region_marker"
        }
        paths <- write_spatial_marker_outputs(out_dir, fig_dir, markers, summary, vars)
        manifest_rows[[length(manifest_rows) + 1L]] <- spatial_manifest_row_05(vars, layer_id, group_by, inference_status = "exploratory_region_marker", status = run$status, reason = run$reason, exploratory_results_tsv = paths$markers, exploratory_summary_tsv = paths$summary, cfg = cfg)
      }
    }
  }
  manifest_df <- if (length(manifest_rows) > 0) do.call(rbind, manifest_rows) else empty_spatial_manifest_05()
}

manifest_tsv <- file.path(cfg$spatial_marker_table_dir, "marker_discovery_manifest.tsv")
spatial_write_tsv(manifest_df, manifest_tsv)
report_path <- file.path(cfg$spatial_marker_de_report_dir, "region_marker_discovery.md")
summarize_spatial_de_manifest(manifest_df, cfg, report_path)

write_manifest_local(
  manifest_path = cfg$module_05_region_marker_manifest_path,
  new_outputs = list(
    marker_discovery_manifest = build_output_entry(manifest_tsv, "tsv", module_name, "spatial region marker discovery manifest", base_dir = cfg$project_root, schema = infer_schema_from_df(manifest_df)),
    report = build_output_entry(report_path, "md", module_name, "spatial marker discovery summary report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(comparisons = cfg$comparison_sheet, panorama = panorama_input),
  version = cfg$module_05_version,
  depends_on = list(spatial_03a_region_annotation_eda = cfg$module_03a_region_annotation_eda_manifest_path, spatial_04c_subcluster_eda = cfg$module_04c_subcluster_eda_manifest_path)
)

message(sprintf("spatial region marker discovery complete: %d manifest rows", nrow(manifest_df)))
