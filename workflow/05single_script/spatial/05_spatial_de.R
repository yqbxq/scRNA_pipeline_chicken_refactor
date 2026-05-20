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
module_name <- "spatial_05_spatial_de"
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

  for (idx in seq_len(nrow(comparisons))) {
    vars <- spatial_comparison_vars_05(comparisons[idx, , drop = FALSE], cfg)
    subset_res <- subset_panorama_for_comparison_05(panorama, vars)
    for (layer_id in spatial_layer_ids_for_comparison(vars)) {
      for (group_by in spatial_group_bys_05(panorama)) {
        obj <- subset_res$object
        values <- if (identical(subset_res$status, "ok")) spatial_valid_group_values_05(obj@meta.data, group_by, include_undetermined) else character(0)
        if (length(values) == 0) {
          values <- ""
        }
        for (region_value in values) {
          out_dir <- file.path(cfg$spatial_de_table_dir, group_by, spatial_safe_id(vars$comparison_id), spatial_safe_id(layer_id), spatial_safe_id(region_value))
          fig_dir <- file.path(cfg$spatial_de_figure_dir, group_by, spatial_safe_id(vars$comparison_id), spatial_safe_id(layer_id), spatial_safe_id(region_value))
          parent_region <- if (nzchar(region_value) && identical(subset_res$status, "ok")) spatial_parent_region_for_value(obj@meta.data, group_by, region_value) else ""
          if (!identical(subset_res$status, "ok")) {
            status <- subset_res$status
            reason <- subset_res$reason
            de_res <- data.frame()
          } else {
            run <- run_spotlevel_findmarkers_st(obj, group_by, region_value, vars)
            status <- run$status
            reason <- run$reason
            de_res <- run$result
            if (nrow(de_res) > 0) {
              de_res$comparison_id <- vars$comparison_id
              de_res$layer_id <- layer_id
              de_res$group_by <- group_by
              de_res$region_value <- region_value
              de_res$parent_region <- parent_region
              de_res$ident_1 <- vars$ident_1
              de_res$ident_2 <- vars$ident_2
            }
          }
          result_path <- write_spatial_de_outputs(out_dir, fig_dir, de_res, vars, status, reason, region_value)
          manifest_rows[[length(manifest_rows) + 1L]] <- spatial_manifest_row_05(vars, layer_id, group_by, region_value, parent_region, inference_status = "exploratory_spotlevel", status = status, reason = reason, exploratory_results_tsv = result_path, cfg = cfg)
        }
      }
    }
  }
  manifest_df <- if (length(manifest_rows) > 0) do.call(rbind, manifest_rows) else empty_spatial_manifest_05()
}

manifest_tsv <- file.path(cfg$spatial_de_table_dir, "spatial_de_manifest.tsv")
spatial_write_tsv(manifest_df, manifest_tsv)
report_path <- file.path(cfg$spatial_marker_de_report_dir, "spatial_de_spotlevel.md")
summarize_spatial_de_manifest(manifest_df, cfg, report_path)

write_manifest_local(
  manifest_path = cfg$module_05_spatial_de_manifest_path,
  new_outputs = list(
    spatial_de_manifest = build_output_entry(manifest_tsv, "tsv", module_name, "spatial spot-level DE manifest", base_dir = cfg$project_root, schema = infer_schema_from_df(manifest_df)),
    report = build_output_entry(report_path, "md", module_name, "spatial marker/DE/composition summary report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(comparisons = cfg$comparison_sheet, panorama = panorama_input),
  version = cfg$module_05_version,
  depends_on = list(spatial_05b_region_composition = cfg$module_05b_region_composition_manifest_path)
)

message(sprintf("spatial spot-level DE complete: %d manifest rows", nrow(manifest_df)))
