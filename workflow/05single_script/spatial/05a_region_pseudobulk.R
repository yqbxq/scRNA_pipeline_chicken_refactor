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
source(file.path(PIPELINE_ROOT, "workflow/05single_script/helpers/evidence_tier_dispatch_utils.R"), encoding = "UTF-8")

load_required_packages(c("Seurat", "Matrix", "jsonlite"))

cfg <- get_spatial_script_config()
module_name <- "spatial_05a_region_pseudobulk"
prepare_dirs_spatial(cfg)

comparisons <- read_spatial_comparisons_05(cfg)
cluster_eligibility_tsv <- manifest_output_optional_05(cfg$module_04c_subcluster_eda_manifest_path, "cluster_eligibility_tsv")
cluster_eligibility <- read_cluster_eligibility_05(cfg$module_04c_subcluster_eda_manifest_path)
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
        gate <- if (identical(subset_res$status, "ok")) replicate_gate_check_st(obj, group_by, vars$group_var, vars) else list(pass = FALSE, status = "skip", reason = subset_res$reason)
        for (region_value in values) {
          out_dir <- file.path(cfg$spatial_pseudobulk_table_dir, group_by, spatial_safe_id(vars$comparison_id), spatial_safe_id(layer_id), spatial_safe_id(region_value))
          fallback_path <- ""
          parent_region <- if (nzchar(region_value) && identical(subset_res$status, "ok")) spatial_parent_region_for_value(obj@meta.data, group_by, region_value) else ""
          evidence_decision <- evidence_decision_05(cluster_eligibility, layer_id, vars$comparison_id, vars$group_var, cluster_id = region_value)
          if (!identical(subset_res$status, "ok")) {
            status <- subset_res$status
            reason <- subset_res$reason
            formal <- data.frame()
            gate_status <- "skip"
            inference_status <- "skipped"
          } else if (!evidence_allows_formal_05(evidence_decision)) {
            status <- evidence_decision$evidence_tier
            reason <- evidence_decision$reason
            formal <- data.frame()
            gate_status <- "evidence_tier"
            inference_status <- evidence_decision$evidence_tier
          } else if (!gate$pass) {
            status <- "exploratory_replicate_gate_blocked"
            reason <- paste(c(evidence_decision$reason, gate$reason)[nzchar(c(evidence_decision$reason, gate$reason))], collapse = "; ")
            formal <- data.frame()
            gate_status <- "skip"
            inference_status <- status
          } else {
            run <- run_spatial_pseudobulk_de(obj, vars, group_by, region_value)
            status <- run$status
            reason <- run$reason
            formal <- run$result
            if (nrow(formal) > 0) {
              formal$evidence_tier <- evidence_decision$evidence_tier
              formal$recommended_action <- evidence_decision$recommended_action
            }
            gate_status <- "pass"
            inference_status <- ifelse(identical(status, "ok"), "formal_pseudobulk", status)
          }
          result_path <- write_spatial_pseudobulk_outputs(out_dir, formal, status, reason, fallback_path)
          manifest_row <- spatial_manifest_row_05(vars, layer_id, group_by, region_value, parent_region, inference_status = inference_status, replicate_gate_status = gate_status, status = status, reason = reason, formal_results_tsv = ifelse(identical(status, "ok"), result_path, ""), exploratory_summary_tsv = ifelse(identical(status, "ok"), "", result_path), exploratory_fallback = "", cfg = cfg)
          manifest_row$evidence_tier <- evidence_decision$evidence_tier
          manifest_row$recommended_action <- evidence_decision$recommended_action
          manifest_row$evidence_tier_summary <- evidence_decision$evidence_tier_summary
          manifest_row$warning_banner <- evidence_decision$warning_banner
          manifest_row$parent_cluster_if_merged <- evidence_decision$parent_cluster_if_merged
          manifest_row$cluster_eligibility_tsv <- cluster_eligibility_tsv
          manifest_rows[[length(manifest_rows) + 1L]] <- manifest_row
        }
      }
    }
  }
  manifest_df <- if (length(manifest_rows) > 0) do.call(rbind, manifest_rows) else empty_spatial_manifest_05()
}

manifest_tsv <- file.path(cfg$spatial_pseudobulk_table_dir, "pseudobulk_de_manifest.tsv")
manifest_df <- evidence_manifest_columns_05(manifest_df)
missing_eligibility_path <- is.na(manifest_df$cluster_eligibility_tsv) | !nzchar(manifest_df$cluster_eligibility_tsv)
manifest_df$cluster_eligibility_tsv[missing_eligibility_path] <- cluster_eligibility_tsv
spatial_write_tsv(manifest_df, manifest_tsv)
report_path <- file.path(cfg$spatial_marker_de_report_dir, "region_pseudobulk_de.md")
summarize_spatial_de_manifest(manifest_df, cfg, report_path)

write_manifest_local(
  manifest_path = cfg$module_05a_region_pseudobulk_manifest_path,
  new_outputs = list(
    pseudobulk_de_manifest = build_output_entry(manifest_tsv, "tsv", module_name, "spatial region pseudobulk DE manifest", base_dir = cfg$project_root, schema = infer_schema_from_df(manifest_df)),
    report = build_output_entry(report_path, "md", module_name, "spatial pseudobulk DE summary report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(comparisons = cfg$comparison_sheet, panorama = panorama_input, cluster_eligibility_tsv = cluster_eligibility_tsv),
  version = cfg$module_05_version,
  depends_on = list(spatial_05_region_marker = cfg$module_05_region_marker_manifest_path, spatial_04c_subcluster_eda = cfg$module_04c_subcluster_eda_manifest_path)
)

message(sprintf("spatial region pseudobulk DE complete: %d manifest rows", nrow(manifest_df)))
