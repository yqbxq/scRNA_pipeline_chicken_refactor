env_or_default_03 <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

env_numeric_03 <- function(name, default) {
  value <- suppressWarnings(as.numeric(env_or_default_03(name, as.character(default))))
  if (is.na(value)) default else value
}

env_integer_03 <- function(name, default) {
  value <- suppressWarnings(as.integer(env_or_default_03(name, as.character(default))))
  if (is.na(value)) default else value
}

env_logical_03 <- function(name, default = FALSE) {
  value <- tolower(env_or_default_03(name, ifelse(isTRUE(default), "yes", "no")))
  value %in% c("yes", "true", "1", "on", "approved")
}

get_single_script_config_03 <- function() {
  base <- get_single_script_config_02()
  layer_id <- env_or_default_03("PANORAMA_LAYER_ID", "panorama")

  base$panorama_layer_id <- layer_id
  base$module_00c_manifest_path <- base$ortholog_manifest_path
  base$object_layer_config_file <- env_or_default_03("OBJECT_LAYER_CONFIG_FILE", file.path(base$config_dir, "object_layers.tsv"))
  base$comparison_sheet <- env_or_default_03("COMPARISON_SHEET", file.path(base$metadata_dir, "comparisons.tsv"))
  base$annotation_marker_targets_sheet <- env_or_default_03("ANNOTATION_MARKER_TARGETS_SHEET", file.path(base$metadata_dir, "annotation_marker_targets.tsv"))
  base$marker_panel_dir <- env_or_default_03("MARKER_PANEL_DIR", file.path(base$config_dir, "marker_panels"))
  base$cluster_selection_file <- env_or_default_03("CLUSTER_SELECTION_FILE", file.path(base$config_dir, "cluster_selection.tsv"))

  base$module_03a1_manifest_path <- file.path(base$manifest_dir, "03a1_normalize_hvg", "_manifest.json")
  base$module_03a2_manifest_path <- file.path(base$manifest_dir, "03a2_reduce_integrate", "_manifest.json")
  base$module_03b_manifest_path <- file.path(base$manifest_dir, "03b_integration_eda", "_manifest.json")
  base$module_03c_manifest_path <- file.path(base$manifest_dir, "03c_cluster", "_manifest.json")
  base$module_03c2_manifest_path <- file.path(base$manifest_dir, "03c2_cluster_marker_audit", "_manifest.json")
  base$module_03c3_manifest_path <- file.path(base$manifest_dir, "03c3_cluster_selection_report", "_manifest.json")
  base$module_03d_manifest_path <- file.path(base$manifest_dir, "03d_annotate", "_manifest.json")
  base$module_03e_manifest_path <- file.path(base$manifest_dir, "03e_annotation_eda", "_manifest.json")
  base$module_03d_marker_risk_manifest_path <- file.path(base$manifest_dir, "03d_marker_risk", "_manifest.json")
  base$module_03e_panel_evidence_manifest_path <- file.path(base$manifest_dir, "03e_panel_evidence", "_manifest.json")
  base$module_03f_manifest_path <- file.path(base$manifest_dir, "03f_apply_manual_annotation", "_manifest.json")
  base$module_03g_manifest_path <- file.path(base$manifest_dir, "03g_annotation_eda", "_manifest.json")

  base$panorama_checkpoint_dir <- file.path(base$checkpoint_dir, "layers", layer_id)
  base$panorama_normalized_dir <- file.path(base$panorama_checkpoint_dir, "normalized")
  base$panorama_reduction_dir <- file.path(base$panorama_checkpoint_dir, "reductions")
  base$panorama_candidate_clustered_rds <- file.path(base$panorama_checkpoint_dir, sprintf("%s_cluster_candidates.rds", layer_id))
  base$panorama_clustered_rds <- file.path(base$panorama_checkpoint_dir, sprintf("%s_after_clustering.rds", layer_id))
  base$panorama_annotated_rds <- file.path(base$panorama_checkpoint_dir, sprintf("%s_after_annotation.rds", layer_id))

  base$compat_clustered_rds <- env_or_default_03("ANNOTATION_HUB_PATH_CLUSTERED", file.path(base$checkpoint_dir, "02_after_clustering.rds"))
  base$compat_annotated_rds <- env_or_default_03("ANNOTATION_HUB_PATH", file.path(base$checkpoint_dir, "03_after_annotation.rds"))

  integration_root <- env_or_default_03("INTEGRATION_REPORT_DIR", file.path(base$eda_report_dir, "integration"))
  annotation_root <- env_or_default_03("ANNOTATION_REPORT_DIR", file.path(base$eda_report_dir, "annotation"))
  base$integration_report_dir_layer <- file.path(integration_root, layer_id)
  base$annotation_report_dir_layer <- file.path(annotation_root, "layers", layer_id)
  base$annotation_table_dir_layer <- file.path(base$table_dir, "annotation", "layers", layer_id)
  base$integration_table_dir_layer <- file.path(base$table_dir, "integration", layer_id)
  base$cluster_marker_table_dir_layer <- env_or_default_03(
    "CLUSTER_MARKER_TABLE_DIR",
    file.path(base$table_dir, "cluster_marker_audit", "layers", layer_id)
  )
  base$cluster_marker_report_dir_layer <- env_or_default_03(
    "CLUSTER_MARKER_REPORT_DIR",
    file.path(base$eda_report_dir, "cluster_marker_audit", "layers", layer_id)
  )
  base$cluster_selection_report_dir_layer <- env_or_default_03(
    "CLUSTER_SELECTION_REPORT_DIR",
    file.path(base$eda_report_dir, "clustering", "layers", layer_id)
  )
  base$marker_risk_table_dir_layer <- env_or_default_03(
    "MARKER_RISK_TABLE_DIR",
    file.path(base$table_dir, "marker_risk", "layers", layer_id)
  )
  base$marker_risk_report_dir_layer <- env_or_default_03(
    "MARKER_RISK_REPORT_DIR",
    file.path(base$eda_report_dir, "marker_risk", "layers", layer_id)
  )
  base$panel_evidence_table_dir_layer <- env_or_default_03(
    "PANEL_EVIDENCE_TABLE_DIR",
    file.path(base$table_dir, "panel_evidence", "layers", layer_id)
  )
  base$panel_evidence_report_dir_layer <- env_or_default_03(
    "PANEL_EVIDENCE_REPORT_DIR",
    file.path(base$eda_report_dir, "panel_evidence", "layers", layer_id)
  )
  base$manual_annotation_file <- env_or_default_03("MANUAL_ANNOTATION_FILE", file.path(base$config_dir, "manual_annotation.tsv"))

  base$layer_status_file <- env_or_default_03("LAYER_STATUS_FILE", file.path(base$table_dir, "layer_status.tsv"))
  base$selected_integration_file <- env_or_default_03("SELECTED_INTEGRATION_FILE", file.path(base$integration_report_dir_layer, "selected_integration.txt"))

  base$normalization_methods_default <- env_or_default_03("NORMALIZATION_METHODS", "lognorm")
  base$integration_modes_default <- env_or_default_03("INTEGRATION_MODES", env_or_default_03("INTEGRATION_MODE", "harmony"))
  base$vars_to_regress_default <- env_or_default_03("VARS_TO_REGRESS_DEFAULT", "")
  base$res_fine_step_default <- env_numeric_03("RES_FINE_STEP", 0.005)
  base$pca_dims_panorama_raw <- env_or_default_03("PCA_DIMS_PANORAMA", env_or_default_03("PCA_DIMS", "1:30"))
  base$pca_dims_subcluster_raw <- env_or_default_03("PCA_DIMS_SUBCLUSTER", "1:20")
  base$cluster_seed_repeats <- env_integer_03("CLUSTER_SEED_REPEATS", 3L)
  base$cluster_subsample_repeats <- env_integer_03("CLUSTER_SUBSAMPLE_REPEATS", 2L)
  base$cluster_subsample_fraction <- env_numeric_03("CLUSTER_SUBSAMPLE_FRACTION", 0.80)
  base$cluster_metric_max_cells <- env_integer_03("CLUSTER_METRIC_MAX_CELLS", 2000L)
  base$cluster_candidate_top_n <- env_integer_03("CLUSTER_CANDIDATE_TOP_N", 5L)
  base$cluster_marker_candidate_top_n <- env_integer_03("CLUSTER_MARKER_CANDIDATE_TOP_N", 3L)
  base$cluster_marker_assay <- env_or_default_03("CLUSTER_MARKER_ASSAY", "RNA")
  base$cluster_marker_min_pct <- env_numeric_03("CLUSTER_MARKER_MIN_PCT", 0.10)
  base$cluster_marker_logfc_threshold <- env_numeric_03("CLUSTER_MARKER_LOGFC_THRESHOLD", 0.10)
  base$cluster_marker_anno_padj <- env_numeric_03("CLUSTER_MARKER_ANNO_PADJ", 0.05)
  base$cluster_marker_anno_logfc <- env_numeric_03("CLUSTER_MARKER_ANNO_LOGFC", 0.25)
  base$cluster_marker_anno_pct1 <- env_numeric_03("CLUSTER_MARKER_ANNO_PCT1", 0.25)
  base$cluster_marker_anno_pct_diff <- env_numeric_03("CLUSTER_MARKER_ANNO_PCT_DIFF", 0.10)
  base$cluster_marker_strict_padj <- env_numeric_03("CLUSTER_MARKER_STRICT_PADJ", 0.01)
  base$cluster_marker_strict_logfc <- env_numeric_03("CLUSTER_MARKER_STRICT_LOGFC", 0.50)
  base$cluster_marker_strict_pct1 <- env_numeric_03("CLUSTER_MARKER_STRICT_PCT1", 0.30)
  base$cluster_marker_strict_pct_diff <- env_numeric_03("CLUSTER_MARKER_STRICT_PCT_DIFF", 0.20)
  base$cluster_marker_min_strict_markers <- env_integer_03("CLUSTER_MARKER_MIN_STRICT_MARKERS", 5L)
  base$panel_evidence_mode <- tolower(env_or_default_03("PANEL_EVIDENCE_MODE", "off"))
  base$annotation_mode <- tolower(env_or_default_03("ANNOTATION_MODE", "manual"))
  base$auto_apply_annotation <- env_logical_03("AUTO_APPLY_ANNOTATION", FALSE)
  base$min_biological_replicates <- env_integer_03("MIN_BIOLOGICAL_REPLICATES", 2L)
  base$module_version <- env_or_default_03("MODULE_03_VERSION", "1.0")
  base
}

prepare_dirs_03 <- function(cfg) {
  ensure_dirs <- c(
    cfg$panorama_checkpoint_dir,
    cfg$panorama_normalized_dir,
    cfg$panorama_reduction_dir,
    cfg$integration_report_dir_layer,
    cfg$annotation_report_dir_layer,
    cfg$annotation_table_dir_layer,
    cfg$integration_table_dir_layer,
    cfg$cluster_marker_table_dir_layer,
    cfg$cluster_marker_report_dir_layer,
    cfg$cluster_selection_report_dir_layer,
    cfg$marker_risk_table_dir_layer,
    cfg$marker_risk_report_dir_layer,
    cfg$panel_evidence_table_dir_layer,
    cfg$panel_evidence_report_dir_layer,
    dirname(cfg$layer_status_file),
    dirname(cfg$selected_integration_file),
    dirname(cfg$manual_annotation_file),
    dirname(cfg$cluster_selection_file),
    dirname(cfg$module_03a1_manifest_path),
    dirname(cfg$module_03a2_manifest_path),
    dirname(cfg$module_03b_manifest_path),
    dirname(cfg$module_03c_manifest_path),
    dirname(cfg$module_03c2_manifest_path),
    dirname(cfg$module_03c3_manifest_path),
    dirname(cfg$module_03d_manifest_path),
    dirname(cfg$module_03e_manifest_path),
    dirname(cfg$module_03d_marker_risk_manifest_path),
    dirname(cfg$module_03e_panel_evidence_manifest_path),
    dirname(cfg$module_03f_manifest_path),
    dirname(cfg$module_03g_manifest_path)
  )
  invisible(lapply(ensure_dirs[nzchar(ensure_dirs)], ensure_dir))
}
