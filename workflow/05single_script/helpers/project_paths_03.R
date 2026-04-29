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

get_single_script_config_03 <- function() {
  base <- get_single_script_config_02()
  layer_id <- env_or_default_03("PANORAMA_LAYER_ID", "panorama")

  base$panorama_layer_id <- layer_id
  base$module_00c_manifest_path <- base$ortholog_manifest_path
  base$object_layer_config_file <- env_or_default_03("OBJECT_LAYER_CONFIG_FILE", file.path(base$config_dir, "object_layers.tsv"))
  base$comparison_sheet <- env_or_default_03("COMPARISON_SHEET", file.path(base$metadata_dir, "comparisons.tsv"))
  base$marker_panel_dir <- env_or_default_03("MARKER_PANEL_DIR", file.path(base$config_dir, "marker_panels"))

  base$module_03a1_manifest_path <- file.path(base$manifest_dir, "03a1_normalize_hvg", "_manifest.json")
  base$module_03a2_manifest_path <- file.path(base$manifest_dir, "03a2_reduce_integrate", "_manifest.json")
  base$module_03b_manifest_path <- file.path(base$manifest_dir, "03b_integration_eda", "_manifest.json")
  base$module_03c_manifest_path <- file.path(base$manifest_dir, "03c_cluster", "_manifest.json")
  base$module_03d_manifest_path <- file.path(base$manifest_dir, "03d_annotate", "_manifest.json")
  base$module_03e_manifest_path <- file.path(base$manifest_dir, "03e_annotation_eda", "_manifest.json")

  base$panorama_checkpoint_dir <- file.path(base$checkpoint_dir, "layers", layer_id)
  base$panorama_normalized_dir <- file.path(base$panorama_checkpoint_dir, "normalized")
  base$panorama_reduction_dir <- file.path(base$panorama_checkpoint_dir, "reductions")
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

  base$layer_status_file <- env_or_default_03("LAYER_STATUS_FILE", file.path(base$table_dir, "layer_status.tsv"))
  base$selected_integration_file <- env_or_default_03("SELECTED_INTEGRATION_FILE", file.path(base$integration_report_dir_layer, "selected_integration.txt"))

  base$normalization_methods_default <- env_or_default_03("NORMALIZATION_METHODS", "lognorm")
  base$integration_modes_default <- env_or_default_03("INTEGRATION_MODES", env_or_default_03("INTEGRATION_MODE", "harmony"))
  base$vars_to_regress_default <- env_or_default_03("VARS_TO_REGRESS_DEFAULT", "")
  base$res_fine_step_default <- env_numeric_03("RES_FINE_STEP", 0.005)
  base$pca_dims_panorama_raw <- env_or_default_03("PCA_DIMS_PANORAMA", env_or_default_03("PCA_DIMS", "1:30"))
  base$pca_dims_subcluster_raw <- env_or_default_03("PCA_DIMS_SUBCLUSTER", "1:20")
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
    dirname(cfg$layer_status_file),
    dirname(cfg$selected_integration_file),
    dirname(cfg$module_03a1_manifest_path),
    dirname(cfg$module_03a2_manifest_path),
    dirname(cfg$module_03b_manifest_path),
    dirname(cfg$module_03c_manifest_path),
    dirname(cfg$module_03d_manifest_path),
    dirname(cfg$module_03e_manifest_path)
  )
  invisible(lapply(ensure_dirs[nzchar(ensure_dirs)], ensure_dir))
}
