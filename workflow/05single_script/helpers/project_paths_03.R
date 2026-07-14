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
  base$cluster_resolution_decision_file <- env_or_default_03("CLUSTER_RESOLUTION_DECISION_FILE", file.path(base$config_dir, "cluster_resolution_decisions.tsv"))
  base$comparison_sheet <- env_or_default_03("COMPARISON_SHEET", file.path(base$metadata_dir, "comparisons.tsv"))
  base$annotation_marker_targets_sheet <- env_or_default_03("ANNOTATION_MARKER_TARGETS_SHEET", file.path(base$metadata_dir, "annotation_marker_targets.tsv"))
  base$marker_panel_dir <- env_or_default_03("MARKER_PANEL_DIR", file.path(base$config_dir, "marker_panels"))

  base$module_03a1_manifest_path <- file.path(base$manifest_dir, "03a1_normalize_hvg", "_manifest.json")
  base$module_03a2_manifest_path <- file.path(base$manifest_dir, "03a2_reduce_integrate", "_manifest.json")
  base$module_03a3_manifest_path <- file.path(base$manifest_dir, "03a3_compute_umap", "_manifest.json")
  base$module_03b_manifest_path <- file.path(base$manifest_dir, "03b_integration_eda", "_manifest.json")
  base$module_03c_manifest_path <- file.path(base$manifest_dir, "03c_cluster", "_manifest.json")
  base$module_03d_manifest_path <- file.path(base$manifest_dir, "03d_annotate", "_manifest.json")
  base$module_03e_manifest_path <- file.path(base$manifest_dir, "03e_annotation_eda", "_manifest.json")

  base$panorama_checkpoint_dir <- file.path(base$checkpoint_dir, "layers", layer_id)
  base$panorama_normalized_dir <- file.path(base$panorama_checkpoint_dir, "normalized")
  base$panorama_reduction_dir <- file.path(base$panorama_checkpoint_dir, "reductions")
  base$panorama_umap_dir <- file.path(base$panorama_checkpoint_dir, "umap")
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
  cluster_resolution_raw <- env_or_default_03("CLUSTER_RESOLUTIONS", paste(base$res_range, collapse = ","))
  parsed_cluster_resolutions <- suppressWarnings(as.numeric(split_csv_02(cluster_resolution_raw)))
  parsed_cluster_resolutions <- parsed_cluster_resolutions[is.finite(parsed_cluster_resolutions)]
  if (length(parsed_cluster_resolutions) > 0) {
    base$res_range <- parsed_cluster_resolutions
  }
  cluster_seeds <- suppressWarnings(as.integer(split_csv_02(env_or_default_03("CLUSTER_SEEDS", "1,11,21,31,41"))))
  cluster_seeds <- cluster_seeds[is.finite(cluster_seeds)]
  if (length(cluster_seeds) < 2L) {
    cluster_seeds <- c(base$random_seed, base$random_seed + 1L)
  }
  base$cluster_selection_mode <- tolower(env_or_default_03("CLUSTER_SELECTION_MODE", "target_clusters"))
  base$cluster_seeds <- unique(cluster_seeds)
  base$cluster_subsample_reps <- env_integer_03("CLUSTER_SUBSAMPLE_REPS", 5L)
  base$cluster_subsample_fraction <- env_numeric_03("CLUSTER_SUBSAMPLE_FRACTION", 0.80)
  base$cluster_silhouette_max_cells <- env_integer_03("CLUSTER_SILHOUETTE_MAX_CELLS", 5000L)
  base$cluster_marker_shortlist_n <- env_integer_03("CLUSTER_MARKER_SHORTLIST_N", 5L)
  base$cluster_min_expected <- env_integer_03("CLUSTER_MIN_EXPECTED", 5L)
  base$cluster_max_expected <- env_integer_03("CLUSTER_MAX_EXPECTED", 30L)
  base$cluster_min_cells_abs <- env_integer_03("CLUSTER_MIN_CELLS_ABS", 30L)
  base$cluster_min_cell_fraction <- env_numeric_03("CLUSTER_MIN_CELL_FRACTION", 0.001)
  base$cluster_min_seed_ari <- env_numeric_03("CLUSTER_MIN_SEED_ARI", 0.80)
  base$cluster_min_subsample_ari <- env_numeric_03("CLUSTER_MIN_SUBSAMPLE_ARI", 0.65)
  base$cluster_stable_local_ari <- env_numeric_03("CLUSTER_STABLE_LOCAL_ARI", 0.90)
  base$cluster_plateau_min_points <- env_integer_03("CLUSTER_PLATEAU_MIN_POINTS", 3L)
  base$cluster_score_tie_delta <- env_numeric_03("CLUSTER_SCORE_TIE_DELTA", 0.02)
  base$cluster_gap_enabled <- tolower(env_or_default_03("CLUSTER_GAP_ENABLED", "no"))
  base$cluster_validation_depth <- tolower(env_or_default_03("CLUSTER_VALIDATION_DEPTH", "standard"))
  base$res_fine_step_default <- env_numeric_03("RES_FINE_STEP", 0.005)
  base$pca_dims_panorama_raw <- env_or_default_03("PCA_DIMS_PANORAMA", env_or_default_03("PCA_DIMS", "1:30"))
  base$pca_dims_subcluster_raw <- env_or_default_03("PCA_DIMS_SUBCLUSTER", "1:20")
  base$min_biological_replicates <- env_integer_03("MIN_BIOLOGICAL_REPLICATES", 2L)
  base$module_version <- env_or_default_03("MODULE_03_VERSION", "1.0")
  base$module_03a3_version <- env_or_default_03("MODULE_03A3_VERSION", base$module_version)
  base
}

prepare_dirs_03 <- function(cfg) {
  ensure_dirs <- c(
    cfg$panorama_checkpoint_dir,
    cfg$panorama_normalized_dir,
    cfg$panorama_reduction_dir,
    cfg$panorama_umap_dir,
    cfg$integration_report_dir_layer,
    cfg$annotation_report_dir_layer,
    cfg$annotation_table_dir_layer,
    cfg$integration_table_dir_layer,
    dirname(cfg$layer_status_file),
    dirname(cfg$selected_integration_file),
    dirname(cfg$cluster_resolution_decision_file),
    dirname(cfg$module_03a1_manifest_path),
    dirname(cfg$module_03a2_manifest_path),
    dirname(cfg$module_03a3_manifest_path),
    dirname(cfg$module_03b_manifest_path),
    dirname(cfg$module_03c_manifest_path),
    dirname(cfg$module_03d_manifest_path),
    dirname(cfg$module_03e_manifest_path)
  )
  invisible(lapply(ensure_dirs[nzchar(ensure_dirs)], ensure_dir))
}
