env_integer_04 <- function(name, default) {
  value <- suppressWarnings(as.integer(env_or_default_03(name, as.character(default))))
  if (is.na(value)) default else value
}

get_single_script_config_04 <- function() {
  base <- get_single_script_config_03()

  base$module_04a_manifest_path <- file.path(base$manifest_dir, "04a_subcluster_build", "_manifest.json")
  base$module_04a_review_manifest_path <- file.path(base$manifest_dir, "04a_review", "_manifest.json")
  base$module_04b_manifest_path <- file.path(base$manifest_dir, "04b_subcluster_annotate", "_manifest.json")
  base$module_04c_manifest_path <- file.path(base$manifest_dir, "04c_subcluster_eda", "_manifest.json")
  base$module_04d_manifest_path <- file.path(base$manifest_dir, "04d_cluster_robustness", "_manifest.json")

  base$subcluster_report_dir <- env_or_default_03("SUBCLUSTER_REPORT_DIR", file.path(base$eda_report_dir, "subcluster"))
  base$subcluster_table_dir <- file.path(base$table_dir, "subcluster")
  base$subcluster_review_summary_file <- env_or_default_03(
    "SUBCLUSTER_REVIEW_SUMMARY_FILE",
    file.path(base$subcluster_table_dir, "subcluster_review_summary.tsv")
  )
  base$cell_count_threshold_file <- env_or_default_03(
    "CELL_COUNT_INVENTORY_THRESHOLD_OVERRIDE_TSV",
    file.path(base$metadata_dir, "cell_count_thresholds.tsv")
  )
  base$subcluster_candidate_layers_count_file <- env_or_default_03(
    "SUBCLUSTER_CANDIDATE_LAYERS_COUNT_FILE",
    file.path(base$subcluster_table_dir, "candidate_layers_count.txt")
  )
  base$max_integration_candidates_per_layer <- env_integer_04("MAX_INTEGRATION_CANDIDATES_PER_LAYER", 4L)
  base$module_version <- env_or_default_03("MODULE_04_VERSION", "1.0")
  base
}

prepare_dirs_04 <- function(cfg) {
  ensure_dirs <- c(
    cfg$subcluster_report_dir,
    cfg$subcluster_table_dir,
    dirname(cfg$subcluster_review_summary_file),
    dirname(cfg$subcluster_candidate_layers_count_file),
    dirname(cfg$module_04a_manifest_path),
    dirname(cfg$module_04a_review_manifest_path),
    dirname(cfg$module_04b_manifest_path),
    dirname(cfg$module_04c_manifest_path),
    dirname(cfg$module_04d_manifest_path)
  )
  invisible(lapply(ensure_dirs[nzchar(ensure_dirs)], ensure_dir))
}

layer_checkpoint_dir_04 <- function(cfg, layer_id) {
  file.path(cfg$checkpoint_dir, "layers", layer_id)
}

layer_integration_report_dir_04 <- function(cfg, layer_id) {
  file.path(env_or_default_03("INTEGRATION_REPORT_DIR", file.path(cfg$eda_report_dir, "integration")), layer_id)
}

layer_integration_table_dir_04 <- function(cfg, layer_id) {
  file.path(cfg$table_dir, "integration", layer_id)
}

layer_subcluster_report_dir_04 <- function(cfg, layer_id) {
  file.path(cfg$subcluster_report_dir, "layers", layer_id)
}

layer_subcluster_table_dir_04 <- function(cfg, layer_id) {
  file.path(cfg$subcluster_table_dir, "layers", layer_id)
}

layer_annotation_report_dir_04 <- function(cfg, layer_id) {
  file.path(env_or_default_03("ANNOTATION_REPORT_DIR", file.path(cfg$eda_report_dir, "annotation")), "layers", layer_id)
}

layer_annotation_table_dir_04 <- function(cfg, layer_id) {
  file.path(cfg$table_dir, "annotation", "layers", layer_id)
}

selected_integration_file_04 <- function(cfg, layer_id) {
  file.path(layer_integration_report_dir_04(cfg, layer_id), "selected_integration.txt")
}
