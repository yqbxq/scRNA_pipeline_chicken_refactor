env_numeric_10 <- function(name, default) {
  value <- suppressWarnings(as.numeric(env_or_default_03(name, as.character(default))))
  if (is.na(value)) default else value
}

env_integer_10 <- function(name, default) {
  value <- suppressWarnings(as.integer(env_or_default_03(name, as.character(default))))
  if (is.na(value)) default else value
}

module_version_10 <- function(cfg, module_name) {
  versions <- cfg$module_versions_10 %||% list()
  normalize_scalar_value(versions[[module_name]], cfg$module_10_version %||% "1.0")
}

get_single_script_config_10 <- function() {
  base <- get_single_script_config_09()
  base$module_10_version <- env_or_default_03("MODULE_10_VERSION", "1.0")
  base$module_versions_10 <- list(
    `10a_run_velocyto` = env_or_default_03("MODULE_10A_VERSION", base$module_10_version),
    `10b_prepare_velocity_reference` = env_or_default_03("MODULE_10B_VERSION", base$module_10_version)
  )

  base$module_10a_manifest_path <- file.path(base$manifest_dir, "10a_run_velocyto", "_manifest.json")
  base$module_10b_manifest_path <- file.path(base$manifest_dir, "10b_prepare_velocity_reference", "_manifest.json")

  base$velocity_dir <- env_or_default_03("VELOCITY_DIR", file.path(base$results_dir, "velocity"))
  base$velocity_input_dir <- env_or_default_03("VELOCITY_INPUT_DIR", file.path(base$velocity_dir, "input"))
  base$velocity_loom_dir <- env_or_default_03("VELOCITY_LOOM_DIR", file.path(base$velocity_dir, "loom"))
  base$velocity_output_dir <- env_or_default_03("VELOCITY_OUTPUT_DIR", file.path(base$velocity_dir, "output"))
  base$velocity_table_dir <- file.path(base$table_dir, "velocity")
  base$velocity_input_table_dir <- file.path(base$velocity_table_dir, "inputs")
  base$velocity_barcode_dir <- file.path(base$velocity_input_dir, "barcodes")

  base$velocity_loom_jobs_tsv <- file.path(base$velocity_input_table_dir, "velocity_loom_jobs.tsv")
  base$velocity_loom_index_tsv <- file.path(base$velocity_input_table_dir, "velocity_loom_index.tsv")
  base$velocity_reference_index_tsv <- file.path(base$velocity_input_table_dir, "velocity_reference_index.tsv")

  base$dnbc4tools_out_dir <- env_or_default_03("DNBC4TOOLS_OUT_DIR", file.path(base$data_dir, "dnbc4tools_out"))
  base$velocity_bam_pattern <- env_or_default_03("VELOCITY_BAM_PATTERN", "anno_decon_sorted.bam")
  base$velocity_h5_pattern <- env_or_default_03("VELOCITY_H5_PATTERN", "filtered_feature_bc_matrix.h5")
  base$velocity_gtf <- env_or_default_03("VELOCITY_GTF", env_or_default_03("CLEAN_GTF", env_or_default_03("REFERENCE_GTF", "")))
  base$velocity_sample_ids <- env_or_default_03("VELOCITY_SAMPLE_IDS", "")
  base$velocyto_repeat_mask_gtf <- env_or_default_03("VELOCYTO_REPEAT_MASK_GTF", "")
  base$velocyto_threads <- env_integer_10("VELOCYTO_THREADS", env_integer_10("MAIN_THREADS", 8L))
  base$velocity_min_reference_cells <- env_integer_10("VELOCITY_MIN_REFERENCE_CELLS", 10L)

  base
}

prepare_dirs_10 <- function(cfg) {
  prepare_dirs_09(cfg)
  ensure_dirs <- c(
    cfg$velocity_dir,
    cfg$velocity_input_dir,
    cfg$velocity_loom_dir,
    cfg$velocity_output_dir,
    cfg$velocity_table_dir,
    cfg$velocity_input_table_dir,
    cfg$velocity_barcode_dir,
    dirname(cfg$velocity_loom_jobs_tsv),
    dirname(cfg$velocity_loom_index_tsv),
    dirname(cfg$velocity_reference_index_tsv),
    dirname(cfg$module_10a_manifest_path),
    dirname(cfg$module_10b_manifest_path)
  )
  invisible(lapply(ensure_dirs[nzchar(ensure_dirs)], ensure_dir))
}

velocity_unit_file_id_10 <- function(pair_id, split_value = "") {
  trajectory_unit_file_id_09(pair_id, split_value)
}

velocity_reference_umap_path_10 <- function(cfg, pair_id, split_value = "") {
  file.path(
    cfg$velocity_input_dir,
    sprintf("velocity_umap_%s.csv", velocity_unit_file_id_10(pair_id, split_value))
  )
}

velocity_reference_metadata_path_10 <- function(cfg, pair_id, split_value = "") {
  file.path(
    cfg$velocity_input_dir,
    sprintf("velocity_metadata_%s.csv", velocity_unit_file_id_10(pair_id, split_value))
  )
}

velocity_reference_legacy_umap_path_10 <- function(cfg) {
  file.path(cfg$velocity_input_dir, "velocity_umap.csv")
}

velocity_reference_legacy_metadata_path_10 <- function(cfg) {
  file.path(cfg$velocity_input_dir, "velocity_metadata.csv")
}
