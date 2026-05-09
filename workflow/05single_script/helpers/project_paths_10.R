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
    `10b_prepare_velocity_reference` = env_or_default_03("MODULE_10B_VERSION", base$module_10_version),
    `10c_scvelo_dynamical` = env_or_default_03("MODULE_10C_VERSION", base$module_10_version),
    `10d_velocyto_steady_state` = env_or_default_03("MODULE_10D_VERSION", base$module_10_version),
    `10e_scvelo_drivers` = env_or_default_03("MODULE_10E_VERSION", base$module_10_version),
    `10f_cellrank_fate` = env_or_default_03("MODULE_10F_VERSION", base$module_10_version),
    `10g_velocity_consistency` = env_or_default_03("MODULE_10G_VERSION", base$module_10_version),
    `10h_velocity_root_terminal` = env_or_default_03("MODULE_10H_VERSION", base$module_10_version),
    `10i_velocity_eda` = env_or_default_03("MODULE_10I_VERSION", base$module_10_version)
  )

  base$module_10a_manifest_path <- file.path(base$manifest_dir, "10a_run_velocyto", "_manifest.json")
  base$module_10b_manifest_path <- file.path(base$manifest_dir, "10b_prepare_velocity_reference", "_manifest.json")
  base$module_10c_manifest_path <- file.path(base$manifest_dir, "10c_scvelo_dynamical", "_manifest.json")
  base$module_10d_manifest_path <- file.path(base$manifest_dir, "10d_velocyto_steady_state", "_manifest.json")
  base$module_10e_manifest_path <- file.path(base$manifest_dir, "10e_scvelo_drivers", "_manifest.json")
  base$module_10f_manifest_path <- file.path(base$manifest_dir, "10f_cellrank_fate", "_manifest.json")
  base$module_10g_manifest_path <- file.path(base$manifest_dir, "10g_velocity_consistency", "_manifest.json")
  base$module_10h_manifest_path <- file.path(base$manifest_dir, "10h_velocity_root_terminal", "_manifest.json")
  base$module_10i_manifest_path <- file.path(base$manifest_dir, "10i_velocity_eda", "_manifest.json")

  base$velocity_dir <- env_or_default_03("VELOCITY_DIR", file.path(base$results_dir, "velocity"))
  base$velocity_input_dir <- env_or_default_03("VELOCITY_INPUT_DIR", file.path(base$velocity_dir, "input"))
  base$velocity_loom_dir <- env_or_default_03("VELOCITY_LOOM_DIR", file.path(base$velocity_dir, "loom"))
  base$velocity_output_dir <- env_or_default_03("VELOCITY_OUTPUT_DIR", file.path(base$velocity_dir, "output"))
  base$velocity_table_dir <- file.path(base$table_dir, "velocity")
  base$velocity_input_table_dir <- file.path(base$velocity_table_dir, "inputs")
  base$velocity_methods_table_dir <- file.path(base$velocity_table_dir, "methods")
  base$velocity_scvelo_table_dir <- file.path(base$velocity_methods_table_dir, "scvelo")
  base$velocity_velocyto_table_dir <- file.path(base$velocity_methods_table_dir, "velocyto")
  base$velocity_driver_table_dir <- file.path(base$velocity_methods_table_dir, "drivers")
  base$velocity_cellrank_table_dir <- file.path(base$velocity_methods_table_dir, "cellrank")
  base$velocity_consistency_table_dir <- file.path(base$velocity_methods_table_dir, "consistency")
  base$velocity_root_terminal_table_dir <- file.path(base$velocity_methods_table_dir, "root_terminal")
  base$velocity_final_table_dir <- file.path(base$velocity_table_dir, "final")
  base$velocity_figure_dir <- file.path(base$figure_dir, "velocity")
  base$velocity_scvelo_figure_dir <- file.path(base$velocity_figure_dir, "scvelo")
  base$velocity_cellrank_figure_dir <- file.path(base$velocity_figure_dir, "cellrank")
  base$velocity_consistency_figure_dir <- file.path(base$velocity_figure_dir, "consistency")
  base$velocity_report_dir <- env_or_default_03("VELOCITY_REPORT_DIR", file.path(base$eda_report_dir, "velocity"))
  base$velocity_barcode_dir <- file.path(base$velocity_input_dir, "barcodes")

  base$velocity_loom_jobs_tsv <- file.path(base$velocity_input_table_dir, "velocity_loom_jobs.tsv")
  base$velocity_loom_index_tsv <- file.path(base$velocity_input_table_dir, "velocity_loom_index.tsv")
  base$velocity_reference_index_tsv <- file.path(base$velocity_input_table_dir, "velocity_reference_index.tsv")
  base$velocity_scvelo_index_tsv <- file.path(base$velocity_scvelo_table_dir, "scvelo_index.tsv")
  base$velocity_scvelo_qc_tsv <- file.path(base$velocity_scvelo_table_dir, "velocity_qc.tsv")
  base$velocity_velocyto_steady_index_tsv <- file.path(base$velocity_velocyto_table_dir, "velocyto_steady_index.tsv")
  base$velocity_driver_index_tsv <- file.path(base$velocity_driver_table_dir, "velocity_driver_index.tsv")
  base$velocity_driver_overlap_tsv <- file.path(base$velocity_driver_table_dir, "velocity_driver_overlap.tsv")
  base$velocity_cellrank_index_tsv <- file.path(base$velocity_cellrank_table_dir, "cellrank_index.tsv")
  base$velocity_consistency_index_tsv <- file.path(base$velocity_consistency_table_dir, "velocity_consistency_index.tsv")
  base$velocity_consistency_summary_tsv <- file.path(base$velocity_consistency_table_dir, "velocity_consistency_summary.tsv")
  base$velocity_split_compare_tsv <- file.path(base$velocity_consistency_table_dir, "velocity_split_compare.tsv")
  base$velocity_root_terminal_index_tsv <- file.path(base$velocity_root_terminal_table_dir, "velocity_root_terminal_index.tsv")
  base$velocity_module_status_tsv <- file.path(base$velocity_final_table_dir, "velocity_module_status.tsv")
  base$velocity_triage_tsv <- file.path(base$velocity_final_table_dir, "velocity_triage.tsv")
  base$velocity_report_md <- file.path(base$velocity_report_dir, "report.md")

  base$dnbc4tools_out_dir <- env_or_default_03("DNBC4TOOLS_OUT_DIR", file.path(base$data_dir, "dnbc4tools_out"))
  base$velocity_bam_pattern <- env_or_default_03("VELOCITY_BAM_PATTERN", "anno_decon_sorted.bam")
  base$velocity_h5_pattern <- env_or_default_03("VELOCITY_H5_PATTERN", "filtered_feature_bc_matrix.h5")
  base$velocity_gtf <- env_or_default_03("VELOCITY_GTF", env_or_default_03("CLEAN_GTF", env_or_default_03("REFERENCE_GTF", "")))
  base$velocity_sample_ids <- env_or_default_03("VELOCITY_SAMPLE_IDS", "")
  base$velocyto_repeat_mask_gtf <- env_or_default_03("VELOCYTO_REPEAT_MASK_GTF", "")
  base$velocyto_threads <- env_integer_10("VELOCYTO_THREADS", env_integer_10("MAIN_THREADS", 8L))
  base$velocity_min_reference_cells <- env_integer_10("VELOCITY_MIN_REFERENCE_CELLS", 10L)
  base$scvelo_min_shared_counts <- env_integer_10("SCVELO_MIN_SHARED_COUNTS", 20L)
  base$scvelo_top_genes <- env_integer_10("SCVELO_TOP_GENES", 2000L)
  base$scvelo_n_pcs <- env_integer_10("SCVELO_N_PCS", 30L)
  base$scvelo_n_neighbors <- env_integer_10("SCVELO_N_NEIGHBORS", 30L)
  base$velocity_driver_top_n <- env_integer_10("VELOCITY_DRIVER_TOP_N", 200L)
  base$cellrank_min_cells <- env_integer_10("CELLRANK_MIN_CELLS", 200L)
  base$cellrank_min_velocity_confidence <- env_numeric_10("CELLRANK_MIN_VELOCITY_CONFIDENCE", 0.05)
  base$cellrank_n_states <- env_integer_10("CELLRANK_N_STATES", 6L)

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
    cfg$velocity_methods_table_dir,
    cfg$velocity_scvelo_table_dir,
    cfg$velocity_velocyto_table_dir,
    cfg$velocity_driver_table_dir,
    cfg$velocity_cellrank_table_dir,
    cfg$velocity_consistency_table_dir,
    cfg$velocity_root_terminal_table_dir,
    cfg$velocity_final_table_dir,
    cfg$velocity_figure_dir,
    cfg$velocity_scvelo_figure_dir,
    cfg$velocity_cellrank_figure_dir,
    cfg$velocity_consistency_figure_dir,
    cfg$velocity_report_dir,
    cfg$velocity_barcode_dir,
    dirname(cfg$velocity_loom_jobs_tsv),
    dirname(cfg$velocity_loom_index_tsv),
    dirname(cfg$velocity_reference_index_tsv),
    dirname(cfg$velocity_scvelo_index_tsv),
    dirname(cfg$velocity_scvelo_qc_tsv),
    dirname(cfg$velocity_velocyto_steady_index_tsv),
    dirname(cfg$velocity_driver_index_tsv),
    dirname(cfg$velocity_driver_overlap_tsv),
    dirname(cfg$velocity_cellrank_index_tsv),
    dirname(cfg$velocity_consistency_index_tsv),
    dirname(cfg$velocity_consistency_summary_tsv),
    dirname(cfg$velocity_split_compare_tsv),
    dirname(cfg$velocity_root_terminal_index_tsv),
    dirname(cfg$velocity_module_status_tsv),
    dirname(cfg$velocity_triage_tsv),
    dirname(cfg$velocity_report_md),
    dirname(cfg$module_10a_manifest_path),
    dirname(cfg$module_10b_manifest_path),
    dirname(cfg$module_10c_manifest_path),
    dirname(cfg$module_10d_manifest_path),
    dirname(cfg$module_10e_manifest_path),
    dirname(cfg$module_10f_manifest_path),
    dirname(cfg$module_10g_manifest_path),
    dirname(cfg$module_10h_manifest_path),
    dirname(cfg$module_10i_manifest_path)
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

velocity_scvelo_output_dir_10 <- function(cfg, pair_id, split_value = "") {
  file.path(cfg$velocity_output_dir, "scvelo", velocity_unit_file_id_10(pair_id, split_value))
}

velocity_scvelo_figure_dir_10 <- function(cfg, pair_id, split_value = "") {
  file.path(cfg$velocity_scvelo_figure_dir, velocity_unit_file_id_10(pair_id, split_value))
}

velocity_scvelo_h5ad_path_10 <- function(cfg, pair_id, split_value = "") {
  file.path(velocity_scvelo_output_dir_10(cfg, pair_id, split_value), sprintf("scvelo_result_%s.h5ad", velocity_unit_file_id_10(pair_id, split_value)))
}

velocity_scvelo_qc_path_10 <- function(cfg, pair_id, split_value = "") {
  file.path(velocity_scvelo_output_dir_10(cfg, pair_id, split_value), sprintf("velocity_qc_%s.tsv", velocity_unit_file_id_10(pair_id, split_value)))
}

velocity_velocyto_steady_dir_10 <- function(cfg, pair_id, split_value = "") {
  file.path(cfg$velocity_output_dir, "velocyto_steady", velocity_unit_file_id_10(pair_id, split_value))
}

velocity_velocyto_steady_rds_path_10 <- function(cfg, pair_id, split_value = "") {
  file.path(velocity_velocyto_steady_dir_10(cfg, pair_id, split_value), sprintf("velocyto_steady_%s.rds", velocity_unit_file_id_10(pair_id, split_value)))
}

velocity_velocyto_direction_path_10 <- function(cfg, pair_id, split_value = "") {
  file.path(velocity_velocyto_steady_dir_10(cfg, pair_id, split_value), sprintf("velocyto_steady_direction_%s.tsv", velocity_unit_file_id_10(pair_id, split_value)))
}

velocity_driver_path_10 <- function(cfg, pair_id, split_value = "") {
  file.path(cfg$velocity_driver_table_dir, sprintf("velocity_drivers_%s_top%s.tsv", velocity_unit_file_id_10(pair_id, split_value), cfg$velocity_driver_top_n))
}

velocity_cellrank_fate_path_10 <- function(cfg, pair_id, split_value = "") {
  file.path(cfg$velocity_cellrank_table_dir, sprintf("cellrank_fate_%s.csv", velocity_unit_file_id_10(pair_id, split_value)))
}

velocity_cellrank_macrostates_path_10 <- function(cfg, pair_id, split_value = "") {
  file.path(cfg$velocity_cellrank_table_dir, sprintf("cellrank_macrostates_%s.tsv", velocity_unit_file_id_10(pair_id, split_value)))
}

velocity_cellrank_figure_path_10 <- function(cfg, pair_id, split_value = "", stem = "Figure_CR_fate") {
  file.path(cfg$velocity_cellrank_figure_dir, velocity_unit_file_id_10(pair_id, split_value), sprintf("%s_%s.png", stem, velocity_unit_file_id_10(pair_id, split_value)))
}

velocity_consistency_path_10 <- function(cfg, pair_id, split_value = "") {
  file.path(cfg$velocity_consistency_table_dir, sprintf("velocity_consistency_%s.tsv", velocity_unit_file_id_10(pair_id, split_value)))
}

velocity_consistency_figure_path_10 <- function(cfg, pair_id, split_value = "") {
  file.path(cfg$velocity_consistency_figure_dir, velocity_unit_file_id_10(pair_id, split_value), sprintf("Figure_VC_%s.png", velocity_unit_file_id_10(pair_id, split_value)))
}

velocity_split_compare_path_10 <- function(cfg, pair_id) {
  file.path(cfg$velocity_consistency_table_dir, sprintf("velocity_split_compare_%s.tsv", safe_id_09(pair_id)))
}

velocity_split_compare_figure_path_10 <- function(cfg, pair_id) {
  file.path(cfg$velocity_consistency_figure_dir, safe_id_09(pair_id), sprintf("Figure_VC_split_%s.png", safe_id_09(pair_id)))
}

velocity_root_terminal_path_10 <- function(cfg, pair_id, split_value = "") {
  file.path(cfg$velocity_output_dir, sprintf("velocity_root_terminal_%s.tsv", velocity_unit_file_id_10(pair_id, split_value)))
}
