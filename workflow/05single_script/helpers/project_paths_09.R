env_numeric_09 <- function(name, default) {
  value <- suppressWarnings(as.numeric(env_or_default_03(name, as.character(default))))
  if (is.na(value)) default else value
}

env_integer_09 <- function(name, default) {
  value <- suppressWarnings(as.integer(env_or_default_03(name, as.character(default))))
  if (is.na(value)) default else value
}

get_single_script_config_09 <- function() {
  base <- get_single_script_config_07()

  base$module_09a_manifest_path <- file.path(base$manifest_dir, "09a_trajectory_inputs", "_manifest.json")
  base$module_09b_manifest_path <- file.path(base$manifest_dir, "09b_trajectory_inputs_eda", "_manifest.json")

  base$trajectory_pairs_sheet <- env_or_default_03("TRAJECTORY_PAIRS_SHEET", file.path(base$metadata_dir, "trajectory_pairs.tsv"))
  base$trajectory_dir <- env_or_default_03("TRAJECTORY_DIR", file.path(base$results_dir, "trajectory"))
  base$trajectory_table_dir <- file.path(base$table_dir, "trajectory")
  base$trajectory_inputs_table_dir <- file.path(base$trajectory_table_dir, "inputs")
  base$trajectory_figure_dir <- file.path(base$figure_dir, "trajectory")
  base$trajectory_inputs_figure_dir <- file.path(base$trajectory_figure_dir, "inputs")
  base$trajectory_report_dir <- env_or_default_03("TRAJECTORY_REPORT_DIR", file.path(base$eda_report_dir, "trajectory"))
  base$trajectory_inputs_report_dir <- file.path(base$trajectory_report_dir, "inputs")

  base$trajectory_input_index_tsv <- file.path(base$trajectory_inputs_table_dir, "trajectory_input_index.tsv")
  base$trajectory_outlier_qc_summary_tsv <- file.path(base$trajectory_inputs_table_dir, "outlier_qc_summary.tsv")
  base$trajectory_outlier_qc_cells_tsv <- file.path(base$trajectory_inputs_table_dir, "outlier_qc_cells.tsv")
  base$trajectory_cc_score_summary_tsv <- file.path(base$trajectory_inputs_table_dir, "cc_score_summary.tsv")
  base$trajectory_pca_summary_tsv <- file.path(base$trajectory_inputs_table_dir, "pca_summary.tsv")
  base$trajectory_split_index_tsv <- file.path(base$trajectory_inputs_table_dir, "split_index.tsv")
  base$trajectory_input_status_tsv <- file.path(base$trajectory_inputs_table_dir, "input_status.tsv")
  base$trajectory_outlier_pre_post_summary_tsv <- file.path(base$trajectory_inputs_table_dir, "outlier_pre_post_summary.tsv")
  base$trajectory_split_balance_summary_tsv <- file.path(base$trajectory_inputs_table_dir, "split_balance_summary.tsv")
  base$trajectory_umap_figure_index_tsv <- file.path(base$trajectory_inputs_table_dir, "umap_figure_index.tsv")
  base$trajectory_inputs_triage_tsv <- file.path(base$trajectory_inputs_table_dir, "triage.tsv")
  base$trajectory_inputs_report_md <- file.path(base$trajectory_inputs_report_dir, "report.md")
  base$trajectory_inputs_qc_violin_png <- file.path(base$trajectory_inputs_figure_dir, "outlier_pre_post_violin.png")

  base$trajectory_hvg_nfeatures <- env_integer_09("TRAJECTORY_HVG_NFEATURES", base$hvg_nfeatures %||% 2000L)
  base$trajectory_pca_dims <- parse_index_spec_local(env_or_default_03("TRAJECTORY_PCA_DIMS", "1:30"), default = 1:30)
  base$trajectory_umap_n_neighbors <- env_integer_09("TRAJECTORY_UMAP_N_NEIGHBORS", 30L)
  base$trajectory_split_min_cells <- env_integer_09("TRAJECTORY_SPLIT_MIN_CELLS", 50L)
  base$trajectory_balance_warn_fraction <- env_numeric_09("TRAJECTORY_BALANCE_WARN_FRACTION", 0.30)
  base$module_version <- env_or_default_03("MODULE_09_VERSION", "1.0")
  base
}

prepare_dirs_09 <- function(cfg) {
  ensure_dirs <- c(
    cfg$trajectory_dir,
    cfg$trajectory_table_dir,
    cfg$trajectory_inputs_table_dir,
    cfg$trajectory_figure_dir,
    cfg$trajectory_inputs_figure_dir,
    cfg$trajectory_report_dir,
    cfg$trajectory_inputs_report_dir,
    dirname(cfg$module_09a_manifest_path),
    dirname(cfg$module_09b_manifest_path)
  )
  invisible(lapply(ensure_dirs[nzchar(ensure_dirs)], ensure_dir))
}

safe_id_09 <- function(value) {
  safe_id_07(value)
}

trajectory_pair_dir_09 <- function(cfg, pair_id) {
  file.path(cfg$trajectory_dir, safe_id_09(pair_id))
}

trajectory_pair_input_dir_09 <- function(cfg, pair_id) {
  file.path(trajectory_pair_dir_09(cfg, pair_id), "inputs")
}

trajectory_input_rds_path_09 <- function(cfg, pair_id, split_value = "") {
  input_dir <- trajectory_pair_input_dir_09(cfg, pair_id)
  if (!nzchar(normalize_scalar_value(split_value))) {
    return(file.path(input_dir, "seurat.rds"))
  }
  file.path(input_dir, sprintf("seurat__%s.rds", safe_id_09(split_value)))
}

trajectory_pair_figure_dir_09 <- function(cfg, pair_id, split_value = "") {
  value_id <- normalize_scalar_value(split_value, "pooled")
  file.path(cfg$trajectory_inputs_figure_dir, safe_id_09(pair_id), safe_id_09(value_id))
}

trajectory_split_balance_path_09 <- function(cfg, pair_id) {
  file.path(cfg$trajectory_inputs_table_dir, sprintf("split_balance_%s.tsv", safe_id_09(pair_id)))
}
