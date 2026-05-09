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
  base$module_version <- env_or_default_03("MODULE_09_VERSION", "1.0")
  base$module_versions_09 <- list(
    `09a_trajectory_inputs` = env_or_default_03("MODULE_09A_VERSION", base$module_version),
    `09b_trajectory_inputs_eda` = env_or_default_03("MODULE_09B_VERSION", base$module_version),
    `09c_trajectory_root` = env_or_default_03("MODULE_09C_VERSION", base$module_version),
    `09d_trajectory_slingshot` = env_or_default_03("MODULE_09D_VERSION", base$module_version),
    `09e_trajectory_monocle3` = env_or_default_03("MODULE_09E_VERSION", base$module_version),
    `09e2_trajectory_monocle2` = env_or_default_03("MODULE_09E2_VERSION", base$module_version),
    `09f_trajectory_paga_dpt` = env_or_default_03("MODULE_09F_VERSION", base$module_version),
    `09g_trajectory_palantir` = env_or_default_03("MODULE_09G_VERSION", base$module_version),
    `09h_trajectory_tradeseq` = env_or_default_03("MODULE_09H_VERSION", base$module_version),
    `09i_trajectory_consensus` = env_or_default_03("MODULE_09I_VERSION", base$module_version),
    `09j_trajectory_velocity_link` = env_or_default_03("MODULE_09J_VERSION", base$module_version),
    `09k_trajectory_figures` = env_or_default_03("MODULE_09K_VERSION", base$module_version),
    `09m_trajectory_split_compare` = env_or_default_03("MODULE_09M_VERSION", base$module_version),
    `09l_trajectory_eda` = env_or_default_03("MODULE_09L_VERSION", base$module_version)
  )

  base$module_09a_manifest_path <- file.path(base$manifest_dir, "09a_trajectory_inputs", "_manifest.json")
  base$module_09b_manifest_path <- file.path(base$manifest_dir, "09b_trajectory_inputs_eda", "_manifest.json")
  base$module_09c_manifest_path <- file.path(base$manifest_dir, "09c_trajectory_root", "_manifest.json")
  base$module_09d_manifest_path <- file.path(base$manifest_dir, "09d_trajectory_slingshot", "_manifest.json")
  base$module_09e_manifest_path <- file.path(base$manifest_dir, "09e_trajectory_monocle3", "_manifest.json")
  base$module_09e2_manifest_path <- file.path(base$manifest_dir, "09e2_trajectory_monocle2", "_manifest.json")
  base$module_09f_manifest_path <- file.path(base$manifest_dir, "09f_trajectory_paga_dpt", "_manifest.json")
  base$module_09g_manifest_path <- file.path(base$manifest_dir, "09g_trajectory_palantir", "_manifest.json")
  base$module_09h_manifest_path <- file.path(base$manifest_dir, "09h_trajectory_tradeseq", "_manifest.json")
  base$module_09i_manifest_path <- file.path(base$manifest_dir, "09i_trajectory_consensus", "_manifest.json")
  base$module_09j_manifest_path <- file.path(base$manifest_dir, "09j_trajectory_velocity_link", "_manifest.json")
  base$module_09k_manifest_path <- file.path(base$manifest_dir, "09k_trajectory_figures", "_manifest.json")
  base$module_09m_manifest_path <- file.path(base$manifest_dir, "09m_trajectory_split_compare", "_manifest.json")
  base$module_09l_manifest_path <- file.path(base$manifest_dir, "09l_trajectory_eda", "_manifest.json")

  base$trajectory_pairs_sheet <- env_or_default_03("TRAJECTORY_PAIRS_SHEET", file.path(base$metadata_dir, "trajectory_pairs.tsv"))
  base$velocity_dir <- env_or_default_03("VELOCITY_DIR", file.path(base$results_dir, "velocity"))
  base$velocity_input_dir <- env_or_default_03("VELOCITY_INPUT_DIR", file.path(base$velocity_dir, "input"))
  base$velocity_loom_dir <- env_or_default_03("VELOCITY_LOOM_DIR", file.path(base$velocity_dir, "loom"))
  base$velocity_output_dir <- env_or_default_03("VELOCITY_OUTPUT_DIR", file.path(base$velocity_dir, "output"))
  base$trajectory_dir <- env_or_default_03("TRAJECTORY_DIR", file.path(base$results_dir, "trajectory"))
  base$trajectory_table_dir <- file.path(base$table_dir, "trajectory")
  base$trajectory_inputs_table_dir <- file.path(base$trajectory_table_dir, "inputs")
  base$trajectory_methods_table_dir <- file.path(base$trajectory_table_dir, "methods")
  base$trajectory_root_table_dir <- file.path(base$trajectory_methods_table_dir, "root")
  base$trajectory_slingshot_table_dir <- file.path(base$trajectory_methods_table_dir, "slingshot")
  base$trajectory_monocle3_table_dir <- file.path(base$trajectory_methods_table_dir, "monocle3")
  base$trajectory_monocle2_table_dir <- file.path(base$trajectory_methods_table_dir, "monocle2")
  base$trajectory_paga_table_dir <- file.path(base$trajectory_methods_table_dir, "paga_dpt")
  base$trajectory_palantir_table_dir <- file.path(base$trajectory_methods_table_dir, "palantir")
  base$trajectory_tradeseq_table_dir <- file.path(base$trajectory_methods_table_dir, "tradeseq")
  base$trajectory_consensus_table_dir <- file.path(base$trajectory_methods_table_dir, "consensus")
  base$trajectory_velocity_link_table_dir <- file.path(base$trajectory_methods_table_dir, "velocity_link")
  base$trajectory_figures_table_dir <- file.path(base$trajectory_table_dir, "figures")
  base$trajectory_split_compare_table_dir <- file.path(base$trajectory_table_dir, "split_compare")
  base$trajectory_final_table_dir <- file.path(base$trajectory_table_dir, "final")
  base$trajectory_figure_dir <- file.path(base$figure_dir, "trajectory")
  base$trajectory_inputs_figure_dir <- file.path(base$trajectory_figure_dir, "inputs")
  base$trajectory_methods_figure_dir <- file.path(base$trajectory_figure_dir, "methods")
  base$trajectory_monocle2_figure_dir <- file.path(base$trajectory_methods_figure_dir, "monocle2")
  base$trajectory_paga_figure_dir <- file.path(base$trajectory_methods_figure_dir, "paga_dpt")
  base$trajectory_palantir_figure_dir <- file.path(base$trajectory_methods_figure_dir, "palantir")
  base$trajectory_tradeseq_figure_dir <- file.path(base$trajectory_methods_figure_dir, "tradeseq")
  base$trajectory_consensus_figure_dir <- file.path(base$trajectory_methods_figure_dir, "consensus")
  base$trajectory_velocity_link_figure_dir <- file.path(base$trajectory_methods_figure_dir, "velocity_link")
  base$trajectory_summary_figure_dir <- file.path(base$trajectory_figure_dir, "summary")
  base$trajectory_split_compare_figure_dir <- file.path(base$trajectory_figure_dir, "split_compare")
  base$trajectory_python_input_dir <- file.path(base$trajectory_dir, "_python_input")
  base$trajectory_report_dir <- env_or_default_03("TRAJECTORY_REPORT_DIR", file.path(base$eda_report_dir, "trajectory"))
  base$trajectory_inputs_report_dir <- file.path(base$trajectory_report_dir, "inputs")
  base$trajectory_split_compare_report_dir <- file.path(base$trajectory_report_dir, "split_compare")

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
  base$trajectory_root_index_tsv <- file.path(base$trajectory_root_table_dir, "root_index.tsv")
  base$trajectory_root_split_agreement_tsv <- file.path(base$trajectory_root_table_dir, "split_root_agreement.tsv")
  base$trajectory_slingshot_index_tsv <- file.path(base$trajectory_slingshot_table_dir, "slingshot_index.tsv")
  base$trajectory_monocle3_index_tsv <- file.path(base$trajectory_monocle3_table_dir, "monocle3_index.tsv")
  base$trajectory_monocle2_index_tsv <- file.path(base$trajectory_monocle2_table_dir, "monocle2_index.tsv")
  base$trajectory_paga_index_tsv <- file.path(base$trajectory_paga_table_dir, "paga_dpt_index.tsv")
  base$trajectory_paga_jobs_tsv <- file.path(base$trajectory_paga_table_dir, "paga_dpt_jobs.tsv")
  base$trajectory_palantir_index_tsv <- file.path(base$trajectory_palantir_table_dir, "palantir_index.tsv")
  base$trajectory_palantir_jobs_tsv <- file.path(base$trajectory_palantir_table_dir, "palantir_jobs.tsv")
  base$trajectory_tradeseq_index_tsv <- file.path(base$trajectory_tradeseq_table_dir, "tradeseq_index.tsv")
  base$trajectory_tradeseq_condition_index_tsv <- file.path(base$trajectory_tradeseq_table_dir, "tradeseq_condition_index.tsv")
  base$trajectory_tradeseq_driver_index_tsv <- file.path(base$trajectory_tradeseq_table_dir, "tradeseq_driver_index.tsv")
  base$trajectory_consensus_index_tsv <- file.path(base$trajectory_consensus_table_dir, "consensus_index.tsv")
  base$trajectory_consensus_correlation_summary_tsv <- file.path(base$trajectory_consensus_table_dir, "consensus_pseudotime_correlation_summary.tsv")
  base$trajectory_consensus_jaccard_summary_tsv <- file.path(base$trajectory_consensus_table_dir, "consensus_jaccard_summary.tsv")
  base$trajectory_consensus_split_agreement_tsv <- file.path(base$trajectory_consensus_table_dir, "consensus_split_agreement.tsv")
  base$trajectory_consensus_triage_tsv <- file.path(base$trajectory_consensus_table_dir, "triage.tsv")
  base$trajectory_velocity_link_index_tsv <- file.path(base$trajectory_velocity_link_table_dir, "velocity_link_index.tsv")
  base$trajectory_figure_index_tsv <- file.path(base$trajectory_figures_table_dir, "figure_index.tsv")
  base$trajectory_split_compare_index_tsv <- file.path(base$trajectory_split_compare_table_dir, "split_compare_index.tsv")
  base$trajectory_pseudotime_distribution_compare_tsv <- file.path(base$trajectory_split_compare_table_dir, "pseudotime_distribution_compare.tsv")
  base$trajectory_lineage_backbone_match_tsv <- file.path(base$trajectory_split_compare_table_dir, "lineage_backbone_match.tsv")
  base$trajectory_driver_overlap_tsv <- file.path(base$trajectory_split_compare_table_dir, "driver_overlap.tsv")
  base$trajectory_option_a_b_consistency_tsv <- file.path(base$trajectory_split_compare_table_dir, "option_a_vs_b_consistency.tsv")
  base$trajectory_entropy_compare_tsv <- file.path(base$trajectory_split_compare_table_dir, "entropy_compare.tsv")
  base$trajectory_split_compare_report_md <- file.path(base$trajectory_split_compare_report_dir, "report.md")
  base$trajectory_module_status_tsv <- file.path(base$trajectory_final_table_dir, "trajectory_module_status.tsv")
  base$trajectory_triage_tsv <- file.path(base$trajectory_final_table_dir, "trajectory_triage.tsv")
  base$trajectory_report_md <- file.path(base$trajectory_report_dir, "report.md")

  base$trajectory_hvg_nfeatures <- env_integer_09("TRAJECTORY_HVG_NFEATURES", base$hvg_nfeatures %||% 2000L)
  base$trajectory_pca_dims <- parse_index_spec_local(env_or_default_03("TRAJECTORY_PCA_DIMS", "1:30"), default = 1:30)
  base$trajectory_umap_n_neighbors <- env_integer_09("TRAJECTORY_UMAP_N_NEIGHBORS", 30L)
  base$trajectory_split_min_cells <- env_integer_09("TRAJECTORY_SPLIT_MIN_CELLS", 50L)
  base$trajectory_balance_warn_fraction <- env_numeric_09("TRAJECTORY_BALANCE_WARN_FRACTION", 0.30)
  base$trajectory_outlier_low_prob <- env_numeric_09("TRAJECTORY_OUTLIER_LOW_PROB", 0.01)
  base$trajectory_outlier_high_prob <- env_numeric_09("TRAJECTORY_OUTLIER_HIGH_PROB", 0.99)
  base$trajectory_outlier_neighbor_cutoff <- env_numeric_09("TRAJECTORY_OUTLIER_NEIGHBOR_CUTOFF", 4)
  base$trajectory_outlier_strict_low_prob <- env_numeric_09("TRAJECTORY_OUTLIER_STRICT_LOW_PROB", 0.02)
  base$trajectory_outlier_strict_high_prob <- env_numeric_09("TRAJECTORY_OUTLIER_STRICT_HIGH_PROB", 0.98)
  base$trajectory_outlier_strict_neighbor_cutoff <- env_numeric_09("TRAJECTORY_OUTLIER_STRICT_NEIGHBOR_CUTOFF", 3)
  base$tradeseq_knots <- env_integer_09("TRADESEQ_KNOTS", 6L)
  base$trajectory_consensus_min_methods <- env_integer_09("TRAJECTORY_CONSENSUS_MIN_METHODS", 2L)
  base$trajectory_consensus_conflict_rho <- env_numeric_09("TRAJECTORY_CONSENSUS_CONFLICT_RHO", 0.30)
  base
}

prepare_dirs_09 <- function(cfg) {
  ensure_dirs <- c(
    cfg$trajectory_dir,
    cfg$trajectory_table_dir,
    cfg$trajectory_inputs_table_dir,
    cfg$trajectory_methods_table_dir,
    cfg$trajectory_root_table_dir,
    cfg$trajectory_slingshot_table_dir,
    cfg$trajectory_monocle3_table_dir,
    cfg$trajectory_monocle2_table_dir,
    cfg$trajectory_paga_table_dir,
    cfg$trajectory_palantir_table_dir,
    cfg$trajectory_tradeseq_table_dir,
    cfg$trajectory_consensus_table_dir,
    cfg$trajectory_velocity_link_table_dir,
    cfg$trajectory_figures_table_dir,
    cfg$trajectory_split_compare_table_dir,
    cfg$trajectory_final_table_dir,
    cfg$trajectory_figure_dir,
    cfg$trajectory_inputs_figure_dir,
    cfg$trajectory_methods_figure_dir,
    cfg$trajectory_monocle2_figure_dir,
    cfg$trajectory_paga_figure_dir,
    cfg$trajectory_palantir_figure_dir,
    cfg$trajectory_tradeseq_figure_dir,
    cfg$trajectory_consensus_figure_dir,
    cfg$trajectory_velocity_link_figure_dir,
    cfg$trajectory_summary_figure_dir,
    cfg$trajectory_split_compare_figure_dir,
    cfg$trajectory_python_input_dir,
    cfg$trajectory_report_dir,
    cfg$trajectory_inputs_report_dir,
    cfg$trajectory_split_compare_report_dir,
    dirname(cfg$module_09a_manifest_path),
    dirname(cfg$module_09b_manifest_path),
    dirname(cfg$module_09c_manifest_path),
    dirname(cfg$module_09d_manifest_path),
    dirname(cfg$module_09e_manifest_path),
    dirname(cfg$module_09e2_manifest_path),
    dirname(cfg$module_09f_manifest_path),
    dirname(cfg$module_09g_manifest_path),
    dirname(cfg$module_09h_manifest_path),
    dirname(cfg$module_09i_manifest_path),
    dirname(cfg$module_09j_manifest_path),
    dirname(cfg$module_09k_manifest_path),
    dirname(cfg$module_09m_manifest_path),
    dirname(cfg$module_09l_manifest_path),
    dirname(cfg$trajectory_root_index_tsv),
    dirname(cfg$trajectory_root_split_agreement_tsv),
    dirname(cfg$trajectory_slingshot_index_tsv),
    dirname(cfg$trajectory_monocle3_index_tsv),
    dirname(cfg$trajectory_monocle2_index_tsv),
    dirname(cfg$trajectory_paga_index_tsv),
    dirname(cfg$trajectory_paga_jobs_tsv),
    dirname(cfg$trajectory_palantir_index_tsv),
    dirname(cfg$trajectory_palantir_jobs_tsv),
    dirname(cfg$trajectory_tradeseq_index_tsv),
    dirname(cfg$trajectory_tradeseq_condition_index_tsv),
    dirname(cfg$trajectory_tradeseq_driver_index_tsv),
    dirname(cfg$trajectory_consensus_index_tsv),
    dirname(cfg$trajectory_consensus_correlation_summary_tsv),
    dirname(cfg$trajectory_consensus_jaccard_summary_tsv),
    dirname(cfg$trajectory_consensus_split_agreement_tsv),
    dirname(cfg$trajectory_consensus_triage_tsv),
    dirname(cfg$trajectory_velocity_link_index_tsv),
    dirname(cfg$trajectory_figure_index_tsv),
    dirname(cfg$trajectory_split_compare_index_tsv),
    dirname(cfg$trajectory_pseudotime_distribution_compare_tsv),
    dirname(cfg$trajectory_lineage_backbone_match_tsv),
    dirname(cfg$trajectory_driver_overlap_tsv),
    dirname(cfg$trajectory_option_a_b_consistency_tsv),
    dirname(cfg$trajectory_entropy_compare_tsv),
    dirname(cfg$trajectory_split_compare_report_md),
    dirname(cfg$trajectory_module_status_tsv),
    dirname(cfg$trajectory_triage_tsv),
    dirname(cfg$trajectory_report_md)
  )
  invisible(lapply(ensure_dirs[nzchar(ensure_dirs)], ensure_dir))
}

safe_id_09 <- function(value) {
  safe_id_07(value)
}

trajectory_method_registry_09 <- function() {
  list(
    root = list(table_dir_key = "trajectory_root_table_dir"),
    slingshot = list(table_dir_key = "trajectory_slingshot_table_dir", figure_dir_key = "trajectory_methods_figure_dir"),
    monocle3 = list(table_dir_key = "trajectory_monocle3_table_dir", figure_dir_key = "trajectory_methods_figure_dir"),
    monocle2 = list(table_dir_key = "trajectory_monocle2_table_dir", figure_dir_key = "trajectory_monocle2_figure_dir"),
    paga_dpt = list(table_dir_key = "trajectory_paga_table_dir", figure_dir_key = "trajectory_paga_figure_dir"),
    palantir = list(table_dir_key = "trajectory_palantir_table_dir", figure_dir_key = "trajectory_palantir_figure_dir"),
    tradeseq = list(table_dir_key = "trajectory_tradeseq_table_dir", figure_dir_key = "trajectory_tradeseq_figure_dir"),
    consensus = list(table_dir_key = "trajectory_consensus_table_dir", figure_dir_key = "trajectory_consensus_figure_dir"),
    velocity_link = list(table_dir_key = "trajectory_velocity_link_table_dir", figure_dir_key = "trajectory_velocity_link_figure_dir"),
    figures = list(table_dir_key = "trajectory_figures_table_dir", figure_dir_key = "trajectory_summary_figure_dir"),
    split_compare = list(table_dir_key = "trajectory_split_compare_table_dir", figure_dir_key = "trajectory_split_compare_figure_dir"),
    final = list(table_dir_key = "trajectory_final_table_dir")
  )
}

trajectory_method_root_dir_09 <- function(cfg, method, kind = c("table", "figure")) {
  kind <- match.arg(kind)
  method <- tolower(normalize_scalar_value(method))
  registry <- trajectory_method_registry_09()
  entry <- registry[[method]]
  key <- if (!is.null(entry)) entry[[paste0(kind, "_dir_key")]] else NULL
  if (!is.null(key) && nzchar(key) && !is.null(cfg[[key]])) {
    return(cfg[[key]])
  }
  if (identical(kind, "table")) {
    return(file.path(cfg$trajectory_methods_table_dir, safe_id_09(method)))
  }
  file.path(cfg$trajectory_methods_figure_dir, safe_id_09(method))
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

trajectory_unit_value_09 <- function(split_value = "") {
  value <- normalize_scalar_value(split_value)
  if (!nzchar(value) || identical(value, "pooled")) "" else value
}

trajectory_unit_file_id_09 <- function(pair_id, split_value = "") {
  value <- trajectory_unit_value_09(split_value)
  if (!nzchar(value)) {
    return(safe_id_09(pair_id))
  }
  sprintf("%s__%s", safe_id_09(pair_id), safe_id_09(value))
}

trajectory_method_table_dir_09 <- function(cfg, method, pair_id = "", split_value = "") {
  root <- trajectory_method_root_dir_09(cfg, method, "table")
  if (!nzchar(normalize_scalar_value(pair_id))) {
    return(root)
  }
  file.path(root, trajectory_unit_file_id_09(pair_id, split_value))
}

trajectory_method_figure_dir_09 <- function(cfg, method, pair_id, split_value = "") {
  root <- trajectory_method_root_dir_09(cfg, method, "figure")
  file.path(root, trajectory_unit_file_id_09(pair_id, split_value))
}

trajectory_root_inference_path_09 <- function(cfg, pair_id, split_value = "") {
  file.path(
    trajectory_method_table_dir_09(cfg, "root", pair_id, split_value),
    sprintf("root_inference_%s.tsv", trajectory_unit_file_id_09(pair_id, split_value))
  )
}

trajectory_method_output_path_09 <- function(cfg, method, pair_id, split_value, stem, ext = "csv") {
  file.path(
    trajectory_method_table_dir_09(cfg, method, pair_id, split_value),
    sprintf("%s_%s.%s", stem, trajectory_unit_file_id_09(pair_id, split_value), ext)
  )
}

trajectory_method_figure_path_09 <- function(cfg, method, pair_id, split_value, stem, ext = "png") {
  file.path(
    trajectory_method_figure_dir_09(cfg, method, pair_id, split_value),
    sprintf("%s_%s.%s", stem, trajectory_unit_file_id_09(pair_id, split_value), ext)
  )
}

trajectory_python_input_dir_09 <- function(cfg, method, pair_id, split_value = "") {
  file.path(cfg$trajectory_python_input_dir, safe_id_09(method), trajectory_unit_file_id_09(pair_id, split_value))
}
