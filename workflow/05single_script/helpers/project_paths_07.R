env_numeric_07 <- function(name, default) {
  value <- suppressWarnings(as.numeric(env_or_default_03(name, as.character(default))))
  if (is.na(value)) default else value
}

env_integer_07 <- function(name, default) {
  value <- suppressWarnings(as.integer(env_or_default_03(name, as.character(default))))
  if (is.na(value)) default else value
}

get_single_script_config_07 <- function() {
  base <- get_single_script_config_06()

  base$module_07a_manifest_path <- file.path(base$manifest_dir, "07a_cellchat", "_manifest.json")
  base$module_07b_manifest_path <- file.path(base$manifest_dir, "07b_liana_consensus", "_manifest.json")
  base$module_07c_manifest_path <- file.path(base$manifest_dir, "07c_nichenet", "_manifest.json")
  base$module_07d_manifest_path <- file.path(base$manifest_dir, "07d_communication_consensus", "_manifest.json")
  base$module_07e_manifest_path <- file.path(base$manifest_dir, "07e_communication_eda", "_manifest.json")
  base$module_07b_liana_manifest_path <- base$module_07b_manifest_path
  base$module_07c_nichenet_manifest_path <- base$module_07c_manifest_path
  base$module_07d_consensus_manifest_path <- base$module_07d_manifest_path
  base$module_07e_eda_manifest_path <- base$module_07e_manifest_path

  base$communication_pairs_sheet <- env_or_default_03("COMMUNICATION_PAIRS_SHEET", file.path(base$metadata_dir, "communication_pairs.tsv"))
  base$communication_report_dir <- env_or_default_03("COMMUNICATION_REPORT_DIR", file.path(base$eda_report_dir, "communication"))
  base$communication_table_dir <- file.path(base$table_dir, "communication")
  base$communication_figure_dir <- file.path(base$figure_dir, "communication")
  base$cellchat_table_dir <- file.path(base$communication_table_dir, "cellchat")
  base$cellchat_figure_dir <- file.path(base$communication_figure_dir, "cellchat")
  base$liana_table_dir <- file.path(base$communication_table_dir, "liana_consensus")
  base$liana_figure_dir <- file.path(base$communication_figure_dir, "liana_consensus")
  base$nichenet_table_dir <- file.path(base$communication_table_dir, "nichenet")
  base$nichenet_figure_dir <- file.path(base$communication_figure_dir, "nichenet")
  base$communication_consensus_table_dir <- file.path(base$communication_table_dir, "consensus")
  base$communication_consensus_figure_dir <- file.path(base$communication_figure_dir, "consensus")

  base$cellchat_index_tsv <- file.path(base$cellchat_table_dir, "cellchat_index.tsv")
  base$cellchat_inventory_gate_log_tsv <- file.path(base$cellchat_table_dir, "cellchat_inventory_gate_log.tsv")
  base$cellchat_mapping_summary_tsv <- file.path(base$cellchat_table_dir, "mapping_summary.tsv")
  base$cellchat_triage_tsv <- file.path(base$cellchat_table_dir, "triage.tsv")
  base$liana_index_tsv <- file.path(base$liana_table_dir, "liana_index.tsv")
  base$liana_consensus_lr_tsv <- file.path(base$liana_table_dir, "liana_consensus_lr.tsv")
  base$liana_triage_tsv <- file.path(base$liana_table_dir, "triage.tsv")
  base$nichenet_index_tsv <- file.path(base$nichenet_table_dir, "nichenet_index.tsv")
  base$nichenet_inventory_gate_log_tsv <- file.path(base$nichenet_table_dir, "nichenet_inventory_gate_log.tsv")
  base$nichenet_roles_resolved_tsv <- file.path(base$nichenet_table_dir, "roles_resolved.tsv")
  base$nichenet_triage_tsv <- file.path(base$nichenet_table_dir, "triage.tsv")
  base$communication_method_consensus_tsv <- file.path(base$communication_consensus_table_dir, "method_consensus.tsv")
  base$communication_evidence_tier_tsv <- file.path(base$communication_consensus_table_dir, "evidence_tiers.tsv")
  base$communication_cross_validation_tsv <- file.path(base$communication_table_dir, "cross_validation.tsv")
  base$communication_consensus_lr_tsv <- file.path(base$communication_table_dir, "consensus_lr.tsv")
  base$communication_triage_tsv <- file.path(base$communication_table_dir, "triage.tsv")
  base$communication_gate_summary_tsv <- file.path(base$communication_table_dir, "communication_gate_summary.tsv")
  base$communication_skipped_low_cells_tsv <- file.path(base$communication_table_dir, "skipped_low_cells.tsv")
  base$communication_fallback_summary_tsv <- file.path(base$communication_table_dir, "communication_fallback_summary.tsv")
  base$derived_communication_eligibility_tsv <- file.path(base$communication_table_dir, "derived_communication_eligibility.tsv")
  base$communication_summary_plot_png <- file.path(base$communication_figure_dir, "summary_jaccard.png")
  base$communication_report_md <- file.path(base$communication_report_dir, "report.md")
  base$communication_report_html <- file.path(base$communication_report_dir, "report.html")
  base$communication_eda_panels_tsv <- file.path(base$communication_consensus_table_dir, "communication_eda_panels.tsv")
  base$communication_eda_filtered_tsv <- file.path(base$communication_consensus_table_dir, "communication_eda_filtered_axes.tsv")
  base$communication_eda_panel_dir <- file.path(base$communication_consensus_figure_dir, "07e_panels")
  base$commot_spatial_summary_tsv <- env_or_default_03("COMMOT_SPATIAL_SUMMARY_TSV", file.path(base$table_dir, "spatial", "08_commot", "commot_spatial_summary.tsv"))

  base$communication_cell_type_col <- env_or_default_03("COMMUNICATION_CELL_TYPE_COL", "cell_subtype")
  base$nichenet_resource_dir <- env_or_default_03("NICHENET_RESOURCE_DIR", file.path(base$project_root, "resources", "nichenet"))
  base$nichenet_lr_network_rds <- env_or_default_03("NICHENET_LR_NETWORK_RDS", file.path(base$nichenet_resource_dir, "lr_network_human_21122021.rds"))
  base$nichenet_ligand_target_matrix_rds <- env_or_default_03("NICHENET_LIGAND_TARGET_MATRIX_RDS", file.path(base$nichenet_resource_dir, "ligand_target_matrix_nsga2r_final.rds"))
  base$nichenet_weighted_networks_rds <- env_or_default_03("NICHENET_WEIGHTED_NETWORKS_RDS", file.path(base$nichenet_resource_dir, "weighted_networks_nsga2r_final.rds"))

  base$communication_min_cells_per_celltype <- env_integer_07("COMMUNICATION_MIN_CELLS_PER_CELLTYPE", 20L)
  base$communication_ortholog_min_coverage <- env_numeric_07("COMMUNICATION_ORTHOLOG_MIN_COVERAGE", 0.30)
  base$cellchat_min_cells_per_group <- env_integer_07("CELLCHAT_MIN_CELLS_PER_GROUP", 20L)
  base$cellchat_prob_cutoff <- env_numeric_07("CELLCHAT_PROB_CUTOFF", 0.05)
  base$nichenet_expression_pct <- env_numeric_07("NICHENET_EXPRESSION_PCT", 0.10)
  base$nichenet_top_ligand_n <- env_integer_07("NICHENET_TOP_LIGAND_N", 20L)
  base$nichenet_top_target_n <- env_integer_07("NICHENET_TOP_TARGET_N", 200L)
  base$communication_require_full_pipeline <- env_or_default_03("COMMUNICATION_REQUIRE_FULL_PIPELINE", "no")
  base$communication_report_filter_tier <- env_or_default_03("COMMUNICATION_REPORT_FILTER_TIER", "auto")
  base$communication_report_top_n_primary <- env_integer_07("COMMUNICATION_REPORT_TOP_N_PRIMARY", 30L)
  base$communication_report_html_enabled <- env_or_default_03("COMMUNICATION_REPORT_HTML", "yes")
  base$consensus_min_methods_for_primary <- env_integer_07("CONSENSUS_MIN_METHODS_FOR_PRIMARY", 2L)
  base$consensus_require_nichenet_for_primary <- env_or_default_03("CONSENSUS_REQUIRE_NICHENET_FOR_PRIMARY", "yes")
  base$liana_consensus_enabled <- env_or_default_03("LIANA_CONSENSUS_ENABLED", "yes")
  base$liana_methods_list <- env_or_default_03("LIANA_METHODS_LIST", "cellphonedb,connectome,sca,natmi,logfc,rank_aggregate")
  base$liana_consensus_aggregate <- env_or_default_03("LIANA_CONSENSUS_AGGREGATE", "rank_aggregate")
  base$liana_min_methods_agreed_inside <- env_integer_07("LIANA_MIN_METHODS_AGREED_INSIDE", 3L)
  base$liana_cellphonedb_pval_threshold <- env_numeric_07("LIANA_CELLPHONEDB_PVAL_THRESHOLD", 0.05)
  base$liana_resource_db <- env_or_default_03("LIANA_RESOURCE_DB", "consensus")
  base$liana_min_cells_per_group <- env_integer_07("LIANA_MIN_CELLS_PER_GROUP", 100L)
  base$liana_condition_col <- env_or_default_03("LIANA_CONDITION_COL", "condition")
  base$ortholog_chicken_human_tsv <- env_or_default_03("ORTHOLOG_CHICKEN_HUMAN_TSV", file.path(base$metadata_dir, "ortholog_chicken_human.tsv"))
  base$multinichenet_enabled <- env_or_default_03("MULTINICHENET_ENABLED", "auto")
  base$nichenet_mode <- env_or_default_03("NICHENET_MODE", "auto")
  base$multinichenet_min_samples_per_group <- env_integer_07("MULTINICHENET_MIN_SAMPLES_PER_GROUP", 2L)
  base$multinichenet_top_n_lr <- env_integer_07("MULTINICHENET_TOP_N_LR", 250L)
  base$multinichenet_top_n_targets <- env_integer_07("MULTINICHENET_TOP_N_TARGETS", 20L)
  base$multinichenet_min_cells <- env_integer_07("MULTINICHENET_MIN_CELLS", 10L)
  base$communication_receiver_de_p_threshold <- env_numeric_07("COMMUNICATION_RECEIVER_DE_P_THRESHOLD", 0.05)
  base$commot_enabled <- env_or_default_03("COMMOT_ENABLED", "auto")
  base$module_07a_cellchat_version <- env_or_default_03(
    "MODULE_07A_CELLCHAT_VERSION",
    env_or_default_03("MODULE_07A_VERSION", env_or_default_03("MODULE_07_VERSION", "1.2"))
  )
  base$module_07b_liana_version <- env_or_default_03("MODULE_07B_LIANA_VERSION", env_or_default_03("MODULE_07_VERSION", "0.1"))
  base$module_07c_nichenet_version <- env_or_default_03(
    "MODULE_07C_NICHENET_VERSION",
    env_or_default_03("MODULE_07_VERSION", "1.1")
  )
  base$module_07d_consensus_version <- env_or_default_03("MODULE_07D_CONSENSUS_VERSION", env_or_default_03("MODULE_07_VERSION", "1.0"))
  base$module_07e_eda_version <- env_or_default_03("MODULE_07E_EDA_VERSION", env_or_default_03("MODULE_07_VERSION", "2.0"))
  base$module_07a_version <- base$module_07a_cellchat_version
  base$module_07b_version <- base$module_07b_liana_version
  base$module_07c_version <- base$module_07c_nichenet_version
  base$module_version <- env_or_default_03("MODULE_07_VERSION", "1.0")
  base
}

prepare_dirs_07 <- function(cfg) {
  ensure_dirs <- c(
    cfg$communication_report_dir,
    cfg$communication_table_dir,
    cfg$communication_figure_dir,
    cfg$cellchat_table_dir,
    cfg$cellchat_figure_dir,
    cfg$liana_table_dir,
    cfg$liana_figure_dir,
    cfg$nichenet_table_dir,
    cfg$nichenet_figure_dir,
    cfg$communication_consensus_table_dir,
    cfg$communication_consensus_figure_dir,
    dirname(cfg$module_07a_manifest_path),
    dirname(cfg$module_07b_manifest_path),
    dirname(cfg$module_07c_manifest_path),
    dirname(cfg$module_07d_manifest_path),
    dirname(cfg$module_07e_manifest_path),
    dirname(cfg$cellchat_index_tsv),
    dirname(cfg$cellchat_inventory_gate_log_tsv),
    dirname(cfg$cellchat_mapping_summary_tsv),
    dirname(cfg$cellchat_triage_tsv),
    dirname(cfg$liana_index_tsv),
    dirname(cfg$liana_consensus_lr_tsv),
    dirname(cfg$liana_triage_tsv),
    dirname(cfg$nichenet_index_tsv),
    dirname(cfg$nichenet_inventory_gate_log_tsv),
    dirname(cfg$nichenet_roles_resolved_tsv),
    dirname(cfg$nichenet_triage_tsv),
    dirname(cfg$communication_method_consensus_tsv),
    dirname(cfg$communication_evidence_tier_tsv),
    dirname(cfg$communication_cross_validation_tsv),
    dirname(cfg$communication_consensus_lr_tsv),
    dirname(cfg$communication_triage_tsv),
    dirname(cfg$communication_gate_summary_tsv),
    dirname(cfg$communication_skipped_low_cells_tsv),
    dirname(cfg$communication_fallback_summary_tsv),
    dirname(cfg$derived_communication_eligibility_tsv),
    dirname(cfg$communication_summary_plot_png),
    dirname(cfg$communication_report_md),
    dirname(cfg$communication_report_html),
    dirname(cfg$communication_eda_panels_tsv),
    dirname(cfg$communication_eda_filtered_tsv),
    cfg$communication_eda_panel_dir
  )
  invisible(lapply(ensure_dirs[nzchar(ensure_dirs)], ensure_dir))
}

safe_id_07 <- function(value) {
  safe_id_06(value)
}

cellchat_paths_07 <- function(cfg, layer_id, comparison_id, group_value) {
  safe_layer <- safe_id_07(layer_id)
  safe_cmp <- safe_id_07(comparison_id)
  safe_group <- safe_id_07(group_value)
  table_dir <- file.path(cfg$cellchat_table_dir, safe_layer, safe_cmp, safe_group)
  figure_dir <- file.path(cfg$cellchat_figure_dir, safe_layer, safe_cmp, safe_group)
  list(
    table_dir = table_dir,
    figure_dir = figure_dir,
    cellchat_rds = file.path(table_dir, "cellchat.rds"),
    lr_table_tsv = file.path(table_dir, "lr_pairs.tsv"),
    pathway_table_tsv = file.path(table_dir, "pathways.tsv"),
    bubble_png = file.path(figure_dir, "bubble.png"),
    network_png = file.path(figure_dir, "network.png"),
    heatmap_png = file.path(figure_dir, "heatmap.png")
  )
}

nichenet_paths_07 <- function(cfg, layer_id, comparison_id, receiver_set) {
  safe_layer <- safe_id_07(layer_id)
  safe_cmp <- safe_id_07(comparison_id)
  safe_receiver <- safe_id_07(receiver_set)
  table_dir <- file.path(cfg$nichenet_table_dir, safe_layer, safe_cmp, safe_receiver)
  figure_dir <- file.path(cfg$nichenet_figure_dir, safe_layer, safe_cmp, safe_receiver)
  list(
    table_dir = table_dir,
    figure_dir = figure_dir,
    ligand_activity_tsv = file.path(table_dir, "ligand_activity.tsv"),
    ligand_target_links_tsv = file.path(table_dir, "ligand_target_links.tsv"),
    nichenet_rds = file.path(table_dir, "nichenet.rds"),
    ligand_activity_heatmap_png = file.path(figure_dir, "ligand_activity_heatmap.png"),
    ligand_target_heatmap_png = file.path(figure_dir, "ligand_target_heatmap.png"),
    circos_png = file.path(figure_dir, "circos.png")
  )
}

communication_layer_report_path_07 <- function(cfg, layer_id) {
  file.path(cfg$communication_report_dir, "layers", safe_id_07(layer_id), "report.md")
}
