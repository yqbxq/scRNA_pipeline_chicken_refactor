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
  base$module_07b_manifest_path <- file.path(base$manifest_dir, "07b_nichenet", "_manifest.json")
  base$module_07c_manifest_path <- file.path(base$manifest_dir, "07c_communication_eda", "_manifest.json")

  base$communication_pairs_sheet <- env_or_default_03("COMMUNICATION_PAIRS_SHEET", file.path(base$metadata_dir, "communication_pairs.tsv"))
  base$communication_report_dir <- env_or_default_03("COMMUNICATION_REPORT_DIR", file.path(base$eda_report_dir, "communication"))
  base$communication_table_dir <- file.path(base$table_dir, "communication")
  base$communication_figure_dir <- file.path(base$figure_dir, "communication")
  base$cellchat_table_dir <- file.path(base$communication_table_dir, "cellchat")
  base$cellchat_figure_dir <- file.path(base$communication_figure_dir, "cellchat")
  base$nichenet_table_dir <- file.path(base$communication_table_dir, "nichenet")
  base$nichenet_figure_dir <- file.path(base$communication_figure_dir, "nichenet")

  base$cellchat_index_tsv <- file.path(base$cellchat_table_dir, "cellchat_index.tsv")
  base$cellchat_mapping_summary_tsv <- file.path(base$cellchat_table_dir, "mapping_summary.tsv")
  base$cellchat_triage_tsv <- file.path(base$cellchat_table_dir, "triage.tsv")
  base$nichenet_index_tsv <- file.path(base$nichenet_table_dir, "nichenet_index.tsv")
  base$nichenet_roles_resolved_tsv <- file.path(base$nichenet_table_dir, "roles_resolved.tsv")
  base$nichenet_triage_tsv <- file.path(base$nichenet_table_dir, "triage.tsv")
  base$communication_cross_validation_tsv <- file.path(base$communication_table_dir, "cross_validation.tsv")
  base$communication_consensus_lr_tsv <- file.path(base$communication_table_dir, "consensus_lr.tsv")
  base$communication_triage_tsv <- file.path(base$communication_table_dir, "triage.tsv")
  base$communication_summary_plot_png <- file.path(base$communication_figure_dir, "summary_jaccard.png")
  base$communication_report_md <- file.path(base$communication_report_dir, "report.md")

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
    cfg$nichenet_table_dir,
    cfg$nichenet_figure_dir,
    dirname(cfg$module_07a_manifest_path),
    dirname(cfg$module_07b_manifest_path),
    dirname(cfg$module_07c_manifest_path),
    dirname(cfg$cellchat_index_tsv),
    dirname(cfg$cellchat_mapping_summary_tsv),
    dirname(cfg$cellchat_triage_tsv),
    dirname(cfg$nichenet_index_tsv),
    dirname(cfg$nichenet_roles_resolved_tsv),
    dirname(cfg$nichenet_triage_tsv),
    dirname(cfg$communication_cross_validation_tsv),
    dirname(cfg$communication_consensus_lr_tsv),
    dirname(cfg$communication_triage_tsv),
    dirname(cfg$communication_summary_plot_png),
    dirname(cfg$communication_report_md)
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
