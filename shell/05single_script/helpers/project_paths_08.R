env_numeric_08 <- function(name, default) {
  value <- suppressWarnings(as.numeric(env_or_default_03(name, as.character(default))))
  if (is.na(value)) default else value
}

env_integer_08 <- function(name, default) {
  value <- suppressWarnings(as.integer(env_or_default_03(name, as.character(default))))
  if (is.na(value)) default else value
}

get_single_script_config_08 <- function() {
  base <- get_single_script_config_07()

  base$module_08a_manifest_path <- file.path(base$manifest_dir, "08a_scenic_export", "_manifest.json")
  base$module_08b_manifest_path <- file.path(base$manifest_dir, "08b_scenic_grn", "_manifest.json")
  base$module_08c_manifest_path <- file.path(base$manifest_dir, "08c_scenic_regulons", "_manifest.json")
  base$module_08d_manifest_path <- file.path(base$manifest_dir, "08d_scenic_downstream", "_manifest.json")

  base$scenic_input_dir <- env_or_default_03("SCENIC_INPUT_DIR", file.path(base$results_dir, "scenic_input"))
  base$scenic_output_dir <- env_or_default_03("SCENIC_OUTPUT_DIR", file.path(base$results_dir, "scenic_output"))
  base$scenic_db_dir <- env_or_default_03("SCENIC_DB_DIR", file.path(base$project_root, "resources", "scenic_db"))
  base$scenic_tf_list <- env_or_default_03("SCENIC_TF_LIST", file.path(base$scenic_db_dir, "hs_hgnc_tfs.txt"))
  base$scenic_motif_ann <- env_or_default_03("SCENIC_MOTIF_ANN", file.path(base$scenic_db_dir, "motifs-v9-nr.hgnc-m0.001-o0.0.tbl"))
  base$scenic_db_500bp <- env_or_default_03("SCENIC_DB_500BP", file.path(base$scenic_db_dir, "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"))
  base$scenic_db_10kb <- env_or_default_03("SCENIC_DB_10KB", file.path(base$scenic_db_dir, "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"))
  base$scenic_resource_manifest <- env_or_default_03("SCENIC_RESOURCE_MANIFEST", file.path(base$scenic_db_dir, "scenic_resources_manifest.json"))

  base$regulation_report_dir <- env_or_default_03("REGULATION_REPORT_DIR", file.path(base$eda_report_dir, "regulation"))
  base$regulation_table_dir <- file.path(base$table_dir, "regulation")
  base$regulation_figure_dir <- file.path(base$figure_dir, "regulation")
  base$regulation_checkpoint_dir <- file.path(base$checkpoint_dir, "regulation")
  base$scenic_table_dir <- file.path(base$regulation_table_dir, "scenic")
  base$scenic_figure_dir <- file.path(base$regulation_figure_dir, "scenic")
  base$scenic_checkpoint_dir <- file.path(base$regulation_checkpoint_dir, "scenic")

  base$scenic_export_index_tsv <- file.path(base$scenic_table_dir, "scenic_export_index.tsv")
  base$scenic_export_mapping_summary_tsv <- file.path(base$scenic_table_dir, "mapping_summary.tsv")
  base$scenic_export_triage_tsv <- file.path(base$scenic_table_dir, "triage.tsv")
  base$scenic_grn_index_tsv <- file.path(base$scenic_table_dir, "scenic_grn_index.tsv")
  base$scenic_regulons_index_tsv <- file.path(base$scenic_table_dir, "scenic_regulons_index.tsv")
  base$scenic_downstream_index_tsv <- file.path(base$scenic_table_dir, "scenic_downstream_index.tsv")

  base$regulation_layers <- split_csv_local(env_or_default_03("REGULATION_LAYERS", base$panorama_layer_id))
  if (length(base$regulation_layers) == 0) {
    base$regulation_layers <- base$panorama_layer_id
  }
  base$scenic_module_top_n <- env_integer_08("SCENIC_MODULE_TOP_N", 50L)
  base$scenic_module_min_genes <- env_integer_08("SCENIC_MODULE_MIN_GENES", 10L)
  base$scenic_regulon_min_targets <- env_integer_08("SCENIC_REGULON_MIN_TARGETS", 5L)
  base$scenic_nes_threshold <- env_numeric_08("SCENIC_NES_THRESHOLD", 3)
  base$scenic_auc_max_rank_fraction <- env_numeric_08("SCENIC_AUC_MAX_RANK_FRACTION", 0.05)
  base$scenic_threads <- env_integer_08("SCENIC_THREADS", 8L)
  base$scenic_ortholog_min_coverage <- env_numeric_08("SCENIC_ORTHOLOG_MIN_COVERAGE", 0.30)
  base$module_version <- env_or_default_03("MODULE_08_VERSION", "1.0")
  base
}

prepare_dirs_08 <- function(cfg) {
  ensure_dirs <- c(
    cfg$scenic_input_dir,
    cfg$scenic_output_dir,
    cfg$scenic_db_dir,
    cfg$regulation_report_dir,
    cfg$regulation_table_dir,
    cfg$regulation_figure_dir,
    cfg$regulation_checkpoint_dir,
    cfg$scenic_table_dir,
    cfg$scenic_figure_dir,
    cfg$scenic_checkpoint_dir,
    dirname(cfg$module_08a_manifest_path),
    dirname(cfg$module_08b_manifest_path),
    dirname(cfg$module_08c_manifest_path),
    dirname(cfg$module_08d_manifest_path),
    dirname(cfg$scenic_export_index_tsv),
    dirname(cfg$scenic_export_mapping_summary_tsv),
    dirname(cfg$scenic_export_triage_tsv),
    dirname(cfg$scenic_grn_index_tsv),
    dirname(cfg$scenic_regulons_index_tsv),
    dirname(cfg$scenic_downstream_index_tsv)
  )
  invisible(lapply(ensure_dirs[nzchar(ensure_dirs)], ensure_dir))
}

safe_id_08 <- function(value) {
  safe_id_07(value)
}

scenic_export_paths_08 <- function(cfg, layer_id) {
  safe_layer <- safe_id_08(layer_id)
  list(
    input_dir = file.path(cfg$scenic_input_dir, safe_layer),
    expr_mat_human_csv = file.path(cfg$scenic_input_dir, safe_layer, "expr_mat_human.csv"),
    ortholog_csv = file.path(cfg$scenic_input_dir, safe_layer, "chicken_to_human_orthologs.csv")
  )
}

scenic_grn_paths_08 <- function(cfg, layer_id) {
  safe_layer <- safe_id_08(layer_id)
  list(
    output_dir = file.path(cfg$scenic_output_dir, safe_layer),
    adjacencies_tsv = file.path(cfg$scenic_output_dir, safe_layer, "adjacencies.tsv")
  )
}

scenic_regulon_paths_08 <- function(cfg, layer_id) {
  safe_layer <- safe_id_08(layer_id)
  list(
    output_dir = file.path(cfg$scenic_output_dir, safe_layer),
    regulons_csv = file.path(cfg$scenic_output_dir, safe_layer, "regulons.csv"),
    regulon_targets_long_csv = file.path(cfg$scenic_output_dir, safe_layer, "regulon_targets_long.csv"),
    regulons_gmt = file.path(cfg$scenic_output_dir, safe_layer, "regulons.gmt"),
    motif_enrichment_tsv = file.path(cfg$scenic_output_dir, safe_layer, "motif_enrichment.tsv"),
    motif_enrichment_rds = file.path(cfg$scenic_output_dir, safe_layer, "motif_enrichment.rds"),
    auc_matrix_csv = file.path(cfg$scenic_output_dir, safe_layer, "auc_matrix.csv"),
    scenic_core_outputs_rds = file.path(cfg$scenic_output_dir, safe_layer, "scenic_core_outputs.rds")
  )
}

scenic_downstream_paths_08 <- function(cfg, layer_id) {
  safe_layer <- safe_id_08(layer_id)
  list(
    table_dir = file.path(cfg$scenic_table_dir, safe_layer),
    figure_dir = file.path(cfg$scenic_figure_dir, safe_layer),
    checkpoint_dir = file.path(cfg$scenic_checkpoint_dir, safe_layer),
    rss_matrix_csv = file.path(cfg$scenic_table_dir, safe_layer, "rss_matrix.csv"),
    top_regulons_by_celltype_csv = file.path(cfg$scenic_table_dir, safe_layer, "top_regulons_by_celltype.csv"),
    csi_matrix_csv = file.path(cfg$scenic_table_dir, safe_layer, "csi_matrix.csv"),
    csi_source_note_txt = file.path(cfg$scenic_table_dir, safe_layer, "csi_original_source_used.txt"),
    scenic_integrated_object_rds = file.path(cfg$scenic_checkpoint_dir, safe_layer, "scenic_integrated_object.rds"),
    figure_6a_png = file.path(cfg$scenic_figure_dir, safe_layer, "Figure_6A_Regulon_Heatmap.png"),
    figure_6b_png = file.path(cfg$scenic_figure_dir, safe_layer, "Figure_6B_RSS_Ranks.png"),
    figure_6b_pdf = file.path(cfg$scenic_figure_dir, safe_layer, "Figure_6B_RSS_Ranks.pdf"),
    figure_6c_png = file.path(cfg$scenic_figure_dir, safe_layer, "Figure_6C_CSI_Clustering_Heatmap.png"),
    figure_6d_png = file.path(cfg$scenic_figure_dir, safe_layer, "Figure_6D_CSI_Module_Activity_Heatmap.png")
  )
}
