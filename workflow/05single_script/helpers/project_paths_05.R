env_numeric_05 <- function(name, default) {
  value <- suppressWarnings(as.numeric(env_or_default_03(name, as.character(default))))
  if (is.na(value)) default else value
}

env_integer_05 <- function(name, default) {
  value <- suppressWarnings(as.integer(env_or_default_03(name, as.character(default))))
  if (is.na(value)) default else value
}

get_single_script_config_05 <- function() {
  base <- get_single_script_config_04()

  base$module_05a_manifest_path <- file.path(base$manifest_dir, "05a_marker_discovery", "_manifest.json")
  base$module_05b_manifest_path <- file.path(base$manifest_dir, "05b_pseudobulk_de", "_manifest.json")
  base$module_05c_manifest_path <- file.path(base$manifest_dir, "05c_composition", "_manifest.json")
  base$module_05d_manifest_path <- file.path(base$manifest_dir, "05d_deg_eda", "_manifest.json")

  base$deg_report_dir <- env_or_default_03("DEG_REPORT_DIR", file.path(base$eda_report_dir, "deg"))
  base$deg_table_dir <- file.path(base$table_dir, "deg")
  base$marker_discovery_table_dir <- file.path(base$table_dir, "marker_discovery")
  base$pseudobulk_table_dir <- file.path(base$table_dir, "pseudobulk_ds")
  base$composition_table_dir <- file.path(base$table_dir, "composition")
  base$marker_discovery_report_dir <- file.path(base$eda_report_dir, "marker_discovery")
  base$pseudobulk_report_dir <- file.path(base$eda_report_dir, "pseudobulk_ds")
  base$composition_report_dir <- file.path(base$eda_report_dir, "composition")

  base$deg_default_min_cells_per_group <- env_integer_05("DEG_MIN_CELLS_PER_GROUP", 3L)
  base$deg_default_logfc_threshold <- env_numeric_05("DEG_LOGFC_THRESHOLD", 0)
  base$deg_ident_1 <- env_or_default_03("DEG_IDENT_1", "")
  base$deg_ident_2 <- env_or_default_03("DEG_IDENT_2", "")
  base$deg_alpha <- env_numeric_05("DEG_ALPHA", 0.05)
  base$gene_program_targets_sheet <- env_or_default_03("GENE_PROGRAM_TARGETS_SHEET", file.path(base$metadata_dir, "gene_program_targets.tsv"))
  base$gene_program_registry_tsv <- env_or_default_03("GENE_PROGRAM_REGISTRY_TSV", file.path(base$deg_table_dir, "gene_program_registry.tsv"))
  base$module_version <- env_or_default_03("MODULE_05_VERSION", "1.0")
  base
}

prepare_dirs_05 <- function(cfg) {
  ensure_dirs <- c(
    cfg$deg_report_dir,
    cfg$deg_table_dir,
    cfg$marker_discovery_table_dir,
    cfg$pseudobulk_table_dir,
    cfg$composition_table_dir,
    cfg$marker_discovery_report_dir,
    cfg$pseudobulk_report_dir,
    cfg$composition_report_dir,
    dirname(cfg$module_05a_manifest_path),
    dirname(cfg$module_05b_manifest_path),
    dirname(cfg$module_05c_manifest_path),
    dirname(cfg$module_05d_manifest_path)
  )
  invisible(lapply(ensure_dirs[nzchar(ensure_dirs)], ensure_dir))
}

safe_id_05 <- function(value) {
  sanitize_layer_id_local(value)
}

marker_discovery_paths_05 <- function(cfg, layer_id, comparison_id = NULL) {
  safe_layer <- safe_id_05(layer_id)
  if (is.null(comparison_id) || !nzchar(comparison_id)) {
    return(list(
      table_dir = file.path(cfg$marker_discovery_table_dir, "layers", safe_layer),
      cluster_markers_tsv = file.path(cfg$marker_discovery_table_dir, "layers", safe_layer, "cluster_markers.tsv")
    ))
  }

  safe_cmp <- safe_id_05(comparison_id)
  list(
    table_dir = file.path(cfg$marker_discovery_table_dir, "comparisons", safe_cmp, safe_layer),
    exploratory_tsv = file.path(cfg$marker_discovery_table_dir, "comparisons", safe_cmp, safe_layer, "exploratory_cell_level_findmarkers.tsv"),
    exploratory_summary_tsv = file.path(cfg$marker_discovery_table_dir, "comparisons", safe_cmp, safe_layer, "exploratory_summary.tsv")
  )
}

pseudobulk_paths_05 <- function(cfg, layer_id, comparison_id = NULL) {
  safe_layer <- safe_id_05(layer_id)
  if (is.null(comparison_id) || !nzchar(comparison_id)) {
    return(list(
      table_dir = file.path(cfg$pseudobulk_table_dir, "layers", safe_layer),
      aggregation_rds = file.path(cfg$pseudobulk_table_dir, "layers", safe_layer, "pseudobulk_counts.rds"),
      aggregation_tsv = file.path(cfg$pseudobulk_table_dir, "layers", safe_layer, "pseudobulk_column_metadata.tsv"),
      aggregation_manifest_tsv = file.path(cfg$pseudobulk_table_dir, "layers", safe_layer, "aggregation_manifest.tsv")
    ))
  }

  safe_cmp <- safe_id_05(comparison_id)
  list(
    table_dir = file.path(cfg$pseudobulk_table_dir, "comparisons", safe_cmp, safe_layer),
    aggregation_rds = file.path(cfg$pseudobulk_table_dir, "comparisons", safe_cmp, safe_layer, "pseudobulk_counts.rds"),
    aggregation_tsv = file.path(cfg$pseudobulk_table_dir, "comparisons", safe_cmp, safe_layer, "pseudobulk_column_metadata.tsv"),
    ds_results_tsv = file.path(cfg$pseudobulk_table_dir, "comparisons", safe_cmp, safe_layer, "formal_pseudobulk_ds.tsv"),
    gate_summary_tsv = file.path(cfg$pseudobulk_table_dir, "comparisons", safe_cmp, safe_layer, "replicate_gate_summary.tsv"),
    status_tsv = file.path(cfg$pseudobulk_table_dir, "comparisons", safe_cmp, safe_layer, "ds_status.tsv"),
    result_rds = file.path(cfg$pseudobulk_table_dir, "comparisons", safe_cmp, safe_layer, "formal_pseudobulk_ds.rds")
  )
}

composition_paths_05 <- function(cfg, layer_id, comparison_id = NULL) {
  safe_layer <- safe_id_05(layer_id)
  if (is.null(comparison_id) || !nzchar(comparison_id)) {
    return(list(
      table_dir = file.path(cfg$composition_table_dir, "layers", safe_layer),
      proportion_tsv = file.path(cfg$composition_table_dir, "layers", safe_layer, "sample_level_proportions.tsv"),
      proportion_plot_png = file.path(cfg$composition_report_dir, "layers", safe_layer, "sample_level_proportions.png")
    ))
  }

  safe_cmp <- safe_id_05(comparison_id)
  list(
    table_dir = file.path(cfg$composition_table_dir, "comparisons", safe_cmp, safe_layer),
    proportion_tsv = file.path(cfg$composition_table_dir, "comparisons", safe_cmp, safe_layer, "sample_level_proportions.tsv"),
    formal_results_tsv = file.path(cfg$composition_table_dir, "comparisons", safe_cmp, safe_layer, "formal_propeller.tsv"),
    gate_summary_tsv = file.path(cfg$composition_table_dir, "comparisons", safe_cmp, safe_layer, "replicate_gate_summary.tsv"),
    status_tsv = file.path(cfg$composition_table_dir, "comparisons", safe_cmp, safe_layer, "composition_status.tsv")
  )
}

deg_report_paths_05 <- function(cfg) {
  list(
    report_md = file.path(cfg$deg_report_dir, "report.md"),
    status_matrix_tsv = file.path(cfg$deg_table_dir, "deg_status_matrix.tsv"),
    summary_tsv = file.path(cfg$deg_table_dir, "deg_summary.tsv"),
    gene_program_registry_tsv = cfg$gene_program_registry_tsv %||% file.path(cfg$deg_table_dir, "gene_program_registry.tsv")
  )
}
