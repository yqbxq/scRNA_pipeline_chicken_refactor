env_numeric_06 <- function(name, default) {
  value <- suppressWarnings(as.numeric(env_or_default_03(name, as.character(default))))
  if (is.na(value)) default else value
}

env_integer_06 <- function(name, default) {
  value <- suppressWarnings(as.integer(env_or_default_03(name, as.character(default))))
  if (is.na(value)) default else value
}

get_single_script_config_06 <- function() {
  base <- get_single_script_config_05()

  base$module_06a_manifest_path <- file.path(base$manifest_dir, "06a_go_enrichment", "_manifest.json")
  base$module_06b_manifest_path <- file.path(base$manifest_dir, "06b_kegg_enrichment", "_manifest.json")
  base$module_06c_manifest_path <- file.path(base$manifest_dir, "06c_enrichment_eda", "_manifest.json")

  base$enrichment_report_dir <- env_or_default_03("ENRICHMENT_REPORT_DIR", file.path(base$eda_report_dir, "enrichment"))
  base$enrichment_table_dir <- file.path(base$table_dir, "enrichment")
  base$go_table_dir <- file.path(base$enrichment_table_dir, "go")
  base$kegg_table_dir <- file.path(base$enrichment_table_dir, "kegg")
  base$enrichment_background_table_dir <- file.path(base$enrichment_table_dir, "background")
  base$enrichment_figure_dir <- file.path(base$figure_dir, "enrichment")
  base$go_figure_dir <- file.path(base$enrichment_figure_dir, "go")
  base$kegg_figure_dir <- file.path(base$enrichment_figure_dir, "kegg")

  base$enrichment_targets_sheet <- env_or_default_03("ENRICHMENT_TARGETS_SHEET", file.path(base$metadata_dir, "enrichment_targets.tsv"))
  base$go_enrichment_manifest_tsv <- file.path(base$go_table_dir, "go_enrichment_manifest.tsv")
  base$kegg_enrichment_manifest_tsv <- file.path(base$kegg_table_dir, "kegg_enrichment_manifest.tsv")
  base$enrichment_summary_tsv <- file.path(base$enrichment_table_dir, "enrichment_summary.tsv")
  base$shared_pathways_tsv <- file.path(base$enrichment_report_dir, "shared_pathways.tsv")
  base$enrichment_report_md <- file.path(base$enrichment_report_dir, "report.md")
  base$cross_cluster_go_heatmap_png <- file.path(base$enrichment_figure_dir, "cross_cluster_go_bp_heatmap.png")
  base$cross_cluster_kegg_heatmap_png <- file.path(base$enrichment_figure_dir, "cross_cluster_kegg_heatmap.png")

  base$enrichment_species_strategy <- env_or_default_03("ENRICHMENT_SPECIES_STRATEGY", "chicken_primary")
  base$enrichment_pvalue_cutoff <- env_numeric_06("ENRICHMENT_PVALUE_CUTOFF", 0.05)
  base$enrichment_qvalue_cutoff <- env_numeric_06("ENRICHMENT_QVALUE_CUTOFF", 0.2)
  base$enrichment_top_n_genes <- env_integer_06("ENRICHMENT_TOP_N_GENES", 100L)
  base$enrichment_min_input_genes <- env_integer_06("ENRICHMENT_MIN_INPUT_GENES", 5L)
  base$enrichment_min_gs_size <- env_integer_06("ENRICHMENT_MIN_GS_SIZE", 10L)
  base$enrichment_max_gs_size <- env_integer_06("ENRICHMENT_MAX_GS_SIZE", 500L)
  base$enrichment_kegg_timeout_sec <- env_integer_06("ENRICHMENT_KEGG_TIMEOUT_SEC", 60L)
  base$module_version <- env_or_default_03("MODULE_06_VERSION", "1.0")
  base
}

prepare_dirs_06 <- function(cfg) {
  ensure_dirs <- c(
    cfg$enrichment_report_dir,
    cfg$enrichment_table_dir,
    cfg$go_table_dir,
    cfg$kegg_table_dir,
    cfg$enrichment_background_table_dir,
    cfg$enrichment_figure_dir,
    cfg$go_figure_dir,
    cfg$kegg_figure_dir,
    dirname(cfg$module_06a_manifest_path),
    dirname(cfg$module_06b_manifest_path),
    dirname(cfg$module_06c_manifest_path),
    dirname(cfg$enrichment_summary_tsv),
    dirname(cfg$shared_pathways_tsv),
    dirname(cfg$cross_cluster_go_heatmap_png),
    dirname(cfg$cross_cluster_kegg_heatmap_png)
  )
  invisible(lapply(ensure_dirs[nzchar(ensure_dirs)], ensure_dir))
}

safe_id_06 <- function(value) {
  safe_id_05(value)
}

enrichment_paths_06 <- function(cfg, layer_id, comparison_id, cluster_id, subtype, source_species, gene_direction = "all") {
  safe_layer <- safe_id_06(layer_id)
  safe_cmp <- safe_id_06(comparison_id)
  safe_cluster <- safe_id_06(cluster_id)
  safe_subtype <- safe_id_06(subtype)
  safe_species <- safe_id_06(source_species)
  safe_direction <- safe_id_06(gene_direction)
  root <- if (startsWith(safe_subtype, "go_")) cfg$go_table_dir else cfg$kegg_table_dir
  fig_root <- if (startsWith(safe_subtype, "go_")) cfg$go_figure_dir else cfg$kegg_figure_dir
  prefix <- sprintf("%s_%s_%s", safe_subtype, safe_direction, safe_species)

  list(
    table_dir = file.path(root, safe_layer, safe_cmp, safe_cluster),
    figure_dir = file.path(fig_root, safe_layer, safe_cmp, safe_cluster),
    enrichment_tsv = file.path(root, safe_layer, safe_cmp, safe_cluster, paste0(prefix, ".tsv")),
    enrichment_significant_tsv = file.path(root, safe_layer, safe_cmp, safe_cluster, paste0(prefix, "_significant.tsv")),
    dotplot_png = file.path(fig_root, safe_layer, safe_cmp, safe_cluster, paste0(prefix, "_dotplot.png")),
    barplot_png = file.path(fig_root, safe_layer, safe_cmp, safe_cluster, paste0(prefix, "_barplot.png"))
  )
}

enrichment_report_paths_06 <- function(cfg) {
  list(
    report_md = cfg$enrichment_report_md,
    summary_tsv = cfg$enrichment_summary_tsv,
    shared_pathways_tsv = cfg$shared_pathways_tsv,
    cross_cluster_go_heatmap_png = cfg$cross_cluster_go_heatmap_png,
    cross_cluster_kegg_heatmap_png = cfg$cross_cluster_kegg_heatmap_png
  )
}
