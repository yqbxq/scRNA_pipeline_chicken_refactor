#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

source_utf8 <- function(path) source(path, encoding = "UTF-8")

source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "config.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_02.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_03.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_04.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_05.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_06.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "deg_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "enrichment_utils.R"))

load_required_packages(c("dplyr", "tibble", "jsonlite", "ggplot2"))

cfg <- get_single_script_config_06()
module_name <- "06c_enrichment_eda"
prepare_dirs_06(cfg)

paths <- enrichment_report_paths_06(cfg)
go_manifest <- read_tsv_optional(cfg$go_enrichment_manifest_tsv)
kegg_manifest <- read_tsv_optional(cfg$kegg_enrichment_manifest_tsv)
enrichment_targets <- read_tsv_optional(cfg$enrichment_targets_sheet)
all_manifest <- dplyr::bind_rows(go_manifest, kegg_manifest)
if (nrow(all_manifest) == 0) {
  all_manifest <- empty_enrichment_manifest_06()
}

go_results <- read_enrichment_tables_from_manifest_06(go_manifest)
kegg_results <- read_enrichment_tables_from_manifest_06(kegg_manifest)
combined_results <- dplyr::bind_rows(go_results, kegg_results)

numeric_col <- function(df, col) {
  if (!col %in% colnames(df)) {
    return(rep(NA_real_, nrow(df)))
  }
  suppressWarnings(as.numeric(df[[col]]))
}

build_summary_06 <- function(manifest_df) {
  if (nrow(manifest_df) == 0) {
    return(empty_df_05(c(
      "layer_id", "comparison_id", "cluster_id", "source_species", "deg_source",
      "go_bp_significant_n", "go_cc_significant_n", "go_mf_significant_n",
      "kegg_significant_n", "total_significant_n", "median_mapping_rate"
    )))
  }
  for (col in c("layer_id", "comparison_id", "cluster_id", "source_species", "deg_source", "analysis_type", "ontology")) {
    if (!col %in% colnames(manifest_df)) {
      manifest_df[[col]] <- ""
    }
  }
  manifest_df$significant_term_n_num <- numeric_col(manifest_df, "significant_term_n")
  manifest_df$mapping_rate_num <- numeric_col(manifest_df, "mapping_rate")
  keys <- unique(manifest_df[, c("layer_id", "comparison_id", "cluster_id", "source_species"), drop = FALSE])
  rows <- lapply(seq_len(nrow(keys)), function(idx) {
    key <- keys[idx, , drop = FALSE]
    hit <- manifest_df[
      manifest_df$layer_id == key$layer_id[[1]] &
        manifest_df$comparison_id == key$comparison_id[[1]] &
        manifest_df$cluster_id == key$cluster_id[[1]] &
        manifest_df$source_species == key$source_species[[1]],
      ,
      drop = FALSE
    ]
    count_terms <- function(analysis_type, ontology) {
      sum(hit$significant_term_n_num[hit$analysis_type == analysis_type & hit$ontology == ontology], na.rm = TRUE)
    }
    mapping <- hit$mapping_rate_num[is.finite(hit$mapping_rate_num)]
    data.frame(
      layer_id = key$layer_id[[1]],
      comparison_id = key$comparison_id[[1]],
      cluster_id = key$cluster_id[[1]],
      source_species = key$source_species[[1]],
      deg_source = paste(sort(unique(hit$deg_source[nzchar(hit$deg_source)])), collapse = ";"),
      go_bp_significant_n = count_terms("go", "BP"),
      go_cc_significant_n = count_terms("go", "CC"),
      go_mf_significant_n = count_terms("go", "MF"),
      kegg_significant_n = count_terms("kegg", "KEGG"),
      total_significant_n = sum(hit$significant_term_n_num, na.rm = TRUE),
      median_mapping_rate = if (length(mapping) == 0) NA_real_ else stats::median(mapping),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

summary_df <- build_summary_06(all_manifest)
write_tsv_local(summary_df, paths$summary_tsv)

if (nrow(combined_results) > 0) {
  combined_results$padj_num <- numeric_col(combined_results, "p.adjust")
  combined_results <- combined_results[is.finite(combined_results$padj_num) & combined_results$padj_num <= cfg$enrichment_pvalue_cutoff, , drop = FALSE]
}

shared_pathways <- build_shared_pathways_06(combined_results, min_cluster_n = 2L)
write_tsv_local(shared_pathways, paths$shared_pathways_tsv)

go_bp_results <- combined_results
if (nrow(go_bp_results) > 0 && all(c("analysis_type", "ontology", "gene_direction") %in% colnames(go_bp_results))) {
  go_bp_results <- go_bp_results[go_bp_results$analysis_type == "go" & go_bp_results$ontology == "BP" & go_bp_results$gene_direction == "all", , drop = FALSE]
}
kegg_heatmap_results <- combined_results
if (nrow(kegg_heatmap_results) > 0 && all(c("analysis_type", "gene_direction") %in% colnames(kegg_heatmap_results))) {
  kegg_heatmap_results <- kegg_heatmap_results[kegg_heatmap_results$analysis_type == "kegg" & kegg_heatmap_results$gene_direction == "all", , drop = FALSE]
}

plot_cross_cluster_heatmap(
  go_bp_results,
  value_col = "p.adjust",
  output_png = paths$cross_cluster_go_heatmap_png,
  title = "GO BP enrichment across clusters",
  top_n = 30L
)
plot_cross_cluster_heatmap(
  kegg_heatmap_results,
  value_col = "p.adjust",
  output_png = paths$cross_cluster_kegg_heatmap_png,
  title = "KEGG enrichment across clusters",
  top_n = 20L
)

low_mapping <- all_manifest[
  is.finite(numeric_col(all_manifest, "mapping_rate")) & numeric_col(all_manifest, "mapping_rate") < 0.3,
  ,
  drop = FALSE
]
no_significant <- summary_df[suppressWarnings(as.numeric(summary_df$total_significant_n)) == 0, , drop = FALSE]
status_summary <- if ("status" %in% colnames(all_manifest)) as.data.frame(table(all_manifest$status), stringsAsFactors = FALSE) else data.frame(stringsAsFactors = FALSE)
colnames(status_summary) <- if (ncol(status_summary) == 2) c("status", "row_n") else colnames(status_summary)
usage_summary <- if (nrow(enrichment_targets) > 0 && all(c("enrichment_usage", "comparison_id") %in% colnames(enrichment_targets))) {
  enrichment_targets %>% dplyr::count(enrichment_usage, name = "target_n")
} else {
  empty_df_05(c("enrichment_usage", "target_n"))
}

deg_status_matrix <- read_tsv_optional(deg_report_paths_05(cfg)$status_matrix_tsv)
exploratory_rows <- deg_status_matrix
if (nrow(exploratory_rows) > 0) {
  for (col in c("pseudobulk_status", "composition_status")) {
    if (!col %in% colnames(exploratory_rows)) {
      exploratory_rows[[col]] <- ""
    }
  }
  exploratory_rows <- exploratory_rows[
    exploratory_rows$pseudobulk_status %in% c("exploratory_only", "exploratory_forced") |
      exploratory_rows$composition_status %in% c("exploratory_only", "exploratory_forced"),
    ,
    drop = FALSE
  ]
}

report_lines <- c(
  "# 06c Enrichment EDA",
  "",
  sprintf("- GO manifest: `%s`", cfg$go_enrichment_manifest_tsv),
  sprintf("- KEGG manifest: `%s`", cfg$kegg_enrichment_manifest_tsv),
  sprintf("- enrichment summary: `%s`", paths$summary_tsv),
  sprintf("- shared pathways: `%s`", paths$shared_pathways_tsv),
  sprintf("- GO BP heatmap: `%s`", paths$cross_cluster_go_heatmap_png),
  sprintf("- KEGG heatmap: `%s`", paths$cross_cluster_kegg_heatmap_png),
  "",
  "## Status Summary",
  render_markdown_table_local(status_summary),
  "",
  "## Enrichment Usage Targets",
  render_markdown_table_local(usage_summary),
  "",
  "## Cluster Summary",
  render_markdown_table_local(head(summary_df, 100)),
  "",
  "## Shared Pathways",
  if (nrow(shared_pathways) == 0) "No pathway was significant in at least two cluster contexts." else render_markdown_table_local(head(shared_pathways, 50)),
  "",
  "## Low Mapping Warnings",
  if (nrow(low_mapping) == 0) "No enrichment rows had ID mapping below 30%." else render_markdown_table_local(head(low_mapping[, intersect(c("layer_id", "comparison_id", "cluster_id", "gene_direction", "analysis_type", "ontology", "source_species", "mapping_rate", "status"), colnames(low_mapping)), drop = FALSE], 100)),
  "",
  "## No Significant Enrichment",
  if (nrow(no_significant) == 0) "Every cluster/species context had at least one significant enrichment row." else render_markdown_table_local(head(no_significant, 100)),
  "",
  "## Exploratory DEG Context",
  if (nrow(exploratory_rows) == 0) "No exploratory-only DEG contexts were recorded in the 05d status matrix." else render_markdown_table_local(head(exploratory_rows[, intersect(c("layer_id", "comparison_id", "pseudobulk_status", "marker_results_tsv"), colnames(exploratory_rows)), drop = FALSE], 100))
)
write_markdown_local(report_lines, paths$report_md)

if (file.exists(cfg$module_06c_manifest_path)) {
  unlink(cfg$module_06c_manifest_path)
}
write_manifest_local(
  manifest_path = cfg$module_06c_manifest_path,
  new_outputs = list(
    report = build_output_entry(paths$report_md, "md", module_name, "enrichment EDA report", base_dir = cfg$project_root),
    enrichment_summary_tsv = build_output_entry(paths$summary_tsv, "tsv", module_name, "one row per layer/comparison/cluster/species enrichment summary", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df)),
    cross_cluster_go_heatmap_png = build_output_entry(paths$cross_cluster_go_heatmap_png, "png", module_name, "GO BP cross-cluster heatmap", base_dir = cfg$project_root),
    cross_cluster_kegg_heatmap_png = build_output_entry(paths$cross_cluster_kegg_heatmap_png, "png", module_name, "KEGG cross-cluster heatmap", base_dir = cfg$project_root),
    shared_pathways_tsv = build_output_entry(paths$shared_pathways_tsv, "tsv", module_name, "pathways significant in at least two cluster contexts", base_dir = cfg$project_root, schema = infer_schema_from_df(shared_pathways))
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    go_enrichment_manifest_tsv = cfg$go_enrichment_manifest_tsv,
    kegg_enrichment_manifest_tsv = cfg$kegg_enrichment_manifest_tsv,
    deg_status_matrix_tsv = deg_report_paths_05(cfg)$status_matrix_tsv
  ),
  version = cfg$module_version,
  depends_on = list(
    module_06a = cfg$module_06a_manifest_path,
    module_06b = cfg$module_06b_manifest_path,
    module_05d = cfg$module_05d_manifest_path
  )
)

message("06c completed. report: ", paths$report_md)
