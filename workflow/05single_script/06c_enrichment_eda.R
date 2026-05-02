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
gene_program_registry <- read_tsv_optional(cfg$gene_program_registry_tsv %||% deg_report_paths_05(cfg)$gene_program_registry_tsv)
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

global_context_warning_06c <- c(
  "D05 is a global panorama stage signature.",
  "It may mix cell-type composition changes, cell-type-specific transcriptional changes, and sampling/capture effects.",
  "Use D05 as contextual background only.",
  "Do not use D05 as receiver DEG for NicheNet or as core cell-type-specific mechanism evidence.",
  "Core mechanism interpretation should prioritize D01/D02/D03 condition_within_type results.",
  "",
  "D05 是 panorama 全局 syf vs f5 背景 signature。",
  "它可能混合细胞组成变化、同一细胞类型内部转录变化以及采样/捕获因素。",
  "D05 只能作为全局背景参考，不应作为 NicheNet receiver DEG，也不应作为核心 cell-type-specific 机制证据。",
  "核心机制解释应优先使用 D01/D02/D03。"
)

ensure_cols_06c <- function(df, cols) {
  for (col in cols) {
    if (!col %in% colnames(df)) {
      df[[col]] <- character(nrow(df))
    }
  }
  df
}

truthy_06c <- function(x) {
  tolower(normalize_scalar_value(x, "no")) %in% c("yes", "true", "1", "on")
}

empty_enrichment_metadata_06c <- function() {
  empty_df_05(c(
    "section", "metadata_join_status", "metadata_join_reason",
    "target_id", "source_question_id", "comparison_id", "layer_id",
    "database", "analysis_mode", "gene_program_role",
    "produces_gene_program", "annotation_only", "qc_only",
    "global_context_only", "nichenet_eligible", "nichenet_usage",
    "enrichment_eligible", "enrichment_usage", "preferred_for_downstream",
    "result_level", "formal_status", "result_status", "ineligible_reason"
  ))
}

assign_enrichment_section_06c <- function(row) {
  analysis_mode <- normalize_scalar_value(row$analysis_mode[[1]])
  gene_program_role <- normalize_scalar_value(row$gene_program_role[[1]])
  enrichment_usage <- normalize_scalar_value(row$enrichment_usage[[1]])
  enrichment_eligible <- normalize_scalar_value(row$enrichment_eligible[[1]])
  if (truthy_06c(row$qc_only[[1]]) ||
      identical(normalize_scalar_value(row$produces_gene_program[[1]]), "no") ||
      identical(enrichment_eligible, "no")) {
    return("QC / Non-gene-program")
  }
  if (identical(analysis_mode, "global_context") ||
      identical(gene_program_role, "global_context") ||
      truthy_06c(row$global_context_only[[1]]) ||
      identical(enrichment_usage, "global_context_enrichment")) {
    return("Global Context Enrichment")
  }
  if (identical(analysis_mode, "condition_within_type") ||
      identical(gene_program_role, "condition_deg") ||
      identical(enrichment_usage, "mechanism_enrichment")) {
    return("Core Mechanism Enrichment")
  }
  if (gene_program_role %in% c("receiver_marker", "subtype_pairwise_deg") ||
      enrichment_usage %in% c("identity_baseline_enrichment", "subtype_pairwise_enrichment")) {
    return("Subtype / Identity Enrichment")
  }
  "Unclassified Enrichment"
}

prepare_enrichment_metadata_06c <- function(targets, registry) {
  target_cols <- c(
    "target_id", "source_question_id", "comparison_id", "layer_scope",
    "analysis_mode", "gene_program_role", "database",
    "enrichment_eligible", "enrichment_usage", "enabled", "notes"
  )
  registry_cols <- c(
    "comparison_id", "source_question_id", "layer_id", "analysis_mode",
    "gene_program_role", "produces_gene_program", "annotation_only",
    "qc_only", "global_context_only", "nichenet_eligible",
    "nichenet_usage", "enrichment_eligible", "enrichment_usage",
    "preferred_for_downstream", "result_level", "formal_status",
    "result_status", "ineligible_reason"
  )
  targets <- ensure_cols_06c(targets, target_cols)
  registry <- ensure_cols_06c(registry, registry_cols)
  if (nrow(targets) == 0) {
    return(empty_enrichment_metadata_06c())
  }
  targets$enabled <- tolower(vapply(targets$enabled, normalize_scalar_value, character(1), default = "yes"))
  targets <- targets[targets$enabled %in% c("yes", "true", "1", "on") & nzchar(targets$comparison_id), , drop = FALSE]
  if (nrow(targets) == 0) {
    return(empty_enrichment_metadata_06c())
  }
  rows <- lapply(seq_len(nrow(targets)), function(idx) {
    target <- targets[idx, , drop = FALSE]
    comparison_id <- normalize_scalar_value(target$comparison_id[[1]])
    layer_id <- normalize_scalar_value(target$layer_scope[[1]])
    hit <- registry[registry$comparison_id == comparison_id, , drop = FALSE]
    if (nrow(hit) > 0) {
      exact <- hit[nzchar(hit$layer_id) & hit$layer_id == layer_id, , drop = FALSE]
      if (nrow(exact) > 0) hit <- exact
      hit <- hit[1, , drop = FALSE]
    }
    registry_hit <- nrow(hit) > 0
    value <- function(col, fallback = "") {
      registry_value <- if (registry_hit && col %in% colnames(hit)) normalize_scalar_value(hit[[col]][[1]]) else ""
      if (nzchar(registry_value)) {
        return(registry_value)
      }
      if (col %in% colnames(target)) {
        return(normalize_scalar_value(target[[col]][[1]], fallback))
      }
      fallback
    }
    metadata_join_status <- if (registry_hit &&
        nzchar(value("analysis_mode")) &&
        nzchar(value("gene_program_role")) &&
        nzchar(value("enrichment_usage"))) "ok" else "incomplete"
    row <- data.frame(
      section = "",
      metadata_join_status = metadata_join_status,
      metadata_join_reason = if (identical(metadata_join_status, "ok")) "" else "missing or incomplete gene_program_registry.tsv metadata",
      target_id = normalize_scalar_value(target$target_id[[1]]),
      source_question_id = value("source_question_id"),
      comparison_id = comparison_id,
      layer_id = layer_id,
      database = normalize_scalar_value(target$database[[1]]),
      analysis_mode = value("analysis_mode"),
      gene_program_role = value("gene_program_role"),
      produces_gene_program = value("produces_gene_program"),
      annotation_only = value("annotation_only"),
      qc_only = value("qc_only"),
      global_context_only = value("global_context_only"),
      nichenet_eligible = value("nichenet_eligible"),
      nichenet_usage = value("nichenet_usage"),
      enrichment_eligible = value("enrichment_eligible"),
      enrichment_usage = value("enrichment_usage"),
      preferred_for_downstream = value("preferred_for_downstream"),
      result_level = value("result_level"),
      formal_status = value("formal_status"),
      result_status = value("result_status"),
      ineligible_reason = value("ineligible_reason"),
      stringsAsFactors = FALSE
    )
    row$section <- assign_enrichment_section_06c(row)
    row
  })
  dplyr::bind_rows(rows)
}

metadata_join_status_06c <- function(section_metadata) {
  if (nrow(section_metadata) == 0) {
    return("no_targets")
  }
  if (all(section_metadata$metadata_join_status == "ok")) {
    return("ok")
  }
  "incomplete"
}

build_section_summary_06c <- function(section_metadata, manifest_df) {
  if (nrow(section_metadata) == 0) {
    return(empty_df_05(c(
      "section", "target_id", "source_question_id", "comparison_id",
      "layer_id", "database", "analysis_mode", "gene_program_role",
      "enrichment_usage", "enrichment_eligible",
      "preferred_for_downstream", "metadata_join_status",
      "manifest_row_n", "significant_term_n", "manifest_status"
    )))
  }
  manifest_df <- ensure_cols_06c(manifest_df, c("target_id", "comparison_id", "layer_id", "status", "significant_term_n"))
  rows <- lapply(seq_len(nrow(section_metadata)), function(idx) {
    row <- section_metadata[idx, , drop = FALSE]
    hit <- manifest_df[manifest_df$target_id == row$target_id[[1]], , drop = FALSE]
    if (nrow(hit) == 0) {
      hit <- manifest_df[
        manifest_df$comparison_id == row$comparison_id[[1]] &
          manifest_df$layer_id == row$layer_id[[1]],
        ,
        drop = FALSE
      ]
    }
    significant_n <- suppressWarnings(as.numeric(hit$significant_term_n))
    data.frame(
      section = row$section[[1]],
      target_id = row$target_id[[1]],
      source_question_id = row$source_question_id[[1]],
      comparison_id = row$comparison_id[[1]],
      layer_id = row$layer_id[[1]],
      database = row$database[[1]],
      analysis_mode = row$analysis_mode[[1]],
      gene_program_role = row$gene_program_role[[1]],
      enrichment_usage = row$enrichment_usage[[1]],
      enrichment_eligible = row$enrichment_eligible[[1]],
      preferred_for_downstream = row$preferred_for_downstream[[1]],
      metadata_join_status = row$metadata_join_status[[1]],
      manifest_row_n = nrow(hit),
      significant_term_n = sum(significant_n, na.rm = TRUE),
      manifest_status = paste(sort(unique(hit$status[nzchar(hit$status)])), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

section_report_lines_06c <- function(section_summary, section_name, include_global_warning = FALSE) {
  section_rows <- section_summary[section_summary$section == section_name, , drop = FALSE]
  lines <- c(sprintf("## %s", section_name))
  if (include_global_warning) {
    lines <- c(lines, "", global_context_warning_06c)
  }
  if (nrow(section_rows) == 0) {
    return(c(lines, "", "No targets in this section."))
  }
  show_cols <- intersect(
    c(
      "source_question_id", "comparison_id", "database", "analysis_mode",
      "gene_program_role", "enrichment_usage", "enrichment_eligible",
      "preferred_for_downstream", "metadata_join_status",
      "manifest_row_n", "significant_term_n", "manifest_status"
    ),
    colnames(section_rows)
  )
  c(lines, "", render_markdown_table_local(section_rows[, show_cols, drop = FALSE]))
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
section_metadata <- prepare_enrichment_metadata_06c(enrichment_targets, gene_program_registry)
metadata_join_status <- metadata_join_status_06c(section_metadata)
section_summary <- build_section_summary_06c(section_metadata, all_manifest)
write_tsv_local(section_summary, paths$section_summary_tsv)

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
  sprintf("- enrichment section summary: `%s`", paths$section_summary_tsv),
  sprintf("- shared pathways: `%s`", paths$shared_pathways_tsv),
  sprintf("- GO BP heatmap: `%s`", paths$cross_cluster_go_heatmap_png),
  sprintf("- KEGG heatmap: `%s`", paths$cross_cluster_kegg_heatmap_png),
  "",
  "## Metadata Join",
  sprintf("- metadata_join_status=%s", metadata_join_status),
  if (identical(metadata_join_status, "ok")) {
    "- metadata_join_reason: all enrichment targets joined to gene_program_registry.tsv."
  } else {
    "- metadata_join_reason: gene_program_registry.tsv is missing or incomplete for at least one enrichment target."
  },
  "",
  section_report_lines_06c(section_summary, "Core Mechanism Enrichment"),
  "",
  section_report_lines_06c(section_summary, "Subtype / Identity Enrichment"),
  "",
  section_report_lines_06c(section_summary, "Global Context Enrichment", include_global_warning = TRUE),
  "",
  section_report_lines_06c(section_summary, "QC / Non-gene-program"),
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
    enrichment_section_summary_tsv = build_output_entry(paths$section_summary_tsv, "tsv", module_name, "one row per enrichment target with section assignment", base_dir = cfg$project_root, schema = infer_schema_from_df(section_summary)),
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
