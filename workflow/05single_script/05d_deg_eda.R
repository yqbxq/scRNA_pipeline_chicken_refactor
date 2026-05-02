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
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "deg_utils.R"))

load_required_packages(c("dplyr", "tibble", "jsonlite"))

cfg <- get_single_script_config_05()
module_name <- "05d_deg_eda"
prepare_dirs_05(cfg)

marker_manifest_tsv <- file.path(cfg$marker_discovery_table_dir, "marker_discovery_manifest.tsv")
pseudobulk_manifest_tsv <- file.path(cfg$pseudobulk_table_dir, "pseudobulk_manifest.tsv")
composition_manifest_tsv <- file.path(cfg$composition_table_dir, "composition_manifest.tsv")

marker_df <- read_tsv_optional(marker_manifest_tsv)
pb_df <- read_tsv_optional(pseudobulk_manifest_tsv)
comp_df <- read_tsv_optional(composition_manifest_tsv)
gene_program_targets <- read_tsv_optional(cfg$gene_program_targets_sheet)

resolve_manifest_output_05d <- function(manifest_path, keys) {
  if (!file.exists(manifest_path)) {
    return("")
  }
  manifest <- read_manifest_local(manifest_path)
  for (key in keys) {
    entry <- manifest$outputs[[key]]
    if (!is.null(entry) && !is.null(entry$path) && nzchar(as.character(entry$path))) {
      return(resolve_output_local(manifest, key))
    }
  }
  ""
}

annotation_cluster_marker_path_05d <- function(cfg, layer_id) {
  layer_id <- normalize_scalar_value(layer_id)
  if (layer_id %in% c("panorama", cfg$panorama_layer_id)) {
    hit <- resolve_manifest_output_05d(cfg$module_03d_manifest_path, c("annotation_cluster_markers_tsv"))
    if (nzchar(hit)) return(hit)
    return(file.path(cfg$annotation_table_dir_layer, "annotation_cluster_markers.tsv"))
  }
  hit <- resolve_manifest_output_05d(cfg$module_04b_manifest_path, c(paste0("annotation_cluster_markers_tsv_", layer_id)))
  if (nzchar(hit)) return(hit)
  file.path(layer_annotation_table_dir_04(cfg, layer_id), sprintf("annotation_cluster_markers_%s.tsv", layer_id))
}

annotation_evidence_path_05d <- function(cfg, layer_id) {
  layer_id <- normalize_scalar_value(layer_id)
  if (layer_id %in% c("panorama", cfg$panorama_layer_id)) {
    hit <- resolve_manifest_output_05d(cfg$module_03d_manifest_path, c("annotation_evidence_tsv"))
    if (nzchar(hit)) return(hit)
    return(file.path(cfg$annotation_table_dir_layer, "annotation_evidence.tsv"))
  }
  hit <- resolve_manifest_output_05d(cfg$module_04b_manifest_path, c(paste0("annotation_evidence_tsv_", layer_id)))
  if (nzchar(hit)) return(hit)
  file.path(layer_annotation_table_dir_04(cfg, layer_id), sprintf("annotation_evidence_%s.tsv", layer_id))
}

key_cols <- c("layer_id", "comparison_id")
empty_keyed <- empty_df_05(key_cols)
key_subset_05d <- function(df) {
  if (!all(key_cols %in% colnames(df))) {
    return(empty_keyed)
  }
  out <- df[, key_cols, drop = FALSE]
  for (col in key_cols) {
    out[[col]] <- as.character(out[[col]])
  }
  out
}

composition_result_path_05d <- function(comp_hit) {
  if (nrow(comp_hit) == 0) {
    return("")
  }
  is_qc <- "is_qc_composition" %in% colnames(comp_hit) &&
    tolower(normalize_scalar_value(comp_hit$is_qc_composition[[1]], "no")) %in% c("yes", "true", "1", "on")
  candidates <- if (is_qc) {
    c("qc_mirror_tsv", "composition_tsv", "proportion_tsv", "formal_results_tsv")
  } else {
    c("formal_results_tsv", "composition_tsv", "proportion_tsv", "qc_mirror_tsv")
  }
  for (col in candidates) {
    if (col %in% colnames(comp_hit)) {
      path <- normalize_scalar_value(comp_hit[[col]][[1]])
      if (nzchar(path)) {
        return(path)
      }
    }
  }
  ""
}

keys <- dplyr::bind_rows(
  key_subset_05d(marker_df),
  key_subset_05d(pb_df),
  key_subset_05d(comp_df)
)
keys <- unique(keys[nzchar(keys$layer_id) & nzchar(keys$comparison_id), , drop = FALSE])

status_rows <- list()
summary_rows <- list()
for (idx in seq_len(nrow(keys))) {
  layer_id <- keys$layer_id[[idx]]
  comparison_id <- keys$comparison_id[[idx]]
  marker_hit <- marker_df[marker_df$layer_id == layer_id & marker_df$comparison_id == comparison_id, , drop = FALSE]
  pb_hit <- pb_df[pb_df$layer_id == layer_id & pb_df$comparison_id == comparison_id, , drop = FALSE]
  comp_hit <- comp_df[comp_df$layer_id == layer_id & comp_df$comparison_id == comparison_id, , drop = FALSE]

  marker_status <- if (nrow(marker_hit) > 0) normalize_scalar_value(marker_hit$status[[1]], "available") else "missing"
  pb_status <- if (nrow(pb_hit) > 0) normalize_scalar_value(pb_hit$inference_status[[1]], "skipped") else "missing"
  comp_status <- if (nrow(comp_hit) > 0) normalize_scalar_value(comp_hit$inference_status[[1]], "skipped") else "missing"
  pb_path <- if (nrow(pb_hit) > 0 && "ds_results_tsv" %in% colnames(pb_hit)) normalize_scalar_value(pb_hit$ds_results_tsv[[1]]) else ""
  comp_path <- composition_result_path_05d(comp_hit)
  marker_path <- if (nrow(marker_hit) > 0 && "exploratory_results_tsv" %in% colnames(marker_hit)) normalize_scalar_value(marker_hit$exploratory_results_tsv[[1]]) else ""

  marker_counts <- count_significant_rows_05(marker_path, alpha = cfg$deg_alpha)
  pb_counts <- count_significant_rows_05(pb_path, alpha = cfg$deg_alpha)
  comp_counts <- count_significant_rows_05(comp_path, alpha = cfg$deg_alpha)

  status_rows[[length(status_rows) + 1]] <- data.frame(
    layer_id = layer_id,
    comparison_id = comparison_id,
    marker_status = marker_status,
    pseudobulk_status = pb_status,
    composition_status = comp_status,
    marker_results_tsv = marker_path,
    pseudobulk_results_tsv = pb_path,
    composition_results_tsv = comp_path,
    stringsAsFactors = FALSE
  )
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    layer_id = layer_id,
    comparison_id = comparison_id,
    exploratory_marker_total_n = marker_counts$total_n,
    exploratory_marker_significant_n = marker_counts$significant_n,
    formal_de_total_n = pb_counts$total_n,
    formal_de_significant_n = pb_counts$significant_n,
    significant_composition_cluster_n = comp_counts$significant_n,
    alpha = cfg$deg_alpha,
    stringsAsFactors = FALSE
  )
}

status_df <- if (length(status_rows) > 0) dplyr::bind_rows(status_rows) else empty_df_05(c(
  "layer_id", "comparison_id", "marker_status", "pseudobulk_status",
  "composition_status", "marker_results_tsv", "pseudobulk_results_tsv", "composition_results_tsv"
))
summary_df <- if (length(summary_rows) > 0) dplyr::bind_rows(summary_rows) else empty_df_05(c(
  "layer_id", "comparison_id", "exploratory_marker_total_n", "exploratory_marker_significant_n",
  "formal_de_total_n", "formal_de_significant_n", "significant_composition_cluster_n", "alpha"
))

paths <- deg_report_paths_05(cfg)
write_tsv_local(status_df, paths$status_matrix_tsv)
write_tsv_local(summary_df, paths$summary_tsv)

registry_cols <- c(
  "comparison_id", "source_question_id", "layer_id", "analysis_mode",
  "gene_program_role", "produces_gene_program", "annotation_only", "qc_only",
  "global_context_only", "nichenet_eligible", "nichenet_usage",
  "enrichment_eligible", "enrichment_usage", "preferred_for_downstream",
  "result_level", "formal_status", "result_status", "skip_reason",
  "eligible_reason", "ineligible_reason", "deg_tsv", "marker_tsv",
  "annotation_evidence_tsv", "composition_tsv", "top_gene_tsv",
  "n_significant", "warning"
)

if (nrow(gene_program_targets) == 0) {
  registry_df <- empty_df_05(registry_cols)
} else {
  for (col in c(
    "comparison_id", "source_question_id", "layer_scope", "analysis_mode",
    "gene_program_role", "produces_gene_program", "annotation_only", "qc_only",
    "global_context_only", "nichenet_eligible", "nichenet_usage",
    "enrichment_eligible", "enrichment_usage", "preferred_for_downstream",
    "formal_preferred", "formal_status", "result_status", "skip_reason",
    "eligible_reason", "ineligible_reason"
  )) {
    if (!col %in% colnames(gene_program_targets)) {
      gene_program_targets[[col]] <- ""
    }
  }
  registry_rows <- list()
  for (idx in seq_len(nrow(gene_program_targets))) {
    target <- gene_program_targets[idx, , drop = FALSE]
    comparison_id <- normalize_scalar_value(target$comparison_id[[1]])
    layer_id <- normalize_scalar_value(target$layer_scope[[1]])
    marker_hit <- marker_df[marker_df$layer_id == layer_id & marker_df$comparison_id == comparison_id, , drop = FALSE]
    pb_hit <- pb_df[pb_df$layer_id == layer_id & pb_df$comparison_id == comparison_id, , drop = FALSE]
    comp_hit <- comp_df[comp_df$layer_id == layer_id & comp_df$comparison_id == comparison_id, , drop = FALSE]

    marker_path <- if (nrow(marker_hit) > 0 && "exploratory_results_tsv" %in% colnames(marker_hit)) normalize_scalar_value(marker_hit$exploratory_results_tsv[[1]]) else ""
    deg_path <- if (nrow(pb_hit) > 0 && "ds_results_tsv" %in% colnames(pb_hit)) normalize_scalar_value(pb_hit$ds_results_tsv[[1]]) else ""
    composition_path <- composition_result_path_05d(comp_hit)
    formal_status <- if (nrow(pb_hit) > 0 && "inference_status" %in% colnames(pb_hit)) normalize_scalar_value(pb_hit$inference_status[[1]], "missing") else normalize_scalar_value(target$formal_status[[1]], "missing")
    marker_status <- if (nrow(marker_hit) > 0 && "status" %in% colnames(marker_hit)) normalize_scalar_value(marker_hit$status[[1]], "missing") else "missing"
    composition_status <- if (nrow(comp_hit) > 0 && "inference_status" %in% colnames(comp_hit)) normalize_scalar_value(comp_hit$inference_status[[1]], "missing") else "missing"
    analysis_mode <- normalize_scalar_value(target$analysis_mode[[1]])
    gene_program_role <- normalize_scalar_value(target$gene_program_role[[1]])
    annotation_evidence_tsv <- ""

    if (identical(analysis_mode, "annotation_cluster_marker") || identical(gene_program_role, "annotation_marker")) {
      marker_path <- annotation_cluster_marker_path_05d(cfg, layer_id)
      annotation_evidence_tsv <- annotation_evidence_path_05d(cfg, layer_id)
      deg_path <- ""
      composition_path <- ""
      formal_status <- "annotation_only"
    }

    formal_preferred <- normalize_scalar_value(target$formal_preferred[[1]], "no")
    produces_gene_program <- normalize_scalar_value(target$produces_gene_program[[1]], "yes")
    preferred_path <- if (identical(produces_gene_program, "yes") && formal_preferred == "yes" && nzchar(deg_path)) {
      deg_path
    } else if (identical(produces_gene_program, "yes") && nzchar(marker_path) && !identical(gene_program_role, "annotation_marker")) {
      marker_path
    } else if (identical(analysis_mode, "global_context") && nzchar(deg_path)) {
      deg_path
    } else {
      ""
    }
    result_level <- if (formal_preferred == "yes" && identical(formal_status, "formal")) {
      "pseudobulk_formal"
    } else if (formal_preferred == "yes" && formal_status %in% c("exploratory_only", "exploratory_forced") && nzchar(deg_path)) {
      "cell_level_exploratory"
    } else if (nzchar(marker_path)) {
      "cell_level_exploratory"
    } else if (identical(gene_program_role, "annotation_marker")) {
      "annotation_cluster_markers"
    } else if (identical(gene_program_role, "qc_only")) {
      "qc_composition"
    } else if (identical(gene_program_role, "none")) {
      "composition"
    } else {
      "unavailable"
    }
    counts <- count_significant_rows_05(preferred_path, alpha = cfg$deg_alpha)
    result_status <- normalize_scalar_value(target$result_status[[1]], "target_defined")
    if (identical(produces_gene_program, "no")) {
      result_status <- normalize_scalar_value(target$skip_reason[[1]], "not_gene_program")
    } else if (identical(gene_program_role, "annotation_marker")) {
      result_status <- if (nzchar(marker_path) && file.exists(marker_path)) "annotation_only" else "missing_annotation_output"
    } else if (nzchar(preferred_path) && file.exists(preferred_path)) {
      result_status <- "available"
    } else {
      result_status <- "unavailable"
    }
    warning <- ""
    if (formal_preferred == "yes" && !formal_status %in% c("formal", "exploratory_only", "exploratory_forced") && identical(produces_gene_program, "yes")) {
      warning <- sprintf("formal preferred but pseudobulk status is %s; using exploratory marker fallback when available", formal_status)
    } else if (!nzchar(preferred_path) && identical(produces_gene_program, "yes") && !identical(gene_program_role, "annotation_marker")) {
      warning <- sprintf("gene program unavailable; marker status is %s, formal status is %s, composition status is %s", marker_status, formal_status, composition_status)
    }

    registry_rows[[length(registry_rows) + 1]] <- data.frame(
      comparison_id = comparison_id,
      source_question_id = normalize_scalar_value(target$source_question_id[[1]]),
      layer_id = layer_id,
      analysis_mode = analysis_mode,
      gene_program_role = gene_program_role,
      produces_gene_program = produces_gene_program,
      annotation_only = normalize_scalar_value(target$annotation_only[[1]], "no"),
      qc_only = normalize_scalar_value(target$qc_only[[1]], "no"),
      global_context_only = normalize_scalar_value(target$global_context_only[[1]], "no"),
      nichenet_eligible = normalize_scalar_value(target$nichenet_eligible[[1]], "no"),
      nichenet_usage = normalize_scalar_value(target$nichenet_usage[[1]], "none"),
      enrichment_eligible = normalize_scalar_value(target$enrichment_eligible[[1]], "no"),
      enrichment_usage = normalize_scalar_value(target$enrichment_usage[[1]], "none"),
      preferred_for_downstream = normalize_scalar_value(target$preferred_for_downstream[[1]]),
      result_level = result_level,
      formal_status = formal_status,
      result_status = result_status,
      skip_reason = normalize_scalar_value(target$skip_reason[[1]]),
      eligible_reason = normalize_scalar_value(target$eligible_reason[[1]]),
      ineligible_reason = normalize_scalar_value(target$ineligible_reason[[1]]),
      deg_tsv = deg_path,
      marker_tsv = marker_path,
      annotation_evidence_tsv = annotation_evidence_tsv,
      composition_tsv = composition_path,
      top_gene_tsv = preferred_path,
      n_significant = counts$significant_n,
      warning = warning,
      stringsAsFactors = FALSE
    )
  }
  registry_df <- dplyr::bind_rows(registry_rows)
}
write_tsv_local(registry_df, paths$gene_program_registry_tsv)

forced_or_exploratory <- status_df[
  status_df$pseudobulk_status %in% c("exploratory_only", "exploratory_forced") |
    status_df$composition_status %in% c("exploratory_only", "exploratory_forced"),
  ,
  drop = FALSE
]

report_lines <- c(
  "# 05d DEG EDA",
  "",
  sprintf("- marker_manifest: `%s`", marker_manifest_tsv),
  sprintf("- pseudobulk_manifest: `%s`", pseudobulk_manifest_tsv),
  sprintf("- composition_manifest: `%s`", composition_manifest_tsv),
  sprintf("- status_matrix: `%s`", paths$status_matrix_tsv),
  sprintf("- summary: `%s`", paths$summary_tsv),
  sprintf("- gene_program_registry: `%s`", paths$gene_program_registry_tsv),
  "",
  "## Analysis Modes",
  render_markdown_table_local(
    if (nrow(status_df) == 0 || !"comparison_id" %in% colnames(status_df)) {
      empty_df_05(c("analysis_mode", "comparison_n"))
    } else {
      mode_df <- status_df
      if (!"marker_status" %in% colnames(mode_df)) mode_df$marker_status <- ""
      target_modes <- if (all(c("comparison_id", "analysis_mode") %in% colnames(gene_program_targets))) {
        gene_program_targets[, c("comparison_id", "analysis_mode"), drop = FALSE]
      } else {
        empty_df_05(c("comparison_id", "analysis_mode"))
      }
      unique(dplyr::left_join(mode_df[, "comparison_id", drop = FALSE], target_modes, by = "comparison_id")) %>%
        dplyr::count(analysis_mode, name = "comparison_n")
    }
  ),
  "",
  "## Gene Program Registry",
  render_markdown_table_local(registry_df),
  "",
  "## Inference Status Matrix",
  render_markdown_table_local(status_df),
  "",
  "## DEG Summary",
  render_markdown_table_local(summary_df),
  "",
  "## Exploratory Warnings",
  if (nrow(forced_or_exploratory) == 0) "No exploratory-only or forced exploratory results recorded." else render_markdown_table_local(forced_or_exploratory[, c("layer_id", "comparison_id", "pseudobulk_status", "composition_status", "marker_results_tsv"), drop = FALSE])
)
write_markdown_local(report_lines, paths$report_md)

if (file.exists(cfg$module_05d_manifest_path)) {
  unlink(cfg$module_05d_manifest_path)
}
write_manifest_local(
  manifest_path = cfg$module_05d_manifest_path,
  new_outputs = list(
    report = build_output_entry(paths$report_md, "md", module_name, "DEG EDA report", base_dir = cfg$project_root),
    deg_status_matrix_tsv = build_output_entry(paths$status_matrix_tsv, "tsv", module_name, "one row per layer/comparison DEG status", base_dir = cfg$project_root, schema = infer_schema_from_df(status_df)),
    deg_summary_tsv = build_output_entry(paths$summary_tsv, "tsv", module_name, "one row per layer/comparison DEG summary", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df)),
    gene_program_registry_tsv = build_output_entry(paths$gene_program_registry_tsv, "tsv", module_name, "gene program availability for downstream 06/07/08", base_dir = cfg$project_root, schema = infer_schema_from_df(registry_df))
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    marker_manifest_tsv = marker_manifest_tsv,
    pseudobulk_manifest_tsv = pseudobulk_manifest_tsv,
    composition_manifest_tsv = composition_manifest_tsv,
    gene_program_targets_tsv = cfg$gene_program_targets_sheet
  ),
  version = cfg$module_version,
  depends_on = list(
    module_05a = cfg$module_05a_manifest_path,
    module_05b = cfg$module_05b_manifest_path,
    module_05c = cfg$module_05c_manifest_path
  )
)

message("05d completed. report: ", paths$report_md)
