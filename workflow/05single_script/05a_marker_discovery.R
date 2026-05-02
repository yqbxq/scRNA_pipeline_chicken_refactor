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
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "ambient_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "comparison_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "deg_utils.R"))

load_required_packages(c("Seurat", "dplyr", "tibble", "jsonlite"))

cfg <- get_single_script_config_05()
module_name <- "05a_marker_discovery"
prepare_dirs_05(cfg)
set.seed(cfg$random_seed)

run_exploratory_findmarkers_05 <- function(seu, vars, cluster_id = "", annotation_label = "") {
  prepared <- prepare_comparison_group_05(seu, vars)
  obj <- maybe_join_layers(prepared$object)
  group_var <- prepared$group_var
  meta <- obj@meta.data
  groups <- as.character(meta[[group_var]])
  n1 <- sum(groups == vars$ident_1, na.rm = TRUE)
  n2 <- sum(groups == vars$ident_2, na.rm = TRUE)
  if (n1 < vars$min_cells_per_group || n2 < vars$min_cells_per_group) {
    return(list(
      result = NULL,
      status = "too_few_cells",
      reason = sprintf("min_cells_per_group=%s; %s=%s; %s=%s", vars$min_cells_per_group, vars$ident_1, n1, vars$ident_2, n2),
      n1 = n1,
      n2 = n2
    ))
  }

  Seurat::Idents(obj) <- group_var
  res <- tryCatch(
    Seurat::FindMarkers(
      obj,
      ident.1 = vars$ident_1,
      ident.2 = vars$ident_2,
      logfc.threshold = vars$logfc_threshold,
      test.use = "wilcox",
      verbose = FALSE
    ),
    error = function(e) {
      attr(e, "deg_status") <- "findmarkers_error"
      e
    }
  )
  if (inherits(res, "error")) {
    return(list(result = NULL, status = "findmarkers_error", reason = conditionMessage(res), n1 = n1, n2 = n2))
  }
  if (is.null(res) || nrow(res) == 0) {
    return(list(result = NULL, status = "empty_result", reason = "FindMarkers returned 0 rows", n1 = n1, n2 = n2))
  }
  out <- tibble::rownames_to_column(as.data.frame(res), "gene")
  out$cluster_id <- cluster_id
  out$annotation_label <- annotation_label
  list(result = out, status = "ok", reason = "", n1 = n1, n2 = n2)
}

layer_status_df <- deg_layer_status(cfg)
comparison_df <- read_deg_comparison_sheet(cfg)

manifest_rows <- list()

for (idx in seq_len(nrow(layer_status_df))) {
  layer_row <- layer_status_df[idx, , drop = FALSE]
  layer_id <- layer_row$layer_id[[1]]
  message("05a marker discovery layer: ", layer_id)
  obj <- load_layer_for_deg(cfg, layer_row)

  for (cmp_idx in seq_len(nrow(comparison_df))) {
    comparison_row <- comparison_df[cmp_idx, , drop = FALSE]
    if (!comparison_applies_to_layer(comparison_row, layer_id)) {
      next
    }

    vars <- resolve_comparison_vars_05(obj, comparison_row)
    out_paths <- marker_discovery_paths_05(cfg, layer_id, vars$comparison_id)
    subset_result <- subset_cells_for_comparison(obj, vars)
    result_rows <- list()
    summary_rows <- list()
    inference_status <- if (vars$force_exploratory) "exploratory_forced" else "exploratory_cell_level"

    if (identical(vars$analysis_mode, "composition")) {
      write_tsv_local(empty_marker_result_05(), out_paths$exploratory_tsv)
      summary_rows[[1]] <- data.frame(
        comparison_id = vars$comparison_id,
        layer_id = layer_id,
        cluster_id = "",
        annotation_label = "",
        analysis_mode = vars$analysis_mode,
        analysis_unit = vars$analysis_unit,
        ident_1_n = NA_integer_,
        ident_2_n = NA_integer_,
        min_cells_per_group = vars$min_cells_per_group,
        result_available = "no",
        status = "skipped_composition",
        reason = "composition rows are handled by 05c",
        stringsAsFactors = FALSE
      )
      write_tsv_local(dplyr::bind_rows(summary_rows), out_paths$exploratory_summary_tsv)
      manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
        layer_id = layer_id,
        comparison_id = vars$comparison_id,
        analysis_mode = vars$analysis_mode,
        analysis_unit = vars$analysis_unit,
        gene_program_role = vars$gene_program_role,
        group_var = vars$group_var,
        subset_column = vars$subset_column,
        subset_value = vars$subset_value,
        force_exploratory = ifelse(vars$force_exploratory, "yes", "no"),
        min_cells_per_group = vars$min_cells_per_group,
        logfc_threshold = vars$logfc_threshold,
        inference_status = "skipped",
        status = "skipped_composition",
        exploratory_results_tsv = normalizePath(out_paths$exploratory_tsv, winslash = "/", mustWork = FALSE),
        exploratory_summary_tsv = normalizePath(out_paths$exploratory_summary_tsv, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
      next
    }

    if (!identical(subset_result$status, "ok")) {
      write_tsv_local(empty_marker_result_05(), out_paths$exploratory_tsv)
      summary_rows[[1]] <- data.frame(
        comparison_id = vars$comparison_id,
        layer_id = layer_id,
        cluster_id = "",
        annotation_label = "",
        analysis_mode = vars$analysis_mode,
        analysis_unit = vars$analysis_unit,
        ident_1_n = NA_integer_,
        ident_2_n = NA_integer_,
        min_cells_per_group = vars$min_cells_per_group,
        result_available = "no",
        status = subset_result$status,
        reason = subset_result$reason,
        stringsAsFactors = FALSE
      )
      write_tsv_local(dplyr::bind_rows(summary_rows), out_paths$exploratory_summary_tsv)
      manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
        layer_id = layer_id,
        comparison_id = vars$comparison_id,
        analysis_mode = vars$analysis_mode,
        analysis_unit = vars$analysis_unit,
        gene_program_role = vars$gene_program_role,
        group_var = vars$group_var,
        subset_column = vars$subset_column,
        subset_value = vars$subset_value,
        force_exploratory = ifelse(vars$force_exploratory, "yes", "no"),
        min_cells_per_group = vars$min_cells_per_group,
        logfc_threshold = vars$logfc_threshold,
        inference_status = "skipped",
        status = subset_result$status,
        exploratory_results_tsv = normalizePath(out_paths$exploratory_tsv, winslash = "/", mustWork = FALSE),
        exploratory_summary_tsv = normalizePath(out_paths$exploratory_summary_tsv, winslash = "/", mustWork = FALSE),
        stringsAsFactors = FALSE
      )
      next
    }

    obj_sub <- subset_result$object
    run_res <- run_exploratory_findmarkers_05(
      obj_sub,
      vars,
      cluster_id = "",
      annotation_label = if (identical(vars$analysis_mode, "subtype_marker")) vars$ident_1 else vars$analysis_mode
    )
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      comparison_id = vars$comparison_id,
      layer_id = layer_id,
      cluster_id = "",
      annotation_label = if (identical(vars$analysis_mode, "subtype_marker")) vars$ident_1 else vars$analysis_mode,
      analysis_mode = vars$analysis_mode,
      analysis_unit = vars$analysis_unit,
      ident_1_n = run_res$n1,
      ident_2_n = run_res$n2,
      min_cells_per_group = vars$min_cells_per_group,
      result_available = ifelse(is.null(run_res$result), "no", "yes"),
      status = run_res$status,
      reason = run_res$reason,
      stringsAsFactors = FALSE
    )
    if (!is.null(run_res$result)) {
      res <- run_res$result
      res$comparison_id <- vars$comparison_id
      res$layer_id <- layer_id
      res$analysis_mode <- vars$analysis_mode
      res$analysis_unit <- vars$analysis_unit
      res$gene_program_role <- vars$gene_program_role
      res$inference_scope <- "exploratory_cell_level"
      res$inference_status <- inference_status
      res$subset_column <- vars$subset_column
      res$subset_value <- vars$subset_value
      result_rows[[length(result_rows) + 1]] <- res
    }

    if (identical(vars$analysis_mode, "annotation_cluster_marker")) {
      cluster_ids <- sort(unique(as.character(obj_sub$cluster_id)))
      cluster_ids <- cluster_ids[nzchar(cluster_ids)]
      result_rows <- list()
      summary_rows <- list()
      for (cluster_id in cluster_ids) {
        cluster_vars <- vars
        cluster_vars$group_var <- "cluster_id"
        cluster_vars$ident_1 <- cluster_id
        cluster_vars$ident_2 <- "__rest__"
        run_res <- run_exploratory_findmarkers_05(
          obj_sub,
          cluster_vars,
          cluster_id = cluster_id,
          annotation_label = cluster_label_for_05(obj_sub, cluster_id)
        )
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        comparison_id = vars$comparison_id,
        layer_id = layer_id,
        cluster_id = cluster_id,
        annotation_label = cluster_label_for_05(obj_sub, cluster_id),
        analysis_mode = vars$analysis_mode,
        analysis_unit = vars$analysis_unit,
        ident_1_n = run_res$n1,
        ident_2_n = run_res$n2,
        min_cells_per_group = vars$min_cells_per_group,
        result_available = ifelse(is.null(run_res$result), "no", "yes"),
        status = run_res$status,
        reason = run_res$reason,
        stringsAsFactors = FALSE
      )
      if (!is.null(run_res$result)) {
        res <- run_res$result
        res$comparison_id <- vars$comparison_id
        res$layer_id <- layer_id
        res$analysis_mode <- vars$analysis_mode
        res$analysis_unit <- vars$analysis_unit
        res$gene_program_role <- vars$gene_program_role
        res$inference_scope <- "exploratory_cell_level"
        res$inference_status <- inference_status
        res$subset_column <- vars$subset_column
        res$subset_value <- vars$subset_value
        result_rows[[length(result_rows) + 1]] <- res
      }
      }
    }

    result_df <- if (length(result_rows) > 0) dplyr::bind_rows(result_rows) else empty_marker_result_05()
    summary_df <- if (length(summary_rows) > 0) dplyr::bind_rows(summary_rows) else empty_marker_summary_05()
    write_tsv_local(result_df, out_paths$exploratory_tsv)
    write_tsv_local(summary_df, out_paths$exploratory_summary_tsv)

    manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
      layer_id = layer_id,
      comparison_id = vars$comparison_id,
      analysis_mode = vars$analysis_mode,
      analysis_unit = vars$analysis_unit,
      gene_program_role = vars$gene_program_role,
      group_var = vars$group_var,
      subset_column = vars$subset_column,
      subset_value = vars$subset_value,
      force_exploratory = ifelse(vars$force_exploratory, "yes", "no"),
      min_cells_per_group = vars$min_cells_per_group,
      logfc_threshold = vars$logfc_threshold,
      inference_status = inference_status,
      status = ifelse(nrow(result_df) > 0, "ok", "empty_result"),
      exploratory_results_tsv = normalizePath(out_paths$exploratory_tsv, winslash = "/", mustWork = FALSE),
      exploratory_summary_tsv = normalizePath(out_paths$exploratory_summary_tsv, winslash = "/", mustWork = FALSE),
      stringsAsFactors = FALSE
    )
  }
}

manifest_df <- if (length(manifest_rows) > 0) dplyr::bind_rows(manifest_rows) else empty_df_05(c(
  "layer_id", "comparison_id", "analysis_mode", "analysis_unit", "gene_program_role",
  "group_var", "subset_column", "subset_value",
  "force_exploratory", "min_cells_per_group", "logfc_threshold",
  "inference_status", "status", "exploratory_results_tsv", "exploratory_summary_tsv"
))
manifest_tsv <- file.path(cfg$marker_discovery_table_dir, "marker_discovery_manifest.tsv")
write_tsv_local(manifest_df, manifest_tsv)

report_path <- file.path(cfg$marker_discovery_report_dir, "report.md")
write_markdown_local(
  c(
    "# 05a Marker Discovery",
    "",
    "- inference_scope: `exploratory_cell_level`",
    sprintf("- manifest: `%s`", manifest_tsv),
    "",
    "## Manifest",
    render_markdown_table_local(head(manifest_df, 100))
  ),
  report_path
)

if (file.exists(cfg$module_05a_manifest_path)) {
  unlink(cfg$module_05a_manifest_path)
}
write_manifest_local(
  manifest_path = cfg$module_05a_manifest_path,
  new_outputs = list(
    marker_discovery_manifest_tsv = build_output_entry(manifest_tsv, "tsv", module_name, "one row per layer/comparison marker discovery result", base_dir = cfg$project_root, schema = infer_schema_from_df(manifest_df)),
    report = build_output_entry(report_path, "md", module_name, "marker discovery report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(layer_status_tsv = cfg$layer_status_file, comparison_sheet = cfg$comparison_sheet),
  version = cfg$module_version
)

message("05a completed. manifest: ", manifest_tsv)
