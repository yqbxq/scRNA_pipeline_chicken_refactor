#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
      getwd()
    }
  }
)

source_utf8 <- function(path) {
  source(path, encoding = "UTF-8")
}

source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "config.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_02.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_03.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_04.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "triage_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "reduction_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "integration_dispatch.R"))
source_utf8(file.path(.script_dir, "helpers", "clustering_utils.R"))

load_required_packages(c("Seurat", "dplyr", "ggplot2", "tibble", "jsonlite"))

cfg <- get_single_script_config_04()
module_name <- "04a_subcluster_build"
prepare_dirs_04(cfg)
set.seed(cfg$random_seed)

empty_candidate_index_04 <- function() {
  data.frame(
    layer_id = character(0),
    normalization = character(0),
    integration = character(0),
    candidate_id = character(0),
    candidate_key = character(0),
    mode = character(0),
    reduction_name = character(0),
    umap_name = character(0),
    runtime_sec = numeric(0),
    downgrade_reason = character(0),
    out_rds = character(0),
    selected_integration_file = character(0),
    stringsAsFactors = FALSE
  )
}

empty_subcluster_index_04 <- function() {
  data.frame(
    layer_id = character(0),
    parent_layer = character(0),
    mode = character(0),
    status = character(0),
    candidate_count = integer(0),
    selected_integration = character(0),
    selected_integration_file = character(0),
    clustered_key = character(0),
    clustered_rds = character(0),
    report_md = character(0),
    reason = character(0),
    stringsAsFactors = FALSE
  )
}

append_triage <- function(rows, layer_id, severity, signal_id, evidence) {
  rows[[length(rows) + 1]] <- make_triage_row(
    sample_id = layer_id,
    severity = severity,
    signal_id = signal_id,
    evidence = evidence
  )
  rows
}

resolve_parent_column <- function(meta_df, requested) {
  requested <- normalize_scalar_value(requested)
  if (!nzchar(requested)) {
    return("")
  }
  if (requested %in% colnames(meta_df)) {
    return(requested)
  }
  aliases <- switch(
    requested,
    parent_cell_type = c("cell_type", "panorama_cell_type"),
    panorama_cell_type = c("cell_type", "panorama_cell_type"),
    parent_cell_type_confidence = c("cell_type_confidence", "panorama_cell_type_confidence"),
    panorama_cell_type_confidence = c("cell_type_confidence", "panorama_cell_type_confidence"),
    parent_annotation_relation = c("annotation_relation", "panorama_annotation_relation"),
    panorama_annotation_relation = c("annotation_relation", "panorama_annotation_relation"),
    character(0)
  )
  hit <- aliases[aliases %in% colnames(meta_df)]
  if (length(hit) > 0) {
    hit[[1]]
  } else {
    ""
  }
}

subset_parent_for_layer <- function(parent_seu, layer_spec) {
  meta_df <- parent_seu@meta.data
  keep <- rep(TRUE, ncol(parent_seu))
  names(keep) <- colnames(parent_seu)
  reason <- character(0)

  sample_col <- if ("sample_id" %in% colnames(meta_df)) "sample_id" else if ("orig.ident" %in% colnames(meta_df)) "orig.ident" else ""
  if (length(layer_spec$sample_include) > 0) {
    if (!nzchar(sample_col)) {
      return(list(object = NULL, reason = "missing_sample_column"))
    }
    keep <- keep & as.character(meta_df[[sample_col]]) %in% layer_spec$sample_include
    reason <- c(reason, sprintf("sample_include=%s", paste(layer_spec$sample_include, collapse = ",")))
  }
  if (length(layer_spec$sample_exclude) > 0) {
    if (!nzchar(sample_col)) {
      return(list(object = NULL, reason = "missing_sample_column"))
    }
    keep <- keep & !as.character(meta_df[[sample_col]]) %in% layer_spec$sample_exclude
    reason <- c(reason, sprintf("sample_exclude=%s", paste(layer_spec$sample_exclude, collapse = ",")))
  }

  if (length(reason) == 0 && nzchar(layer_spec$selection_column) && length(layer_spec$selection_values) > 0) {
    selection_col <- resolve_parent_column(meta_df, layer_spec$selection_column)
    if (!nzchar(selection_col)) {
      return(list(object = NULL, reason = sprintf("missing_selection_column:%s", layer_spec$selection_column)))
    }
    keep <- keep & as.character(meta_df[[selection_col]]) %in% layer_spec$selection_values
    reason <- c(reason, sprintf("%s in %s", selection_col, paste(layer_spec$selection_values, collapse = ",")))
  }

  if (!any(keep)) {
    return(list(object = NULL, reason = paste(c(reason, "no_matching_cells"), collapse = "; ")))
  }
  list(object = subset(parent_seu, cells = names(keep)[keep]), reason = paste(reason, collapse = "; "))
}

prepare_candidate_object <- function(parent_seu, layer_spec, norm_method, integration_mode) {
  seu <- maybe_join_layers(parent_seu)
  seu <- strip_reduction_state_preserving_parent_meta(seu, parent_prefix = "panorama_")
  seu <- normalize_layer(seu, norm_method, layer_spec, cfg)
  assay_name <- if (identical(norm_method, "sct") && "SCT" %in% Assays(seu)) "SCT" else "RNA"
  seu <- reduce_pca_umap(
    seu,
    layer_spec,
    assay = assay_name,
    reduction_key_prefix = sprintf("PCA%s_", toupper(norm_method)),
    umap_name = sprintf("umap_rna_%s", norm_method)
  )
  dispatch_result <- dispatch_integration(
    seu,
    integration_mode,
    group_var = "orig.ident",
    layer_spec = layer_spec,
    normalization_method = norm_method
  )
  seu <- dispatch_result$object
  dims <- usable_reduction_dims(seu, layer_spec$pca_dims, dispatch_result$reduction_name)
  if (length(dims) >= 2 && !dispatch_result$umap_name %in% Reductions(seu)) {
    seu <- RunUMAP(
      seu,
      reduction = dispatch_result$reduction_name,
      dims = dims,
      reduction.name = dispatch_result$umap_name,
      reduction.key = paste0(gsub("[^A-Za-z0-9]", "", toupper(dispatch_result$umap_name)), "_"),
      verbose = FALSE
    )
  }
  seu@misc$layer_spec <- layer_spec
  seu@misc$normalization_method <- norm_method
  seu@misc$integration_mode <- integration_mode
  list(object = seu, dispatch = dispatch_result)
}

finalize_subcluster_candidate <- function(seu, layer_spec, candidate_row, selected_value) {
  layer_id <- layer_spec$layer_id
  layer_dir <- layer_checkpoint_dir_04(cfg, layer_id)
  table_dir <- layer_subcluster_table_dir_04(cfg, layer_id)
  report_dir <- layer_subcluster_report_dir_04(cfg, layer_id)
  ensure_dir(layer_dir)
  ensure_dir(table_dir)
  ensure_dir(report_dir)

  cluster_col <- paste0(layer_id, "_cluster")
  resolution_search_tsv <- file.path(table_dir, "resolution_search.tsv")
  cluster_result <- run_resolution_search(
    seu,
    reduction_name = candidate_row$reduction_name[[1]],
    layer_spec = layer_spec,
    seed = cfg$random_seed,
    log_path = resolution_search_tsv
  )
  seu <- finalize_layer_object(
    cluster_result$seu,
    reduction_name = candidate_row$reduction_name[[1]],
    umap_name = candidate_row$umap_name[[1]],
    cluster_col = cluster_col
  )
  seu$cell_type <- ""
  seu$cell_type_confidence <- "未定"
  seu$annotation_relation <- "无关"
  seu@misc$selected_integration <- selected_value
  seu@misc$selected_resolution <- cluster_result$selected_resolution

  clustered_rds <- file.path(layer_dir, sprintf("%s_after_clustering.rds", layer_id))
  saveRDS(seu, clustered_rds)

  cluster_summary <- build_cluster_summary_local(seu, cluster_col)
  cluster_summary_csv <- file.path(table_dir, "cluster_summary.csv")
  write_csv_local(cluster_summary, cluster_summary_csv)

  selected_resolution_txt <- file.path(table_dir, "selected_resolution.txt")
  writeLines(
    c(
      sprintf("selected_resolution=%s", cluster_result$selected_resolution),
      sprintf("selected_cluster_count=%s", cluster_result$selected_cluster_count),
      sprintf("target_clusters=%s", layer_spec$target_clusters),
      sprintf("fallback_used=%s", ifelse(cluster_result$fallback_used, "true", "false")),
      sprintf("selected_integration=%s", selected_value)
    ),
    selected_resolution_txt,
    useBytes = TRUE
  )

  report_path <- file.path(report_dir, "report.md")
  report_lines <- c(
    sprintf("# 04a Subcluster Build Report: %s", layer_id),
    "",
    sprintf("- parent_layer: `%s`", layer_spec$parent_layer),
    sprintf("- selected_integration: `%s`", selected_value),
    sprintf("- cluster_column: `%s`", cluster_col),
    sprintf("- selected_resolution: `%s`", cluster_result$selected_resolution),
    sprintf("- selected_cluster_count: `%s`", cluster_result$selected_cluster_count),
    sprintf("- target_clusters: `%s`", layer_spec$target_clusters),
    sprintf("- fallback_used: `%s`", ifelse(cluster_result$fallback_used, "true", "false")),
    sprintf("- clustered_rds: `%s`", clustered_rds),
    "",
    "## Cluster Summary",
    render_markdown_table_local(head(cluster_summary, 50))
  )
  write_markdown_local(report_lines, report_path)

  status_row <- data.frame(
    layer_id = layer_id,
    status = "built",
    parent_layer = layer_spec$parent_layer,
    cluster_column = cluster_col,
    selected_normalization = candidate_row$normalization[[1]],
    selected_integration = candidate_row$integration[[1]],
    selected_reduction = candidate_row$reduction_name[[1]],
    selected_umap = candidate_row$umap_name[[1]],
    selected_resolution = as.character(cluster_result$selected_resolution),
    selected_cluster_count = as.character(cluster_result$selected_cluster_count),
    fallback_used = ifelse(cluster_result$fallback_used, "true", "false"),
    clustered_rds = normalizePath(clustered_rds, winslash = "/", mustWork = FALSE),
    annotated_rds = "",
    stringsAsFactors = FALSE
  )
  invisible(upsert_layer_status(cfg$layer_status_file, status_row))

  list(
    object = seu,
    clustered_rds = clustered_rds,
    resolution_search_tsv = resolution_search_tsv,
    cluster_summary_csv = cluster_summary_csv,
    selected_resolution_txt = selected_resolution_txt,
    report_path = report_path,
    cluster_result = cluster_result
  )
}

manifest_03d <- read_manifest_local(cfg$module_03d_manifest_path)
panorama_annotated_rds <- resolve_output_local(manifest_03d, "annotated_object")
layer_df <- validate_layer_config(read_object_layer_config(cfg))
enabled_layers <- topo_order_layers(filter_enabled_layers(layer_df))
sub_rows <- enabled_layers[enabled_layers$layer_role == "subcluster", , drop = FALSE]

parent_objects <- list()
selected_panorama <- NULL
if (nrow(sub_rows) > 0) {
  panorama_seu <- readRDS(panorama_annotated_rds)
  panorama_seu <- maybe_join_layers(panorama_seu)
  parent_objects[[cfg$panorama_layer_id]] <- panorama_seu

  selected_panorama_value <- read_selected_integration_value_local(cfg$selected_integration_file)
  if (!nzchar(selected_panorama_value)) {
    stop(sprintf("缺少 panorama selected integration: %s", cfg$selected_integration_file), call. = FALSE)
  }
  selected_panorama <- parse_selected_integration_local(selected_panorama_value)
}

candidate_rows <- list()
subcluster_rows <- list()
triage_rows <- list()
output_entries <- list()

if (nrow(sub_rows) > 0) {
  for (idx in seq_len(nrow(sub_rows))) {
    layer_spec <- layer_config_row_to_spec(sub_rows[idx, , drop = FALSE])
    layer_id <- layer_spec$layer_id
    selected_file <- selected_integration_file_04(cfg, layer_id)
    clustered_key <- paste0("clustered_", layer_id)
    ensure_dir(layer_checkpoint_dir_04(cfg, layer_id))
    ensure_dir(file.path(layer_checkpoint_dir_04(cfg, layer_id), "candidates"))
    ensure_dir(layer_integration_report_dir_04(cfg, layer_id))
    ensure_dir(layer_integration_table_dir_04(cfg, layer_id))
    ensure_dir(layer_subcluster_report_dir_04(cfg, layer_id))
    ensure_dir(layer_subcluster_table_dir_04(cfg, layer_id))

    if (!layer_spec_has_filter(layer_spec)) {
      triage_rows <- append_triage(triage_rows, layer_id, "info", "subcluster_skipped_no_filter", "missing sample_include/sample_exclude and selection_column/selection_values")
      status_row <- data.frame(
        layer_id = layer_id,
        status = "skipped_no_filter",
        parent_layer = layer_spec$parent_layer,
        cluster_column = "",
        selected_normalization = "",
        selected_integration = "",
        selected_reduction = "",
        selected_umap = "",
        selected_resolution = "",
        selected_cluster_count = "",
        fallback_used = "",
        clustered_rds = "",
        annotated_rds = "",
        stringsAsFactors = FALSE
      )
      invisible(upsert_layer_status(cfg$layer_status_file, status_row))
      subcluster_rows[[length(subcluster_rows) + 1]] <- data.frame(
        layer_id = layer_id,
        parent_layer = layer_spec$parent_layer,
        mode = "skipped",
        status = "skipped_no_filter",
        candidate_count = 0L,
        selected_integration = "",
        selected_integration_file = selected_file,
        clustered_key = clustered_key,
        clustered_rds = "",
        report_md = "",
        reason = "missing_filter",
        stringsAsFactors = FALSE
      )
      next
    }

    if (!layer_spec$parent_layer %in% names(parent_objects)) {
      subcluster_rows[[length(subcluster_rows) + 1]] <- data.frame(
        layer_id = layer_id,
        parent_layer = layer_spec$parent_layer,
        mode = "skipped",
        status = "skipped_parent_unavailable",
        candidate_count = 0L,
        selected_integration = "",
        selected_integration_file = selected_file,
        clustered_key = clustered_key,
        clustered_rds = "",
        report_md = "",
        reason = "parent_unavailable",
        stringsAsFactors = FALSE
      )
      next
    }

    subset_payload <- subset_parent_for_layer(parent_objects[[layer_spec$parent_layer]], layer_spec)
    if (is.null(subset_payload$object) || ncol(subset_payload$object) < 3) {
      reason <- subset_payload$reason %||% "too_few_cells"
      subcluster_rows[[length(subcluster_rows) + 1]] <- data.frame(
        layer_id = layer_id,
        parent_layer = layer_spec$parent_layer,
        mode = "skipped",
        status = "skipped_empty_or_too_small",
        candidate_count = 0L,
        selected_integration = "",
        selected_integration_file = selected_file,
        clustered_key = clustered_key,
        clustered_rds = "",
        report_md = "",
        reason = reason,
        stringsAsFactors = FALSE
      )
      next
    }

    candidate_grid <- cap_layer_candidates(
      layer_spec,
      selected_normalization = selected_panorama$normalization,
      selected_integration = selected_panorama$integration,
      cap = cfg$max_integration_candidates_per_layer
    )
    if (nrow(candidate_grid) == 0) {
      stop(sprintf("layer `%s` 没有可运行的 normalization/integration 候选。", layer_id), call. = FALSE)
    }
    if (any(candidate_grid$capped == "true")) {
      triage_rows <- append_triage(
        triage_rows,
        layer_id,
        "low",
        "subcluster_candidate_capped",
        sprintf("original_candidate_count=%s; cap=%s", candidate_grid$original_candidate_count[[1]], cfg$max_integration_candidates_per_layer)
      )
    }

    selected_value <- read_selected_integration_value_local(selected_file)
    selected_from_file <- NULL
    if (nzchar(selected_value)) {
      selected_from_file <- parse_selected_integration_local(selected_value)
      if (!selected_from_file$value %in% candidate_grid$candidate_id) {
        stop(sprintf("layer `%s` selected_integration.txt 选择不合法: %s", layer_id, selected_from_file$value), call. = FALSE)
      }
    }

    should_finalize <- nrow(candidate_grid) == 1 || !is.null(selected_from_file)
    selected_candidate_id <- if (!is.null(selected_from_file)) selected_from_file$value else candidate_grid$candidate_id[[1]]
    candidates_to_run <- if (should_finalize) {
      candidate_grid[candidate_grid$candidate_id == selected_candidate_id, , drop = FALSE]
    } else {
      candidate_grid
    }

    selected_candidate_row <- NULL
    selected_candidate_obj <- NULL
    for (candidate_idx in seq_len(nrow(candidates_to_run))) {
      c_row <- candidates_to_run[candidate_idx, , drop = FALSE]
      candidate_key <- sprintf("candidate_%s_%s_%s", layer_id, c_row$normalization[[1]], c_row$integration[[1]])
      candidate_rds <- file.path(layer_checkpoint_dir_04(cfg, layer_id), "candidates", sprintf("%s__%s__%s.rds", layer_id, c_row$normalization[[1]], c_row$integration[[1]]))

      if (file.exists(candidate_rds) && should_finalize) {
        seu_candidate <- readRDS(candidate_rds)
        dispatch <- list(
          reduction_name = normalize_scalar_value(seu_candidate@misc$selected_reduction_candidate, "pca"),
          umap_name = normalize_scalar_value(seu_candidate@misc$selected_umap_candidate, sprintf("umap_%s_%s", c_row$integration[[1]], c_row$normalization[[1]])),
          diagnostics = list(runtime_sec = 0, downgrade_reason = "")
        )
      } else {
        message("04a candidate: ", layer_id, " / ", c_row$candidate_id[[1]])
        prepared <- prepare_candidate_object(subset_payload$object, layer_spec, c_row$normalization[[1]], c_row$integration[[1]])
        seu_candidate <- prepared$object
        dispatch <- prepared$dispatch
        seu_candidate@misc$selected_reduction_candidate <- dispatch$reduction_name
        seu_candidate@misc$selected_umap_candidate <- dispatch$umap_name
        saveRDS(seu_candidate, candidate_rds)
      }

      candidate_row <- data.frame(
        layer_id = layer_id,
        normalization = c_row$normalization[[1]],
        integration = c_row$integration[[1]],
        candidate_id = c_row$candidate_id[[1]],
        candidate_key = candidate_key,
        mode = c_row$mode[[1]],
        reduction_name = dispatch$reduction_name,
        umap_name = dispatch$umap_name,
        runtime_sec = as.numeric(dispatch$diagnostics$runtime_sec %||% 0),
        downgrade_reason = normalize_scalar_value(dispatch$diagnostics$downgrade_reason),
        out_rds = normalizePath(candidate_rds, winslash = "/", mustWork = FALSE),
        selected_integration_file = selected_file,
        stringsAsFactors = FALSE
      )
      candidate_rows[[length(candidate_rows) + 1]] <- candidate_row
      output_entries[[candidate_key]] <- build_output_entry(
        candidate_rds,
        "rds",
        module_name,
        sprintf("subcluster candidate layer=%s normalization=%s integration=%s", layer_id, c_row$normalization[[1]], c_row$integration[[1]]),
        base_dir = cfg$project_root
      )
      if (identical(c_row$candidate_id[[1]], selected_candidate_id)) {
        selected_candidate_row <- candidate_row
        selected_candidate_obj <- seu_candidate
      }
    }

    if (should_finalize) {
      if (is.null(selected_candidate_row) || is.null(selected_candidate_obj)) {
        stop(sprintf("layer `%s` 未能准备 selected candidate: %s", layer_id, selected_candidate_id), call. = FALSE)
      }
      finalized <- finalize_subcluster_candidate(selected_candidate_obj, layer_spec, selected_candidate_row, selected_candidate_id)
      output_entries[[clustered_key]] <- build_output_entry(
        finalized$clustered_rds,
        "rds",
        module_name,
        sprintf("clustered subcluster layer %s", layer_id),
        base_dir = cfg$project_root
      )
      output_entries[[paste0("resolution_search_tsv_", layer_id)]] <- build_output_entry(finalized$resolution_search_tsv, "tsv", module_name, sprintf("resolution search for %s", layer_id), base_dir = cfg$project_root)
      output_entries[[paste0("cluster_summary_csv_", layer_id)]] <- build_output_entry(finalized$cluster_summary_csv, "csv", module_name, sprintf("cluster summary for %s", layer_id), base_dir = cfg$project_root)
      output_entries[[paste0("selected_resolution_txt_", layer_id)]] <- build_output_entry(finalized$selected_resolution_txt, "txt", module_name, sprintf("selected resolution for %s", layer_id), base_dir = cfg$project_root)
      output_entries[[paste0("report_md_", layer_id)]] <- build_output_entry(finalized$report_path, "md", module_name, sprintf("04a report for %s", layer_id), base_dir = cfg$project_root)

      if (isTRUE(finalized$cluster_result$fallback_used)) {
        triage_rows <- append_triage(
          triage_rows,
          layer_id,
          "medium",
          "subcluster_resolution_fallback_used",
          sprintf("target_clusters=%s; selected_cluster_count=%s; selected_resolution=%s", layer_spec$target_clusters, finalized$cluster_result$selected_cluster_count, finalized$cluster_result$selected_resolution)
        )
      }
      parent_objects[[layer_id]] <- finalized$object
      subcluster_rows[[length(subcluster_rows) + 1]] <- data.frame(
        layer_id = layer_id,
        parent_layer = layer_spec$parent_layer,
        mode = if (nrow(candidate_grid) == 1) candidate_grid$mode[[1]] else "candidate",
        status = "built",
        candidate_count = nrow(candidate_grid),
        selected_integration = selected_candidate_id,
        selected_integration_file = selected_file,
        clustered_key = clustered_key,
        clustered_rds = normalizePath(finalized$clustered_rds, winslash = "/", mustWork = FALSE),
        report_md = normalizePath(finalized$report_path, winslash = "/", mustWork = FALSE),
        reason = subset_payload$reason,
        stringsAsFactors = FALSE
      )
    } else {
      subcluster_rows[[length(subcluster_rows) + 1]] <- data.frame(
        layer_id = layer_id,
        parent_layer = layer_spec$parent_layer,
        mode = "candidate",
        status = "pending_integration_review",
        candidate_count = nrow(candidate_grid),
        selected_integration = "",
        selected_integration_file = selected_file,
        clustered_key = clustered_key,
        clustered_rds = "",
        report_md = "",
        reason = subset_payload$reason,
        stringsAsFactors = FALSE
      )
    }
  }
}

candidate_index <- if (length(candidate_rows) > 0) dplyr::bind_rows(candidate_rows) else empty_candidate_index_04()
subcluster_index <- if (length(subcluster_rows) > 0) dplyr::bind_rows(subcluster_rows) else empty_subcluster_index_04()
triage_df <- if (length(triage_rows) > 0) dplyr::bind_rows(triage_rows) else empty_triage_df(include_sample = TRUE)

candidate_index_tsv <- file.path(cfg$subcluster_table_dir, "subcluster_candidate_index.tsv")
subcluster_index_tsv <- file.path(cfg$subcluster_table_dir, "subcluster_index.tsv")
triage_tsv <- file.path(cfg$subcluster_table_dir, "subcluster_triage.tsv")
write_tsv_local(candidate_index, candidate_index_tsv)
write_tsv_local(subcluster_index, subcluster_index_tsv)
write_tsv_local(triage_df, triage_tsv)

output_entries$candidate_index_tsv <- build_output_entry(candidate_index_tsv, "tsv", module_name, "one row per subcluster integration candidate", base_dir = cfg$project_root, schema = infer_schema_from_df(candidate_index))
output_entries$subcluster_index_tsv <- build_output_entry(subcluster_index_tsv, "tsv", module_name, "one row per enabled subcluster layer", base_dir = cfg$project_root, schema = infer_schema_from_df(subcluster_index))
output_entries$subcluster_triage_tsv <- build_output_entry(triage_tsv, "tsv", module_name, "one row per subcluster build triage signal", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df))
output_entries$layer_status_tsv <- build_output_entry(cfg$layer_status_file, "tsv", module_name, "one row per built/annotated object layer", base_dir = cfg$project_root)

if (file.exists(cfg$module_04a_manifest_path)) {
  unlink(cfg$module_04a_manifest_path)
}
write_manifest_local(
  manifest_path = cfg$module_04a_manifest_path,
  new_outputs = output_entries,
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    module_03d_manifest = cfg$module_03d_manifest_path,
    panorama_annotated_object = panorama_annotated_rds,
    panorama_selected_integration_file = cfg$selected_integration_file,
    object_layer_config = cfg$object_layer_config_file
  ),
  version = cfg$module_version,
  depends_on = list(module_03d = cfg$module_03d_manifest_path)
)

message("04a 完成。subcluster_index: ", subcluster_index_tsv)
