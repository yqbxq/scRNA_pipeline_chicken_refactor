velocity_empty_df_10 <- function(cols) {
  trajectory_empty_df_09(cols)
}

velocity_loom_index_cols_10 <- c(
  "sample_id", "bam_path", "matrix_h5", "whitelist_path", "loom_path",
  "status", "reason", "runtime_s", "command"
)

velocity_loom_jobs_cols_10 <- c(
  "sample_id", "bam_path", "matrix_h5", "whitelist_path", "loom_path", "command"
)

velocity_reference_index_cols_10 <- c(
  "pair_id", "source_question_id", "layer_scope", "split_mode",
  "split_var", "split_value", "input_rds", "umap_csv", "metadata_csv",
  "n_cells", "coarse_label_var", "fine_label_var", "status", "reason",
  "produced_at"
)

velocity_method_index_cols_10 <- c(
  "pair_id", "split_value", "method", "methods_enabled", "input_path",
  "output_path", "extra_path", "figure_path", "n_cells", "status",
  "reason", "runtime_s"
)

velocity_scvelo_index_cols_10 <- c(
  velocity_method_index_cols_10,
  "qc_tsv", "vector_tsv", "stochastic_status"
)

velocity_scvelo_qc_cols_10 <- c(
  "pair_id", "split_value", "n_cells", "n_genes", "spliced_total",
  "unspliced_total", "unspliced_spliced_ratio", "fit_likelihood_mean",
  "fit_likelihood_median", "fit_likelihood_max", "velocity_gene_n",
  "velocity_confidence_mean", "velocity_confidence_median",
  "velocity_confidence_q05", "velocity_confidence_q95",
  "velocity_length_mean", "velocity_confidence_stochastic_mean",
  "velocity_length_stochastic_mean", "status", "reason"
)

velocity_driver_index_cols_10 <- c(
  velocity_method_index_cols_10,
  "driver_tsv", "driver_n", "velocity_gene_n"
)

velocity_driver_overlap_cols_10 <- c(
  "pair_id", "split_a", "split_b", "driver_n_a", "driver_n_b",
  "overlap_n", "jaccard", "overlap_genes"
)

velocity_velocyto_steady_index_cols_10 <- c(
  velocity_method_index_cols_10,
  "direction_tsv"
)

velocity_cellrank_index_cols_10 <- c(
  velocity_method_index_cols_10,
  "fate_csv", "macrostates_tsv", "terminal_figure_path"
)

velocity_consistency_index_cols_10 <- c(
  velocity_method_index_cols_10,
  "split_compare_path"
)

velocity_consistency_summary_cols_10 <- c(
  "pair_id", "split_value", "metric_id", "group_a", "group_b",
  "n", "value", "p_value", "status", "reason"
)

velocity_split_compare_cols_10 <- c(
  "pair_id", "comparison_id", "metric_id", "group_a", "group_b",
  "n_a", "n_b", "value", "p_value", "status", "reason"
)

velocity_root_terminal_index_cols_10 <- c(
  velocity_method_index_cols_10,
  "root_terminal_tsv"
)

velocity_root_terminal_cols_10 <- c(
  "pair_id", "split_value", "cluster", "cell_n", "initial_score",
  "root_score", "terminal_score", "initial_probability",
  "terminal_probability", "latent_time_mean", "status", "reason"
)

velocity_module_status_cols_10 <- c(
  "pair_id", "value", "method", "methods_enabled", "status",
  "n_cells_used", "runtime_s", "metric_value", "notes"
)

velocity_read_pairs_10 <- function(cfg) {
  pairs <- trajectory_read_pairs_09(cfg)
  if (nrow(pairs) == 0) {
    return(pairs)
  }
  pairs[pairs$enabled != "no" & pairs$method == "velocity", , drop = FALSE]
}

velocity_manifest_output_optional_10 <- function(manifest_path, key) {
  trajectory_manifest_output_optional_09(manifest_path, key)
}

velocity_read_loom_index_10 <- function(cfg) {
  path <- velocity_manifest_output_optional_10(cfg$module_10a_manifest_path, "velocity_loom_index")
  if (!nzchar(path)) {
    path <- cfg$velocity_loom_index_tsv
  }
  idx <- read_tsv_optional(path)
  if (nrow(idx) == 0) {
    return(velocity_empty_df_10(velocity_loom_index_cols_10))
  }
  for (col in velocity_loom_index_cols_10) {
    if (!col %in% colnames(idx)) {
      idx[[col]] <- ""
    }
  }
  idx
}

velocity_read_reference_index_10 <- function(cfg, include_not_ok = FALSE) {
  path <- velocity_manifest_output_optional_10(cfg$module_10b_manifest_path, "velocity_reference_index")
  if (!nzchar(path)) {
    path <- cfg$velocity_reference_index_tsv
  }
  idx <- read_tsv_optional(path)
  if (nrow(idx) == 0) {
    return(velocity_empty_df_10(velocity_reference_index_cols_10))
  }
  for (col in velocity_reference_index_cols_10) {
    if (!col %in% colnames(idx)) {
      idx[[col]] <- ""
    }
  }
  if (!isTRUE(include_not_ok)) {
    idx <- idx[idx$status == "ok", , drop = FALSE]
  }
  idx
}

velocity_read_scvelo_index_10 <- function(cfg, include_not_ok = FALSE) {
  path <- velocity_manifest_output_optional_10(cfg$module_10c_manifest_path, "scvelo_index_tsv")
  if (!nzchar(path)) {
    path <- cfg$velocity_scvelo_index_tsv
  }
  idx <- read_tsv_optional(path)
  if (nrow(idx) == 0) {
    return(velocity_empty_df_10(velocity_scvelo_index_cols_10))
  }
  for (col in velocity_scvelo_index_cols_10) {
    if (!col %in% colnames(idx)) {
      idx[[col]] <- ""
    }
  }
  if (!isTRUE(include_not_ok)) {
    idx <- idx[idx$status == "ok", , drop = FALSE]
  }
  idx
}

velocity_read_cellrank_index_10 <- function(cfg, include_not_ok = FALSE) {
  path <- velocity_manifest_output_optional_10(cfg$module_10f_manifest_path, "cellrank_index_tsv")
  if (!nzchar(path)) {
    path <- cfg$velocity_cellrank_index_tsv
  }
  idx <- read_tsv_optional(path)
  if (nrow(idx) == 0) {
    return(velocity_empty_df_10(velocity_cellrank_index_cols_10))
  }
  for (col in velocity_cellrank_index_cols_10) {
    if (!col %in% colnames(idx)) {
      idx[[col]] <- ""
    }
  }
  if (!isTRUE(include_not_ok)) {
    idx <- idx[idx$status == "ok", , drop = FALSE]
  }
  idx
}

velocity_execution_units_10 <- function(cfg, include_not_ok = FALSE) {
  pairs <- velocity_read_pairs_10(cfg)
  ref <- velocity_read_reference_index_10(cfg, include_not_ok = include_not_ok)
  if (nrow(pairs) == 0 || nrow(ref) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  colnames(pairs)[colnames(pairs) == "trajectory_id"] <- "pair_id"
  merged <- merge(
    ref,
    pairs,
    by = "pair_id",
    all.x = TRUE,
    suffixes = c("", ".pair"),
    sort = FALSE
  )
  for (base_col in c("source_question_id", "layer_scope", "coarse_label_var", "fine_label_var", "split_mode", "tools_to_run", "methods_extra")) {
    pair_col <- paste0(base_col, ".pair")
    if (!base_col %in% colnames(merged)) {
      merged[[base_col]] <- ""
    }
    if (pair_col %in% colnames(merged)) {
      empty <- !nzchar(normalize_flag(merged[[base_col]], ""))
      merged[[base_col]][empty] <- merged[[pair_col]][empty]
    }
  }
  pair_cols <- grep("\\.pair$", colnames(merged), value = TRUE)
  merged[, !colnames(merged) %in% pair_cols, drop = FALSE]
}

velocity_method_enabled_10 <- function(pair_row, method, default = TRUE) {
  trajectory_method_enabled_09(pair_row, method, default = default)
}

velocity_method_disabled_reason_10 <- function(pair_row, method) {
  trajectory_method_disabled_reason_09(pair_row, method)
}

velocity_pair_id_10 <- function(pair_row) {
  normalize_scalar_value(pair_row$trajectory_id[[1]])
}

velocity_display_split_10 <- function(split_value) {
  display_scalar_value(split_value, "pooled")
}

velocity_split_values_for_pair_10 <- function(seu, pair_row) {
  split_mode <- tolower(normalize_scalar_value(pair_row$split_mode[[1]], "auto"))
  split_var <- normalize_scalar_value(pair_row$condition_split_var[[1]])
  if (split_mode %in% c("pooled", "force_pooled") || !nzchar(split_var)) {
    return(list(split_var = "", values = ""))
  }
  if (!split_var %in% colnames(seu@meta.data)) {
    stop(sprintf("condition_split_var=%s is absent from object metadata", split_var), call. = FALSE)
  }
  values <- split_csv_local(pair_row$condition_split_values[[1]])
  if (length(values) == 0) {
    values <- sort(unique(trimws(as.character(seu@meta.data[[split_var]]))))
    values <- values[nzchar(values)]
  }
  if (length(values) == 0) {
    return(list(split_var = "", values = ""))
  }
  list(split_var = split_var, values = values)
}

velocity_resolve_layer_10 <- function(cfg, layer_scope) {
  layer_scope <- normalize_scalar_value(layer_scope, cfg$panorama_layer_id)
  layers <- communication_layer_status_07(cfg)
  if (nrow(layers) == 0) {
    stop("no usable 03d/04b annotation objects found; complete annotation/subcluster outputs first", call. = FALSE)
  }

  aliases <- unique(c(layer_scope, safe_id_09(layer_scope)))
  if (tolower(layer_scope) %in% c("panorama", tolower(cfg$panorama_layer_id))) {
    aliases <- unique(c(aliases, cfg$panorama_layer_id, "panorama"))
  }

  hit <- layers[layers$layer_id %in% aliases, , drop = FALSE]
  if (nrow(hit) == 0) {
    hit <- layers[tolower(layers$layer_id) %in% tolower(aliases), , drop = FALSE]
  }
  if (nrow(hit) == 0) {
    stop(
      sprintf(
        "layer_scope=%s did not match any usable annotation object; available layers: %s",
        layer_scope,
        paste(layers$layer_id, collapse = ",")
      ),
      call. = FALSE
    )
  }
  hit[1, , drop = FALSE]
}

velocity_manifest_10 <- function(cfg, manifest_path, module_name, outputs, inputs, depends_on = list()) {
  if (file.exists(manifest_path)) {
    unlink(manifest_path)
  }
  write_manifest_local(
    manifest_path = manifest_path,
    new_outputs = outputs,
    module_name = module_name,
    base_dir = cfg$project_root,
    inputs = inputs,
    version = module_version_10(cfg, module_name),
    depends_on = depends_on
  )
}

velocity_reference_output_key_10 <- function(kind, pair_id, split_value = "") {
  sprintf("%s__%s", kind, velocity_unit_file_id_10(pair_id, split_value))
}
