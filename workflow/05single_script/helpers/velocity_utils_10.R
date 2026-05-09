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

velocity_read_pairs_10 <- function(cfg) {
  pairs <- trajectory_read_pairs_09(cfg)
  if (nrow(pairs) == 0) {
    return(pairs)
  }
  pairs[pairs$enabled != "no" & pairs$method == "velocity", , drop = FALSE]
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
