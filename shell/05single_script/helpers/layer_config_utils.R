LAYER_NORMALIZATION_METHODS <- c("lognorm", "sct")
LAYER_INTEGRATION_METHODS <- c("none", "harmony", "seurat_cca", "seurat_rpca", "scanorama")

sanitize_layer_id_local <- function(x) {
  value <- normalize_scalar_value(x)
  if (!nzchar(value)) {
    return("")
  }
  gsub("[^A-Za-z0-9._-]+", "_", value)
}

split_csv_local <- function(value) {
  if (length(value) == 0) {
    return(character(0))
  }
  value <- normalize_scalar_value(value[[1]])
  if (!nzchar(value)) {
    return(character(0))
  }
  parts <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  parts[nzchar(parts)]
}

parse_index_spec_local <- function(value, default = integer(0)) {
  value <- normalize_scalar_value(value)
  if (!nzchar(value)) {
    return(as.integer(default))
  }
  if (grepl("^[0-9]+:[0-9]+$", value)) {
    parts <- as.integer(strsplit(value, ":", fixed = TRUE)[[1]])
    return(seq.int(parts[1], parts[2]))
  }
  out <- suppressWarnings(as.integer(split_csv_local(value)))
  out <- out[is.finite(out)]
  if (length(out) == 0) {
    return(as.integer(default))
  }
  out
}

parse_numeric_csv_local <- function(value, default = numeric(0)) {
  out <- suppressWarnings(as.numeric(split_csv_local(value)))
  out <- out[is.finite(out)]
  if (length(out) == 0) {
    return(as.numeric(default))
  }
  out
}

parse_int_scalar_local <- function(value, default) {
  out <- suppressWarnings(as.integer(normalize_scalar_value(value)))
  if (length(out) == 0 || is.na(out)) {
    return(as.integer(default))
  }
  out
}

parse_num_scalar_local <- function(value, default) {
  out <- suppressWarnings(as.numeric(normalize_scalar_value(value)))
  if (length(out) == 0 || is.na(out)) {
    return(as.numeric(default))
  }
  out
}

parse_yes_no_local <- function(value, default = "yes") {
  out <- tolower(normalize_scalar_value(value, default))
  if (!out %in% c("yes", "no")) {
    out <- default
  }
  out
}

default_object_layer_config_df <- function(cfg) {
  data.frame(
    layer_id = c(cfg$panorama_layer_id, "subcluster_1", "subcluster_2"),
    layer_role = c("panorama", "subcluster", "subcluster"),
    enabled = c("yes", "no", "no"),
    parent_layer = c("", cfg$panorama_layer_id, cfg$panorama_layer_id),
    sample_include = c("", "", ""),
    sample_exclude = c("", "", ""),
    selection_column = c("", "", ""),
    selection_values = c("", "", ""),
    rebuild_normalization = c("yes", "yes", "yes"),
    hvg_nfeatures = c(cfg$hvg_nfeatures, cfg$hvg_nfeatures, cfg$hvg_nfeatures),
    pca_dims = c(cfg$pca_dims_panorama_raw, cfg$pca_dims_subcluster_raw, cfg$pca_dims_subcluster_raw),
    target_clusters = c(cfg$target_clusters, 8L, 6L),
    res_range = c(
      paste(cfg$res_range, collapse = ","),
      "0.10,0.15,0.20,0.25,0.30,0.35,0.40",
      "0.10,0.15,0.20,0.25,0.30,0.35,0.40"
    ),
    res_fine_step = c(cfg$res_fine_step_default, cfg$res_fine_step_default, cfg$res_fine_step_default),
    normalization_methods = c(cfg$normalization_methods_default, "lognorm", "lognorm"),
    integration_mode = c(cfg$integration_modes_default, cfg$integration_modes_default, cfg$integration_modes_default),
    vars_to_regress = c(cfg$vars_to_regress_default, "", ""),
    description = c(
      "Root panorama object built from all post-QC cells.",
      "Fill sample_include or selection_column/selection_values from panorama results before enabling this subcluster.",
      "Fill sample_include or selection_column/selection_values from panorama results before enabling this subcluster."
    ),
    stringsAsFactors = FALSE
  )
}

read_object_layer_config <- function(cfg) {
  df <- read_tsv_optional(cfg$object_layer_config_file)
  if (nrow(df) == 0) {
    df <- default_object_layer_config_df(cfg)
  }

  expected_cols <- c(
    "layer_id", "layer_role", "enabled", "parent_layer",
    "sample_include", "sample_exclude",
    "selection_column", "selection_values",
    "rebuild_normalization",
    "hvg_nfeatures", "pca_dims", "target_clusters",
    "res_range", "res_fine_step",
    "normalization_methods", "integration_mode", "vars_to_regress",
    "description"
  )
  for (col in expected_cols) {
    if (!col %in% colnames(df)) {
      df[[col]] <- ""
    }
  }

  # Backward-compatible defaults for older object_layers.tsv files.
  df$normalization_methods[!nzchar(trimws(df$normalization_methods))] <- cfg$normalization_methods_default
  df$integration_mode[!nzchar(trimws(df$integration_mode))] <- cfg$integration_modes_default
  df$vars_to_regress[!nzchar(trimws(df$vars_to_regress))] <- cfg$vars_to_regress_default
  df$res_fine_step[!nzchar(trimws(df$res_fine_step))] <- as.character(cfg$res_fine_step_default)

  df <- df[, expected_cols, drop = FALSE]
  df$layer_id <- vapply(df$layer_id, sanitize_layer_id_local, character(1))
  df$parent_layer <- vapply(df$parent_layer, sanitize_layer_id_local, character(1))
  df$layer_role <- tolower(vapply(df$layer_role, normalize_scalar_value, character(1), default = "subcluster"))
  df$enabled <- vapply(df$enabled, parse_yes_no_local, character(1), default = "yes")
  df$rebuild_normalization <- vapply(df$rebuild_normalization, parse_yes_no_local, character(1), default = "yes")
  df$description <- vapply(df$description, normalize_scalar_value, character(1))
  df <- df[nzchar(df$layer_id), , drop = FALSE]

  if (nrow(df) == 0) {
    stop("object_layers.tsv 没有有效 layer_id。", call. = FALSE)
  }

  df$sample_include <- I(lapply(df$sample_include, split_csv_local))
  df$sample_exclude <- I(lapply(df$sample_exclude, split_csv_local))
  df$selection_values <- I(lapply(df$selection_values, split_csv_local))
  df$pca_dims <- I(lapply(df$pca_dims, parse_index_spec_local, default = cfg$pca_dims))
  df$res_range <- I(lapply(df$res_range, parse_numeric_csv_local, default = cfg$res_range))
  df$normalization_methods <- I(lapply(df$normalization_methods, function(x) tolower(split_csv_local(x))))
  df$integration_mode <- I(lapply(df$integration_mode, function(x) tolower(split_csv_local(x))))
  df$vars_to_regress <- I(lapply(df$vars_to_regress, split_csv_local))
  df$hvg_nfeatures <- vapply(df$hvg_nfeatures, parse_int_scalar_local, integer(1), default = cfg$hvg_nfeatures)
  df$target_clusters <- vapply(df$target_clusters, parse_int_scalar_local, integer(1), default = cfg$target_clusters)
  df$res_fine_step <- vapply(df$res_fine_step, parse_num_scalar_local, numeric(1), default = cfg$res_fine_step_default)

  df
}

filter_enabled_layers <- function(df) {
  df[df$enabled == "yes", , drop = FALSE]
}

validate_layer_config <- function(df) {
  panorama <- df[df$layer_role == "panorama" & df$enabled == "yes", , drop = FALSE]
  if (nrow(panorama) != 1) {
    stop("object_layers.tsv 必须且只能启用一行 panorama。", call. = FALSE)
  }

  for (i in seq_len(nrow(df))) {
    row <- df[i, , drop = FALSE]
    bad_norm <- setdiff(row$normalization_methods[[1]], LAYER_NORMALIZATION_METHODS)
    if (length(bad_norm) > 0) {
      stop(sprintf("layer `%s` normalization_methods 非法: %s", row$layer_id, paste(bad_norm, collapse = ",")), call. = FALSE)
    }
    bad_integration <- setdiff(row$integration_mode[[1]], LAYER_INTEGRATION_METHODS)
    if (length(bad_integration) > 0) {
      stop(sprintf("layer `%s` integration_mode 非法: %s", row$layer_id, paste(bad_integration, collapse = ",")), call. = FALSE)
    }
    if (row$enabled == "yes" && row$layer_role == "subcluster") {
      if (!nzchar(row$parent_layer)) {
        stop(sprintf("subcluster layer `%s` 必须设置 parent_layer。", row$layer_id), call. = FALSE)
      }
      if (length(row$sample_include[[1]]) == 0 && !nzchar(row$selection_column)) {
        stop(sprintf("subcluster layer `%s` 至少需要 sample_include 或 selection_column。", row$layer_id), call. = FALSE)
      }
    }
  }
  df
}

topo_order_layers <- function(df) {
  df <- df[order(ifelse(df$layer_role == "panorama", 0L, 1L)), , drop = FALSE]
  ordered <- df[df$layer_role == "panorama", , drop = FALSE]
  remaining <- df[df$layer_role != "panorama", , drop = FALSE]
  while (nrow(remaining) > 0) {
    ready <- remaining$parent_layer %in% ordered$layer_id
    if (!any(ready)) {
      stop("object_layers.tsv parent_layer 形成循环或引用了不存在的 layer。", call. = FALSE)
    }
    ordered <- rbind(ordered, remaining[ready, , drop = FALSE])
    remaining <- remaining[!ready, , drop = FALSE]
  }
  ordered
}

layer_config_row_to_spec <- function(row) {
  if (nrow(row) != 1) {
    stop("layer_config_row_to_spec 需要单行 data.frame。", call. = FALSE)
  }
  list(
    layer_id = row$layer_id[[1]],
    layer_role = row$layer_role[[1]],
    enabled = identical(row$enabled[[1]], "yes"),
    parent_layer = row$parent_layer[[1]],
    sample_include = row$sample_include[[1]],
    sample_exclude = row$sample_exclude[[1]],
    selection_column = normalize_scalar_value(row$selection_column[[1]]),
    selection_values = row$selection_values[[1]],
    rebuild_normalization = identical(row$rebuild_normalization[[1]], "yes"),
    hvg_nfeatures = row$hvg_nfeatures[[1]],
    pca_dims = row$pca_dims[[1]],
    target_clusters = row$target_clusters[[1]],
    res_range = row$res_range[[1]],
    res_fine_step = row$res_fine_step[[1]],
    integration_mode = row$integration_mode[[1]],
    normalization_methods = row$normalization_methods[[1]],
    vars_to_regress = row$vars_to_regress[[1]],
    description = row$description[[1]]
  )
}

panorama_layer_spec <- function(cfg) {
  df <- validate_layer_config(read_object_layer_config(cfg))
  enabled <- filter_enabled_layers(df)
  row <- enabled[enabled$layer_role == "panorama", , drop = FALSE]
  layer_config_row_to_spec(row)
}
