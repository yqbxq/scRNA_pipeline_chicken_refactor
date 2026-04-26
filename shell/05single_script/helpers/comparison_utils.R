read_comparisons_with_subset <- function(cfg) {
  df <- read_tsv_optional(cfg$comparison_sheet)
  expected_cols <- c(
    "comparison_id", "ident_1", "ident_2", "enabled", "group_var", "batch_var",
    "layer_scope", "min_biological_replicates", "subset_column", "subset_value"
  )
  for (col in expected_cols) {
    if (!col %in% colnames(df)) {
      df[[col]] <- ""
    }
  }
  if (nrow(df) == 0) {
    return(df[, expected_cols, drop = FALSE])
  }

  df <- df[, expected_cols, drop = FALSE]
  for (col in setdiff(expected_cols, "min_biological_replicates")) {
    df[[col]] <- vapply(df[[col]], normalize_scalar_value, character(1))
  }
  df$enabled <- normalize_flag(df$enabled, "yes")
  df$min_biological_replicates <- suppressWarnings(as.integer(df$min_biological_replicates))
  df$min_biological_replicates[is.na(df$min_biological_replicates)] <- cfg$min_biological_replicates %||% 2L
  df$subset_column[is.na(df$subset_column)] <- ""
  df$subset_value[is.na(df$subset_value)] <- ""
  df
}

read_comparison_sheet_local <- function(cfg) {
  read_comparisons_with_subset(cfg)
}

comparison_applies_to_layer <- function(comparison_row, layer_id) {
  scope <- normalize_scalar_value(comparison_row$layer_scope[[1]], "*")
  if (!nzchar(scope) || identical(scope, "*")) {
    return(TRUE)
  }
  layer_id %in% split_csv_local(scope)
}

apply_subset_filter <- function(seu, comparison_row) {
  subset_column <- normalize_scalar_value(comparison_row$subset_column[[1]])
  subset_value <- normalize_scalar_value(comparison_row$subset_value[[1]])
  if (!nzchar(subset_column) || !nzchar(subset_value)) {
    return(seu)
  }
  if (!subset_column %in% colnames(seu@meta.data)) {
    warning(sprintf("comparison subset_column 不存在，跳过 subset: %s", subset_column), call. = FALSE)
    return(seu)
  }
  keep_values <- split_csv_local(subset_value)
  keep <- as.character(seu@meta.data[[subset_column]]) %in% keep_values
  if (!any(keep)) {
    warning(sprintf("comparison subset 没有匹配细胞: %s in %s", subset_value, subset_column), call. = FALSE)
    return(seu[, FALSE])
  }
  subset(seu, cells = colnames(seu)[keep])
}

apply_comparison_subset <- function(seu, comparison_row) {
  apply_subset_filter(seu, comparison_row)
}
